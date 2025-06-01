# Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
# Copyright (C) 2012-2025 team free-astro (see more in AUTHORS file)
# Reference site is https://siril.org
# SPDX-License-Identifier: GPL-3.0-or-later

"""
GPU helper module for Siril Python interface providing helper functions for detection,
installation and testing of GPU-related modules.
Initial scope is ONNX, torch and jax
"""

import re
import os
import sys
import json
import time
import platform
import tempfile
import importlib
import subprocess
from typing import Optional, Dict, Any
from packaging import version as pkg_version
import requests
import numpy as np
from .utility import ensure_installed, _check_package_installed, _install_package, \
                     SuppressedStderr

def _detect_cuda_version(system) -> Optional[str]:
    """
    Detects the CUDA version by parsing the output of 'nvcc -v'.

    Returns:
        Optional[str]: The CUDA version as a string (e.g., '11.7') if detected,
                        or "0.0" if nvcc is not installed or version cannot be determined.
    """
    try:
        nvcc_command = 'nvcc.exe' if system == 'windows' else 'nvcc'

        # Run nvcc -V and capture the output
        result = subprocess.run([nvcc_command, '-V'],
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                text=True,
                                check=False)

        output = result.stderr if result.stderr else result.stdout

        version_match = re.search(r'release (\d+\.\d+)', output)
        if version_match:
            return version_match.group(1)

        alt_match = re.search(r'V(\d+\.\d+\.\d+)', output)
        if alt_match:
            # Return just the major.minor part (e.g., 11.7 from 11.7.0)
            version_parts = alt_match.group(1).split('.')
            return f"{version_parts[0]}.{version_parts[1]}"

        return "0.0"
    except (subprocess.SubprocessError, FileNotFoundError):
        # nvcc is not installed or an error occurred
        try:
            # possibly we hve torch installed (this is a bit slow to import
            # which is why we don't try it first...)
            torch_helper = TorchHelper()
            if torch_helper.is_torch_installed():
                import torch
                cuda_version = torch.version.cuda  # e.g., '12.1'
                return cuda_version
            return "0.0"
        except Exception:
            return "0.0"

def _detect_nvidia_gpu(system):
    """
    Detect if NVIDIA GPU is available. (Required to check Windows and Linux, not required on MacOS.)

    Returns:
        bool: True if NVIDIA GPU is detected, False otherwise.
    """
    try:
        if system == "windows":
            # Windows detection using PowerShell
            output = subprocess.check_output(
                ["powershell", "-Command", "Get-WmiObject -Class Win32_VideoController"],
                text=True
            )
            return "NVIDIA" in output
        # Linux and macOS detection using lspci or similar
        if system == "linux":
            try:
                output = subprocess.check_output(["nvidia-smi"], stderr=subprocess.DEVNULL, text=True)
                return True
            except (subprocess.SubprocessError, FileNotFoundError):
                pass

            try:
                output = subprocess.check_output(["lspci"], text=True)
                return "NVIDIA" in output
            except (subprocess.SubprocessError, FileNotFoundError):
                return False
    except Exception:
        pass
    return False

def _detect_amd_gpu():
    """
    Detect if AMD GPU is available (Only needed on Linux, as on Windows the directml
    package is used.)

    Returns:
        bool: True if AMD GPU is detected, False otherwise.
    """
    try:
        # Try rocm-smi first
        try:
            output = subprocess.check_output(["rocm-smi"], stderr=subprocess.DEVNULL, text=True)
            return True
        except (subprocess.SubprocessError, FileNotFoundError):
            pass

        # Fallback to lspci
        try:
            output = subprocess.check_output(["lspci"], text=True)
            return any(gpu in output for gpu in ["AMD", "ATI", "Radeon"])
        except (subprocess.SubprocessError, FileNotFoundError):
            return False
    except Exception:
        pass
    return False

def _detect_intel_gpu_for_openvino():
    """
    Detect if an Intel GPU compatible with OpenVINO is available.

    Returns:
        bool: True if compatible Intel GPU is detected
    """
    # First try using OpenVINO's native detection if available
    try:
        import openvino as ov
        core = ov.Core()
        available_devices = core.available_devices
        return any(device.startswith("GPU") for device in available_devices)
    except ImportError:
        print("OpenVINO not installed, falling back to hardware detection")
    except Exception as e:
        print(f"Error using OpenVINO device detection: {e}")

    # Fall back to lspci detection if OpenVINO is not available
    try:
        output = subprocess.check_output(["lspci"], text=True).lower()

        # Intel integrated and discrete GPU keywords
        intel_keywords = ["intel", "graphics", "iris", "uhd", "hd graphics", "arc",
                        "a3", "a5", "a7", "a370m", "a730m", "a770", "a750",
                        "a580", "a380", "battlemage"]

        # Check for any Intel GPU
        return "intel" in output and any(keyword in output for keyword in intel_keywords)

    except (subprocess.SubprocessError, FileNotFoundError):
        return False
    except Exception as e:
        print(f"Error detecting Intel GPU: {e}")
        return False

class ONNXHelper:
    """
    A class to handle detection and installation of the appropriate ONNX Runtime
    package based on the system hardware and configuration.

    Example usage (this should be used instead of
    ``sirilpy.ensure_installed("onnxruntime")`` to install the correct package for
    the user's system.)

    .. code-block:: python

       installer = sirilpy.ONNXHelper()
       installer.install_onnxruntime()

    """

    def __init__(self):
        """Initialize the ONNXHelper."""
        ensure_installed("platformdirs")
        from platformdirs import user_data_dir
        self.system = platform.system().lower()
        self.providers = None
        self.config_file = os.path.join(user_data_dir(appname="siril"), "siril_onnx.conf")

    def create_simple_onnx_model(self):
        """Create a simple ONNX model with matrix multiplication and ReLU."""
        import onnx
        from onnx import helper, TensorProto

        input_shape = [1, 128, 256]
        weight_shape = [256, 512]

        input_tensor = helper.make_tensor_value_info('input', TensorProto.FLOAT, input_shape)
        output_tensor = helper.make_tensor_value_info('output', TensorProto.FLOAT, [1, 128, 512])

        weight_data = np.random.randn(*weight_shape).astype(np.float32)
        weight_tensor = helper.make_tensor('weight', TensorProto.FLOAT, weight_shape, weight_data.flatten())

        matmul_node = helper.make_node('MatMul', inputs=['input', 'weight'], outputs=['matmul_output'])
        relu_node = helper.make_node('Relu', inputs=['matmul_output'], outputs=['output'])

        graph = helper.make_graph(
            nodes=[matmul_node, relu_node],
            name='SimpleGraph',
            inputs=[input_tensor],
            outputs=[output_tensor],
            initializer=[weight_tensor]
        )

        model = helper.make_model(graph, producer_name='ep_test')
        model.opset_import[0].version = 11
        model.ir_version = 8
        onnx.checker.check_model(model)
        return model

    def try_provider(self, ort, model_path, input_data, provider, reference_output=None):
        """Try executing the model with a specific provider (no fallback)."""
        try:
            sess_options = ort.SessionOptions()
            session = ort.InferenceSession(
                model_path,
                sess_options=sess_options,
                providers=[provider]
            )

            actual_provider = session.get_providers()[0]
            if actual_provider != provider:
                print(f"✗ {provider}: fallback occurred (used {actual_provider})")
                return False

            output = session.run(None, {'input': input_data})

            print(f"✓ {provider} ran successfully")
            if reference_output is not None:
                if not np.allclose(reference_output, output[0], rtol=1e-3, atol=1e-5):
                    print("  ⚠ Output mismatch with CPU")
                    return False
            return True

        except Exception as e:
            print(f"✗ {provider} failed: {e}")
            return False

    def onnx_test(self, ort):
        """
        Test an imported onnxruntime.
        Args:
            ort: The ONNX runtime module to test
        Returns:
            list: a list of confirmed working ONNXRuntime ExecutionProviders in priority order
        """
        import onnx
        import os
        import tempfile

        print("=== ONNX Execution Provider Tester ===")
        print("Creating ONNX model...")
        model = self.create_simple_onnx_model()

        # Create temporary file for the model - Windows compatible approach
        temp_file = tempfile.NamedTemporaryFile(suffix='.onnx', delete=False)
        model_path = temp_file.name
        temp_file.close()  # Close the file handle immediately

        try:
            onnx.save(model, model_path)
            print(f"Model saved to temporary file: {model_path}")

            input_data = np.random.randn(1, 128, 256).astype(np.float32)
            print("\nRunning reference on CPU...")
            cpu_output = None

            try:
                cpu_session = ort.InferenceSession(model_path, providers=['CPUExecutionProvider'])
                cpu_output = cpu_session.run(None, {'input': input_data})[0]
                print("✓ CPU output computed.")
            except Exception as e:
                print(f"✗ Failed to run on CPU: {e}")
                return []

            all_providers = ort.get_available_providers()
            print("\nAvailable execution providers:")
            for p in all_providers:
                print(f"  - {p}")

            print("\nTesting each provider without fallback...")
            working_providers = []

            for provider in all_providers:
                print(f"\nTesting {provider}...")
                with SuppressedStderr():
                    if self.try_provider(ort, model_path, input_data, provider, reference_output=cpu_output):
                        working_providers.append(provider)

            print("\n=== Summary ===")
            if working_providers:
                print("✓ Working providers (in priority order):")
                for p in working_providers:
                    print(f"  - {p}")
                print(f"\n→ Best available provider: {working_providers[0]}")
                # Save to cache
                self._save_providers_to_cache(working_providers)
                self.providers = working_providers
                return working_providers

            print("✗ No execution providers were able to run the model.")
            return []

        finally:
            # Clean up the temporary file
            try:
                os.unlink(model_path)
            except (OSError, FileNotFoundError):
                # File might already be deleted or inaccessible, ignore
                pass

    def run(self, session, model_path, output_names, input_feed, run_options=None, return_first_output=False):
        """
        Run inference with automatic CPU fallback if the session fails.

        Args:
            session: The ONNX runtime inference session
            model_path (str): Path to the ONNX model file (needed for CPU fallback)
            output_names: Names of the outputs to compute, or None for all outputs
            input_feed: Dictionary mapping input names to input tensors
            run_options: Optional run options for the inference session
            return_first_output (bool): If True, return only the first output instead of the full list

        Returns:
            tuple: (result, session) where result is the inference output (or first output if return_first_output=True) and
                   session is the (potentially updated) inference session
        """
        import onnxruntime

        try:
            # Try running with the provided session
            result = session.run(output_names, input_feed, run_options)
            if return_first_output:
                result = result[0]
            return result, session
        except Exception:
            print("Warning: falling back to CPU.")

            # Create a new CPU-only session
            providers = ['CPUExecutionProvider']
            cpu_session = onnxruntime.InferenceSession(model_path, providers=providers)

            # Run with the CPU session
            result = cpu_session.run(output_names, input_feed, run_options)
            if return_first_output:
                result = result[0]
            return result, cpu_session

    def install_onnxruntime(self):
        """
        Detect system configuration and install the appropriate ONNX Runtime package.

        Returns:
            bool: True if installation was successful or already installed, False otherwise.

        Raises:
            TimooutError: if a TimeoutError occurs in ensure_installed() - this avoids falling
                        back to the CPU-only package purely because of network issues
        """
        # First check if any onnxruntime is already installed
        is_installed, package_name = self.check_onnxruntime_installed()

        if is_installed:
            print(f"ONNX Runtime is already installed: {package_name}")
            return True

        # If not installed, get recommended package
        onnxruntime_pkg, from_url, index_url = self.get_onnxruntime_package()

        # Check if the package exists
        if not self._check_onnxruntime_availability(onnxruntime_pkg):
            print(f"Package {onnxruntime_pkg} not found. Falling back to default onnxruntime.")
            onnxruntime_pkg = "onnxruntime"

        # Install the package
        try:
            if not _check_package_installed(onnxruntime_pkg):
                _install_package(onnxruntime_pkg, None, from_url=from_url, index_url=index_url)
                # For openvino we need to set the runtime environment:
                # https://github.com/intel/onnxruntime/releases/tag/v5.6
                if self.system == 'windows' and onnxruntime_pkg == 'onnxruntime-openvino':
                    import onnxruntime.tools.add_openvino_win_libs as utils
                    utils.add_openvino_libs_to_path()
        except TimeoutError as e:
            print(f"Failed to install {onnxruntime_pkg}: timeout error {str(e)}")
            raise TimeoutError("Error: timeout in install_onnxruntime()") from e

        except Exception as e:
            print(f"Failed to install {onnxruntime_pkg}: {str(e)}")
            print("Falling back to default onnxruntime package.")
            try:
                ensure_installed("onnxruntime")
            except Exception as err:
                print(f"Failed to install default onnxruntime: {str(err)}")
                return False
        return True

    def check_onnxruntime_installed(self):
        """
        Check if any onnxruntime package is already installed and usable.

        Returns:
            tuple: (is_installed, package_name) where package_name could be
                  'onnxruntime', 'onnxruntime-gpu', 'onnxruntime-silicon', etc.
        """
        try:
            # Try importing onnxruntime - if this fails, the package is not usable
            import onnxruntime

        except ImportError as e:
            # If import fails, package is not usable regardless of pip list
            # One of the error messages is a clue that MSVC runtime need updating
            if platform.system().lower() == "windows" and \
                                "DLL load failed" in str(e):
                print("DLL load failed. This means you need to update system "
                    "libraries. Usually updating Microsoft Visual C++ Runtime "
                    "will solve the issue.", file=sys.stderr)
            return False, None
        except Exception as e:
            print(f"Error checking for installed onnxruntime: {e}")
            return False, None

        # If we get here, some version of onnxruntime is installed and working
        package_name = "onnxruntime"  # Default assumption

        # Check provider information to determine specific package variant
        providers = onnxruntime.get_available_providers()
        print(f"Detected ONNX Runtime with providers: {providers}")

        # Check for specific provider patterns
        if any(p for p in providers if "CUDA" in p or "GPU" in p):
            package_name = "onnxruntime-gpu"
        elif any(p for p in providers if "DirectML" in p):
            package_name = "onnxruntime-directml"
        elif any(p for p in providers if "ROCm" in p):
            package_name = "onnxruntime-rocm"
        elif any(p for p in providers if "OpenVINO" in p or "DML" in p):
            package_name = "onnxruntime-intel"

        return True, package_name

    def get_onnxruntime_package(self):
        """
        Determine the appropriate ONNX Runtime package based on system and available hardware.

        Returns:
            tuple: (package_name, url) where url is None except for special cases like ROCm
        """
        from_url = None
        index_url = None
        onnxruntime_pkg = "onnxruntime"
        cuda_version = _detect_cuda_version(self.system)
        if self.system == "windows":
                # NVidia GPU runtime that supports CUDA
            if _detect_nvidia_gpu(self.system) and pkg_version.Version(cuda_version).major >= 11:
                onnxruntime_pkg = "onnxruntime-gpu"
                if pkg_version.Version(cuda_version).major == 11:
                    index_url = "https://aiinfra.pkgs.visualstudio.com/PublicPackages/_packaging/onnxruntime-cuda-11/pypi/simple/"
            else:
                # DirectML provides hardware acceleration for various GPUs on Windows
                onnxruntime_pkg = "onnxruntime-directml"

        elif self.system == "darwin":  # macOS
            onnxruntime_pkg = "onnxruntime"

        elif self.system == "linux":
            if _detect_nvidia_gpu(self.system) and pkg_version.Version(cuda_version).major >= 11:
                onnxruntime_pkg = "onnxruntime-gpu"
                if pkg_version.Version(cuda_version).major == 11:
                    index_url = "https://aiinfra.pkgs.visualstudio.com/PublicPackages/_packaging/onnxruntime-cuda-11/pypi/simple/"
            elif _detect_amd_gpu():
                # ROCm support for AMD GPUs on Linux
                onnxruntime_pkg = "onnxruntime-rocm"
                from_url = "https://repo.radeon.com/rocm/manylinux/rocm-rel-6.4/"
            elif _detect_intel_gpu_for_openvino():
                onnxruntime_pkg = "onnxruntime-openvino"

        return onnxruntime_pkg, from_url, index_url

    def _check_onnxruntime_availability(self, package_name):
        """
        Check if the specified ONNX Runtime package is available on PyPI.

        Args:
            package_name (str): Package name to check

        Returns:
            bool: True if the package is available, False otherwise.
        """
        try:
            output = subprocess.check_output(
                [sys.executable, "-m", "pip", "index", "versions", package_name],
                stderr=subprocess.DEVNULL,
                text=True
            )
            return "No matching distribution found" not in output
        except subprocess.SubprocessError:
            # If the command fails, check directly from PyPI
            try:
                url = f"https://pypi.org/pypi/{package_name}/json"
                response = requests.get(url)
                return response.status_code == 200
            except Exception:
                return False

    def get_execution_providers_ordered(self, ai_gpu_acceleration=True):
        """
        Get execution providers ordered by priority.
        This function returns a list of available ONNX Runtime execution providers
        in a reasonable order of priority, covering major GPU platforms:
        The CPU provider is always included as the final fallback option.
        Args:
            ai_gpu_acceleration (bool): Whether to include GPU acceleration providers.
                                        Defaults to True.
        Returns:
            list: Ordered list of available execution providers.
        """
        if ai_gpu_acceleration is False:
            return ["CPUExecutionProvider"]

        if self.providers is not None:
            return self.providers

        import onnxruntime as ort

        # Try to load cached providers first
        cached_providers = self._load_cached_providers(ort)
        if cached_providers:
            self.providers = cached_providers
            return self.providers

        # If no valid cache, run the test
        return self.onnx_test(ort)

    def _load_cached_providers(self, ort):
        """
        Load cached execution providers from config file if they're still valid.
        Args:
            ort: The ONNX runtime module
        Returns:
            list or None: List of cached providers if valid, None otherwise
        """
        if not os.path.exists(self.config_file):
            return None

        try:
            with open(self.config_file, 'r', encoding='utf-8') as f:
                cached_data = json.load(f)

            cached_providers = cached_data.get('execution_providers', [])
            if not cached_providers:
                return None

            # Check if all cached providers are still available
            available_providers = set(ort.get_available_providers())
            valid_providers = [p for p in cached_providers if p in available_providers]

            # Only use cache if all providers are still available
            if len(valid_providers) == len(cached_providers):
                print(f"Using cached execution providers from {self.config_file}")
                return valid_providers
            print("Cached providers outdated, will re-test")
            return None

        except (json.JSONDecodeError, IOError, KeyError) as e:
            print(f"Failed to load cached providers: {e}", file=sys.stderr)
            return None

    def _save_providers_to_cache(self, providers):
        """
        Save execution providers to config file.
        Args:
            providers (list): List of working execution providers
        """
        try:
            # Ensure directory exists
            os.makedirs(os.path.dirname(self.config_file), exist_ok=True)

            cache_data = {
                'execution_providers': providers,
                'system': self.system,
                'cached_at': platform.platform()
            }

            with open(self.config_file, 'w', encoding='utf-8') as f:
                json.dump(cache_data, f, indent=2)

            print(f"Cached execution providers to {self.config_file}")

        except IOError as e:
            print(f"Failed to save providers cache: {e}", file=sys.stderr)

    def uninstall_onnxruntime(self):
        """
        Detects and uninstalls all variants of onnxruntime packages.
        Checks for any package starting with 'onnxruntime'.

        Returns:
            list: A list of uninstalled packages
        """
        # Get all installed packages
        try:
            result = subprocess.run(
                [sys.executable, "-m", "pip", "list"],
                capture_output=True,
                text=True,
                check=True
            )
            installed_packages = result.stdout.splitlines()
        except subprocess.CalledProcessError as e:
            print(f"Error getting installed packages: {e}")
            return []

        # Find all packages that start with 'onnxruntime'
        onnx_packages = []
        for line in installed_packages:
            parts = line.split()
            if parts and parts[0].lower().startswith('onnxruntime'):
                onnx_packages.append(parts[0])

        # Uninstall found packages
        if not onnx_packages:
            print("No onnxruntime packages found.")
            return []

        print(f"Found onnxruntime packages: {', '.join(onnx_packages)}")
        uninstalled = []

        for package in onnx_packages:
            print(f"Uninstalling {package}...")
            try:
                subprocess.run(
                    [sys.executable, "-m", "pip", "uninstall", "-y", package],
                    check=True
                )
                uninstalled.append(package)
                print(f"Successfully uninstalled {package}")
            except subprocess.CalledProcessError:
                print(f"Failed to uninstall {package}")

        return uninstalled

class TorchHelper:
    """Helper class for PyTorch GPU testing and benchmarking with optional installation."""

    def __init__(self):
        """Initialize TorchHelper without importing torch at module level."""
        self.torch: Optional[Any] = None
        self.nn: Optional[Any] = None
        self.F: Optional[Any] = None
        self.np: Optional[Any] = None
        self.device_info: Optional[dict] = None
        self._torch_installed = False

    def is_torch_installed(self) -> bool:
        """Check if PyTorch is installed without importing it."""
        if self._torch_installed:
            return True

        # Check if torch is available
        torch_spec = importlib.util.find_spec("torch")
        if torch_spec is not None:
            self._torch_installed = True
            return True
        return False

    def install_torch(self, cuda_version: str = "cu126", force_reinstall: bool = False) -> bool:
        """
        Install PyTorch with CUDA support.

        Args:
            cuda_version: compute platform to install (e.g., 'cu118', 'cu126', 'rocm', 'cpu')
            force_reinstall: Whether to reinstall even if already installed

        Returns:
            bool: True if installation successful, False otherwise
        """
        if self.is_torch_installed() and not force_reinstall:
            print("PyTorch is already installed. Use force_reinstall=True to reinstall.")
            return self._import_torch()

        print(f"Installing PyTorch for compute platform {cuda_version}...")
        host = platform.system().lower()
        try:
            if host == 'darwin' or (host == 'windows' and cuda_version == 'cpu'):
                _install_package("torch")
            elif cuda_version.lower() == 'cpu':
                _install_package("torch", index_url="https://download.pytorch.org/whl/cpu")
            else:
                _install_package("torch", index_url=f"https://download.pytorch.org/whl/{cuda_version}")

            # Try to import after installation
            return self._import_torch()

        except subprocess.CalledProcessError as e:
            print(f"✗ PyTorch installation failed: {e}")
            if e.stderr:
                print("Error output:", e.stderr[-500:])  # Show last 500 chars
            return False
        except Exception as e:
            print(f"✗ Unexpected error during installation: {e}")
            return False

    def _import_torch(self) -> bool:
        """Import torch modules and set instance variables."""
        try:
            import torch
            from torch import nn
            from torch.nn import functional as F

            self.torch = torch
            self.nn = nn
            self.F = F
            self.np = np
            self._torch_installed = True

            # Get device info after successful import
            self.device_info = self._get_device_info()

            print(f"✓ PyTorch {torch.__version__} imported successfully!")
            return True

        except ImportError as e:
            print(f"✗ Failed to import PyTorch: {e}")
            return False
        except Exception as e:
            print(f"✗ Unexpected error importing PyTorch: {e}")
            return False

    def ensure_torch(self, cuda_version: str = "cu126") -> bool:
        """
        Ensure PyTorch is available. Optionally install if not found.

        Args:
            auto_install: Whether to automatically install PyTorch if not found
            cuda_version: CUDA version to install if auto_install is True

        Returns:
            bool: True if PyTorch is available, False otherwise
        """
        if self.torch is not None:
            return True

        # check it imports correctly
        if self.is_torch_installed() and self._import_torch():
            return True

        print("PyTorch not found. Attempting automatic installation...")
        return self.install_torch(cuda_version=cuda_version)

    def _get_device_info(self) -> dict:
        """Get information about available devices."""
        if self.torch is None:
            return {}

        info = {
            'pytorch_version': self.torch.__version__,
            'cuda_available': self.torch.cuda.is_available(),
            'cuda_version': self.torch.version.cuda if self.torch.cuda.is_available() else None,
            'gpu_count': self.torch.cuda.device_count() if self.torch.cuda.is_available() else 0,
            'gpu_names': []
        }

        if self.torch.cuda.is_available():
            for i in range(self.torch.cuda.device_count()):
                info['gpu_names'].append(self.torch.cuda.get_device_name(i))

        return info

    def _create_simple_model(self, input_dim=256, hidden_dim=512, output_dim=128):
        """Create a simple neural network model."""
        if self.torch is None or self.nn is None:
            raise RuntimeError("PyTorch not available. Call ensure_torch() first.")

        nn = self.nn
        F = self.torch.nn.functional

        class SimpleModel(nn.Module):
            """ Simple model to test torch can utilise a host """
            def __init__(self, input_dim, hidden_dim, output_dim):
                super().__init__()
                self.linear1 = nn.Linear(input_dim, hidden_dim)
                self.linear2 = nn.Linear(hidden_dim, hidden_dim)
                self.linear3 = nn.Linear(hidden_dim, output_dim)
                self.dropout = nn.Dropout(0.1)

            def forward(self, x):
                """ Simple function within the test model """
                x = F.relu(self.linear1(x))
                x = self.dropout(x)
                x = F.relu(self.linear2(x))
                x = self.linear3(x)
                return x

        return SimpleModel(input_dim, hidden_dim, output_dim)

    def create_test_data(self, batch_size=32, input_dim=256):
        """Create random test input data."""
        if self.torch is None:
            raise RuntimeError("PyTorch not available. Call ensure_torch() first.")
        return self.torch.randn(batch_size, input_dim)

    def benchmark_model(self, model, input_data, device, num_runs=10):
        """Benchmark model execution on specified device."""
        if self.torch is None:
            raise RuntimeError("PyTorch not available. Call ensure_torch() first.")

        # Set model to evaluation mode to disable dropout
        model.eval()
        model = model.to(device)
        input_data = input_data.to(device)

        # Warm up
        for _ in range(3):
            with self.torch.no_grad():
                _ = model(input_data)

        # Benchmark
        if device.type == 'cuda':
            self.torch.cuda.synchronize()
        start_time = time.time()

        for _ in range(num_runs):
            with self.torch.no_grad():
                output = model(input_data)

        if device.type == 'cuda':
            self.torch.cuda.synchronize()
        end_time = time.time()

        avg_time = (end_time - start_time) / num_runs
        return output, avg_time

    def print_device_info(self):
        """Print device information."""
        if not self.ensure_torch():
            print("PyTorch not available - cannot display device info")
            return

        print("=== Device Information ===")
        print(f"PyTorch version: {self.device_info['pytorch_version']}")
        print(f"CUDA available: {self.device_info['cuda_available']}")

        if self.device_info['cuda_available']:
            print(f"CUDA version: {self.device_info['cuda_version']}")
            print(f"GPU count: {self.device_info['gpu_count']}")
            for i, name in enumerate(self.device_info['gpu_names']):
                print(f"  GPU {i}: {name}")

    def test_torch_gpu(self):
        """Test PyTorch model execution on GPU."""
        print("=== PyTorch GPU Test ===")

        if not self.ensure_torch():
            print("PyTorch not available. Please install it first using install_torch()")
            return False

        # Print device info
        self.print_device_info()

        # Create model and test data with fixed seeds for reproducible results
        print("\nCreating model and test data...")
        self.torch.manual_seed(42)  # Set seed for reproducible results
        model = self._create_simple_model()
        self.torch.manual_seed(42)  # Reset seed for consistent input data
        input_data = self.create_test_data(batch_size=64, input_dim=256)

        print(f"Model parameters: {sum(p.numel() for p in model.parameters()):,}")
        print(f"Input shape: {input_data.shape}")

        # Test CPU execution
        print("\nTesting CPU execution (for accuracy reference)...")
        cpu_device = self.torch.device('cpu')
        model_cpu = self._create_simple_model()  # Create fresh model
        self.torch.manual_seed(42)  # Ensure same initialization
        model_cpu.load_state_dict(model.state_dict())  # Copy weights
        cpu_output, _ = self.benchmark_model(model_cpu, input_data, cpu_device)

        # Test GPU execution
        if self.torch.cuda.is_available():
            print("\nTesting GPU execution...")
            gpu_device = self.torch.device('cuda:0')

            try:
                model_gpu = self._create_simple_model()  # Create fresh model for GPU
                model_gpu.load_state_dict(model.state_dict())  # Copy same weights
                gpu_output, _ = self.benchmark_model(model_gpu, input_data, gpu_device)
                print("✓ GPU execution successful!")

                # Compare outputs using appropriate tolerance for CPU vs GPU
                cpu_output_np = cpu_output.numpy()
                gpu_output_np = gpu_output.cpu().numpy()

                # Check if outputs are close (accounting for floating point differences)
                are_close = self.np.allclose(cpu_output_np, gpu_output_np, rtol=1e-3, atol=1e-4)

                if are_close:
                    print("✓ CPU and GPU outputs match within tolerance!")
                else:
                    print("⚠ CPU and GPU outputs differ more than expected")

            except Exception as e:
                print(f"✗ GPU execution failed: {e}")
        else:
            print("\n⚠ CUDA not available - cannot test GPU execution")
            print("To enable GPU support:")
            print("  1. Install CUDA toolkit")
            print("  2. Use install_torch() with appropriate CUDA version")

        print("\n" + "="*50)
        return True

    def test_tensor_operations(self):
        """Test basic tensor operations on GPU."""
        print("=== Tensor Operations Test ===")

        if not self.ensure_torch():
            print("PyTorch not available. Please install it first using install_torch()")
            return False

        if not self.torch.cuda.is_available():
            print("⚠ CUDA not available - skipping tensor operations test")
            return False

        # Create large tensors for meaningful GPU test with fixed seed
        size = 2048
        print(f"Creating {size}x{size} tensors...")

        # Set seed for reproducible random tensors
        self.torch.manual_seed(123)
        a_cpu = self.torch.randn(size, size)
        b_cpu = self.torch.randn(size, size)

        # GPU tensors (copy the same data)
        a_gpu = a_cpu.clone().cuda()
        b_gpu = b_cpu.clone().cuda()

        # CPU matrix multiplication
        start_time = time.time()
        c_cpu = self.torch.mm(a_cpu, b_cpu)
        cpu_time = time.time() - start_time

        # GPU matrix multiplication
        self.torch.cuda.synchronize()
        start_time = time.time()
        c_gpu = self.torch.mm(a_gpu, b_gpu)
        self.torch.cuda.synchronize()
        gpu_time = time.time() - start_time

        # Verify results using appropriate tolerance for large matrix operations
        cpu_np = c_cpu.numpy()
        gpu_np = c_gpu.cpu().numpy()

        # Use more appropriate tolerances for large matrix multiplication
        are_close = self.np.allclose(cpu_np, gpu_np, rtol=1e-3, atol=1e-3)

        if are_close:
            print("✓ Results match within tolerance!")
        else:
            print("⚠ Results differ more than expected")

        print(f"CPU time: {cpu_time:.4f}s")
        print(f"GPU time: {gpu_time:.4f}s")
        print(f"Speedup: {cpu_time/gpu_time:.2f}x")

        print("\n" + "="*50)
        return True

    def run_all_tests(self, cuda_version: str = "cu126"):
        """
        Run all available tests.

        Args:
            auto_install: Whether to automatically install PyTorch if not found
            cuda_version: CUDA version to install if auto_install is True
        """
        if not self.ensure_torch(cuda_version=cuda_version):
            print("Cannot run tests - PyTorch not available")
            print("Use helper.install_torch() or helper.run_all_tests(auto_install=True)")
            return False

        self.test_torch_gpu()
        self.test_tensor_operations()
        return True

class JaxHelper:
    """
    A helper class for detecting, installing, and testing JAX with appropriate hardware acceleration.

    This class automatically detects the system configuration and installs the correct JAX variant
    (CPU, CUDA, ROCm, etc.) based on available hardware.
    """

    def __init__(self):
        self.system = platform.system().lower()
        self.jax_installed = False
        self.detected_config = None

    def detect_hardware_config(self) -> Dict[str, Any]:
        """
        Detect the hardware configuration and determine the appropriate JAX variant.

        Returns:
            Dict containing detected hardware info and recommended JAX installation.
        """
        config = {
            'system': self.system,
            'has_nvidia_gpu': False,
            'has_amd_gpu': False,
            'has_intel_gpu': False,
            'cuda_version': None,
            'recommended_jax_variant': 'jax[cpu]',
            'install_url': None,
            'index_url': None
        }

        # Detect GPUs
        config['has_nvidia_gpu'] = _detect_nvidia_gpu(self.system)
        config['has_amd_gpu'] = _detect_amd_gpu()
        config['has_intel_gpu'] = _detect_intel_gpu_for_openvino()

        # Detect CUDA version if NVIDIA GPU is present
        if config['has_nvidia_gpu']:
            config['cuda_version'] = _detect_cuda_version(self.system)

        # Determine JAX variant based on hardware
        config = self._determine_jax_variant(config)

        self.detected_config = config
        return config

    def _determine_jax_variant(self, config: Dict[str, Any], force_cpu=False) -> Dict[str, Any]:
        """
        Determine the appropriate JAX variant based on detected hardware.

        Args:
            config: Hardware configuration dictionary
            force_cpu: forces CPU-only installation

        Returns:
            Updated configuration with JAX variant recommendation
        """
        if force_cpu:
            config['recommended_jax_variant'] = 'jax[cpu]'
            print("⚠ Warning: performance of the CPU-only jax variant "
                "is slow and is intended for development use only. It "
                "is generally recommended that you do NOT enable jax "
                "optimisation in scripts that offer it!")
            return config

        if config['has_nvidia_gpu']:
            if config['system'] == 'linux':
                # NVIDIA GPU with CUDA
                config['recommended_jax_variant'] = 'jax[cuda12]'
            elif config['system'] == 'windows':
                print("⚠ Windows CUDA support for jax is experimental: "
                    "if you have problems, reinstall forcing CPU support")

        elif config['has_amd_gpu'] and config['system'] == 'linux':
            # AMD GPU with ROCm (Linux only)
            config['recommended_jax_variant'] = 'jax[rocm]'

        elif config['has_intel_gpu']:
            # Current pypi package for intel-extension-for-openxla
            # has dependency clashes (still requires numpy 1.x)
            print("⚠ Intel GPU detected: unfortunately the experimental "
                "jax plugin for Intel GPUs has dependency clashes and "
                "still requires numpy 1.x therefore only the CPU variant "
                "can be installed. Hopefully this will change in the "
                "future.")
            config['recommended_jax_variant'] = 'jax[cpu]'

        elif config['system'] == 'darwin':
            # macOS - use Metal wheel and print an experimental warning
            print("⚠ jax Metal support on MacOS is still experimental")
            config['recommended_jax_variant'] = 'jax-metal'

        else:
            # Default to CPU
            config['recommended_jax_variant'] = 'jax[cpu]'

        return config

    def install_jax(self, force_variant: Optional[str] = None,
                   version_constraint: Optional[str] = None) -> bool:
        """
        Install JAX with the appropriate variant for the detected hardware.

        Args:
            force_variant: Override auto-detection with specific variant (e.g., 'jax[cpu]')
            version_constraint: Version constraint string (e.g., '>=0.4.0')

        Returns:
            bool: True if installation succeeded, False otherwise
        """
        if not self.detected_config:
            self.detect_hardware_config()

        variant = force_variant or self.detected_config['recommended_jax_variant']

        try:
            print(f"Installing {variant}...")

            # Use the provided install_package function
            _install_package(
                package_name=variant,
                version_constraint=version_constraint,
                from_url=self.detected_config.get('install_url'),
                index_url=self.detected_config.get('index_url')
            )

            self.jax_installed = True
            print(f"Successfully installed {variant}")
            return True

        except Exception as e:
            print(f"Failed to install {variant}: {e}")
            return False

    def test_jax_installation(self) -> str:
        """
        Test JAX functionality and return execution provider.

        Returns:
            str: "gpu" if JAX is using GPU, "cpu" if using CPU

        Raises:
            RuntimeError: If JAX computation fails or accuracy check fails
            ImportError: If JAX is not installed
        """
        try:
            import jax
            import jax.numpy as jnp
        except ImportError as e:
            raise ImportError(f"JAX is not installed or not importable: {e}") from e

        try:
            # Get the default device
            default_device = jax.devices()[0]

            # Create test data with fixed seed for reproducible results
            np.random.seed(42)
            test_data_np = np.random.randn(100, 100).astype(np.float32)

            # Define a JIT-compiled function to test compilation
            @jax.jit
            def test_computation_jax(x):
                # Test matrix operations that benefit from GPU acceleration
                y = jnp.dot(x, x.T)
                z = jnp.sin(y) + jnp.cos(y)
                return jnp.sum(z)

            # Define equivalent numpy computation for cross-check
            def test_computation_numpy(x):
                y = np.dot(x, x.T)
                z = np.sin(y) + np.cos(y)
                return np.sum(z)

            # Convert to JAX array
            test_data_jax = jnp.array(test_data_np)

            # Execute JAX computation (JIT-compiled)
            jax_result = test_computation_jax(test_data_jax)
            jax_result.block_until_ready()

            # Execute numpy computation for cross-check
            numpy_result = test_computation_numpy(test_data_np)

            # Cross-check accuracy with appropriate tolerance
            # Using rtol=1e-5, atol=1e-6 for float32 precision
            if not np.allclose(jax_result, numpy_result, rtol=1e-5, atol=1e-6):
                raise RuntimeError(f"JAX result {jax_result} does not match numpy result {numpy_result} within tolerance")
            print("Accuracy cross-check between jax.numpy and numpy succeeded")

            # Check the platform of the device that was actually used
            if default_device.platform.lower() == 'gpu':
                return "gpu"
            return "cpu"

        except Exception as e:
            # If main test fails, try fallback to basic CPU operations
            if "does not match numpy" in str(e):
                # Re-raise accuracy errors immediately
                raise

            try:
                # Try basic CPU operations with accuracy check
                test_array_np = np.array([1.0, 2.0, 3.0], dtype=np.float32)
                test_array_jax = jnp.array(test_array_np)

                jax_sum = jnp.sum(test_array_jax)
                jax_sum.block_until_ready()
                numpy_sum = np.sum(test_array_np)

                if not np.allclose(jax_sum, numpy_sum, rtol=1e-7, atol=1e-8):
                    raise RuntimeError(f"Basic JAX operation failed accuracy check: {jax_sum} vs {numpy_sum}") from e

                return "cpu"
            except Exception:
                # If even CPU fails, something is seriously wrong
                raise RuntimeError("JAX is not functioning properly") from e

    def get_jax_info(self) -> Dict[str, Any]:
        """
        Get information about the current JAX installation.

        Returns:
            Dict containing JAX version, devices, and backend info
        """
        try:
            import jax
            return {
                'version': jax.__version__,
                'devices': [str(device) for device in jax.devices()],
                'default_backend': jax.default_backend(),
                'available_backends': list(jax.lib.xla_bridge.get_backend_names())
            }
        except ImportError:
            return {'error': 'JAX not installed'}
        except Exception as e:
            return {'error': f'Error getting JAX info: {e}'}

    def setup_jax(self, force_variant: Optional[str] = None,
                  version_constraint: Optional[str] = None,
                  test_after_install: bool = True) -> Dict[str, Any]:
        """
        Complete setup: detect hardware, install JAX, and optionally test.

        Args:
            force_variant: Override auto-detection with specific variant
            version_constraint: Version constraint for JAX installation
            test_after_install: Whether to test JAX after installation

        Returns:
            Dict containing setup results and information
        """
        results = {
            'detection_successful': False,
            'installation_successful': False,
            'test_successful': False,
            'execution_provider': None,
            'detected_config': None,
            'jax_info': None,
            'errors': []
        }

        try:
            # Step 1: Detect hardware configuration
            print("Detecting hardware configuration...")
            config = self.detect_hardware_config()
            results['detected_config'] = config
            results['detection_successful'] = True
            print(f"Detected configuration: {config['recommended_jax_variant']}")

            # Step 2: Install JAX
            print("Installing JAX...")
            install_success = self.install_jax(force_variant, version_constraint)
            results['installation_successful'] = install_success

            if install_success:
                # Step 3: Get JAX info
                jax_info = self.get_jax_info()
                results['jax_info'] = jax_info

                # Step 4: Test if requested
                if test_after_install:
                    print("Testing JAX installation...")
                    try:
                        execution_provider = self.test_jax_installation()
                        results['test_successful'] = True
                        results['execution_provider'] = execution_provider
                        print(f"JAX test successful! Using: {execution_provider}")
                    except Exception as e:
                        results['errors'].append(f"Test failed: {e}")
                        print(f"JAX test failed: {e}")
            else:
                results['errors'].append("JAX installation failed")

        except Exception as e:
            results['errors'].append(f"Setup error: {e}")
            print(f"Setup error: {e}")

        return results

    def detect_and_uninstall_jax(self, dry_run: bool = False) -> Dict[str, Any]:
        """
        Detect and uninstall any existing JAX-related packages.

        This is useful when you need to clean up a problematic JAX installation
        before installing a different variant (e.g., falling back from GPU to CPU).

        Args:
            dry_run: If True, only detect packages without uninstalling them

        Returns:
            Dict containing information about detected and uninstalled packages
        """
        results = {
            'detected_packages': [],
            'uninstalled_packages': [],
            'errors': [],
            'dry_run': dry_run
        }

        # Common JAX-related package patterns
        jax_packages = [
            'jax',
            'jaxlib',
            'jax-cuda',
            'jax-cuda11-local',
            'jax-cuda12-local',
            'jax-rocm',
            'jax-metal',
            'intel-extension-for-openxla'
        ]

        try:
            # Get list of installed packages
            result = subprocess.run(
                [sys.executable, '-m', 'pip', 'list', '--format=freeze'],
                capture_output=True,
                text=True,
                check=True
            )

            installed_packages = result.stdout.strip().split('\n')
            installed_dict = {}

            for package_line in installed_packages:
                if '==' in package_line:
                    name, version = package_line.split('==', 1)
                    installed_dict[name.lower()] = version

            # Find JAX-related packages
            detected = []
            for pkg_name in jax_packages:
                if pkg_name.lower() in installed_dict:
                    detected.append({
                        'name': pkg_name,
                        'version': installed_dict[pkg_name.lower()],
                        'installed_name': pkg_name.lower()
                    })

            # Also check for packages that start with 'jax'
            for installed_name in installed_dict.keys():
                if installed_name.startswith('jax') and installed_name not in [p['installed_name'] for p in detected]:
                    detected.append({
                        'name': installed_name,
                        'version': installed_dict[installed_name],
                        'installed_name': installed_name
                    })

            results['detected_packages'] = detected

            if detected:
                print(f"Found {len(detected)} JAX-related packages:")
                for pkg in detected:
                    print(f"  - {pkg['name']}=={pkg['version']}")
            else:
                print("No JAX-related packages found.")
                return results

            # Uninstall packages if not dry run
            if not dry_run and detected:
                print("Uninstalling JAX packages...")

                # Uninstall in reverse dependency order - start with main packages
                uninstall_order = ['jax'] + [pkg['name'] for pkg in detected if pkg['name'] != 'jax']

                for pkg_name in uninstall_order:
                    # Find the package info
                    pkg_info = next((p for p in detected if p['name'] == pkg_name), None)
                    if not pkg_info:
                        continue

                    try:
                        print(f"Uninstalling {pkg_name}...")
                        subprocess.run(
                            [sys.executable, '-m', 'pip', 'uninstall', pkg_name, '-y'],
                            check=True,
                            capture_output=True,
                            text=True
                        )
                        results['uninstalled_packages'].append(pkg_info)
                        print(f"Successfully uninstalled {pkg_name}")

                    except subprocess.CalledProcessError as e:
                        error_msg = f"Failed to uninstall {pkg_name}: {e}"
                        results['errors'].append(error_msg)
                        print(error_msg)
                        # Continue with other packages

                # Reset installation status
                self.jax_installed = False
                self.detected_config = None

                print(f"Uninstallation complete. Removed {len(results['uninstalled_packages'])} packages.")

            elif detected:
                print("Dry run - no packages were uninstalled.")

        except subprocess.CalledProcessError as e:
            error_msg = f"Error running pip list: {e}"
            results['errors'].append(error_msg)
            print(error_msg)
        except Exception as e:
            error_msg = f"Unexpected error during JAX detection/uninstallation: {e}"
            results['errors'].append(error_msg)
            print(error_msg)

        return results

    def clean_install_jax(self, force_variant: Optional[str] = None,
                         version_constraint: Optional[str] = None) -> bool:
        """
        Clean install JAX by first uninstalling any existing versions.

        This is useful when switching between JAX variants or fixing problematic installations.

        Args:
            force_variant: JAX variant to install (e.g., 'jax[cpu]', 'jax[cuda12]')
            version_constraint: Version constraint for installation

        Returns:
            bool: True if clean installation succeeded, False otherwise
        """
        print("Performing clean JAX installation...")

        # Step 1: Detect and uninstall existing JAX
        uninstall_results = self.detect_and_uninstall_jax(dry_run=False)

        if uninstall_results['errors']:
            print("Warning: Some errors occurred during uninstallation:")
            for error in uninstall_results['errors']:
                print(f"  - {error}")

        # Step 2: Install fresh JAX
        success = self.install_jax(force_variant, version_constraint)

        if success:
            print("Clean installation completed successfully.")
        else:
            print("Clean installation failed.")

        return success

    def __repr__(self) -> str:
        return f"JaxHelper(system={self.system}, jax_installed={self.jax_installed})"
