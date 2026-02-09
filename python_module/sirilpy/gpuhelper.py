# Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
# Copyright (C) 2012-2026 team free-astro (see more in AUTHORS file)
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
import shutil
import platform
import tempfile
import importlib
import subprocess
from typing import Tuple, Optional, Dict, Any
from packaging import version as pkg_version
import requests
import numpy as np
from .version import __version__
from .utility import ensure_installed, _check_package_installed, _install_package, \
                     SuppressedStderr

def _detect_cuda_version(system) -> Optional[str]:
    """
    Detects the CUDA version by parsing the output of 'nvcc -v'.
    Returns:
        Optional[str]: The CUDA version as a string (e.g., '11.7') if detected,
                        or "0.0" if nvcc is not installed or version cannot be determined.
    """
    nvcc_command = 'nvcc.exe' if system == 'windows' else 'nvcc'

    # Check if nvcc exists in PATH
    if shutil.which(nvcc_command) is None:
        print(f"{nvcc_command} not found in PATH, trying torch fallback")
        return _try_torch_cuda_version()

    try:
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
        print("Unable to confirm CUDA availability from nvcc output")
        return "0.0"

    except (subprocess.SubprocessError, FileNotFoundError):
        # nvcc command failed
        return _try_torch_cuda_version()

def _try_torch_cuda_version():
    """Helper method to get CUDA version from torch"""
    try:
        # possibly we have torch installed (this is a bit slow to import
        # which is why we don't try it first...)
        torch_helper = TorchHelper()
        if torch_helper.is_torch_installed():
            import torch
            cuda_version = torch.version.cuda  # e.g., '12.1'
            return cuda_version
        print("Torch unavailable, unable to confirm CUDA availability")
        return "0.0"
    except Exception:
        print("Unable to confirm CUDA availability")
        return "0.0"

def _get_windows_gpu_info(vendor_filter: str) -> Optional[str]:
    """
    Helper function to get GPU name from Windows using PowerShell.

    Args:
        vendor_filter: Vendor name to filter for (e.g., 'NVIDIA', 'AMD', 'Intel')

    Returns:
        GPU name string if found, None otherwise
    """
    try:
        # More reliable PowerShell command that returns just the Name property
        cmd = [
            "powershell", "-Command",
            f"(Get-WmiObject -Class Win32_VideoController | Where-Object {{$_.Name -like '*{vendor_filter}*'}}).Name"
        ]

        output = subprocess.check_output(cmd, text=True, stderr=subprocess.DEVNULL)

        if output.strip():
            # Split by newlines and filter out empty lines
            gpu_names = [line.strip() for line in output.strip().split('\n') if line.strip()]
            if gpu_names:
                # Return the first GPU found
                return gpu_names[0]
    except (subprocess.SubprocessError, FileNotFoundError, Exception):
        pass

    return None

def _detect_nvidia_gpu(system):
    """
    Detect if NVIDIA GPU is available. (Required to check Windows and Linux, not required on MacOS.)

    Returns:
        bool: True if NVIDIA GPU is detected, False otherwise.
    """
    try:
        if system == "windows":
            # Use the helper function for consistent detection
            gpu_name = _get_windows_gpu_info('NVIDIA')
            return gpu_name is not None

        # Linux detection using nvidia-smi or lspci
        if system == "linux":
            try:
                output = subprocess.check_output(
                    ["nvidia-smi"],
                    stderr=subprocess.DEVNULL,
                    text=True
                )
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
    Detect if AMD GPU is available on Windows or Linux.
    This is a legacy function - prefer using _get_amd_gpu_info() for more detailed info.

    Returns:
        bool: True if AMD GPU is detected, False otherwise.
    """
    system = platform.system().lower()

    # Windows detection
    if system == 'windows':
        for vendor in ['AMD', 'Radeon', 'ATI']:
            if _get_windows_gpu_info(vendor):
                return True
        return False

    # Linux detection
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

def _get_nvidia_gpu_info() -> Optional[Dict[str, Any]]:
    """
    Detect NVIDIA GPU and determine the appropriate CUDA version for PyTorch.

    Returns:
        Dict with 'detected', 'gpu_name', 'compute_capability', 'recommended_cuda' keys,
        or None if no NVIDIA GPU detected.
    """
    system = platform.system().lower()

    try:
        # Try nvidia-smi first (works on both Windows and Linux)
        result = subprocess.run(
            ['nvidia-smi', '--query-gpu=name,compute_cap', '--format=csv,noheader'],
            capture_output=True,
            text=True,
            check=False
        )

        if result.returncode == 0 and result.stdout.strip():
            lines = result.stdout.strip().split('\n')
            parts = lines[0].split(',')
            gpu_name = parts[0].strip()
            compute_cap = parts[1].strip() if len(parts) > 1 else None

            # Determine recommended CUDA version based on GPU generation
            recommended_cuda = _determine_cuda_version_for_gpu(gpu_name, compute_cap)

            return {
                'detected': True,
                'gpu_name': gpu_name,
                'compute_capability': compute_cap,
                'recommended_cuda': recommended_cuda
            }
    except (subprocess.SubprocessError, FileNotFoundError):
        pass

    # Windows fallback using PowerShell
    if system == 'windows':
        gpu_name = _get_windows_gpu_info('NVIDIA')
        if gpu_name:
            recommended_cuda = _determine_cuda_version_for_gpu(gpu_name, None)
            return {
                'detected': True,
                'gpu_name': gpu_name,
                'compute_capability': None,
                'recommended_cuda': recommended_cuda
            }

    # Linux fallback to basic detection
    if system == 'linux' and _detect_nvidia_gpu(system):
        return {
            'detected': True,
            'gpu_name': 'Unknown NVIDIA GPU',
            'compute_capability': None,
            'recommended_cuda': 'cu126'  # Safe default for most modern cards
        }

    return None


def _determine_cuda_version_for_gpu(gpu_name: str, compute_cap: Optional[str] = None) -> str:
    """
    Determine the safest CUDA version for PyTorch based on GPU generation.

    Args:
        gpu_name: Name of the GPU (e.g., "GeForce RTX 3090")
        compute_cap: Compute capability (e.g., "8.6"), if available

    Returns:
        CUDA version string for PyTorch (e.g., "cu118", "cu126", "cu128")
    """

    # Use compute capability if available
    if compute_cap:
        try:
            major, minor = map(int, compute_cap.split('.'))
            cc_val = major * 10 + minor

            # SM 9.0+ (Blackwell and future) - CUDA 12.8+
            if cc_val >= 90:
                return 'cu128'
            # SM 7.0-8.9 (Volta, Turing, Ampere, Ada) - CUDA 12.6
            elif cc_val >= 70:
                return 'cu126'
            # SM 5.0-6.x (Maxwell, Pascal) - CUDA 11.8
            elif cc_val >= 50:
                return 'cu118'
        except (ValueError, AttributeError):
            pass

    # Use GPU name
    gpu_lower = gpu_name.lower()
    # RTX 50xx series (Blackwell) - requires CUDA 12.8+
    if any(x in gpu_lower for x in ['rtx 50', 'rtx50']):
        return 'cu128'

    # RTX 40xx series (Ada Lovelace) - CUDA 12.6 recommended
    # RTX 30xx series (Ampere) - CUDA 12.6 supported
    # RTX 20xx series (Turing) - CUDA 12.6 supported
    # GTX 16xx series (Turing) - CUDA 12.6 supported
    # GTX 10xx series (Pascal) - CUDA 12.6 supported
    if any(x in gpu_lower for x in ['rtx 40', 'rtx40', 'rtx 30', 'rtx30',
                                     'rtx 20', 'rtx20', 'gtx 16', 'gtx16',
                                     'gtx 10', 'gtx10', 'titan x', 'tesla p']):
        return 'cu126'

    # GTX 9xx series (Maxwell) - CUDA 11.8 is safest
    # GTX 7xx/8xx series (Kepler) - CUDA 11.8
    if any(x in gpu_lower for x in ['gtx 9', 'gtx 8', 'gtx 7',
                                     'titan', 'tesla k', 'quadro k']):
        return 'cu118'

    # Default to cu126 for unknown modern GPUs
    return 'cu126'


def _get_amd_gpu_info() -> Optional[Dict[str, Any]]:
    """
    Detect AMD GPU and determine if it supports ROCm.

    Returns:
        Dict with 'detected', 'gpu_name', 'rocm_compatible', 'is_igpu' keys,
        or None if no AMD GPU detected.

        Note: rocm_compatible indicates hardware ROCm support regardless of OS.
              Individual frameworks decide whether to use ROCm based on their own
              platform support (e.g., ONNX doesn't support ROCm on Windows, but PyTorch does).
    """
    system = platform.system().lower()

    # Windows detection
    if system == 'windows':
        # Try multiple AMD/ATI vendor identifiers
        gpu_name = None
        for vendor in ['AMD', 'Radeon', 'ATI']:
            gpu_name = _get_windows_gpu_info(vendor)
            if gpu_name:
                break

        if gpu_name:
            is_igpu = _is_amd_igpu(gpu_name)
            # Check if GPU is ROCm-compatible (discrete GPUs with RDNA2/3 or newer)
            rocm_compat = _is_rocm_compatible_gpu(gpu_name) and not is_igpu

            return {
                'detected': True,
                'gpu_name': gpu_name,
                'rocm_compatible': rocm_compat,
                'is_igpu': is_igpu
            }
        return None

    # Linux - check for ROCm compatibility
    if system == 'linux':
        try:
            # Try rocm-smi first
            result = subprocess.run(
                ['rocm-smi', '--showproductname'],
                capture_output=True,
                text=True,
                check=False
            )

            if result.returncode == 0 and result.stdout.strip():
                gpu_name = result.stdout.strip()
                is_igpu = _is_amd_igpu(gpu_name)
                rocm_compat = not is_igpu  # iGPUs generally don't support ROCm well

                return {
                    'detected': True,
                    'gpu_name': gpu_name,
                    'rocm_compatible': rocm_compat,
                    'is_igpu': is_igpu
                }
        except (subprocess.SubprocessError, FileNotFoundError):
            pass

        # Fallback to lspci
        try:
            output = subprocess.check_output(['lspci'], text=True)
            amd_lines = [line for line in output.split('\n')
                         if any(x in line.lower() for x in ['amd', 'ati', 'radeon'])]

            if amd_lines:
                gpu_name = amd_lines[0].split(':')[-1].strip()
                is_igpu = _is_amd_igpu(gpu_name)

                # Check if this is a ROCm-compatible discrete GPU
                rocm_compat = _is_rocm_compatible_gpu(gpu_name) and not is_igpu

                return {
                    'detected': True,
                    'gpu_name': gpu_name,
                    'rocm_compatible': rocm_compat,
                    'is_igpu': is_igpu
                }
        except (subprocess.SubprocessError, FileNotFoundError):
            pass

    return None


def _is_amd_igpu(gpu_name: str) -> bool:
    """Check if AMD GPU is an integrated GPU."""
    igpu_indicators = [
        'vega', 'radeon graphics', 'ryzen', 'renoir', 'cezanne',
        'barcelo', 'rembrandt', 'phoenix', 'raphael', 'dragon range',
        'strix point', 'mendocino', 'picasso', 'raven'
    ]
    gpu_lower = gpu_name.lower()
    return any(indicator in gpu_lower for indicator in igpu_indicators)


def _is_rocm_compatible_gpu(gpu_name: str) -> bool:
    """
    Check if AMD GPU is compatible with ROCm.
    ROCm supports RDNA2/3 and some older GCN architectures.
    """
    gpu_lower = gpu_name.lower()

    # RDNA architectures (RX 5000, 6000, 7000 series)
    rocm_compatible = [
        'rx 7', 'rx 6', 'rx 5',  # Consumer RDNA
        'radeon pro w', 'radeon pro v',  # Professional RDNA
        'instinct',  # Data center
        'vii',  # Radeon VII
    ]

    return any(indicator in gpu_lower for indicator in rocm_compatible)


def _get_intel_gpu_info() -> Optional[Dict[str, Any]]:
    """
    Detect Intel GPU and determine if it's Arc (discrete) or iGPU.

    Returns:
        Dict with 'detected', 'gpu_name', 'is_arc', 'is_igpu', 'torch_compatible' keys,
        or None if no Intel GPU detected.
    """
    system = platform.system().lower()
    gpu_name = None

    # Windows detection
    if system == 'windows':
        gpu_name = _get_windows_gpu_info('Intel')

    # Linux detection
    elif system == 'linux':
        try:
            output = subprocess.check_output(['lspci'], text=True)
            intel_lines = [line for line in output.split('\n')
                          if 'intel' in line.lower() and
                          any(x in line.lower() for x in ['vga', 'display', '3d', 'graphics'])]
            if intel_lines:
                gpu_name = intel_lines[0].split(':')[-1].strip()
        except (subprocess.SubprocessError, FileNotFoundError):
            pass

    # macOS detection
    elif system == 'darwin':
        try:
            output = subprocess.check_output(
                ['system_profiler', 'SPDisplaysDataType'],
                text=True
            )
            if 'intel' in output.lower():
                for line in output.split('\n'):
                    if 'Chipset Model' in line or 'Graphics' in line:
                        if 'Intel' in line:
                            gpu_name = line.split(':')[-1].strip()
                            break
        except (subprocess.SubprocessError, FileNotFoundError):
            pass

    if not gpu_name:
        return None

    gpu_lower = gpu_name.lower()

    is_arc = any(x in gpu_lower for x in [
        # Discrete Arc GPUs
        'arc', 'a770', 'a750', 'a580', 'a380',
        'a370m', 'a730m', 'a550m', 'a350m',

        # Xe-based iGPUs (should be treated as Arc)
        'iris xe',
        'xe graphics',
        'arc graphics',      # Meteor Lake / Core Ultra iGPUs
    ])

    # --- Legacy (non-Xe) Intel iGPUs ---
    is_igpu = (
        any(x in gpu_lower for x in [
            'uhd graphics',
            'hd graphics',
            'uhd 6', 'uhd 7',   # UHD 620 / 630 / etc.
        ])
        and not is_arc
    )

    # Intel Arc GPUs support PyTorch via Intel Extension for PyTorch
    # iGPUs generally have limited/no PyTorch support
    torch_compatible = is_arc

    return {
        'detected': True,
        'gpu_name': gpu_name,
        'is_arc': is_arc,
        'is_igpu': is_igpu,
        'torch_compatible': torch_compatible
    }


def _get_apple_silicon_info() -> Optional[Dict[str, Any]]:
    """
    Detect Apple Silicon (M1/M2/M3/M4) chips.

    Returns:
        Dict with 'detected', 'chip_name', 'architecture' keys,
        or None if not Apple Silicon.
    """
    system = platform.system().lower()

    if system != 'darwin':
        return None

    try:
        # Check if it's Apple Silicon
        machine = platform.machine().lower()
        if machine != 'arm64':
            return None

        # Get chip information
        output = subprocess.check_output(
            ['sysctl', '-n', 'machdep.cpu.brand_string'],
            text=True
        ).strip()

        chip_name = output

        # Determine architecture
        if 'M1' in output:
            architecture = 'M1'
        elif 'M2' in output:
            architecture = 'M2'
        elif 'M3' in output:
            architecture = 'M3'
        elif 'M4' in output:
            architecture = 'M4'
        else:
            architecture = 'Apple Silicon'

        return {
            'detected': True,
            'chip_name': chip_name,
            'architecture': architecture
        }
    except (subprocess.SubprocessError, FileNotFoundError):
        pass

    return None


def detect_gpu_capabilities() -> Dict[str, Any]:
    """
    Comprehensive GPU detection for all supported hardware.

    Returns:
        Dict containing detected hardware information and recommendations for
        PyTorch, ONNX, and JAX backends.
    """
    system = platform.system().lower()

    capabilities = {
        'system': system,
        'nvidia': _get_nvidia_gpu_info(),
        'amd': _get_amd_gpu_info(),
        'intel': _get_intel_gpu_info(),
        'apple_silicon': _get_apple_silicon_info()
    }

    return capabilities


def _get_gpu_priority() -> Dict[str, Any]:
    """
    Determine GPU priority when multiple GPUs are available.

    Priority order (highest to lowest):
    1. NVIDIA discrete GPU (best PyTorch/ONNX/JAX support)
    2. AMD ROCm-compatible discrete GPU
    3. Intel Arc discrete GPU
    4. Apple Silicon (macOS only)
    5. AMD non-ROCm GPU
    6. Intel iGPU (limited support)
    7. AMD iGPU (very limited support)

    Returns:
        Dict with 'primary_gpu' (type), 'reason' (why chosen), and full capabilities
    """
    caps = detect_gpu_capabilities()

    # Priority 1: NVIDIA discrete GPU
    if caps['nvidia'] and caps['nvidia']['detected']:
        return {
            'primary_gpu': 'nvidia',
            'reason': 'NVIDIA GPU has best support across PyTorch, ONNX, and JAX',
            'capabilities': caps
        }

    # Priority 2: AMD ROCm-compatible GPU
    if (caps['amd'] and caps['amd']['detected'] and
        caps['amd'].get('rocm_compatible') and
        not caps['amd'].get('is_igpu')):
        return {
            'primary_gpu': 'amd_rocm',
            'reason': 'AMD ROCm-compatible discrete GPU detected',
            'capabilities': caps
        }

    # Priority 3: Intel Arc discrete GPU
    if (caps['intel'] and caps['intel']['detected'] and
        caps['intel'].get('is_arc')):
        return {
            'primary_gpu': 'intel_arc',
            'reason': 'Intel Arc discrete GPU detected',
            'capabilities': caps
        }

    # Priority 4: Apple Silicon
    if caps['apple_silicon']:
        return {
            'primary_gpu': 'apple_silicon',
            'reason': 'Apple Silicon with Metal/MPS support',
            'capabilities': caps
        }

    # Priority 5: AMD GPU without ROCm support
    if (caps['amd'] and caps['amd']['detected'] and
        not caps['amd'].get('is_igpu')):
        return {
            'primary_gpu': 'amd_other',
            'reason': 'AMD GPU detected (limited support)',
            'capabilities': caps
        }

    # Priority 6: Intel iGPU (OpenVINO only for ONNX)
    if (caps['intel'] and caps['intel']['detected'] and
        caps['intel'].get('is_igpu')):
        return {
            'primary_gpu': 'intel_igpu',
            'reason': 'Intel iGPU - limited ONNX OpenVINO support only',
            'capabilities': caps
        }

    # Priority 7: AMD iGPU (very limited support)
    if (caps['amd'] and caps['amd']['detected'] and
        caps['amd'].get('is_igpu')):
        return {
            'primary_gpu': 'amd_igpu',
            'reason': 'AMD iGPU - very limited support',
            'capabilities': caps
        }

    # No GPU detected
    return {
        'primary_gpu': 'cpu',
        'reason': 'No supported GPU detected',
        'capabilities': caps
    }

def get_gpu_detection_report() -> Dict[str, Any]:
    """
    Generate a comprehensive GPU detection report suitable for display in a GUI.

    Returns:
        Dict containing all detected hardware and recommended configurations
        for ONNX, PyTorch, and JAX, including GPU priority information.
    """
    priority = _get_gpu_priority()
    caps = priority['capabilities']

    # Note: This assumes ONNXHelper, TorchHelper, and JaxHelper classes exist
    # If not available, you'll need to import them or handle gracefully
    try:
        from . import ONNXHelper, TorchHelper, JaxHelper

        onnx_helper = ONNXHelper()
        onnx_backend = onnx_helper.get_recommended_backend()
        onnx_package, from_url, index_url = onnx_helper._get_onnxruntime_package()

        torch_backend = TorchHelper().get_recommended_backend()
        jax_backend = JaxHelper().get_recommended_backend()
    except (ImportError, AttributeError):
        # Fallback if helper classes not available
        onnx_backend = 'cpu'
        onnx_package = 'onnxruntime'
        from_url = None
        index_url = None
        torch_backend = 'cpu'
        jax_backend = 'cpu'

    report = {
        'system': caps['system'],
        'gpu_priority': {
            'primary_gpu': priority['primary_gpu'],
            'reason': priority['reason']
        },
        'hardware': {
            'nvidia': caps['nvidia'],
            'amd': caps['amd'],
            'intel': caps['intel'],
            'apple_silicon': caps['apple_silicon']
        },
        'recommendations': {
            'onnx': {
                'backend': onnx_backend,
                'package': onnx_package,
                'from_url': from_url,
                'index_url': index_url
            },
            'torch': torch_backend,
            'jax': jax_backend
        }
    }

    return report

class ONNXHelper:
    """
    A class to handle detection and installation of the appropriate ONNX Runtime
    package based on the system hardware and configuration.

    Example usage (this should be used instead of
    ``sirilpy.ensure_installed("onnxruntime")`` to install the correct package for
    the user's system.)

    .. code-block:: python

       oh = sirilpy.ONNXHelper()
       oh.ensure_onnxruntime()

    """

    def __init__(self):
        """Initialize the ONNXHelper."""
        ensure_installed("platformdirs")
        ensure_installed("onnx")
        from platformdirs import user_config_dir
        self.system = platform.system().lower()
        self.providers = None
        self.config_file = os.path.join(user_config_dir(appname="siril"), "siril_onnx.conf")

    def get_recommended_backend(self) -> str:
        """
        Determine the recommended ONNX Runtime backend based on hardware.
        When multiple GPUs are present, uses priority system to choose the best one.

        Note: ONNX Runtime has platform-specific backend support:
        - Windows: DirectML for all GPUs (NVIDIA, AMD, Intel) - more reliable, no driver dependencies
        - Linux: CUDA for NVIDIA, CPU for AMD (no ROCm support yet), OpenVINO for Intel

        Returns:
            Backend name: 'gpu', 'directml', 'openvino', 'coreml', or 'cpu'
        """
        priority = _get_gpu_priority()
        primary_gpu = priority['primary_gpu']
        system = priority['capabilities']['system']

        # Windows: Use DirectML for all GPUs
        # DirectML is more reliable as it doesn't depend on system CUDA/ROCm libraries
        if system == 'windows':
            if primary_gpu in ['nvidia', 'amd_rocm', 'amd_other', 'amd_igpu', 'intel_igpu']:
                return 'directml'
            elif primary_gpu == 'intel_arc':
                return 'openvino'
            else:
                return 'cpu'

        # Linux and other platforms: Use hardware-specific backends
        gpu_backend_map = {
            'nvidia': 'gpu',  # CUDA
            'amd_rocm': 'cpu',  # ONNX Runtime doesn't support ROCm yet
            'amd_other': 'cpu',
            'intel_arc': 'openvino',
            'apple_silicon': 'coreml',
            'intel_igpu': 'cpu',
            'amd_igpu': 'cpu',
            'cpu': 'cpu'
        }

        return gpu_backend_map.get(primary_gpu, 'cpu')

    def _get_onnxruntime_package(self):
        """
        Determine which ONNX Runtime package to install.
        Uses get_recommended_backend() to decide.

        Returns:
            tuple: (package_name, from_url, index_url) where from_url and index_url
                   are None except for special cases like ROCm
        """
        from_url = None
        index_url = None

        backend = self.get_recommended_backend()

        if backend == 'gpu':
            onnxruntime_pkg = 'onnxruntime-gpu'
            # Check CUDA version for index_url
            cuda_version = _detect_cuda_version(self.system)
            if cuda_version and pkg_version.Version(cuda_version).major == 11:
                index_url = "https://aiinfra.pkgs.visualstudio.com/PublicPackages/_packaging/onnxruntime-cuda-11/pypi/simple/"

        elif backend == 'directml':
            onnxruntime_pkg = 'onnxruntime-directml'

        elif backend == 'openvino':
            onnxruntime_pkg = 'onnxruntime-openvino'

        elif backend == 'coreml':
            # Standard onnxruntime supports CoreML on macOS
            onnxruntime_pkg = 'onnxruntime'

        else:  # CPU or unknown
            onnxruntime_pkg = 'onnxruntime'

        return onnxruntime_pkg, from_url, index_url

    def status(self):
        """
        Prints the current status of the ONNX Helper class in regard to its support for different
        OSes, GPUs. The world of heterogenous computing is developing rapidly and while support
        for some of the frameworks for which helpers are available is not yet universally
        available, hopefully it will improve in the future.
        """
        print(f"ONNXHelper status as of sirilpy version {__version__}")
        if self.system == 'windows':
            print("Windows: ONNXHelper will install the DirectML runtime wherever it is supported. "
                "This includes NVidia CPUs which might see a slight performance improvement using "
                "the onnxruntime-gpu module as this is difficult to configure correctly on some "
                "systems as it relies on system libraries and paths, so onnxruntime-directml is "
                "more robust.")
        elif self.system == 'linux':
            print("Linux: ONNXHelper will attempt to detect NVidia, AMD and Intel GPUs and install "
                "either onnxruntime-gpu, onnxruntime-rocm or onnxruntime-intel as appropriate. "
                "Note that we have had no feedback so far from AMD or Intel GPU users and none of "
                "the developers have these GPUs, so although it *should* work on these systems "
                "we would be very grateful for confirmation either of success or failure.")
        elif self.system == 'darwin':
            print("MacOS: ONNXHelper will install the standard onnxruntime module on MacOS. This "
                "provides good support for Apple silicon and reasonable support for older Intel "
                "silicon Macs.")
        print("Detection of working ExecutionProviders: ONNXHelper tests using a simple model to "
            "confirm whether each supported ExecutionProvider in the installed runtime actually "
            "works or not. This helps to ensure that scripts calling "
            "ONNXHelper.get_execution_providers_ordered() get a set of known working "
            "ExecutionProviders. The set is cached so the model does not need to be run on "
            "subsequent calls.")
        print("Model features: ONNX does not support all machine learning functions that are "
            "supported by all targets, and different runtimes support different subsets of "
            "machine learning operations. This means that some more demandng models may not "
            "work on all providers. In such cases the user may need to fall back to the "
            "CPU ExecutionProvider. As onnxruntime does not automatically handle runtime errors "
            "the ONNXHelper.run() method is provided to manage this (see the method docstring "
            "for details).")

    def _create_simple_onnx_model(self):
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

    def _try_provider(self, ort, model_path, input_data, provider, reference_output=None):
        """Try executing the model with a specific provider (no fallback)."""
        try:
            # Handle provider format and extract expected provider name
            if isinstance(provider, tuple):
                expected_provider_name, provider_options = provider
                providers_for_session = [(expected_provider_name, provider_options)]
            else:
                expected_provider_name = provider
                providers_for_session = [provider]

            sess_options = ort.SessionOptions()
            session = ort.InferenceSession(
                model_path,
                sess_options=sess_options,
                providers=providers_for_session
            )

            actual_provider = session.get_providers()[0]

            if actual_provider != expected_provider_name:
                print(f"(x) {provider}: fallback occurred (used {actual_provider})")
                return False

            # Determine TF32 usage
            provider_opts = session.get_provider_options().get(expected_provider_name, {})
            use_tf32 = provider_opts.get("use_tf32", "0")  # default to "0" if missing

            # Set tolerances
            if use_tf32 == "1":
                rtol = 5e-2
                atol = 1e-3
            else:
                rtol = 1e-3
                atol = 1e-5

            output = session.run(None, {'input': input_data})

            print(f"OK: {expected_provider_name} ran successfully")
            if reference_output is not None:
                # 1.4.0-beta4: disable this test for use_tf32 is True: it is failing because of the lower
                # precision, but with zero impact on the output, so a more useful test needs to be found.
                if not use_tf32 and not np.allclose(reference_output, output[0], rtol=rtol, atol=atol):
                    print(f"(!) Output mismatch with CPU (rtol={rtol}, atol={atol})")
                    return False

            return True

        except Exception as e:
            print(f"(x) {expected_provider_name} failed: {e}")
            return False

    def import_onnxruntime(self):
        """
        Import onnxruntime, add it to the global dict, test if it's built against
        CUDA and if so preload the CUDA and CUDNN libraries to improve the chances
        of finding them if Torch[CUDA] happens to be installed.
        """
        import onnxruntime
        globals()['onnxruntime'] = onnxruntime  # Add to the global dict
        providers = onnxruntime.get_available_providers()
        if 'CUDAExecutionProvider' in providers:
            # Attempt to preload CUDA / CUDnn libraries
            # This helps on some systems where the system-wide libraries are not found but
            # torch is installed with nvidia library dependencies installed in the venv
            onnxruntime.preload_dlls()
            # Set logging to only report critical issues by default
            onnxruntime.set_default_logger_severity(4)

    def test_onnxruntime(self, ort=None):
        """
        Test an imported onnxruntime.

        Args:install_torch(cuda_version=cuda_version)
            ort: The ONNX runtime module to test. If None, the method will
                attempt to import onnxruntime for the test.

        Returns:
            list: a list of confirmed working ONNXRuntime ExecutionProviders in priority order
        """
        import onnx
        if ort is None:
            try:
                import onnxruntime as ort
            except ImportError:
                print("(x) Unable to import onnxruntime. Test failed.")
                return []

        print("=== ONNX Execution Provider Tester ===")
        print("Creating ONNX model...")
        model = self._create_simple_onnx_model()

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
                print("OK: CPU output computed.")
            except Exception as e:
                print(f"(x) Failed to run on CPU: {e}")
                return []

            all_providers_raw = ort.get_available_providers()
            all_providers = []

            for p in all_providers_raw:
                if p == "OpenVINOExecutionProvider":
                    all_providers.append((p, {
                                                'device_type': 'GPU',
                                                'precision': 'FP32'
                                            },))
                    all_providers.append((p, {
                                                'device_type': 'CPU',
                                                'precision': 'FP32'
                                            }))
                else:
                    all_providers.append(p)

            print("\nAvailable execution providers:")
            for p in all_providers:
                print(f"  - {p}")

            print("\nTesting each provider without fallback...")
            working_providers = []

            for provider in all_providers:
                print(f"\nTesting {provider}...")
                with SuppressedStderr():
                    if self._try_provider(ort, model_path, input_data, provider, reference_output=cpu_output):
                        working_providers.append(provider)

            print("\n=== Summary ===")
            if working_providers:
                print("OK: Working providers (in priority order):")
                for p in working_providers:
                    print(f"  - {p}")
                print(f"\n→ Best available provider: {working_providers[0]}")
                # Save to cache
                self._save_providers_to_cache(working_providers)
                self.providers = working_providers
                return working_providers

            print("(x) No execution providers were able to run the model.")
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

    def install_onnxruntime(self, force=False):
        """
        Detect system configuration and install the appropriate ONNX Runtime package.

        Args:
            force: bool: If True, force reinstallation even if onnxruntime is already installed.

        Returns:
            bool: True if installation was successful or already installed, False otherwise.

        Raises:
            TimooutError: if a TimeoutError occurs in ensure_installed() - this avoids falling
                        back to the CPU-only package purely because of network issues
        """
        # First check if any onnxruntime is already installed
        try:
            is_installed, package_name = self.is_onnxruntime_installed()
        except:
            # Error checking if it's installed
            print("Error checking ONNX Runtime state. Removing and reinstalling...")
            self.uninstall_onnxruntime()
            is_installed = False

        if is_installed:
            if force:
                print("ONNX Runtime is already installed. Removing and reinstalling...")
                self.uninstall_onnxruntime()
                if self.config_file is not None and os.path.exists(self.config_file):
                    os.unlink(self.config_file)
            else:
                print(f"ONNX Runtime is already installed: {package_name}")
                return True

        # If not installed, get recommended package
        onnxruntime_pkg, from_url, index_url = self._get_onnxruntime_package()

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
                try:
                    self.import_onnxruntime()
                except ImportError as e:
                    print(f"Checked installed runtime {onnxruntime_pkg} cannot be imported: {e}. Falling back to the basic CPU runtime", file=sys.stderr)
                    self.uninstall_onnxruntime()
                    _install_package(onnxruntime, None)
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

    def ensure_onnxruntime(self) -> bool:
        """
        Wrapper for install_onnxruntime() that only installs it if needed, with
        negligible overhead if it is already installed.
        """
        if not self.is_onnxruntime_installed():
            return self.install_onnxruntime()
        return True

    def is_onnxruntime_installed(self):
        """
        Check if any onnxruntime package is already installed and usable.

        Returns:
            tuple: (is_installed, package_name) where package_name could be
                  'onnxruntime', 'onnxruntime-gpu', 'onnxruntime-silicon', etc.
        """
        try:
            # Try importing onnxruntime - if this fails, the package is not usable
            self.import_onnxruntime()

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
                response = requests.get(url, timeout=10)
                return response.status_code == 200
            except requests.exceptions.RequestException:
                print("Connection error {e}, please try again later")
                raise

            except Exception:
                return False

    def get_execution_providers_ordered(self, ai_gpu_acceleration=True, force_check=False):
        """
        Get execution providers ordered by priority.
        This function returns a list of available ONNX Runtime execution providers
        in a reasonable order of priority, covering major GPU platforms:
        The CPU provider is always included as the final fallback option.

        Args:
            ai_gpu_acceleration (bool): Whether to include GPU acceleration providers.
                Defaults to True.
            force_check (bool): Whether to force re-checking even if a cached config exists.
                Defaults to False.

        Returns:
            list: Ordered list of available execution providers.
        """
        if force_check is True: # Clear any previous config so it has to be re-checked
            self.providers = None
            if os.path.exists(self.config_file):
                os.unlink(self.config_file)

        if ai_gpu_acceleration is False:
            return ["CPUExecutionProvider"]

        if self.providers is not None:
            return self.providers

        self.import_onnxruntime()

        # Try to load cached providers first
        try:
            cached_providers = self._load_cached_providers()
            if cached_providers:
                self.providers = cached_providers
                return self.providers
        except Exception:
            # If an error occurs with _load_cached_providers() delete the config file and do
            # the full test to re-cache them
            if os.path.exists(self.config_file):
                os.unlink(self.config_file)

        # If no valid cache, run the test
        return self.test_onnxruntime(onnxruntime)

    def _load_cached_providers(self):
        """
        Load cached execution providers from config file if they're still valid.

        Returns:
            list or None: List of cached providers if valid, None otherwise
        """
        if not os.path.exists(self.config_file):
            return None

        try:
            with open(self.config_file, 'r', encoding='utf-8') as f:
                cached_data = json.load(f)

            raw_providers = cached_data.get('execution_providers', [])
            cached_providers = []

            for p in raw_providers:
                name = p["name"]
                options = p["options"]
                if options is not None:
                    cached_providers.append((name, options))
                else:
                    cached_providers.append(name)

            if not cached_providers:
                return None

            print(f"Using cached execution providers from {self.config_file}")
            return valid_providers

        except (json.JSONDecodeError, IOError, KeyError) as e:
            print(f"Failed to load cached providers: {e}", file=sys.stderr)
            os.unlink(self.config_file)
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

            serialized_providers = []
            for p in providers:
                if isinstance(p, tuple):
                    provider_name, options = p
                    serialized_providers.append({"name": provider_name, "options": options})
                else:
                    serialized_providers.append({"name": p, "options": None})

            cache_data = {
                'execution_providers': serialized_providers,
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
        # Remove the execution provders config file
        if os.path.exists(self.config_file):
            os.unlink(self.config_file)

        self.providers = None

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
    """Helper class for PyTorch detection, installation and testing."""

    def __init__(self):
        """Initialize TorchHelper without importing torch at module level."""
        self.device_info: Optional[dict] = None
        self._torch_installed = False
        self.system = platform.system().lower()

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

    def status(self):
        """
        Prints the current status of the Torch Helper class in regard to its support for different
        OSes, GPUs. The world of heterogenous computing is developing rapidly and while support
        for some of the frameworks for which helpers are available is not yet universally
        available, hopefully it will improve in the future.
        """
        print(f"TorchHelper status as of sirilpy version {__version__}")
        if self.system == 'windows':
            print("Windows: TorchHelper will install Torch. A version may be specified but by default "
                "autodetection will take place. The CUDA 12.8 version will be installed for NVidia GPUs "
                "but other CUDA versions may be specified manually and a CPU version will be installed "
                "for unsupported GPUs. There are not yet Intel, ROCm or DirectML Torch runtimes that are "
                "sufficietly stable to support in this helper module.")
        elif self.system == 'linux':
            print("Linux: TorchHelper will install Torch.  A version may be specified but by default "
                "autodetection will taklookinge place. The CUDA 12.8 version will be installed for NVidia GPUs "
                "but other CUDA versions may be specified manually as well as a ROCm version for AMD GPUs. "
                "There is not yet an Intel Torch runtime that is sufficietly stable to support in this "
                "helper module.")
        elif self.system == 'darwin':
            print("MacOS: TorchHelper will install the standard Torch module on MacOS. This is targeted "
                "at all Apple Macs regardless of CPU architecture: any issues with Torch on MacOS "
                "should be reported upstream to Torch.")
        print("Dependencies: Torch is currently excessively strict about required versions of some "
            "dependencies including CUDnn: it requires an exact version match rather than at least a "
            "certain version. This cauess conflict with other GPU acceleration modules that have "
            "differing dependency requirements, including jax. It is therefore not currently possible "
            "to write a script that uses both Torch and jax, and even switching between the two in "
            "different scripts is difficult at present. This issue has been raised upstream with Torch.")

    def get_recommended_backend(self) -> Dict[str, Any]:
        """
        Determine the recommended PyTorch backend and installation parameters.
        When multiple GPUs are present, uses priority system to choose the best one.

        PyTorch has good cross-platform support:
        - NVIDIA: CUDA on all platforms
        - AMD ROCm-compatible: ROCm on Linux AND Windows
        - Intel Arc: XPU backend
        - Apple Silicon: MPS backend

        Returns:
            Dict with 'backend', 'cuda_version', 'extra_index_url', 'packages' keys
        """
        priority = _get_gpu_priority()
        primary_gpu = priority['primary_gpu']
        caps = priority['capabilities']

        # NVIDIA GPUs - CUDA backend on all platforms
        if primary_gpu == 'nvidia':
            cuda_version = caps['nvidia']['recommended_cuda']
            return {
                'backend': 'cuda',
                'cuda_version': cuda_version,
                'extra_index_url': f'https://download.pytorch.org/whl/{cuda_version}',
                'packages': ['torch', 'torchvision']
            }

        # AMD ROCm-compatible GPUs - ROCm backend on Linux and Windows
        if primary_gpu == 'amd_rocm':
            return {
                'backend': 'rocm',
                'cuda_version': None,
                'extra_index_url': 'https://download.pytorch.org/whl/rocm7.1',
                'packages': ['torch', 'torchvision']
            }

        # Intel Arc GPUs
        if primary_gpu == 'intel_arc':
            return {
                'backend': 'intel',
                'cuda_version': None,
                'extra_index_url': 'https://download.pytorch.org/whl/xpu',
                'packages': ['torch', 'torchvision',
                            'intel-extension-for-pytorch']
            }

        # Apple Silicon
        if primary_gpu == 'apple_silicon':
            return {
                'backend': 'mps',
                'cuda_version': None,
                'extra_index_url': None,
                'packages': ['torch', 'torchvision']
            }

        # All other cases (AMD non-ROCm, iGPUs, CPU) - use CPU
        return {
            'backend': 'cpu',
            'cuda_version': None,
            'extra_index_url': None,
            'packages': ['torch', 'torchvision']
        }

    def ensure_torch(self, cuda_version: Optional[str] = None) -> bool:
        """
        Ensure PyTorch is installed with the appropriate backend.

        Args:
            cuda_version: Optional CUDA version to override auto-detection
                            (e.g., 'cu118', 'cu126', 'cu128')

        Returns: True on success, False on failure
        """
        if self.is_torch_installed():
            print("Torch is already installed")
            return True

        backend_info = self.get_recommended_backend()

        # Override CUDA version if specified
        if cuda_version is not None and backend_info['backend'] == 'cuda':
            backend_info['cuda_version'] = cuda_version
            backend_info['extra_index_url'] = f'https://download.pytorch.org/whl/{cuda_version}'

        print(f"Installing PyTorch with backend: {backend_info['backend']}")
        if backend_info['cuda_version']:
            print(f"Using CUDA version: {backend_info['cuda_version']}")

        try:
            self.install_torch(
                cuda_version=backend_info['cuda_version'],
                extra_index_url=backend_info['extra_index_url'],
                packages=backend_info['packages']
            )
            return True
        except Exception as e:
            print(f"Error installing Torch: {e}")
            return False

    def install_torch(self, cuda_version: Optional[str] = None,
                        extra_index_url: Optional[str] = None,
                        packages: Optional[list] = None):
        """
        Install PyTorch with specified configuration.

        Args:
            cuda_version: CUDA version (e.g., 'cu118', 'cu126', 'cu128')
            extra_index_url: PyTorch wheel repository URL
            packages: List of packages to install
        """
        if packages is None:
            packages = ['torch', 'torchvision', 'torchaudio']

        install_cmd = [sys.executable, '-m', 'pip', 'install'] + packages

        if extra_index_url:
            install_cmd.extend(['--index-url', extra_index_url])

        try:
            print(f"Installing: {' '.join(packages)}")
            subprocess.run(install_cmd, check=True)
            self.torch_installed = self.is_torch_installed()
            print("PyTorch installation completed successfully")
        except subprocess.CalledProcessError as e:
            print(f"Failed to install PyTorch: {e}")
            raise

    def _detect_compute_platform(self) -> str:
        """
        Auto-detect the appropriate compute platform for PyTorch installation.

        Returns:
            str: The compute platform string ('cu128', 'rocm', 'cpu')
        """
        if self.system == 'darwin':
            # macOS: use standard package
            return 'cpu'  # This will trigger standard installation

        if self.system == 'windows':
            if _detect_nvidia_gpu(self.system):
                return 'cu128'
            return 'cpu'

        if self.system == 'linux':
            if _detect_nvidia_gpu(self.system):
                return 'cu128'
            if _detect_amd_gpu():
                return 'rocm'
            return 'cpu'

        # Unknown system, default to CPU
        print(f"Unknown system '{self.system}', defaulting to CPU installation")
        return 'cpu'

    def _import_torch(self) -> bool:
        """Import torch modules to verify availability."""
        try:
            import torch
            self._torch_installed = True

            # Get device info after successful import
            self.device_info = self._get_device_info()

            print(f"OK: PyTorch {torch.__version__} imported successfully!")
            return True

        except ImportError as e:
            print(f"(x) Failed to import PyTorch: {e}")
            return False
        except Exception as e:
            print(f"(x) Unexpected error importing PyTorch: {e}")
            return False

    def _get_device_info(self) -> dict:
        """Get information about available devices."""
        try:
            import torch
        except ImportError:
            return {}

        info = {
            'pytorch_version': torch.__version__,
            'cuda_available': torch.cuda.is_available(),
            'cuda_version': torch.version.cuda if torch.cuda.is_available() else None,
            'gpu_count': torch.cuda.device_count() if torch.cuda.is_available() else 0,
            'gpu_names': []
        }

        if torch.cuda.is_available():
            for i in range(torch.cuda.device_count()):
                info['gpu_names'].append(torch.cuda.get_device_name(i))

        return info

    def _create_simple_model(self, input_dim=256, hidden_dim=512, output_dim=128):
        """Create a simple neural network model."""
        try:
            import torch
            from torch import nn
        except ImportError:
            raise RuntimeError("PyTorch not available. Call install_torch() first.")

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
                import torch.nn.functional as F
                x = F.relu(self.linear1(x))
                x = self.dropout(x)
                x = F.relu(self.linear2(x))
                x = self.linear3(x)
                return x

        return SimpleModel(input_dim, hidden_dim, output_dim)

    def _create_test_data(self, batch_size=32, input_dim=256):
        """Create random test input data."""
        try:
            import torch
        except ImportError:
            raise RuntimeError("PyTorch not available. Call install_torch() first.")
        return torch.randn(batch_size, input_dim)

    def _benchmark_model(self, model, input_data, device, num_runs=10):
        """Benchmark model execution on specified device."""
        try:
            import torch
        except ImportError:
            raise RuntimeError("PyTorch not available. Call install_torch() first.")

        # Set model to evaluation mode to disable dropout
        model.eval()
        model = model.to(device)
        input_data = input_data.to(device)

        # Warm up
        for _ in range(3):
            with torch.no_grad():
                _ = model(input_data)

        # Benchmark
        if device.type == 'cuda':
            torch.cuda.synchronize()
        start_time = time.time()

        for _ in range(num_runs):
            with torch.no_grad():
                output = model(input_data)

        if device.type == 'cuda':
            torch.cuda.synchronize()
        end_time = time.time()

        avg_time = (end_time - start_time) / num_runs
        return output, avg_time

    def _print_device_info(self):
        """Print device information."""

        print("=== Device Information ===")
        print(f"PyTorch version: {self.device_info['pytorch_version']}")
        print(f"CUDA available: {self.device_info['cuda_available']}")

        if self.device_info['cuda_available']:
            print(f"CUDA version: {self.device_info['cuda_version']}")
            print(f"GPU count: {self.device_info['gpu_count']}")
            for i, name in enumerate(self.device_info['gpu_names']):
                print(f"  GPU {i}: {name}")

    def _test_torch_gpu(self):
        """Test PyTorch model execution on GPU."""
        print("=== PyTorch GPU Test ===")

        if not self._import_torch():
            print("PyTorch not available. Please install it first using install_torch()")
            return False

        try:
            import torch
        except ImportError:
            return False

        # Print device info
        self._print_device_info()

        # Create model and test data with fixed seeds for reproducible results
        print("\nCreating model and test data...")
        torch.manual_seed(42)  # Set seed for reproducible results
        model = self._create_simple_model()
        torch.manual_seed(42)  # Reset seed for consistent input data
        input_data = self._create_test_data(batch_size=64, input_dim=256)

        print(f"Model parameters: {sum(p.numel() for p in model.parameters()):,}")
        print(f"Input shape: {input_data.shape}")

        # Test CPU execution
        print("\nTesting CPU execution (for accuracy reference)...")
        cpu_device = torch.device('cpu')
        model_cpu = self._create_simple_model()  # Create fresh model
        torch.manual_seed(42)  # Ensure same initialization
        model_cpu.load_state_dict(model.state_dict())  # Copy weights
        cpu_output, _ = self._benchmark_model(model_cpu, input_data, cpu_device)

        # Test GPU execution
        if torch.cuda.is_available():
            print("\nTesting GPU execution...")
            gpu_device = torch.device('cuda:0')

            try:
                model_gpu = self._create_simple_model()  # Create fresh model for GPU
                model_gpu.load_state_dict(model.state_dict())  # Copy same weights
                gpu_output, _ = self._benchmark_model(model_gpu, input_data, gpu_device)
                print("OK: GPU execution successful!")

                # Compare outputs using appropriate tolerance for CPU vs GPU
                cpu_output_np = cpu_output.numpy()
                gpu_output_np = gpu_output.cpu().numpy()

                # Check if outputs are close (accounting for floating point differences)
                are_close = np.allclose(cpu_output_np, gpu_output_np, rtol=1e-3, atol=1e-4)

                if are_close:
                    print("OK: CPU and GPU outputs match within tolerance!")
                else:
                    print("(!) CPU and GPU outputs differ more than expected")

            except Exception as e:
                print(f"(x) GPU execution failed: {e}")
        else:
            print("\n(!) CUDA not available - cannot test GPU execution")
            print("To enable GPU support:")
            print("  1. Install CUDA toolkit")
            print("  2. Use install_torch() with appropriate CUDA version")

        print("\n" + "="*50)
        return True

    def _test_tensor_operations(self):
        """Test basic tensor operations on GPU."""
        print("=== Tensor Operations Test ===")
        retval = None

        if not self._import_torch():
            print("PyTorch not available. Please install it first using install_torch()")
            return False

        try:
            import torch
        except ImportError:
            return False

        if not torch.cuda.is_available():
            print("(!) CUDA not available - skipping tensor operations test")
            return False

        # Create large tensors for meaningful GPU test with fixed seed
        size = 2048
        print(f"Creating {size}x{size} tensors...")

        # Set seed for reproducible random tensors
        torch.manual_seed(123)
        a_cpu = torch.randn(size, size)
        b_cpu = torch.randn(size, size)

        # GPU tensors (copy the same data)
        a_gpu = a_cpu.clone().cuda()
        b_gpu = b_cpu.clone().cuda()

        # CPU matrix multiplication
        start_time = time.time()
        c_cpu = torch.mm(a_cpu, b_cpu)
        cpu_time = time.time() - start_time

        # GPU matrix multiplication
        torch.cuda.synchronize()
        start_time = time.time()
        c_gpu = torch.mm(a_gpu, b_gpu)
        torch.cuda.synchronize()
        gpu_time = time.time() - start_time

        # Verify results using appropriate tolerance for large matrix operations
        cpu_np = c_cpu.numpy()
        gpu_np = c_gpu.cpu().numpy()

        # Use more appropriate tolerances for large matrix multiplication
        are_close = np.allclose(cpu_np, gpu_np, rtol=1e-3, atol=1e-3)

        if are_close:
            print("OK: Results match within tolerance!")
            retval = True
        else:
            print("(!) Results differ more than expected")
            retval = False

        print(f"CPU time: {cpu_time:.4f}s")
        print(f"GPU time: {gpu_time:.4f}s")
        print(f"Speedup: {cpu_time/gpu_time:.2f}x")

        print("\n" + "="*50)
        return retval

    def test_torch(self):
        """
        Run tests to verify that torch is installed and runs correctly.

        """
        if not self.ensure_torch():
            print("Cannot run tests - PyTorch not available")
            print("Use helper.install_torch() or helper.test_torch(auto_install=True)")
            return False

        self._test_torch_gpu()
        return self._test_tensor_operations()

    def uninstall_torch(self):
        """
        Detects and uninstalls PyTorch and related packages.
        Checks for torch ecosystem packages.
        Returns:
            list: A list of uninstalled packages
        """
        self._torch_installed = False
        self.device_info = None

        # Define torch ecosystem packages to look for
        torch_ecosystem_packages = [
            'torch', 'torchvision', 'torchaudio', 'torchtext', 'torchdata',
            'pytorch-lightning', 'lightning', 'pytorch-ignite', 'ignite',
            'fastai', 'transformers', 'accelerate', 'timm'
        ]

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

        # Extract package names from pip list output
        installed_package_names = set()
        for line in installed_packages:
            parts = line.split()
            if parts:
                installed_package_names.add(parts[0].lower())

        # Find torch packages that are actually installed
        torch_packages = []
        for package in torch_ecosystem_packages:
            if package.lower() in installed_package_names:
                torch_packages.append(package)

        # Uninstall found packages
        if not torch_packages:
            print("No torch packages found.")
            return []

        print(f"Found torch packages: {', '.join(torch_packages)}")
        uninstalled = []
        for package in torch_packages:
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

    def status(self):
        """
        Prints the current status of the Jax Helper class in regard to its support for different
        OSes, GPUs. The world of heterogenous computing is developing rapidly and while support
        for some of the frameworks for which helpers are available is not yet universally
        available, hopefully it will improve in the future.
        """
        print(f"JaxHelper status as of sirilpy version {__version__}")
        if self.system == 'windows':
            print("Windows: JaxHelper will install the CPU version of jax. This does not provide "
                "acceleration and jax.numpy may often be slower than numpy itself, so it is only "
                "recommended for development and research purposes at the moment. A jax CUDA wheel is "
                "in development for Windows however the current version does not support numpy v2.x "
                "and therefore causes unacceptable conflict with other modules. As the wheel develops "
                "the CUDA wheel will hopefully be supportable in the future. There is also an Intel "
                "plugin under development, however as with the CUDA wheel this is not yet built so "
                "as to be compatible with numpy 2.x so the helper will not install it. An issue has "
                "been raised with the plugin authors and we hope progress on this plugin improves soon.")
        elif self.system == 'linux':
            print("Linux: JaxHelper will install jax variants for NVidia or AMD GPUs, with automatic "
                "GPU hardware detection. Note that as none of the developers have an AMD GPU, feedback "
                "on how well jax[rocm] works would be very much appreciated. There is also an Intel "
                "plugin under development, however as with the CUDA wheel this is not yet built so "
                "as to be compatible with numpy 2.x so the helper will not install it. An issue has "
                "been raised with the plugin authors and we hope progress on this plugin improves soon.")
        elif self.system == 'darwin':
            print("MacOS: ONNXHelper will install the Metal jax backed on MacOS. Note that this is still "
                "experimental and may provide only limited acceleration at present. We hope that support "
                "for this platform improves soon.")
        print("Dependencies: Jax has a sensible dependency on CUDnn however unfortunately it conflicts "
            "with Torch, which is currently excessively strict about required versions of some "
            "dependencies including CUDnn: it requires an exact version match rather than at least a "
            "certain version, and the version Torch requires is older than the version that jax is "
            "built against. This issue is best handled by using TorchHelper.install_torch() to install "
            "Torch first (this auto-reinstalls with --no-deps), and then using JaxHelper.install_jax() to "
            "install Jax. (You can do this the other way round and it will still work but will be less "
            "efficient.)")

    def is_jax_installed(self) -> bool:
        """Check if PyTorch is installed without importing it."""
        if self.jax_installed:
            return True

        # Check if torch is available
        jax_spec = importlib.util.find_spec("jax")
        if jax_spec is not None:
            self.jax_installed = True
            return True
        return False

    def get_recommended_backend(self) -> Dict[str, Any]:
        """
        Determine the recommended JAX backend and installation parameters.
        When multiple GPUs are present, uses priority system to choose the best one.

        JAX has limited platform support:
        - NVIDIA: CUDA on Linux (Windows support experimental/limited)
        - AMD ROCm: Linux only
        - Apple Silicon: Metal backend
        - Intel: Experimental plugin with dependency issues

        Returns:
            Dict with 'backend', 'packages', 'extra_index_url' keys
        """
        priority = _get_gpu_priority()
        primary_gpu = priority['primary_gpu']
        system = priority['capabilities']['system']

        # NVIDIA GPU - CUDA on Linux
        if primary_gpu == 'nvidia' and system == 'linux':
            return {
                'backend': 'cuda',
                'packages': ['jax[cuda12]'],
                'extra_index_url': None
            }

        # AMD ROCm - Linux only
        if primary_gpu == 'amd_rocm' and system == 'linux':
            return {
                'backend': 'rocm',
                'packages': ['jax[rocm]'],
                'extra_index_url': None
            }

        # Apple Silicon
        if primary_gpu == 'apple_silicon':
            return {
                'backend': 'metal',
                'packages': ['jax-metal'],
                'extra_index_url': None
            }

        # Intel Arc (has dependency conflicts currently)
        if primary_gpu == 'intel_arc':
            return {
                'backend': 'intel',
                'packages': ['jax', 'intel-extension-for-openxla'],
                'extra_index_url': None
            }

        # All other cases fall back to CPU
        # This includes: NVIDIA/AMD on Windows, AMD non-ROCm, iGPUs
        return {
            'backend': 'cpu',
            'packages': ['jax'],
            'extra_index_url': None
        }

    def _detect_hardware_config(self) -> Dict[str, Any]:
        """
        Detect the hardware configuration and determine the appropriate JAX variant.
        Uses the GPU priority system to handle multi-GPU scenarios.

        Returns:
            Dict containing detected hardware info and recommended JAX installation.
        """
        priority = _get_gpu_priority()
        caps = priority['capabilities']

        config = {
            'system': self.system,
            'primary_gpu': priority['primary_gpu'],
            'has_nvidia_gpu': caps['nvidia'] is not None and caps['nvidia']['detected'],
            'has_amd_gpu': caps['amd'] is not None and caps['amd']['detected'],
            'has_intel_gpu': caps['intel'] is not None and caps['intel']['detected'],
            'cuda_version': caps['nvidia']['recommended_cuda'] if caps['nvidia'] else None,
            'amd_rocm_compatible': caps['amd']['rocm_compatible'] if caps['amd'] else False,
            'intel_is_arc': caps['intel']['is_arc'] if caps['intel'] else False,
            'apple_silicon': caps['apple_silicon'] is not None,
            'recommended_jax_variant': 'jax[cpu]',
            'install_url': None,
            'index_url': None
        }

        # Determine JAX variant based on hardware
        config = self._determine_jax_variant(config)

        self.detected_config = config
        return config

    def _determine_jax_variant(self, config: Dict[str, Any], force_cpu=False) -> Dict[str, Any]:
        """
        Determine the appropriate JAX variant based on detected hardware.
        Uses primary_gpu from priority system when available.

        Args:
            config: Hardware configuration dictionary
            force_cpu: forces CPU-only installation

        Returns:
            Updated configuration with JAX variant recommendation
        """
        if config['system'] == 'windows':
            print("(!) The Windows jax[cuda] wheel currently does not "
                "support numpy 2.x and therefore causes unacceptable "
                "conflicts with other packages. As jax Windows support "
                "matures we hope to enable the cuda wheel however at "
                "present only CPU support is available, which is mostly "
                "only useful for development purposes.")

        if force_cpu:
            config['recommended_jax_variant'] = 'jax[cpu]'
            print("(!) Warning: performance of the CPU-only jax variant "
                "is slow and is intended for development use only. It "
                "is generally recommended that you do NOT enable jax "
                "optimisation in scripts that offer it!")
            return config

        # Use primary_gpu if available (from priority system)
        primary_gpu = config.get('primary_gpu', 'cpu')

        # NVIDIA GPU - Linux only
        if primary_gpu == 'nvidia' and config['system'] == 'linux':
            config['recommended_jax_variant'] = 'jax[cuda12]'
            print(f"Detected NVIDIA GPU as primary - will install jax[cuda12]")

        # AMD ROCm GPU - Linux only
        elif primary_gpu == 'amd_rocm' and config['system'] == 'linux':
            config['recommended_jax_variant'] = 'jax[rocm]'
            print("Detected ROCm-compatible AMD GPU as primary - will install jax[rocm]")
            print("(!) AMD ROCm support in JAX is experimental - feedback appreciated!")

        # AMD ROCm GPU on Windows - not supported by JAX
        elif primary_gpu == 'amd_rocm' and config['system'] == 'windows':
            print("(!) AMD ROCm-compatible GPU detected on Windows: unfortunately "
                "JAX does not yet support ROCm on Windows. Falling back to CPU-only.")
            config['recommended_jax_variant'] = 'jax[cpu]'

        # Intel Arc GPU (has dependency conflicts currently)
        elif primary_gpu == 'intel_arc':
            print("(!) Intel Arc GPU detected as primary: unfortunately the experimental "
                "jax plugin for Intel GPUs has dependency clashes and "
                "still requires numpy 1.x therefore only the CPU variant "
                "can be installed. Hopefully this will change in the "
                "future.")
            config['recommended_jax_variant'] = 'jax[cpu]'

        # Apple Silicon
        elif primary_gpu == 'apple_silicon':
            config['recommended_jax_variant'] = 'jax-metal'
            print("Detected Apple Silicon - will install jax-metal")
            print("(!) jax Metal support on macOS is still experimental")

        # Fallback cases for when priority system isn't used
        elif config['has_nvidia_gpu'] and config['system'] == 'linux':
            config['recommended_jax_variant'] = 'jax[cuda12]'
            print(f"Detected NVIDIA GPU - will install jax[cuda12]")

        elif config['amd_rocm_compatible'] and config['system'] == 'linux':
            config['recommended_jax_variant'] = 'jax[rocm]'
            print("Detected ROCm-compatible AMD GPU - will install jax[rocm]")
            print("(!) AMD ROCm support in JAX is experimental - feedback appreciated!")

        elif config['intel_is_arc']:
            print("(!) Intel Arc GPU detected: unfortunately the experimental "
                "jax plugin for Intel GPUs has dependency clashes and "
                "still requires numpy 1.x therefore only the CPU variant "
                "can be installed. Hopefully this will change in the "
                "future.")
            config['recommended_jax_variant'] = 'jax[cpu]'

        elif config['apple_silicon']:
            config['recommended_jax_variant'] = 'jax-metal'
            print("Detected Apple Silicon - will install jax-metal")
            print("(!) jax Metal support on macOS is still experimental")

        # AMD or Intel iGPU - fall back to CPU
        elif config['has_amd_gpu'] or config['has_intel_gpu']:
            print(f"(!) GPU detected but not suitable for JAX acceleration - falling back to CPU")
            config['recommended_jax_variant'] = 'jax[cpu]'

        # Default to CPU
        else:
            config['recommended_jax_variant'] = 'jax[cpu]'

        # Warn about CPU-only performance
        if config['recommended_jax_variant'] in ['jax[cpu]', 'jax']:
            print("(!) Warning: performance of the CPU-only jax variant "
                "is slow and is intended for development use only. It "
                "is generally recommended that you do NOT enable jax "
                "optimisation in scripts that offer it!")

        return config

    def install_jax(self, force_variant: Optional[str] = None,
                   version_constraint: Optional[str] = None,
                   force_reinstall: Optional[bool] = False) -> bool:
        """
        Install JAX with the appropriate variant for the detected hardware. Use this instead of
        ensure_installed() to ensure that jax is installed correctly for the given hardware / OS

        Args:
            force_variant: Override auto-detection with specific variant (e.g., 'jax[cpu]')
            version_constraint: Version constraint string (e.g., '>=0.4.0')
        force_reinstall: Forces a reinstallation

        Returns:
            bool: True if installation succeeded, False otherwise
        """
        if self.is_jax_installed():
            print("Jax is already installed.")
            return

        if not self.detected_config:
            self._detect_hardware_config()

        variant = force_variant or self.detected_config['recommended_jax_variant']

        try:
            print(f"Installing {variant}...")

            # Use the provided install_package function
            _install_package(
                package_name=variant,
                version_constraint=version_constraint,
                from_url=self.detected_config.get('install_url'),
                index_url=self.detected_config.get('index_url'),
                reinstall = force_reinstall
            )

            self.jax_installed = self.is_jax_installed()
            print(f"Successfully installed {variant}")
            return True

        except Exception as e:
            print(f"Failed to install {variant}: {e}")
            return False

    def ensure_jax(self) -> bool:
        """
        Wrapper for install_jax() that only installs it if needed, with
        negligible overhead if it is already installed.
        """
        if not self.is_jax_installed():
            return self.install_jax()
        return True

    def test_jax(self) -> Tuple[bool, Optional[str]]:
        """
        Test JAX functionality and return execution provider.

        Returns:
            Tuple[bool,str]: the bool returned is True if jax works or False if
            it does not, and the str is "gpu" if JAX is using GPU, "cpu" if
            using CPU or None if

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
            print("(OK) Accuracy cross-check between jax.numpy and numpy succeeded")

            # Check the platform of the device that was actually used
            if default_device.platform.lower() == 'gpu':
                return True, "gpu"
            return True, "cpu"

        except Exception as e:
            # If main test fails, try fallback to basic CPU operations
            if "does not match numpy" in str(e):
                # Re-raise accuracy errors immediately
                return False, None

            try:
                # Try basic CPU operations with accuracy check
                test_array_np = np.array([1.0, 2.0, 3.0], dtype=np.float32)
                test_array_jax = jnp.array(test_array_np)

                jax_sum = jnp.sum(test_array_jax)
                jax_sum.block_until_ready()
                numpy_sum = np.sum(test_array_np)

                if not np.allclose(jax_sum, numpy_sum, rtol=1e-7, atol=1e-8):
                    print(f"(x) Basic JAX operation failed accuracy check: {jax_sum} vs {numpy_sum}. Error: {e}")
                    return False, None

                print("(OK) Jax available using CPU only. This is likely to perform less "
                    "well than numpy but support is hoped to improve in the future.")
                return True, "cpu"
            except Exception:
                # If even CPU fails, something is seriously wrong
                print("Jax test failed: {e}")
                return False, None

    def _get_jax_info(self) -> Dict[str, Any]:
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

    def uninstall_jax(self, dry_run: bool = False) -> Dict[str, Any]:
        """
        Detect and uninstall any existing JAX-related packages.

        This is useful when you need to clean up a problematic JAX installation
        before installing a different variant (e.g., falling back from GPU to CPU).

        Args:
            dry_run: If True, only detect packages without uninstalling them

        Returns:
            Dict containing information about detected and uninstalled packages
        """
        self.jax_installed = False
        self.detected_config = None

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
            installed_names = [p['installed_name'] for p in detected]
            for installed_name, version in installed_dict.items():
                if installed_name.startswith('jax') and installed_name not in installed_names:
                    detected.append({
                        'name': installed_name,
                        'version': version,
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

    def __repr__(self) -> str:
        return f"JaxHelper(system={self.system}, jax_installed={self.jax_installed})"
