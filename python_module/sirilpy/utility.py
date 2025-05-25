# Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
# Copyright (C) 2012-2025 team free-astro (see more in AUTHORS file)
# Reference site is https://siril.org
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Utility module for Siril Python interface providing helper functions for file operations,
package management, and standard I/O control to support Siril's scripting capabilities.
"""




import os
import sys
import io
import re
import time
import platform
import threading
import subprocess
from importlib import metadata, util
from typing import Union, List, Optional, TYPE_CHECKING, Tuple
import requests
from packaging import version as pkg_version
from packaging.specifiers import SpecifierSet
from packaging.requirements import Requirement
from .exceptions import SirilError

if TYPE_CHECKING:
    from .connection import SirilInterface

def truncate_utf8(data, max_bytes):
    """
    Truncates utf8 input. Accepts either bytes or str as input and
    returns data in the same format as the input.

    Args:
        data (bytes or str): The data to be truncated

    Returns:
        bytes or str: The truncated data

    Raises:
        TypeError: if the input was not bytes or str
    """
    if isinstance(data, bytes):
        # Truncate raw bytes without breaking UTF-8
        truncated = data[:max_bytes]
        # Try decoding, backtrack if needed to avoid decoding errors
        while True:
            try:
                truncated.decode('utf-8')
                return truncated
            except UnicodeDecodeError:
                truncated = truncated[:-1]
                if not truncated:
                    return b''
    elif isinstance(data, str):
        encoded = data.encode('utf-8')
        if len(encoded) <= max_bytes:
            return data
        # Truncate the string by character until it encodes within the limit
        for i in range(len(data) - 1, -1, -1):
            sub = data[:i]
            if len(sub.encode('utf-8')) <= max_bytes:
                return sub
        return ''
    else:
        raise TypeError("Expected str or bytes, got " + type(data).__name__)

def human_readable_size(bytes_size: int) -> str:
    """
    Convert bytes to human-readable format.

    Args:
        bytes_size (int): Size in bytes

    Returns:
        str: Formatted size with appropriate unit (B, KB, MB, GB, TB)

    Raises:
        TypeError: on incorrect input type
    """
    if not isinstance(bytes_size, int):
        raise TypeError("bytes_size must be an int")

    units = [' B', ' KB', ' MB', ' GB', ' TB']
    size = float(bytes_size)
    unit_index = 0

    while size >= 1024 and unit_index < len(units) - 1:
        size /= 1024
        unit_index += 1

    # Round to 2 decimal places, remove trailing zeros
    return f"{size:.2f}".rstrip('0').rstrip('.') + units[unit_index]

def download_with_progress(
    siril: 'SirilInterface',
    url: str,
    file_path: str,
    max_retries: int = 3,
    retry_delay: int = 5,
    resume: bool = True
    ) -> bool:
    """
    Robust file download method with native Siril progress tracking
    and error handling using retries and a resume mechanism.

    Args:
        siril (SirilInterface): SirilInterface to use to update the progress bar
        url (str): URL of the file to download
        file_path (str): Local path to save the downloaded file
        max_retries (int): Number of download retry attempts
        retry_delay (int): Delay between retry attempts in seconds
        resume (bool): Whether or not to resume a partially downloaded file or start again

    Returns:
        bool: True if download successful, False otherwise

    Raises:
        SirilError: On unhandled errors
    """
    temp_file_path = file_path + '.part'

    # If resume is False and the temporary file exists, delete it
    if not resume and os.path.exists(temp_file_path):
        os.remove(temp_file_path)

    def get_file_size_and_resume_point() -> Tuple[int, int]:
        """Determine the current file size for resuming download."""
        if os.path.exists(temp_file_path):
            return os.path.getsize(temp_file_path), 1
        return 0, 0

    for attempt in range(max_retries):
        try:
            # Get initial file size and determine if resuming
            initial_size, resume_attempt = get_file_size_and_resume_point()

            # Prepare headers for partial content
            headers = {}
            if initial_size > 0:
                headers['Range'] = f'bytes={initial_size}-'

            # Establish connection with timeout
            response = requests.get(url, stream=True, headers=headers, timeout=30)
            response.raise_for_status()

            # Determine total file size and content range
            total_size = int(response.headers.get('content-length', 0))
            if headers.get('Range'):
                # If resuming, adjust total size
                content_range = response.headers.get('Content-Range', '')
                if content_range:
                    total_size = int(content_range.split('/')[-1])

            downloaded_size = initial_size

            # Progress update rate limiting
            max_update_frequency = 5.0
            last_update_time = 0
            min_update_interval = 1 / max_update_frequency

            # Open file in append mode or write mode
            mode = 'ab' if initial_size > 0 else 'wb'
            with open(temp_file_path, mode) as f:
                for chunk in response.iter_content(chunk_size=8192):
                    if not chunk:
                        continue

                    f.write(chunk)
                    downloaded_size += len(chunk)

                    current_time = time.time()

                    # Update progress
                    if total_size > 0 and current_time - last_update_time >= min_update_interval:
                        progress = downloaded_size / total_size
                        # Line shortened to comply with line length limit
                        status = (f"Downloading... (Attempt {resume_attempt}, "
                                 f"{human_readable_size(downloaded_size)}/{human_readable_size(total_size)})")
                        siril.update_progress(status, progress)
                        last_update_time = current_time

            # Verify download completeness
            if downloaded_size >= total_size:
                # Rename temp file to final file
                os.replace(temp_file_path, file_path)
                return True

            # If download is incomplete, will retry
            time.sleep(retry_delay)

        except requests.exceptions.RequestException as e:
            # Comprehensive error handling for network-related issues
            error_message = f"Download error (Attempt {attempt + 1}/{max_retries}): {str(e)}"

            # Log or print error
            print(error_message)

            # Provide progress update for error state
            siril.update_progress(error_message, 0.0)

            # Wait before retrying, with exponential backoff
            time.sleep(retry_delay * (attempt + 1))

        except Exception as e:
            # Catch any unexpected errors
            error_message = f"Unexpected error during download: {str(e)}"
            print(error_message)

            siril.update_progress(error_message, 0.0)

            raise SirilError(error_message) from e

    # All retry attempts failed
    raise SirilError(f"Failed to download file from {url} after {max_retries} attempts")

def ensure_installed(*packages: Union[str, List[str]],
                     version_constraints: Optional[Union[str, List[str]]] = None):
    """
    Ensures that the specified package(s) are installed and meet optional version constraints.

    Args:
        *packages (str or List[str]): Name(s) of the package(s) to ensure are installed.
        version_constraints (str or List[str], optional): Version constraint string(s)
            (e.g. ">=1.5", "==2.0"). Can be a single constraint or a list matching packages.

    Returns:
        bool: True if all packages are successfully installed or already meet constraints.

    Raises:
        SirilError: If package installation fails,
        ValueError: If a different number of constraints is provided to the number
                    of packages to be installed.
        TimeoutError: If pip fails with an apparent timeout.
    """
    # Normalize inputs to lists
    if isinstance(packages[0], list):
        packages = packages[0]

    # Handle version constraints
    if version_constraints is None:
        version_constraints = [None] * len(packages)
    elif isinstance(version_constraints, str):
        version_constraints = [version_constraints] * len(packages)

    # Ensure length consistency
    if len(version_constraints) != len(packages):
        raise ValueError("Number of packages must match number of version constraints")

    # Track installation results
    all_installed = True

    for package, constraint in zip(packages, version_constraints):
        # Special handling for core/builtin modules
        if util.find_spec(package) is not None:
            continue

        try:
            # Check if package is installed and meets version constraint
            if _check_package_installed(package, constraint):
                print(f"{package} {'is installed' if constraint is None else f'meets version {constraint}'}")
                continue

            # Attempt installation
            _install_package(package, constraint)

        except TimeoutError as e:
            all_installed = False
            print(f"Timeout error processing {package}: {e}")
            raise

        except Exception as e:
            all_installed = False
            print(f"Error processing {package}: {e}")
            raise SirilError(f"Failed to install or verify package {package}") from e

    return all_installed

def _check_package_installed(package_name: str, version_constraint: Optional[str] = None) -> bool:
    """
    Check if a package is installed and meets version constraint.

    Args:
        package_name (str): Name of the package to check.
        version_constraint (str, optional): Version constraint to validate.

    Returns:
        bool: True if package is installed and meets version constraint.
    """
    try:
        # Check package existence
        installed_version = metadata.version(package_name)

        # If no version constraint, any version is fine
        if version_constraint is None:
            return True

        # Validate version constraint
        try:
            req_string = f"{package_name}{version_constraint}"
            requirement = Requirement(req_string)
            return pkg_version.parse(installed_version) in requirement.specifier

        except ImportError:
            # Fallback if packaging is not available
            print("Warning: packaging library not found. Skipping precise version check.")
            return True

    except metadata.PackageNotFoundError:
        return False

def _stream_output(process):
    """
    Helper function to stream subprocess output to stdout
    """
    for line in io.TextIOWrapper(process.stdout, encoding='utf-8', errors='replace'):
        print(line.rstrip())

def _install_package(package_name: str, version_constraint: Optional[str] = None, from_url: Optional[str] = None, index_url: Optional[str] = None):
    """
    Install a package with optional version constraint, streaming pip output to stdout.

    Args:
        package_name (str): Name of the package to install.
        version_constraint (str, optional): Version constraint for installation.
        from_url (str, optional): URL to find packages at, passed as "-f URL" to pip.
        index_url (str, optional): repository URL, passed as "--index-url URL" to pip.

    Note: the from_url and index_url parameters are not for general use and are only
    required for certain very specific circustances.

    Raises:
        subprocess.CalledProcessError: If pip installation fails.
        TimeoutError: if pip appears to have encountered a TimeOutError internally
    """
    print(f"Installing {package_name}. This may take a few seconds...")

    # Construct installation target
    install_target = f"{package_name}{version_constraint}" if version_constraint else package_name

    try:
        # Build pip command
        pip_command = [sys.executable, "-m", "pip", "install"]

        # Add index-url option if index_url is provided
        if index_url:
            pip_command.append(["--index-url", index_url])

        # Add find-links option if from_url is provided
        if from_url:
            pip_command.extend(["-f", from_url])

        # Add the package to install
        pip_command.append(install_target)

        # Start pip install process with pipe for stdout
        with subprocess.Popen(
            pip_command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            bufsize=-1,
            universal_newlines=False
        ) as process:
            # Create and start output streaming thread
            output_thread = threading.Thread(target=_stream_output, args=(process,))
            output_thread.start()

            # Wait for process to complete
            return_code = process.wait()
            output_thread.join()

            if return_code == 0:
                print(f"Successfully installed {install_target}")
            else:
                raise subprocess.CalledProcessError(return_code, process.args)

    except subprocess.CalledProcessError as e:
        print(f"Failed to install {install_target}")
        if e.stderr and "timed out" in e.stderr.lower():
            raise TimeoutError(f"Likely timeout error in pip: {e}") from e
        raise

def check_module_version(requires=None):
    """
    Check the version of the Siril module is sufficient to support the
    script. This is not mandatory if you are only using classes,
    methods etc. that are provided in the initial public release, but
    if you rely on methods that are noted int he API documentation as
    having been added at a particular version of the module then you
    must check the running sirilpy module supports your script by
    calling this function.

    Args:
        requires (str): A version format specifier string following the
                        same format used by pip, i.e. it may contain
                        '==1.2', '!=3.4', '>5.6', '>=7.8', or a
                        combination such as '>=1.2,<3.4'

    Returns:
        True if requires = None or if the available sirilpy module version
        satisfies the version specifier, otherwise False

    Raises:
        ValueError: if requires is an invalid version specifier.
    """
    import sirilpy  # pylint: disable=import-outside-toplevel

    if requires is None:
        return True  # No version requirement

    try:
        # Create a SpecifierSet from the `requires` string
        specifiers = SpecifierSet(requires)
        # Use pkg_version from top-level import
        return pkg_version.parse(sirilpy.__version__) in specifiers
    except (pkg_version.InvalidVersion, ValueError) as exc:
        raise ValueError(f"Invalid version specifier: {requires}") from exc

class SuppressedStdout:
    """
    This context manager allows suppression of the script's stdout,
    which can be useful to avoid flooding the log with stdout messages
    from an excessively verbose module used in the script.

    Example:
        .. code-block:: python

            import sirilpy as s
            siril = s.SirilInterface()
            print("This message will appear in the Siril log")
            with s.SuppressedStdout():
                print("This message will not appear")
            print("This message will appear again")

    """
    def __init__(self):
        """Initialize attributes that will be used in context management."""
        self.devnull = None
        self.original_stdout_fd = None
        self.original_stdout = None

    def __enter__(self):
        self.devnull = open(os.devnull, 'w', encoding='utf-8')
        self.original_stdout_fd = os.dup(1)  # Duplicate stdout (fd 1)
        os.dup2(self.devnull.fileno(), 1)  # Redirect stdout to devnull
        self.original_stdout = sys.stdout
        sys.stdout = self.devnull  # Also redirect Python stdout
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        os.dup2(self.original_stdout_fd, 1)  # Restore stdout
        os.close(self.original_stdout_fd)
        sys.stdout = self.original_stdout  # Restore Python stdout
        self.devnull.close()


class SuppressedStderr:
    """
    This context manager allows suppression of the script's stderr, which
    can be useful if you are using module functions that are known to
    produce warnings that you want to avoid distracting the user with,
    such as FutureWarnings of features that have become deprecated but
    are in a dependency rather than your own code. The class should
    be used sparingly and should **not** be used to hide evidence of
    bad code.
    """
    def __init__(self):
        """Initialize attributes that will be used in context management."""
        self.devnull = None
        self.original_stderr_fd = None
        self.original_stderr = None

    def __enter__(self):
        self.devnull = open(os.devnull, 'w', encoding='utf-8')
        self.original_stderr_fd = os.dup(2)  # Duplicate stderr (fd 2)
        os.dup2(self.devnull.fileno(), 2)  # Redirect stderr to devnull
        self.original_stderr = sys.stderr
        sys.stderr = self.devnull  # Also redirect Python stderr
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        os.dup2(self.original_stderr_fd, 2)  # Restore stderr
        os.close(self.original_stderr_fd)
        sys.stderr = self.original_stderr  # Restore Python stderr
        self.devnull.close()

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
        self.system = platform.system().lower()

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
        except ImportError:
            # If import fails, package is not usable regardless of pip list
            return False, None
        except Exception as e:
            print(f"Error checking for installed onnxruntime: {e}")
            return False, None

    def get_onnxruntime_package(self):
        """
        Determine the appropriate ONNX Runtime package based on system and available hardware.

        Returns:
            tuple: (package_name, url) where url is None except for special cases like ROCm
        """
        from_url = None
        index_url = None
        onnxruntime_pkg = "onnxruntime"
        cuda_version = self._detect_cuda_version()
        if self.system == "windows":
                # NVidia GPU runtime that supports CUDA
            if self._detect_nvidia_gpu() and pkg_version.Version(cuda_version).major >= 11:
                onnxruntime_pkg = "onnxruntime-gpu"
                if pkg_version.Version(cuda_version).major == 11:
                    index_url = "https://aiinfra.pkgs.visualstudio.com/PublicPackages/_packaging/onnxruntime-cuda-11/pypi/simple/"
            else:
                # DirectML provides hardware acceleration for various GPUs on Windows
                onnxruntime_pkg = "onnxruntime-directml"

        elif self.system == "darwin":  # macOS
            onnxruntime_pkg = "onnxruntime"

        elif self.system == "linux":
            if self._detect_nvidia_gpu() and pkg_version.Version(cuda_version).major >= 11:
                onnxruntime_pkg = "onnxruntime-gpu"
                if pkg_version.Version(cuda_version).major == 11:
                    index_url = "https://aiinfra.pkgs.visualstudio.com/PublicPackages/_packaging/onnxruntime-cuda-11/pypi/simple/"
            elif self._detect_amd_gpu():
                # ROCm support for AMD GPUs on Linux
                onnxruntime_pkg = "onnxruntime-rocm"
                from_url = "https://repo.radeon.com/rocm/manylinux/rocm-rel-6.4/"
            elif self._detect_intel_gpu_for_openvino():
                onnxruntime_pkg = "onnxruntime-openvino"

        return onnxruntime_pkg, from_url, index_url

    def _detect_cuda_version(self) -> Optional[str]:
        """
        Detects the CUDA version by parsing the output of 'nvcc -v'.

        Returns:
            Optional[str]: The CUDA version as a string (e.g., '11.7') if detected,
                          or "0.0" if nvcc is not installed or version cannot be determined.
        """
        try:
            nvcc_command = 'nvcc.exe' if self.system == 'windows' else 'nvcc'

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
            return "0.0"

    def _detect_nvidia_gpu(self):
        """
        Detect if NVIDIA GPU is available. (Required to check Windows and Linux, not required on MacOS.)

        Returns:
            bool: True if NVIDIA GPU is detected, False otherwise.
        """
        try:
            if self.system == "windows":
                # Windows detection using PowerShell
                output = subprocess.check_output(
                    ["powershell", "-Command", "Get-WmiObject -Class Win32_VideoController"],
                    text=True
                )
                return "NVIDIA" in output
            # Linux and macOS detection using lspci or similar
            if self.system == "linux":
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

    def _detect_amd_gpu(self):
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

    def _detect_intel_gpu_for_openvino(self):
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
        def has_tensorrt_library():
            """
            Checks if TensorRT libraries are available and can be loaded.
            """
            import ctypes.util
            import ctypes

            lib_names = {
                "linux": "nvinfer",
                "windows": "nvinfer",
                "darwin": None  # TensorRT isn't available on macOS
            }
            lib_name = lib_names.get(platform.system().lower())
            if not lib_name:
                return False

            path = ctypes.util.find_library(lib_name)
            if path:
                try:
                    ctypes.CDLL(path)
                    return True
                except OSError:
                    return False
            return False

        import onnxruntime
        providers = []
        available_providers = onnxruntime.get_available_providers()

        if ai_gpu_acceleration:
            if self.system == "darwin":
#                is_apple_silicon = (
#                    platform.machine().startswith(("arm", "aarch"))
#                )
                # Apple Silicon / Neural Engine (faster if available and supports the model features)
#                if is_apple_silicon:
#                    if "CoreMLExecutionProvider" in available_providers:
                        # On Apple silicon it's safe to create a MLPROGRAM instead of the older MLMODEL
#                        providers.append(
#                            (
#                                "CoreMLExecutionProvider",
#                                {
#                                    "flags": "COREML_FLAG_CREATE_MLPROGRAM",
#                                },
#                            )
#                        )
#                else:
#                    # On Intel silicon we omit the MLPROGRAM flag
#                    providers.append("CoreMLExecutionProvider")
                providers.append("CoreMLExecutionProvider")

            elif self.system == "windows":
                # NVIDIA TensorRT - best performance if the model features are supported
                if "TensorrtExecutionProvider" in available_providers and has_tensorrt_library():
                    providers.append("TensorrtExecutionProvider")

                # NVIDIA GPU - good fallback option, still much faster than CPU
                if "CUDAExecutionProvider" in available_providers:
                    providers.append("CUDAExecutionProvider")

                # DirectML for Windows GPUs (supports various vendors)
                if "DmlExecutionProvider" in available_providers:
                    providers.append("DmlExecutionProvider")

            elif self.system == "linux":
                # NVIDIA TensorRT - best performance if the model features are supported
                if "TensorrtExecutionProvider" in available_providers and has_tensorrt_library():
                    providers.append("TensorrtExecutionProvider")

                # NVIDIA GPU - good fallback option, still much faster than CPU
                if "CUDAExecutionProvider" in available_providers:
                    providers.append("CUDAExecutionProvider")

                # AMD GPU - fastest possible but commented out until not marked as experimental
                #if "MIGraphXExecutionProvider" in available_providers:
                #    providers.append("MIGraphXExecutionProvider")

                # AMD GPU
                if "ROCmExecutionProvider" in available_providers:
                    providers.append("ROCmExecutionProvider")

                # Intel GPU via OpenVINO
                if "OpenVINOExecutionProvider" in available_providers:
                    providers.append("OpenVINOExecutionProvider")

        # CPU is always the fallback option
        providers.append("CPUExecutionProvider")

        return providers

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

def parse_fits_header(header_text: str) -> dict:
    """
    Parse FITS header from text content into a dictionary
    
    Parameters:
    header_text (str): Content of the FITS header text file
    
    Returns:
    dict: Dictionary containing all header keywords and values
    """
    header_dict = {}
    
    for line in header_text.split('\n'):
        # Skip empty lines, COMMENT, HISTORY, and END
        if not line.strip() or line.startswith('COMMENT') or line.startswith('HISTORY') or line.startswith('END'):
            continue
            
        # Split the line into key and value parts
        parts = line.split('=')
        if len(parts) != 2:
            continue
            
        key = parts[0].strip()
        value_part = parts[1].strip()
        
        # Handle the value part (removing comments after /)
        if '/' in value_part:
            value_part = value_part.split('/')[0].strip()
            
        # Convert value to appropriate type
        try:
            # Try converting to float first
            if value_part.startswith("'") and value_part.endswith("'"):
                # String value
                value = value_part.strip("'").strip()
            elif value_part == 'T':
                value = True
            elif value_part == 'F':
                value = False
            else:
                value = float(value_part)
        except ValueError:
            # If conversion fails, keep as string
            value = value_part.strip()
            
        header_dict[key] = value
        
    return header_dict
