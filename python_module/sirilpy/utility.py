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
            response = requests.get(url, stream=True, headers=headers, timeout=(10, 30))
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
                     version_constraints: Optional[Union[str, List[str]]] = None,
                     reinstall: Optional[bool] = False):
    """
    Ensures that the specified package(s) are installed and meet optional version constraints.

    Args:
        *packages (str or List[str]): Name(s) of the package(s) to ensure are installed.
        version_constraints (str or List[str], optional): Version constraint string(s)
            (e.g. ">=1.5", "==2.0"). Can be a single constraint or a list matching packages.
        reinstall (bool, optional): Forces reinstallation. Defaults to False.

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
            _install_package(package, constraint, reinstall=reinstall)

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

def _install_package(package_name: str, version_constraint: Optional[str] = None, from_url: Optional[str] = None,
                     index_url: Optional[str] = None, reinstall: Optional[bool] = False, nodeps: Optional[bool] = False):
    """
    Install a package with optional version constraint, streaming pip output to stdout.

    Args:
        package_name (str): Name of the package to install.
        version_constraint (str, optional): Version constraint for installation.
        from_url (str, optional): URL to find packages at, passed as "-f URL" to pip.
        index_url (str, optional): repository URL, passed as "--index-url URL" to pip.
        reinstall (bool optional): whether to force reinstallation. Defaults to False.
        nodeps (bool, optional): whether to install without dependencies. Defaults to
            False, but occasionally useful together with reinstall=True for problematic
            packages such as torch.

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
            pip_command.extend(["--index-url", index_url])

        # Add find-links option if from_url is provided
        if from_url:
            pip_command.extend(["-f", from_url])

        # If required, add the --force-reinstall flag
        if reinstall:
            pip_command.append("--force-reinstall")

        # If required, add the --no-deps flag
        if nodeps:
            pip_command.append("--no-deps")

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
