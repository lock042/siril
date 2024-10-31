import os
import sys
import struct
import socket
import ctypes
import time
import mmap
from datetime import datetime
import calendar
from pathlib import Path
from enum import IntEnum
from typing import Tuple, Optional, List, Union
import numpy as np
from .translations import _, N_
from .shm import SharedMemoryWrapper
from .exceptions import SirilError, ConnectionError, CommandError, DataError, NoImageError
from .models import DataType, ImageStats, FKeywords, FFit, Homography, StarProfile, PSFStar, RegData, ImgData, Sequence, SequenceType

if os.name == 'nt':
    import win32pipe
    import win32file
    import win32event
    import pywintypes
    import winerror

# Developer-facing strings such as enums are not given translatable
# doctext as they are not intended for end users. All user-facing
# doctext strings should be made translatable.
# "debug" strings should be handled using print() and directed to
# stderr, so that they can be presented to the Siril log using the
# configurable stderr handler (which defaults to disabled but is
# intended as a debugging aid)

class _Status(IntEnum):
    """ Returns the status of a command. NONE is for commands that
    may legitimately fail to return data but which should not be
    regarded as an error, instead this triggers the command processor
    to return the special python value None"""

    OK = 0
    NONE = 1
    ERROR = 0xFF

class _Command(IntEnum):
    """ Enumerates the commands. This enum MUST match the one in
    siril_pythonmodule.h"""
    SEND_COMMAND = 1
    LOG_MESSAGE = 2
    UPDATE_PROGRESS = 3
    GET_WORKING_DIRECTORY = 4
    GET_FILENAME = 5
    GET_DIMENSIONS = 6
    GET_PIXELDATA = 7
    GET_PIXELDATA_REGION = 8
    RELEASE_SHM = 9
    SET_PIXELDATA = 10
    GET_IMAGE_STATS = 11
    GET_KEYWORDS = 12
    GET_ICC_PROFILE = 13
    GET_FITS_HEADER = 14
    GET_FITS_HISTORY = 15
    GET_FITS_UNKNOWN_KEYS = 16
    GET_IMAGE = 17
    GET_PSFSTARS = 18
    GET_SEQ_STATS = 19
    GET_SEQ_REGDATA = 20
    GET_SEQ_IMGDATA = 21
    GET_SEQ_PIXELDATA = 22
    GET_SEQ_IMAGE = 23
    GET_SEQ = 24
    GET_CONFIG = 25
    ERROR = 0xFF

class ConfigType(IntEnum):
    """Enumerates config variable types for use with the
    get_config() method"""
    BOOL = 0
    INT = 1
    DOUBLE = 2
    STR = 3
    STRDIR = 4
    STRLIST = 5

class SharedMemoryInfo(ctypes.Structure):
    """Structure matching the C-side shared memory info"""
    _fields_ = [
        ("size", ctypes.c_size_t),
        ("data_type", ctypes.c_int),  # 0 for WORD, 1 for float
        ("width", ctypes.c_int),
        ("height", ctypes.c_int),
        ("channels", ctypes.c_int),
        ("shm_name", ctypes.c_char * 256)
    ]

class SirilInterface:
    def __init__(self):
        __doc__ = N_("Initialize the SirilInterface, automatically determining the "
            "correct pipe or socket path based on the environment variable and "
            "operating system.")
        if os.name == 'nt':  # Windows
            self.pipe_path = os.getenv('MY_PIPE')
            if not self.pipe_path:
                raise ConnectionError(_("Environment variable MY_PIPE not set"))
            self.event_pipe_path = self.pipe_path  # Assuming event pipe is the same path
        else:  # POSIX (Linux/macOS)
            self.socket_path = os.getenv('MY_SOCKET')
            if not self.socket_path:
                raise ConnectionError(_("Environment variable MY_SOCKET not set"))
            self.event_pipe_path = self.socket_path  # Assuming event socket is the same path

    def connect(self):
        __doc__ = N_("Establish a connection to Siril based on the pipe or "
            "socket path. Returns True if the connection is successful, otherwise "
            "False.")
        try:
            if os.name == 'nt':
                print(f"Connecting to Windows pipe: {self.pipe_path}")
                try:
                    self.pipe_handle = win32file.CreateFile(
                        self.pipe_path,
                        win32file.GENERIC_READ | win32file.GENERIC_WRITE,
                        0,  # No sharing
                        None,  # Default security
                        win32file.OPEN_EXISTING,
                        win32file.FILE_FLAG_OVERLAPPED,  # Use overlapped I/O
                        None
                    )
                    # Create event for overlapped I/O
                    self.overlap_read = pywintypes.OVERLAPPED()
                    self.overlap_read.hEvent = win32event.CreateEvent(None, True, False, None)
                    self.overlap_write = pywintypes.OVERLAPPED()
                    self.overlap_write.hEvent = win32event.CreateEvent(None, True, False, None)
                    return True
                except pywintypes.error as e:
                    if e.winerror == winerror.ERROR_PIPE_BUSY:
                        raise ConnectionError(_("Pipe is busy"))
                    raise ConnectionError(_("Failed to connect to pipe: {}").format(e))
            else:
                self.sock = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
                self.sock.connect(self.socket_path)
                return True

        except Exception as e:
            raise ConnectionError(_("Failed to connect: {}").format(e))

    def close_connection(self):
        """
        Closes the established socket or pipe connection.
        """
        if os.name == 'nt':
            if hasattr(self, 'pipe_handle'):
                # Close the pipe handle
                win32file.CloseHandle(self.pipe_handle)
                # Close the event handles
                if hasattr(self, 'overlap_read'):
                    win32event.CloseHandle(self.overlap_read.hEvent)
                if hasattr(self, 'overlap_write'):
                    win32event.CloseHandle(self.overlap_write.hEvent)
                return True
            else:
                raise SirilError(_("No pipe connection to close"))
        else:
            if hasattr(self, 'sock'):
                self.sock.close()
                return True
            else:
                raise SirilError(_("No socket connection to close"))

    def _recv_exact(self, n: int, timeout: float = 5.0) -> Optional[bytes]:
        """
        Helper method to receive exactly n bytes from the socket or pipe.
        """
        if n < 0:
            raise ValueError(_("Cannot receive negative number of bytes"))

        if os.name == 'nt':
            try:
                data = bytearray()
                timeout_ms = int(timeout * 1000)

                while len(data) < n:
                    # Calculate remaining bytes to read
                    to_read = n - len(data)

                    # Prepare buffer for reading
                    buf = win32file.AllocateReadBuffer(to_read)

                    # Start overlapped read
                    win32file.ReadFile(self.pipe_handle, buf, self.overlap_read)

                    # Wait for completion or timeout
                    rc = win32event.WaitForSingleObject(self.overlap_read.hEvent, timeout_ms)

                    if rc == win32event.WAIT_TIMEOUT:
                        # Cancel the I/O operation
                        win32file.CancelIo(self.pipe_handle)
                        raise ConnectionError(_("Timeout while receiving data"))

                    if rc != win32event.WAIT_OBJECT_0:
                        raise ConnectionError(_("Error waiting for pipe read completion"))

                    # Get results of the operation
                    err, bytes_read = win32file.GetOverlappedResult(self.pipe_handle, self.overlap_read, False)
                    if bytes_read == 0:
                        raise ConnectionError(_("Pipe closed during read"))

                    # Extend our data buffer
                    data.extend(buf[:bytes_read])

                return bytes(data)

            except pywintypes.error as e:
                raise ConnectionError(_("Windows pipe error during receive: {}").format(e))

        else:
            # Existing socket implementation
            original_timeout = self.sock.gettimeout()
            self.sock.settimeout(timeout)

            try:
                data = bytearray()
                while len(data) < n:
                    try:
                        packet = self.sock.recv(n - len(data))
                        if not packet:
                            raise ConnectionError(_("Connection closed during data transfer"))
                        data.extend(packet)
                    except socket.timeout:
                        raise ConnectionError(_("Timeout while receiving data"))
                    except Exception as e:
                        raise ConnectionError(_("Error receiving data: {}").format(e))
                return bytes(data)
            finally:
                self.sock.settimeout(original_timeout)

    def _map_shared_memory(self, name: str, size: int) -> SharedMemoryWrapper:
        """
        Create or open a shared memory mapping using SharedMemoryWrapper.

        Args:
            name: Name of the shared memory segment
            size: Size of the shared memory segment in bytes

        Returns:
            SharedMemoryWrapper: A wrapper object for the shared memory segment

        Raises:
            RuntimeError: If the shared memory mapping fails
        """
        try:
            return SharedMemoryWrapper(name=name, size=size)
        except Exception as e:
            raise RuntimeError(_("Failed to create shared memory mapping: {}").format(e))

    def _send_command(self, command: _Command, data: Optional[bytes] = None) -> Tuple[Optional[int], Optional[bytes]]:
        """
        Send a command to Siril and receive the response.
        """
        try:
            data_length = len(data) if data else -1

            if data_length > 65529:
                raise RuntimeError(_("Command data too long. Maximum command data 65529 bytes"))

            header = struct.pack('!Bi', command, data_length)

            if os.name == 'nt':
                try:
                    # Send header
                    err, bytes_written = win32file.WriteFile(self.pipe_handle, header, self.overlap_write)
                    rc = win32event.WaitForSingleObject(self.overlap_write.hEvent, 3000)  # 3 second timeout
                    if rc != win32event.WAIT_OBJECT_0:
                        raise ConnectionError(_("Timeout or error while sending header"))

                    # Send data if present
                    if data and data_length > 0:
                        err, bytes_written = win32file.WriteFile(self.pipe_handle, data, self.overlap_write)
                        rc = win32event.WaitForSingleObject(self.overlap_write.hEvent, 3000)
                        if rc != win32event.WAIT_OBJECT_0:
                            raise ConnectionError(_("Timeout or error while sending data"))

                    # Receive response header
                    response_header = self._recv_exact(5)  # 1 byte status + 4 bytes length
                    if not response_header:
                        return None, None

                    status, response_length = struct.unpack('!BI', response_header)

                    # Receive response data if any
                    response_data = None
                    if response_length > 0:
                        response_data = self._recv_exact(response_length)
                        if not response_data:
                            return None, None

                    return status, response_data

                except pywintypes.error as e:
                    raise CommandError(_("Windows pipe error during command: {}").format(e))
            else:
                # Existing socket implementation
                self.sock.sendall(header)
                if data and data_length > 0:
                    self.sock.sendall(data)

                response_header = self._recv_exact(5)
                if not response_header:
                    return None, None

                status, response_length = struct.unpack('!BI', response_header)
                response_data = None
                if response_length > 0:
                    response_data = self._recv_exact(response_length)
                    if not response_data:
                        return None, None

                return status, response_data

        except Exception as e:
            raise CommandError(_("Error sending command: {}").format(e))

# execute_command and request_data are not intended to be used directly in scripts
# They are used to implement more user-friendly commands (see below)

    def execute_command(self, command: _Command, payload: Optional[bytes] = None) -> bool:
        """
        High-level method to execute a command and handle the response.

        Args:
            command: The command to execute
            payload: Optional command payload

        Returns:
            True if command was successful, False otherwise
        """
        status, response = self._send_command(command, payload)

        if status is None:
            return False

        if status == _Status.NONE:
            # This indicates "allowed failure" - no data to return but not an error
            return None

        if status == _Status.ERROR:
            error_msg = response.decode('utf-8') if response else _("Unknown error")
            print(f"Command failed: {error_msg}", file=sys.stderr)
            return False

        return True

    def request_data(self, command: _Command, payload: Optional[bytes] = None) -> Optional[bytes]:
        """
        High-level method to request small-volume data from Siril. The
        payload limit is 63336 bytes. For commands expected to return
        larger volumes of data, SHM should be used.

        Args:
            command: The data request command
            payload: Optional request parameters

        Returns:
            Requested data or None if error
        """
        status, response = self._send_command(command, payload)

        if status is None:
            return None

        if status == _Status.NONE:
            # This indicates "allowed failure" - no data to return but not an error
            return None

        if status is None or status == _Command.ERROR:
            error_msg = response.decode('utf-8') if response else _("Unknown error")
            print(f"Data request failed: {error_msg}", file=sys.stderr)
            return None

        return response

    # Specific commands follow below here

    def log(self, my_string: str) -> bool:
        """
        Send a log message to Siril. The maximum message length is 1022 bytes:
        longer messages will be truncated.

        Args:
            my_string: The message to log

        Returns:
            bool: True if the message was successfully logged, False otherwise
        """
        try:
            # Append a newline character to the string
            truncated_string = my_string[:1021] + '\n'
            # Convert string to bytes using UTF-8 encoding
            message_bytes = truncated_string.encode('utf-8')
            return self.execute_command(_Command.LOG_MESSAGE, message_bytes)

        except Exception as e:
            print(f"Error sending log message: {e}", file=sys.stderr)
            return False

    def update_progress(self, message: str, progress: float) -> bool:
        """
        Send a progress update to Siril with a message and completion percentage.

        Args:
            message: Status message to display
            progress: Progress value between 0.0 and 1.0

        Returns:
            bool: True if the progress update was successfully sent, False otherwise
        """
        try:
            # Validate progress value
            if not 0.0 <= progress <= 1.0:
                print(_("Progress value must be between 0.0 and 1.0"), file=sys.stderr)
                return False

            # Convert string to UTF-8 bytes
            message_bytes = message.encode('utf-8')

            # Create payload: network-order float followed by string
            # '!f' for network byte order 32-bit float
            float_bytes = struct.pack('!f', progress)

            # Combine float and string bytes
            payload = float_bytes + message_bytes

            return self.execute_command(_Command.UPDATE_PROGRESS, payload)

        except Exception as e:
            print(f"Error sending progress update: {e}", file=sys.stderr)
            return False

    def reset_progress(self) -> bool:
        """
        Send a progress update to Siril resetting the progress bar.

        Args:
            none

        Returns:
            bool: True if the progress update was successfully sent, False otherwise
        """
        return self.update_progress("", 0.0)

    def cmd(self, *args: str) -> bool:
        """
        Send a command to Siril by combining multiple arguments into a single string.

        Args:
            *args: Variable number of string arguments to be combined into a command

        Returns:
            bool: True if the command was successfully executed, False otherwise

        Example:
            siril.cmd("ght", "-D=0.5", "-b=2.0")
        """
        try:
            # Join arguments with spaces between them
            command_string = " ".join(str(arg) for arg in args)

            # Convert to bytes for transmission
            command_bytes = command_string.encode('utf-8')

            return self.execute_command(_Command.SEND_COMMAND, command_bytes)

        except Exception as e:
            print(f"Error sending command: {e}", file=sys.stderr)
            return False

    def get_width(self) -> Optional[int]:
        """
        Request the shape of the image from Siril.

        Returns:
            The width of the image as an integer, or None if an error occurred.
        """
        response = self.request_data(_Command.GET_DIMENSIONS)

        if response is None:
            return None

        try:
            # Assuming the response is in the format: width (4 bytes), height (4 bytes), nb_channels (4 bytes)
            width, height, channels = struct.unpack('!III', response)
            return width
        except struct.error as e:
            print(f"Error unpacking image dimensions: {e}", file=sys.stderr)
            return None

    def get_height(self) -> Optional[int]:
        """
        Request the shape of the image from Siril.

        Returns:
            The width of the image as an integer, or None if an error occurred.
        """
        response = self.request_data(_Command.GET_DIMENSIONS)

        if response is None:
            return None

        try:
            # Assuming the response is in the format: width (4 bytes), height (4 bytes), nb_channels (4 bytes)
            width, height, channels = struct.unpack('!III', response)
            return width
        except struct.error as e:
            print(f"Error unpacking image dimensions: {e}", file=sys.stderr)
            return None

    def get_channels(self) -> Optional[int]:
        """
        Request the shape of the image from Siril.

        Returns:
            The width of the image as an integer, or None if an error occurred.
        """
        response = self.request_data(_Command.GET_DIMENSIONS)

        if response is None:
            return None

        try:
            # Assuming the response is in the format: width (4 bytes), height (4 bytes), nb_channels (4 bytes)
            width, height, channels = struct.unpack('!III', response)
            return channels
        except struct.error as e:
            print(f"Error unpacking image dimensions: {e}", file=sys.stderr)
            return None

    def get_shape(self) -> Optional[Tuple[int, int, int]]:
        """
        Request the shape of the image from Siril.

        Returns:
            A tuple (height, width, channels) representing the shape of the image,
            or None if an error occurred.
        """
        response = self.request_data(_Command.GET_DIMENSIONS)

        if response is None:
            return None

        try:
            # Assuming the response is in the format: width (4 bytes), height (4 bytes), nb_channels (4 bytes)
            width, height, channels = struct.unpack('!III', response)
            return height, width, channels  # Returning as (height, width, channels)
        except struct.error as e:
            print(f"Error unpacking image dimensions: {e}", file=sys.stderr)
            return None

    def get_pixel_data(self, shape: Optional[list[int]] = None) -> Optional[np.ndarray]:
        """
        Retrieve the pixel data using shared memory.

        Args:
            shape: Optional list of [x, y, w, h] specifying the region to retrieve.
                If provided, gets data for just that region.
                If None, gets data for the entire image.

        Returns:
            numpy.ndarray: The image data as a numpy array

        Raises:
            NoImageError: If no image is currently loaded
            RuntimeError: For other errors during pixel data retrieval
            ValueError: If the received data format is invalid or shape is invalid
        """
        shm = None
        try:
            # Validate shape if provided
            if shape is not None:
                if len(shape) != 4:
                    raise ValueError("Shape must be a list of [x, y, w, h]")
                if any(not isinstance(v, int) for v in shape):
                    raise ValueError("All shape values must be integers")
                if any(v < 1 for v in shape):
                    raise ValueError("All shape values must be non-negative")

                # Pack shape data for the command
                shape_data = struct.pack('!IIII', *shape)
                command = _Command.GET_PIXELDATA_REGION
            else:
                shape_data = None
                command = _Command.GET_PIXELDATA

            # Request shared memory setup
            status, response = self._send_command(command, shape_data)

            # Handle error responses
            if status == _Command.ERROR:
                if response:
                    error_msg = response.decode('utf-8', errors='replace')
                    if "no image loaded" in error_msg.lower():
                        raise NoImageError("No image is currently loaded in Siril")
                    else:
                        raise RuntimeError(_("Server error: {}").formt(error_msg))
                else:
                    raise RuntimeError("Failed to initiate shared memory transfer: Empty response")

            if not response:
                raise RuntimeError("Failed to initiate shared memory transfer: No data received")

            try:
                # Parse the shared memory information
                shm_info = SharedMemoryInfo.from_buffer_copy(response)
            except (AttributeError, BufferError, ValueError) as e:
                raise ValueError(_("Invalid shared memory information received: {}").format(e))

            # Validate dimensions
            if any(dim <= 0 for dim in (shm_info.width, shm_info.height)):
                raise ValueError(_("Invalid image dimensions: {}x{}").format(shm_info.width, shm_info.height))

            if shm_info.channels <= 0 or shm_info.channels > 3:
                raise ValueError(_("Invalid number of channels: {}").format(shm_info.channels))

            # Map the shared memory
            try:
                shm = self._map_shared_memory(
                    shm_info.shm_name.decode('utf-8'),
                    shm_info.size
                )
            except (OSError, ValueError) as e:
                raise RuntimeError(_("Failed to map shared memory: {}").format(e))

            # Create numpy array from shared memory
            dtype = np.float32 if shm_info.data_type == 1 else np.uint16
            try:
                buffer = memoryview(shm.buf).cast('B')  # Use shm.buf instead of shm directly
                arr = np.frombuffer(buffer, dtype=dtype)
            except (BufferError, ValueError, TypeError) as e:
                raise RuntimeError(_("Failed to create array from shared memory: {}").format(e))

            # Validate array size matches expected dimensions
            expected_size = shm_info.width * shm_info.height * shm_info.channels
            if arr.size != expected_size:
                raise ValueError(
                    f"Data size mismatch: got {arr.size} elements, "
                    f"expected {expected_size} for dimensions "
                    f"{shm_info.width}x{shm_info.height}x{shm_info.channels}"
                )

            # Reshape the array according to the image dimensions
            try:
                if shm_info.channels > 1:
                    arr = arr.reshape((shm_info.height, shm_info.width, shm_info.channels))
                else:
                    arr = arr.reshape((shm_info.height, shm_info.width))
            except ValueError as e:
                raise ValueError(_("Failed to reshape array to image dimensions: {}").format(e))

            # Make a copy of the data since we'll be releasing the shared memory
            result = np.copy(arr)

            return result

        except NoImageError:
            # Re-raise NoImageError without wrapping
            raise
        except Exception as e:
            # Wrap all other exceptions with context
            raise RuntimeError(_("Error retrieving pixel data: {}").format(e)) from e
        finally:
            # Clean up shared memory using the wrapper's methods
            if shm is not None:
                try:
                    # Signal that Python is done with the shared memory and wait for C to finish
                    finish_info = struct.pack('256s', shm_info.shm_name)
                    if not self.execute_command(_Command.RELEASE_SHM, finish_info):
                        raise RuntimeError(_("Failed to cleanup shared memory"))

                    shm.close()  # First close the memory mapping
                    shm.unlink()  # Then unlink/remove the shared memory segment

                except Exception:
                    pass

    def set_pixel_data(self, image_data: np.ndarray) -> bool:
        """
        Send image data to Siril using shared memory.

        Args:
            image_data: numpy.ndarray containing the image data.
                    Must be 2D (single channel) or 3D (multi-channel) array
                    with dtype either np.float32 or np.uint16.

        Returns:
            bool: True if successful, False otherwise
        """
        shm = None
        try:
            # Validate input array
            if not isinstance(image_data, np.ndarray):
                raise ValueError(_("Image data must be a numpy array"))

            if image_data.ndim not in (2, 3):
                raise ValueError(_("Image must be 2D or 3D array"))

            if image_data.dtype not in (np.float32, np.uint16):
                raise ValueError(_("Image data must be float32 or uint16"))

            # Get dimensions
            if image_data.ndim == 2:
                height, width = image_data.shape
                channels = 1
                image_data = image_data.reshape(height, width, 1)
            else:
                height, width, channels = image_data.shape

            if channels > 3:
                raise ValueError(_("Image cannot have more than 3 channels"))

            if any(dim <= 0 for dim in (width, height)):
                raise ValueError(_("Invalid image dimensions: {}x{}").format(width, height))

            # Calculate total size
            element_size = 4 if image_data.dtype == np.float32 else 2
            total_bytes = width * height * channels * element_size

            # Generate unique name for shared memory
            timestamp = int(time.time() * 1000)  # Millisecond precision
            shm_name = f"siril_shm_{os.getpid()}_{timestamp}"
            if sys.platform == 'win32':
                shm_name = shm_name[1:]  # Remove leading slash on Windows

            # Create shared memory using our wrapper
            try:
                shm = SharedMemoryWrapper(shm_name, total_bytes)
            except Exception as e:
                raise RuntimeError(_("Failed to create shared memory: {}").format(e))

            # Copy data to shared memory
            try:
                buffer = memoryview(shm.buf).cast('B')
                shared_array = np.frombuffer(buffer, dtype=image_data.dtype).reshape(image_data.shape)
                np.copyto(shared_array, image_data)
            except Exception as e:
                raise RuntimeError(_("Failed to copy data to shared memory: {}").format(e))

            # Pack the image info structure
            info = struct.pack(
                '!IIIIQ256s',
                width,
                height,
                channels,
                1 if image_data.dtype == np.float32 else 0,
                total_bytes,
                shm_name.encode('utf-8').ljust(256, b'\x00')
            )

            # Send command using the existing execute_command method
            if not self.execute_command(_Command.SET_PIXELDATA, info):
                raise RuntimeError("Failed to send pixel data command")

            return True

        except Exception as e:
            print("Error sending pixel data: {e}", file=sys.stderr)
            return False

        finally:
            if shm is not None:
                try:
                    shm.close()
                    shm.unlink()
                except:
                    pass

    def get_icc_profile(self) -> Optional[bytes]:
        """
        Retrieve the ICC profile of the current Siril image using shared memory.

        Args:
            none.

        Returns:
            bytes: The image ICC profile as a byte array, or None if the current
            image has no ICC profile.

        Raises:
            NoImageError: If no image is currently loaded
            RuntimeError: For other errors during  data retrieval
            ValueError: If the received data format is invalid or shape is invalid
        """
        shm = None
        try:
            # Request shared memory setup
            status, response = self._send_command(_Command.GET_ICC_PROFILE)

            # Handle error responses
            if status == _Status.ERROR:
                if response:
                    error_msg = response.decode('utf-8', errors='replace')
                    if "no image loaded" in error_msg.lower():
                        raise NoImageError(_("No image is currently loaded in Siril"))
                    else:
                        raise RuntimeError(_("Server error: {}").format(error_msg))
                else:
                    raise RuntimeError("Failed to initiate shared memory transfer: Empty response")

            if not response:
                raise RuntimeError("Failed to initiate shared memory transfer: No data received")

            if len(response) < 25: # No payload
                return None

            try:
                # Parse the shared memory information
                shm_info = SharedMemoryInfo.from_buffer_copy(response)
            except (AttributeError, BufferError, ValueError) as e:
                raise ValueError(_("Invalid shared memory information received: {}").format(e))

            # Map the shared memory
            try:
                shm = self._map_shared_memory(
                    shm_info.shm_name.decode('utf-8'),
                    shm_info.size
                )
            except (OSError, ValueError) as e:
                raise RuntimeError(_("Failed to map shared memory: {}").format(e))

            try:
                result = bytes(memoryview(shm.buf).cast('B'))
            except (BufferError, ValueError, TypeError) as e:
                raise RuntimeError(_("Failed to create bytes from shared memory: {}").format(e))

            return result

        except NoImageError:
            # Re-raise NoImageError without wrapping
            raise
        except Exception as e:
            # Wrap all other exceptions with context
            raise RuntimeError(_("Error retrieving ICC data: {}").format(e)) from e
        finally:
            if shm is not None:
                try:
                    # Signal that Python is done with the shared memory and wait for C to finish
                    finish_info = struct.pack('256s', shm_info.shm_name)
                    if not self.execute_command(_Command.RELEASE_SHM, finish_info):
                        raise RuntimeError(_("Failed to cleanup shared memory"))
                    shm.close()
                    shm.unlink()
                except BufferError:
                    pass

    def get_fits_header(self) -> Optional[str]:
        """
        Retrieve the full FITS header of the current Siril image using shared memory.

        Args:
            none.

        Returns:
            bytes: The image FITS header as a string.

        Raises:
            NoImageError: If no image is currently loaded
            RuntimeError: For other errors during  data retrieval
            ValueError: If the received data format is invalid or shape is invalid
        """
        shm = None
        try:
            # Request shared memory setup
            status, response = self._send_command(_Command.GET_FITS_HEADER)

            # Handle error responses
            if status == _Status.ERROR:
                if response:
                    error_msg = response.decode('utf-8', errors='replace')
                    if "no image loaded" in error_msg.lower():
                        raise NoImageError("No image is currently loaded in Siril")
                    else:
                        raise RuntimeError(_("Server error: {}").format(error_msg))
                else:
                    raise RuntimeError(_("Failed to initiate shared memory transfer: Empty response"))

            if not response:
                raise RuntimeError(_("Failed to initiate shared memory transfer: No data received"))

            if len(response) < 25: # No payload
                return None

            try:
                # Parse the shared memory information
                shm_info = SharedMemoryInfo.from_buffer_copy(response)
            except (AttributeError, BufferError, ValueError) as e:
                raise ValueError(_("Invalid shared memory information received: {}").format(e))

            # Map the shared memory
            try:
                shm = self._map_shared_memory(
                    shm_info.shm_name.decode('utf-8'),
                    shm_info.size
                )
            except (OSError, ValueError) as e:
                raise RuntimeError(_("Failed to map shared memory: {}").format(e))

            try:
                # Read entire buffer at once using memoryview
                buffer = memoryview(shm.buf).cast('B')
                result = buffer.tobytes().decode('utf-8', errors='ignore')
            except (BufferError, ValueError, TypeError) as e:
                raise RuntimeError(_("Failed to create string from shared memory: {}").format(e))

            return result

        except NoImageError:
            # Re-raise NoImageError without wrapping
            raise
        except Exception as e:
            # Wrap all other exceptions with context
            raise RuntimeError(_("Error retrieving FITS header: {}").format(e)) from e
        finally:
            if shm is not None:
                try:
                    shm.close()
                    shm.unlink()
                except BufferError:
                    pass

    def get_unknown_keys(self) -> Optional[str]:
        """
        Retrieve unknown keys in the FITS header of the current Siril image using shared memory.

        Args:
            none.

        Returns:
            string: The unknown keys in the FITS header as a string.

        Raises:
            NoImageError: If no image is currently loaded
            RuntimeError: For other errors during  data retrieval
            ValueError: If the received data format is invalid or shape is invalid
        """
        shm = None
        try:
            # Request shared memory setup
            status, response = self._send_command(_Command.GET_FITS_UNKNOWN_KEYS)

            # Handle error responses
            if status == _Status.ERROR:
                if response:
                    error_msg = response.decode('utf-8', errors='replace')
                    if "no image loaded" in error_msg.lower():
                        raise NoImageError("No image is currently loaded in Siril")
                    else:
                        raise RuntimeError(_("Server error: {}").format(error_msg))
                else:
                    raise RuntimeError(_("Failed to initiate shared memory transfer: Empty response"))

            if status == _Status.NONE:
                return None

            if not response:
                raise RuntimeError(_("Failed to initiate shared memory transfer: No data received"))
            if len(response) < 25: # No payload
                return None

            try:
                # Parse the shared memory information
                shm_info = SharedMemoryInfo.from_buffer_copy(response)
            except (AttributeError, BufferError, ValueError) as e:
                raise ValueError(_("Invalid shared memory information received: {}").format(e))

            # Map the shared memory
            try:
                shm = self._map_shared_memory(
                    shm_info.shm_name.decode('utf-8'),
                    shm_info.size
                )
            except (OSError, ValueError) as e:
                raise RuntimeError(_("Failed to map shared memory: {}").format(e))

            try:
                # Read entire buffer at once using memoryview
                buffer = memoryview(shm.buf).cast('B')
                result = buffer.tobytes().decode('utf-8', errors='ignore')
            except (BufferError, ValueError, TypeError) as e:
                raise RuntimeError(_("Failed to create string from shared memory: {}").format(e))

            return result

        except NoImageError:
            # Re-raise NoImageError without wrapping
            raise
        except Exception as e:
            # Wrap all other exceptions with context
            raise RuntimeError(_("Error retrieving FITS unknown keys: {}").format(e)) from e
        finally:
            if shm is not None:
                try:
                    shm.close()
                    shm.unlink()
                except BufferError:
                    pass

    def get_history(self) -> Optional[list[str]]:
        """
        Retrieve history entries in the FITS header of the current Siril image using shared memory.
        Args:
            none.
        Returns:
            list: The history entries in the FITS header as a list of strings.
        Raises:
            NoImageError: If no image is currently loaded
            RuntimeError: For other errors during data retrieval
            ValueError: If the received data format is invalid or shape is invalid
        """
        shm = None
        try:
            # Request shared memory setup
            status, response = self._send_command(_Command.GET_FITS_HISTORY)

            # Handle error responses
            if status == _Status.ERROR:
                if response:
                    error_msg = response.decode('utf-8', errors='replace')
                    if "no image loaded" in error_msg.lower():
                        raise NoImageError("No image is currently loaded in Siril")
                    else:
                        raise RuntimeError(_("Server error: {}").format(error_msg))
                else:
                    raise RuntimeError("Failed to initiate shared memory transfer: Empty response")

            if status == _Status.NONE:
                return None

            if not response:
                raise RuntimeError("Failed to initiate shared memory transfer: No data received")

            if len(response) < 29:  # No payload
                return None

            try:
                # Parse the shared memory information
                shm_info = SharedMemoryInfo.from_buffer_copy(response)
            except (AttributeError, BufferError, ValueError) as e:
                raise ValueError(_("Invalid shared memory information received: {}").format(e))

            # Map the shared memory
            try:
                shm = self._map_shared_memory(
                    shm_info.shm_name.decode('utf-8'),
                    shm_info.size
                )
            except (OSError, ValueError) as e:
                raise RuntimeError(_("Failed to map shared memory: {}").format(e))

            try:
                # Read entire buffer at once using memoryview
                buffer = memoryview(shm.buf).cast('B')
                string_data = buffer.tobytes().decode('utf-8', errors='ignore')
                string_list = string_data.split('\x00')
            except (BufferError, ValueError, TypeError) as e:
                raise RuntimeError(_("Failed to create string from shared memory: {}").format(e))

            return [s for s in string_list if s]

        except NoImageError:
            # Re-raise NoImageError without wrapping
            raise
        except Exception as e:
            # Wrap all other exceptions with context
            raise RuntimeError(_("Error retrieving FITS history: {}").format(e)) from e
        finally:
            if shm is not None:
                try:
                    shm.close()
                    shm.unlink()
                except BufferError:
                    pass

    def get_wd(self) -> Optional[str]:
        """
        Request the working directory from Siril.

        Returns:
            The working directory as a string, or None if an error occurred.
        """
        response = self.request_data(_Command.GET_WORKING_DIRECTORY)

        if response is None:
            return None

        try:
            # Assuming the response is a null-terminated UTF-8 encoded string
            wd = response.decode('utf-8').rstrip('\x00')
            return wd
        except UnicodeDecodeError as e:
            print("Error decoding working directory: {e}", file=sys.stderr)
            return None

    def get_filename(self) -> Optional[str]:
        """
        Request the working directory from Siril.

        Returns:
            The working directory as a string, or None if an error occurred.
        """
        response = self.request_data(_Command.GET_FILENAME)

        if response is None:
            return None

        try:
            # Assuming the response is a null-terminated UTF-8 encoded string
            wd = response.decode('utf-8').rstrip('\x00')
            return wd
        except UnicodeDecodeError as e:
            print(f"Error decoding loaded image filename: {e}", file=sys.stderr)
            return None

    def get_image_stats(self, channel: int) -> Optional[ImageStats]:
        """
        Request image statistics from Siril for a specific channel.

        Args:
            channel: Integer specifying which channel to get statistics for (typically 0, 1, or 2)

        Returns:
            ImageStats object containing the statistics, or None if an error occurred
        """
        # Convert channel number to network byte order bytes
        channel_payload = struct.pack('!I', channel)  # '!I' for network byte order uint32_t

        # Request data with the channel number as payload
        response = self.request_data(_Command.GET_IMAGE_STATS, payload=channel_payload)

        if response is None:
            return None

        try:
            # Define the format string for unpacking the C struct
            # '!' for network byte order (big-endian)
            # 'q' for long (total, ngoodpix) - using 64-bit integers to match gint64
            # 'd' for double (all floating point values)
            # We don't include the gint *nbrefs as it's not needed in Python
            format_string = '!2q12d'  # '!' ensures network byte order

            # Calculate expected size
            expected_size = struct.calcsize(format_string)

            # Verify we got the expected amount of data
            if len(response) != expected_size:
                print(f"Received stats data size {len(response)} doesn't match expected size {expected_size}",
                    file=sys.stderr)
                return None

            # Unpack the binary data
            values = struct.unpack(format_string, response)

            # Create and return an ImageStats object with the unpacked values
            return ImageStats(
                total=values[0],
                ngoodpix=values[1],
                mean=values[2],
                median=values[3],
                sigma=values[4],
                avgDev=values[5],
                mad=values[6],
                sqrtbwmv=values[7],
                location=values[8],
                scale=values[9],
                min=values[10],
                max=values[11],
                normValue=values[12],
                bgnoise=values[13]
            )

        except struct.error as e:
            print(f"Error unpacking image statistics data: {e}", file=sys.stderr)
            return None

    def get_seq_regdata(self, frame: int, channel: int) -> Optional[RegData]:
        """
        Request sequence frame registration data from Siril.

        Returns:
            RegData object containing the registration data, or None if an error occurred
        """
        data_payload = struct.pack('!II', frame, channel)  # '!I' for network byte order uint32_t

        # Request data with the channel number as payload
        response = self.request_data(_Command.GET_SEQ_REGDATA, payload=data_payload)

        if response is None:
            return None

        try:
            format_string = '!5dQ9d2Q'

            values = struct.unpack(format_string, response)

            return RegData (
                fwhm = values[0],
                weighted_fwhm = values[1],
                roundness = values[2],
                quality = values[3],
                background_lvl = values[4],
                number_of_stars = values[5],
                H = Homography (
                    h00=values[6],
                    h01=values[7],
                    h02=values[8],
                    h10=values[9],
                    h11=values[10],
                    h12=values[11],
                    h20=values[12],
                    h21=values[13],
                    h22=values[14],
                    pair_matched=values[15],
                    Inliers=values[16]
                )
            )
        except struct.error as e:
            print(f"Error unpacking frame registration data: {e}", file=sys.stderr)
            return None

    def get_seq_imstats(self, frame: int, channel: int) -> Optional[RegData]:
        """
        Request sequence frame registration data from Siril.

        Returns:
            RegData object containing the registration data, or None if an error occurred
        """
        data_payload = struct.pack('!II', frame, channel)  # '!I' for network byte order uint32_t

        # Request data with the channel number as payload
        response = self.request_data(_Command.GET_SEQ_STATS, payload=data_payload)

        if response is None:
            return None
        try:
            format_string = '!2q12d'

            values = struct.unpack(format_string, response)

            return ImageStats (
                total = values[0],
                ngoodpix = values[1],
                mean = values[2],
                median = values[3],
                sigma = values[4],
                avgDev = values[5],
                mad = values[6],
                sqrtbwmv = values[7],
                location = values[8],
                scale = values[9],
                min = values[10],
                max = values[11],
                normValue = values[12],
                bgnoise = values[13]
            )
        except struct.error as e:
            print(f"Error unpacking frame statistics: {e}", file=sys.stderr)
            return None

    def get_seq_imgdata(self, frame: int) -> Optional[ImgData]:
        """
        Request sequence frame registration data from Siril.

        Returns:
            ImgData object containing the image data, or None if an error occurred
        """
        data_payload = struct.pack('!I', frame)  # '!I' for network byte order uint32_t

        # Request data with the channel number as payload
        response = self.request_data(_Command.GET_SEQ_IMGDATA, payload=data_payload)
        if response is None:
            return None
        try:
            format_string = '!3qd2q'

            values = struct.unpack(format_string, response)

            return ImgData (
                filenum = values[0],
                incl = values[1],
                date_obs = datetime.fromtimestamp(values[2]) if values[2] != 0 else None,
                airmass = values[3],
                rx = values[4],
                ry = values[5]
            )
        except struct.error as e:
            print(f"Error unpacking frame image data: {e}", file=sys.stderr)
            return None

    def get_seq(self) -> Optional[ImgData]:
        """
        Request metadata for the current sequence loaded in Siril.

        Returns:
            Sequence object containing the current sequence metadata, or None
            if an error occurred
        """
        # Request data with the channel number as payload
        response = self.request_data(_Command.GET_SEQ)

        if response is None:
            return None

        try:
            format_string = '!4q3Q4qdQqQq'
            fixed_length = struct.calcsize(format_string)

            values = struct.unpack(format_string, response[:fixed_length])
            # Extract remaining bytes for the null-terminated string
            remaining_data = response[fixed_length:]
            seqname_string = remaining_data.decode('utf-8')

            number = values[0]
            nb_layers = values[3]

            imgparam_list = [self.get_seq_imgdata(frame)
                             for frame in range(number)]

            regdata_list = [[self.get_seq_regdata(frame, channel)
                             for channel in range(nb_layers)]
                            for frame in range(number)]

            stats_list = [[self.get_seq_imstats(frame, channel)
                           for channel in range(nb_layers)]
                          for frame in range(number)]

            return Sequence (
                number = values[0],
                selnum = values[1],
                fixed = values[2],
                nb_layers = values[3],
                rx = values[4],
                ry = values[5],
                is_variable = bool(values[6]),
                bitpix = values[7],
                reference_image = values[8],
                imgparam = imgparam_list,
                regparam = regdata_list,
                stats = stats_list,
                beg = values[9],
                end = values[10],
                exposure = values[11],
                fz = bool(values[12]),
                type = SequenceType(values[13]),
                cfa_opened_monochrome = bool(values[14]),
                current = values[15],
                seqname = seqname_string
            )
        except struct.error as e:
            print(f"Error unpacking sequence data: {e}", file=sys.stderr)
            return None

    def get_keywords(self) -> Optional[FKeywords]:
        """
        Request FITS keywords data from Siril.

        Returns:
            FKeywords object containing the FITS keywords, or None if an error occurred
        """
        response = self.request_data(_Command.GET_KEYWORDS)
        if response is None:
            return None

        try:
            # Constants matching C implementation
            FLEN_VALUE = 71  # Standard FITS keyword length

            # Build format string for struct unpacking
            # Network byte order for all values
            format_parts = [
                f'{FLEN_VALUE}s',  # program
                f'{FLEN_VALUE}s',  # filename
                f'{FLEN_VALUE}s',  # row_order
                f'{FLEN_VALUE}s',  # filter
                f'{FLEN_VALUE}s',  # image_type
                f'{FLEN_VALUE}s',  # object
                f'{FLEN_VALUE}s',  # instrume
                f'{FLEN_VALUE}s',  # telescop
                f'{FLEN_VALUE}s',  # observer
                f'{FLEN_VALUE}s',  # sitelat_str
                f'{FLEN_VALUE}s',  # sitelong_str
                f'{FLEN_VALUE}s',  # bayer_pattern
                f'{FLEN_VALUE}s',  # focname
                'd',  # bscale
                'd',  # bzero
                'Q',  # lo padded to 64bit
                'Q',  # hi padded to 64bit
                'd',  # flo padded to 64bit
                'd',  # fhi padded to 64bit
                'd',  # data_max
                'd',  # data_min
                'd',  # pixel_size_x
                'd',  # pixel_size_y
                'Q',  # binning_x (padded to uint64_t)
                'Q',  # binning_y (padded to uint64_t)
                'd',  # expstart
                'd',  # expend
                'd',  # centalt
                'd',  # centaz
                'd',  # sitelat
                'd',  # sitelong
                'd',  # siteelev
                'q',  # bayer_xoffset
                'q',  # bayer_yoffset
                'd',  # airmass
                'd',  # focal_length
                'd',  # flength
                'd',  # iso_speed
                'd',  # exposure
                'd',  # aperture
                'd',  # ccd_temp
                'd',  # set_temp
                'd',  # livetime
                'Q',  # stackcnt
                'd',  # cvf
                'q',  # key_gain
                'q',  # key_offset
                'q',  # focuspos
                'q',  # focussz
                'd',  # foctemp
                'q',  # date (int64 unix timestamp)
                'q'  # date_obs (int64 unix timestamp)
            ]

            format_string = '!' + ''.join(format_parts)

            # Verify data size
            expected_size = struct.calcsize(format_string)
            if len(response) != expected_size:
                print(f"Received keyword data size {len(response)} doesn't match expected size {expected_size}",
                    file=sys.stderr)
                return None

            # Unpack the binary data
            values = struct.unpack(format_string, response)

            # Helper function to decode and strip null-terminated strings
            def decode_string(s: bytes) -> str:
                return s.decode('utf-8').rstrip('\x00')

            # Helper function to convert timestamp to datetime
            def timestamp_to_datetime(timestamp: int) -> Optional[datetime]:
                return datetime.fromtimestamp(timestamp) if timestamp != 0 else None

            # Create FKeywords object
            return FKeywords(
                program=decode_string(values[0]),
                filename=decode_string(values[1]),
                row_order=decode_string(values[2]),
                filter=decode_string(values[3]),
                image_type=decode_string(values[4]),
                object=decode_string(values[5]),
                instrume=decode_string(values[6]),
                telescop=decode_string(values[7]),
                observer=decode_string(values[8]),
                sitelat_str=decode_string(values[9]),
                sitelong_str=decode_string(values[10]),
                bayer_pattern=decode_string(values[11]),
                focname=decode_string(values[12]),
                bscale=values[13],
                bzero=values[14],
                lo=values[15],
                hi=values[16],
                flo=values[17],
                fhi=values[18],
                data_max=values[19],
                data_min=values[20],
                pixel_size_x=values[21],
                pixel_size_y=values[22],
                binning_x=values[23],
                binning_y=values[24],
                expstart=values[25],
                expend=values[26],
                centalt=values[27],
                centaz=values[28],
                sitelat=values[29],
                sitelong=values[30],
                siteelev=values[31],
                bayer_xoffset=values[32],
                bayer_yoffset=values[33],
                airmass=values[34],
                focal_length=values[35],
                flength=values[36],
                iso_speed=values[37],
                exposure=values[38],
                aperture=values[39],
                ccd_temp=values[40],
                set_temp=values[41],
                livetime=values[42],
                stackcnt=values[43],
                cvf=values[44],
                key_gain=values[45],
                key_offset=values[46],
                focuspos=values[47],
                focussz=values[48],
                foctemp=values[49],
                date=timestamp_to_datetime(values[50]),
                date_obs=timestamp_to_datetime(values[51])
            )

        except struct.error as e:
            print(f"Error unpacking FITS keywords data: {e}", file=sys.stderr)
            return None
        except Exception as e:
            print(f"Error processing FITS keywords data: {e}", file=sys.stderr)
            return None

    def get_image(self, get_pixels: Optional[bool] = True) -> Optional[FFit]:
        """
        Request a copy of the current image open in Siril.

        Args:
            none

        Returns:
            FFit object containing the image pixel data, statistics and metadata, or None
            if an error occurred
        """

        # Request data with the channel number as payload
        response = self.request_data(_Command.GET_IMAGE)
        if response is None:
            return None

        try:
            shape = self.get_shape()
            # Build format string for struct unpacking
            # Network byte order for all values
            format_parts = [
                'q',  # bitpix padded to 64bit
                'q',  # orig_bitpix padded to 64bit
                'Q',  # gboolean checksum padded to 64bit
                'd',  # mini
                'd',  # maxi
                'd',  # neg_ratio padded to 64bit
                'Q',  # data_type (padded to uint64_t)
                'Q',  # gboolean top_down (padded to uint64_t)
                'Q',  # gboolean focalkey (padded to uint64_t)
                'Q',  # gboolean pixelkey (padded to uint64_t)
                'Q',  # gboolean color_managed (padded to uint64_t)
            ]

            format_string = '!' + ''.join(format_parts)

            # Verify data size
            expected_size = struct.calcsize(format_string)
            if len(response) != expected_size:
                print(f"Received image data size {len(response)} doesn't match expected size {expected_size}",
                    file=sys.stderr)
                return None

            # Unpack the binary data
            values = struct.unpack(format_string, response)

            # Helper function to decode and strip null-terminated strings
            def decode_string(s: bytes) -> str:
                return s.decode('utf-8').rstrip('\x00')

            # Create FFit object
            try:
                img_header = self.get_fits_header()
            except Exception as e:
                img_header = ""

            try:
                img_unknown_keys = self.get_unknown_keys()
            except Exception as e:
                img_unknown_keys = ""

            try:
                img_history = self.get_history()
            except Exception as e:
                img_history = []

            return FFit(
                _naxes = (shape[1], shape[0], shape[2]),
                naxis = 2 if shape[2] == 1 else 3,
                bitpix=values[0],
                checksum=True if values[1] else False,
                mini=values[2],
                maxi=values[3],
                neg_ratio=values[4],
                top_down=True if values[6] else False,
                focalkey=True if values[7] else False,
                pixelkey=True if values[8] else False,
                color_managed=True if values[9] else False,
                _data = self.get_pixel_data() if get_pixels == True else None,
                stats = (
                    self.get_image_stats(0),
                    self.get_image_stats(1) if shape[2] > 1 else None,
                    self.get_image_stats(2) if shape[2] > 1 else None
                ),
                keywords = self.get_keywords(),
                _icc_profile = self.get_icc_profile(),
                header = img_header,
                unknown_keys = img_unknown_keys,
                history = img_history
            )

        except struct.error as e:
            print(f"Error unpacking FITS metadata: {e}", file=sys.stderr)
            return None
        except Exception as e:
            print(f"Error processing FITS metadata: {e}", file=sys.stderr)
            return None

    def get_stars(self) -> List[PSFStar]:
        """
        Request PSF star data from Siril, iterating through all available stars.
        Returns:
            List of PSFStar objects containing the star data.
        """
        stars = []
        shm = None

        try:
            # Request shared memory setup
            status, response = self._send_command(_Command.GET_PSFSTARS)

            # Handle error responses
            if status == _Status.ERROR:
                if response:
                    error_msg = response.decode('utf-8', errors='replace')
                    if "no image loaded" in error_msg.lower():
                        raise NoImageError("No image is currently loaded in Siril")
                    else:
                        raise RuntimeError(_("Server error: {}").format(error_msg))
                else:
                    raise RuntimeError("Failed to initiate shared memory transfer: Empty response")

            if status == _Status.NONE:
                return None

            if not response:
                raise RuntimeError(_("Failed to initiate shared memory transfer: No data received"))

            try:
                # Parse the shared memory information
                shm_info = SharedMemoryInfo.from_buffer_copy(response)
            except (AttributeError, BufferError, ValueError) as e:
                raise ValueError(_("Invalid shared memory information received: {}").format(e))

            # Map the shared memory
            try:
                shm = self._map_shared_memory(
                    shm_info.shm_name.decode('utf-8'),
                    shm_info.size
                )
            except (OSError, ValueError) as e:
                raise RuntimeError(_("Failed to map shared memory: {}").format(e))

            format_string = '!13d2qdq16dqdd'  # Define the format string based on PSFStar structure
            fixed_size = struct.calcsize(format_string)

            # Read entire buffer at once using memoryview
            buffer = memoryview(shm.buf).cast('B')

            # Validate buffer size
            if len(buffer) % fixed_size != 0:
                raise ValueError(_("Buffer size {} is not a multiple "
                                   "of struct size {}").format(len(buffer), fixed_size))

            num_stars = len(buffer) // fixed_size

            # Sanity check for number of stars
            if num_stars <= 0:  # adjust max limit as needed
                raise ValueError(_("Invalid number of stars: {}").format(num_stars))

            if num_stars > 200000: # to match the #define MAX_STARS
                num_stars = 200000
                self.log(_("Limiting stars to max 200000"))

            for i in range(num_stars):
                # Calculate start and end positions for each struct
                start = i * fixed_size
                end = start + fixed_size

                try:
                    # Extract the bytes for this struct and unpack
                    values = struct.unpack(format_string, buffer[start:end].tobytes())

                    star = PSFStar(
                        B=values[0], A=values[1], x0=values[2], y0=values[3],
                        sx=values[4], sy=values[5], fwhmx=values[6], fwhmy=values[7],
                        fwhmx_arcsec=values[8], fwhmy_arcsec=values[9], angle=values[10],
                        rmse=values[11], sat=values[12], R=values[13],
                        has_saturated=bool(values[14]), beta=values[15],
                        profile=values[16], xpos=values[17], ypos=values[18],
                        mag=values[19], Bmag=values[20], s_mag=values[21],
                        s_Bmag=values[22], SNR=values[23], BV=values[24],
                        B_err=values[25], A_err=values[26], x_err=values[27],
                        y_err=values[28], sx_err=values[29], sy_err=values[30],
                        ang_err=values[31], beta_err=values[32], layer=values[33],
                        ra=values[34], dec=values[35]
                    )
                    stars.append(star)

                except struct.error as e:
                    print(f"Error unpacking star data for index {i}: {e}", file=sys.stderr)
                    break

            return stars

        except Exception as e:
            raise RuntimeError(_("Error processing star data: {}").format(e))

        finally:
            if shm is not None:
                try:
                    shm.close()
                    shm.unlink()
                except Exception as e:
                    print(f"Error closing shared memory: {e}", file=sys.stderr)

    def get_config(self, group: str, key: str) -> Optional[Union[bool, int, float, str, List[str]]]:
        """
        Request a configuration value from Siril.

        Args:
            group: Configuration group name
            key: Configuration key name within the group

        Returns:
            The configuration value with appropriate Python type, or None if an error occurred
        """
        try:
            # Encode both group and key with null terminator
            group_bytes = group.encode('utf-8') + b'\0'
            key_bytes = key.encode('utf-8') + b'\0'
            payload = group_bytes + key_bytes

            # Request the config value
            response = self.request_data(_Command.GET_CONFIG, payload=payload)
            if response is None:
                return None

            # First byte should be the type
            config_type = ConfigType(response[0])
            value_data = response[1:]

            # Parse based on type
            if config_type == ConfigType.BOOL:
                return bool(struct.unpack('!I', value_data)[0])
            elif config_type == ConfigType.INT:
                return struct.unpack('!i', value_data)[0]
            elif config_type == ConfigType.DOUBLE:
                return struct.unpack('!d', value_data)[0]
            elif config_type in (ConfigType.STR, ConfigType.STRDIR):
                # Assume null-terminated string
                string_value = value_data.split(b'\0')[0].decode('utf-8')
                return string_value
            elif config_type == ConfigType.STRLIST:
                # Split on null bytes and decode each string
                strings = value_data.split(b'\0')
                # Remove empty strings at the end
                while strings and not strings[-1]:
                    strings.pop()
                return [s.decode('utf-8') for s in strings]

        except Exception as e:
            raise RuntimeError(_("Error getting config value: {}").format(e))
            return None