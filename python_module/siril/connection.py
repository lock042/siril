import os
import sys
import struct
import socket
import ctypes
import time
import mmap
from pathlib import Path
from enum import IntEnum
from typing import Tuple, Optional
import numpy as np
from .shm import SharedMemoryWrapper
from .exceptions import SirilError, ConnectionError, CommandError, DataError, NoImageError
if os.name == 'nt':
    import win32pipe
    import win32file
    import win32event
    import pywintypes
    import winerror

class _Command(IntEnum):
    GET_DIMENSIONS = 1
    GET_PIXELDATA = 2
    GET_PIXELDATA_REGION = 3
    NOTIFY_PROGRESS = 4
    SEND_COMMAND = 5
    LOG_MESSAGE = 6
    GET_WORKING_DIRECTORY = 7
    GET_FILENAME = 8
    SET_PIXELDATA = 9
    ERROR = 0xFF

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
        """
        Initialize the SirilInterface, automatically determining the correct pipe or socket path
        based on the environment variable and operating system.
        """
        if os.name == 'nt':  # Windows
            self.pipe_path = os.getenv('MY_PIPE')
            if not self.pipe_path:
                raise ConnectionError("Environment variable MY_PIPE not set")
            self.event_pipe_path = self.pipe_path  # Assuming event pipe is the same path
        else:  # POSIX (Linux/macOS)
            self.socket_path = os.getenv('MY_SOCKET')
            if not self.socket_path:
                raise ConnectionError("Environment variable MY_SOCKET not set")
            self.event_pipe_path = self.socket_path  # Assuming event socket is the same path

    def connect(self):
        """
        Establish a connection to Siril based on the pipe or socket path.
        Returns True if the connection is successful, otherwise False.
        """
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
                        raise ConnectionError("Pipe is busy")
                    raise ConnectionError(f"Failed to connect to pipe: {e}")
            else:
                self.sock = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
                self.sock.connect(self.socket_path)
                return True

        except Exception as e:
            raise ConnectionError(f"Failed to connect: {e}")

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
                raise SirilError("No pipe connection to close")
        else:
            if hasattr(self, 'sock'):
                self.sock.close()
                return True
            else:
                raise SirilError("No socket connection to close")

    def _recv_exact(self, n: int, timeout: float = 3.0) -> Optional[bytes]:
        """
        Helper method to receive exactly n bytes from the socket or pipe.
        """
        if n < 0:
            raise ValueError("Cannot receive negative number of bytes")

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
                        raise ConnectionError("Timeout while receiving data")

                    if rc != win32event.WAIT_OBJECT_0:
                        raise ConnectionError("Error waiting for pipe read completion")

                    # Get results of the operation
                    err, bytes_read = win32file.GetOverlappedResult(self.pipe_handle, self.overlap_read, False)
                    if bytes_read == 0:
                        raise ConnectionError("Pipe closed during read")

                    # Extend our data buffer
                    data.extend(buf[:bytes_read])

                return bytes(data)

            except pywintypes.error as e:
                raise ConnectionError(f"Windows pipe error during receive: {e}")

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
                            raise ConnectionError("Connection closed during data transfer")
                        data.extend(packet)
                    except socket.timeout:
                        raise ConnectionError("Timeout while receiving data")
                    except Exception as e:
                        raise ConnectionError(f"Error receiving data: {e}")
                return bytes(data)
            finally:
                self.sock.settimeout(original_timeout)

    def _map_shared_memory(self, name: str, size: int) -> mmap.mmap:
        """Create or open a shared memory mapping"""
        if os.name == 'nt':
            try:
                # Create a file mapping
                mapping_handle = win32file.CreateFileMapping(
                    win32file.INVALID_HANDLE_VALUE,  # Use paging file
                    None,  # Default security
                    win32file.PAGE_READWRITE,  # Read/write access
                    0,  # Maximum size (high-order DWORD)
                    size,  # Maximum size (low-order DWORD)
                    name  # Mapping name
                )

                try:
                    # Create a memory map view
                    map_view = win32file.MapViewOfFile(
                        mapping_handle,
                        win32file.FILE_MAP_READ | win32file.FILE_MAP_WRITE,
                        0,  # Offset high
                        0,  # Offset low
                        size
                    )

                    # Create mmap object from the handle
                    return mmap.mmap(-1, size, name)

                finally:
                    win32file.CloseHandle(mapping_handle)

            except pywintypes.error as e:
                raise RuntimeError(f"Failed to create shared memory mapping: {e}")
        else:
            fd = os.open(f"/dev/shm/{name}", os.O_RDONLY)
            try:
                return mmap.mmap(fd, size, mmap.MAP_SHARED, mmap.PROT_READ)
            finally:
                os.close(fd)

    def _send_command(self, command: _Command, data: Optional[bytes] = None) -> Tuple[Optional[int], Optional[bytes]]:
        """
        Send a command to Siril and receive the response.
        """
        try:
            data_length = len(data) if data else 0
            header = struct.pack('!BI', command, data_length)

            if os.name == 'nt':
                try:
                    # Send header
                    err, bytes_written = win32file.WriteFile(self.pipe_handle, header, self.overlap_write)
                    rc = win32event.WaitForSingleObject(self.overlap_write.hEvent, 3000)  # 3 second timeout
                    if rc != win32event.WAIT_OBJECT_0:
                        raise ConnectionError("Timeout or error while sending header")

                    # Send data if present
                    if data and data_length > 0:
                        err, bytes_written = win32file.WriteFile(self.pipe_handle, data, self.overlap_write)
                        rc = win32event.WaitForSingleObject(self.overlap_write.hEvent, 3000)
                        if rc != win32event.WAIT_OBJECT_0:
                            raise ConnectionError("Timeout or error while sending data")

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
                    raise CommandError(f"Windows pipe error during command: {e}")
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
            raise CommandError(f"Error sending command: {e}")
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

        if status == _Command.ERROR:
            error_msg = response.decode('utf-8') if response else "Unknown error"
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

        if status is None or status == _Command.ERROR:
            error_msg = response.decode('utf-8') if response else "Unknown error"
            print(f"Data request failed: {error_msg}", file=sys.stderr)
            return None

        return response

    # Specific commands follow below here

    def log(self, my_string: str) -> bool:
        """
        Send a log message to Siril.

        Args:
            my_string: The message to log

        Returns:
            bool: True if the message was successfully logged, False otherwise
        """
        try:
            # Append a newline character to the string
            my_string += '\n'  # Add newline character
            # Convert string to bytes using UTF-8 encoding
            message_bytes = my_string.encode('utf-8')
            return self.execute_command(_Command.LOG_MESSAGE, message_bytes)

        except Exception as e:
            print(f"Error sending log message: {e}", file=sys.stderr)
            return False

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
                        raise RuntimeError(f"Server error: {error_msg}")
                else:
                    raise RuntimeError("Failed to initiate shared memory transfer: Empty response")

            if not response:
                raise RuntimeError("Failed to initiate shared memory transfer: No data received")

            try:
                # Parse the shared memory information
                shm_info = SharedMemoryInfo.from_buffer_copy(response)
            except (AttributeError, BufferError, ValueError) as e:
                raise ValueError(f"Invalid shared memory information received: {e}")

            # Validate dimensions
            if any(dim <= 0 for dim in (shm_info.width, shm_info.height)):
                raise ValueError(f"Invalid image dimensions: {shm_info.width}x{shm_info.height}")

            if shm_info.channels <= 0 or shm_info.channels > 3:
                raise ValueError(f"Invalid number of channels: {shm_info.channels}")

            # Map the shared memory
            try:
                shm = self._map_shared_memory(
                    shm_info.shm_name.decode('utf-8'),
                    shm_info.size
                )
            except (OSError, ValueError) as e:
                raise RuntimeError(f"Failed to map shared memory: {e}")

            # Create numpy array from shared memory
            dtype = np.float32 if shm_info.data_type == 1 else np.uint16
            try:
                buffer = memoryview(shm).cast('B')
                arr = np.frombuffer(buffer, dtype=dtype)
            except (BufferError, ValueError, TypeError) as e:
                raise RuntimeError(f"Failed to create array from shared memory: {e}")

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
                raise ValueError(f"Failed to reshape array to image dimensions: {e}")

            # Make a copy of the data since we'll be releasing the shared memory
            result = np.copy(arr)

            # Notify C side to clean up shared memory
            try:
                # Ensure that shm_info.shm_name is exactly 256 bytes (trimmed or padded with null bytes)
                shm_name = shm_info.shm_name[:256].ljust(256, b'\x00')

                # Pack the header and payload
                cmd_header = struct.pack('!BL', _Command.NOTIFY_PROGRESS, 256)  # 'B' for the command and 'L' for length
                self.sock.sendall(cmd_header + shm_name)
            except socket.error as e:
                raise RuntimeError(f"Failed to send completion notification: {e}")

            return result

        except NoImageError:
            # Re-raise NoImageError without wrapping
            raise
        except Exception as e:
            # Wrap all other exceptions with context
            raise RuntimeError(f"Error retrieving pixel data: {e}") from e
        finally:
            if shm is not None:
                try:
                    shm.close()
                except BufferError:
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
                raise ValueError("Image data must be a numpy array")

            if image_data.ndim not in (2, 3):
                raise ValueError("Image must be 2D or 3D array")

            if image_data.dtype not in (np.float32, np.uint16):
                raise ValueError("Image data must be float32 or uint16")

            # Get dimensions
            if image_data.ndim == 2:
                height, width = image_data.shape
                channels = 1
                image_data = image_data.reshape(height, width, 1)
            else:
                height, width, channels = image_data.shape

            if channels > 3:
                raise ValueError("Image cannot have more than 3 channels")

            if any(dim <= 0 for dim in (width, height)):
                raise ValueError(f"Invalid image dimensions: {width}x{height}")

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
                raise RuntimeError(f"Failed to create shared memory: {e}")

            # Copy data to shared memory
            try:
                buffer = memoryview(shm.buf).cast('B')
                shared_array = np.frombuffer(buffer, dtype=image_data.dtype).reshape(image_data.shape)
                np.copyto(shared_array, image_data)
            except Exception as e:
                raise RuntimeError(f"Failed to copy data to shared memory: {e}")

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
            print(f"Error sending pixel data: {e}", file=sys.stderr)
            return False

        finally:
            if shm is not None:
                try:
                    shm.close()
                    shm.unlink()
                except:
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
            print(f"Error decoding working directory: {e}", file=sys.stderr)
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
