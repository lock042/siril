# Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
# Copyright (C) 2012-2026 team free-astro (see more in AUTHORS file)
# Reference site is https://siril.org
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Connection module for Siril, providing the ability to connect to a running
Siril instance and communicate with it. Includes an extensive range
of methods that can be used to get and set data from / to Siril.
"""

import os
import sys
import re
import atexit
import socket
import struct
import threading
from datetime import datetime
from typing import Tuple, Optional, List, Union
import numpy as np
from .translations import _
from .shm import SharedMemoryWrapper, _SharedMemoryInfo
from .plot import PlotData
from .exceptions import SirilError, DataError, SirilConnectionError, CommandError, \
        NoImageError, NoSequenceError, SharedMemoryError, ProcessingThreadBusyError, \
        ImageDialogOpenError, MouseModeError
from .models import ImageStats, FKeywords, FFit, PSFStar, BGSample, RegData, ImgData, \
        DistoData, Sequence, SequenceType, Polygon, ImageAnalysis
from .enums import _Command, _Status, CommandStatus, _ConfigType, LogColor, SirilVport, \
        STFType, SlidersMode
from .utility import truncate_utf8, parse_fits_header

DEFAULT_TIMEOUT = 10.

if os.name == 'nt':
    import win32file
    import win32event
    import pywintypes
    import winerror

class SirilInterface:
    """
    SirilInterface is the main class providing an interface to a running
    Siril instance and access to methods to interact with it through
    Siril's inbuilt command system and accessing image and sequence data.
    """

    _connected = False

    def __init__(self):
        """
        Initialize the SirilInterface, automatically determining the
        correct pipe or socket path based on the environment variable and
        operating system. Internal method.
        """
        self.connected = False

        if os.name == 'nt':
            self.pipe_path = os.getenv('MY_PIPE')
            if not self.pipe_path:
                raise SirilConnectionError(_("Environment variable MY_PIPE not set"))
            self.event_pipe_path = self.pipe_path

            # Windows-specific attributes
            self.pipe_handle = None
            self.overlap_read = None
            self.overlap_write = None
        else:
            self.socket_path = os.getenv('MY_SOCKET')
            if not self.socket_path:
                raise SirilConnectionError(_("Environment variable MY_SOCKET not set"))
            self.event_pipe_path = self.socket_path

            # Unix socket attribute
            self.sock = None

        # Add synchronization lock
        if os.name == 'nt':
            self.command_lock = win32event.CreateMutex(None, False, None)
        else:
            self.command_lock = threading.Lock()

        self.debug = bool(os.getenv('SIRIL_PYTHON_DEBUG') is not None)
        self._is_cli = bool(os.getenv('SIRIL_PYTHON_CLI') is not None)

    def connect(self) -> bool:
        """
        Establish a connection to Siril based on the pipe or socket path.

        Returns:
            True on success

        Raises:
            SirilConnectionError: if a connection error occurred
        """

        if self.connected is True:
            print(_("Already connected"))
            return True

        if SirilInterface._connected:
            raise SirilConnectionError(_("Error: a SirilInterface is already connected to Siril"))

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
                    if self.debug:
                        current_pid = os.getpid()
                        print(f'Current ProcessID is {current_pid}')
                        self.info_messagebox(f'Current ProcessID is {current_pid}', True)
                    SirilInterface._connected = True
                    self.connected = True
                    atexit.register(self._cleanup)
                    return True
                except pywintypes.error as e:
                    if e.winerror == winerror.ERROR_PIPE_BUSY:
                        raise SirilConnectionError(_("Pipe is busy")) from e
                    raise SirilConnectionError(_("Failed to connect to pipe: {}").format(e)) from e
            else:
                self.sock = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
                self.sock.connect(self.socket_path)
                if self.debug:
                    current_pid = os.getpid()
                    print(f'Current ProcessID is {current_pid}')
                    self.info_messagebox(f'Current ProcessID is {current_pid}', False)
                SirilInterface._connected = True
                self.connected = True
                atexit.register(self._cleanup)
                return True

        except Exception as e:
            raise SirilConnectionError(_("Failed to connect: {}").format(e)) from e

    def disconnect(self):
        """
        Closes the established socket or pipe connection. Note there is
        not usually any need to close this unless for some reason you wish
        to close a connection and subsequently reopen another one. This
        method is automatically called at script termination using an ``atexit``
        handler so there is no need to do so manually.
        Calling this method will reset the progress bar.

        Raises:
            SirilConnectionError: if the connection cannot be closed.
        """

        # Attempt to reset the progress bar in case the script exits while
        # disconnected and the progress bar has been left in a bad state.
        if self.connected is False:
            print("Not connected")
            return

        try:
            self.reset_progress()
        except Exception:
            print("Warning: unable to reset progress bar in disconnect()", file=sys.stderr)
            pass

        atexit.unregister(self._cleanup)
        if os.name == 'nt':
            if hasattr(self, 'pipe_handle'):
                # Close the pipe handle
                win32file.CloseHandle(self.pipe_handle)
                # Close the event handles
                if hasattr(self, 'overlap_read'):
                    win32file.CloseHandle(self.overlap_read.hEvent)
                if hasattr(self, 'overlap_write'):
                    win32file.CloseHandle(self.overlap_write.hEvent)
                self.connected = False
                SirilInterface._connected = False
                return
            raise SirilError(_("No pipe connection to close"))
        if hasattr(self, 'sock'):
            self.sock.close()
            self.connected = False
            SirilInterface._connected = False
            return
        raise SirilConnectionError(_("No socket connection to close"))

    def _cleanup(self):
        """
        Cleanup method that will run on exit. This is registered once,
        automatically, when the first SirilInterface class is created.
        If the connection persists it attempts to reset the progress bar;
        any other "at exit" methods may be added here as required.
        """
        try:
            self._release_thread()
            self.disconnect()
        except Exception:
            print("Warning: failed to clean up python module state")

    def _recv_exact(self, n: int, timeout: Optional[float] = DEFAULT_TIMEOUT) -> Optional[bytes]:
        """
        Helper method to receive exactly n bytes from the socket or pipe.
        Private method, not for direct use in scripts.

        Args:
            n: Number of bytes to receive
            timeout: Timeout in seconds. None for indefinite timeout

        Returns: bytes

        Raises:
            ValueError: if negative number of bytes specified
            SirilConnectionError: if there is an error related to the
                                  pipe / socket connection
        """
        if n < 0:
            raise ValueError(_("Cannot receive negative number of bytes"))

        if self.debug:
            timeout = None

        if os.name == 'nt':
            # Pipe implementation
            try:
                data = bytearray()
                # Convert None timeout to effectively infinite wait
                timeout_ms = int(timeout * 1000) if timeout is not None else win32event.INFINITE

                while len(data) < n:
                    # Calculate remaining bytes to read
                    to_read = n - len(data)

                    # Prepare buffer for reading
                    buf = win32file.AllocateReadBuffer(to_read)

                    # Start overlapped read
                    win32file.ReadFile(self.pipe_handle, buf, self.overlap_read)

                    # Wait for completion or timeout
                    rc = win32event.WaitForSingleObject(self.overlap_read.hEvent, timeout_ms)

                    if timeout is not None and rc == win32event.WAIT_TIMEOUT:
                        # Cancel the I/O operation
                        win32file.CancelIo(self.pipe_handle)
                        raise SirilConnectionError(_("Timeout while receiving data"))

                    if rc != win32event.WAIT_OBJECT_0:
                        raise SirilConnectionError(_("Error waiting for pipe read completion"))

                    # Get results of the operation
                    bytes_read = win32file.GetOverlappedResult(self.pipe_handle, self.overlap_read, False)
                    if bytes_read == 0:
                        raise SirilConnectionError(_("Pipe closed during read"))

                    # Extend our data buffer
                    data.extend(buf[:bytes_read])

                return bytes(data)

            except pywintypes.error as e:
                raise SirilConnectionError(_("Windows pipe error during receive: {}").format(e)) from e
            except Exception as ee:
                raise SirilConnectionError(_("Connection error: {}").format(ee)) from ee

        else:
            # Socket implementation
            original_timeout = self.sock.gettimeout()
            # Set to None if timeout is None, otherwise to the specified timeout
            self.sock.settimeout(timeout)

            try:
                data = bytearray()
                while len(data) < n:
                    try:
                        packet = self.sock.recv(n - len(data))
                        if not packet:
                            raise SirilConnectionError(_("Connection closed during data transfer"))
                        data.extend(packet)
                    except socket.timeout as exc:
                        raise SirilConnectionError(_("Timeout while receiving data")) from exc
                    except Exception as e:
                        raise SirilConnectionError(_("Error receiving data: {}").format(e)) from e
                return bytes(data)
            finally:
                self.sock.settimeout(original_timeout)

    def _send_command(self, command: _Command, data: Optional[bytes] = None, timeout: Optional[float] = DEFAULT_TIMEOUT) -> Tuple[Optional[int], Optional[bytes]]:
        """
        Send a command and receive response with optional timeout. Private method, not
        for use in scripts.

        Args:
            command: (_Command) Command to send
            data: (bytes) Optional data payload
            timeout: (int) Timeout in seconds for receive operations. None for
                     indefinite timeout

        Returns: Tuple (status, response_data)

        Raises: SirilConnectionError: if there is an error with the pipe / socket
                connection
        """

        if self.debug:
            timeout = None

        try:
            data_length = len(data) if data else 0
            if data_length > 65529:
                raise RuntimeError(_("Command data too long. Maximum command data 65529 bytes"))

            # Acquire lock before sending command
            if os.name == 'nt':
                win32event.WaitForSingleObject(self.command_lock, win32event.INFINITE)
            else:
                self.command_lock.acquire()

            try:
                # Pack command and length into fixed-size header
                header = struct.pack('!Bi', command, data_length)

                if os.name == 'nt':
                    # Create event for this write operation
                    event_handle = win32event.CreateEvent(None, True, False, None)
                    if not event_handle:
                        raise SirilConnectionError(_("Failed to create event for write operation"))

                    try:
                        # Create OVERLAPPED structure for this write
                        write_overlapped = pywintypes.OVERLAPPED()
                        write_overlapped.hEvent = event_handle

                        # Combine header and data into single buffer for atomic write
                        complete_message = header
                        if data and data_length > 0:
                            complete_message += data

                        # Send everything in one write operation
                        err, bytes_written = win32file.WriteFile(self.pipe_handle, complete_message, write_overlapped)
                        rc = win32event.WaitForSingleObject(event_handle, 3000)
                        if rc != win32event.WAIT_OBJECT_0:
                            raise SirilConnectionError(_("Timeout while sending message"))

                        # Ensure all bytes were written
                        bytes_transferred = win32file.GetOverlappedResult(self.pipe_handle, write_overlapped, True)
                        if bytes_transferred != len(complete_message):
                            raise SirilConnectionError(_("Incomplete write operation"))

                        # Wait for and receive complete response
                        response_header = self._recv_exact(5, timeout)  # Pass timeout
                        if not response_header:
                            return None, None
                        status, response_length = struct.unpack('!BI', response_header)

                        response_data = None
                        if response_length > 0:
                            response_data = self._recv_exact(response_length, timeout)  # Pass timeout
                            if not response_data:
                                return None, None

                        return status, response_data

                    finally:
                        # Clean up event handle
                        win32file.CloseHandle(event_handle)

                else:
                    # Socket implementation
                    msg = header
                    if data and data_length > 0:
                        msg += data

                    self.sock.sendall(msg)

                    response_header = self._recv_exact(5, timeout)  # Pass timeout
                    if not response_header:
                        return None, None

                    status, response_length = struct.unpack('!BI', response_header)
                    response_data = None
                    if response_length > 0:
                        response_data = self._recv_exact(response_length, timeout)  # Pass timeout
                        if not response_data:
                            return None, None

                    return status, response_data

            finally:
                # Always release the lock
                if os.name == 'nt':
                    win32event.ReleaseMutex(self.command_lock)
                else:
                    self.command_lock.release()

        except SirilError:
            raise # Simply raise any exceptions that are already subtypes of SirilError
        except Exception as e:
            raise SirilConnectionError(_("Error in _send_command(): {}").format(e)) from e

    def _execute_command(self, command: _Command, payload: Optional[bytes] = None, timeout: Optional[float] = DEFAULT_TIMEOUT) -> bool:
        """
        High-level method to execute a command and handle the response.
        Private method, not for end-user use.

        Args:
            command: The command to execute
            payload: Optional command payload
            timeout: Timeout for command execution. None for indefinite timeout.

        Returns:
            True if command was successful, False otherwise

        Raises:
            NoImageError: if no image is loaded when one should be,
            NoSequenceError: if no sequence is loaded when one should be,
            SirilError: for unknown errors.
        """
        if self.debug:
            timeout = None

        try:
            status, response = self._send_command(command, payload, timeout)

            if status is None:
                error_msg = response.decode('utf-8') if response else _("Error: no status returned")
                raise SirilError(error_msg)

            if status == _Status.NONE:
                # This indicates "allowed failure" - no data to return but not an error
                return None

            # Handle error responses
            if status == _Command.ERROR:
                if response:
                    error_msg = response.decode('utf-8', errors='replace').lower()
                    if "no image loaded" in error_msg:
                        raise NoImageError(_("No image is currently loaded in Siril"))
                    if "no seqeunce loaded" in error_msg:
                        raise NoSequenceError(_("No sequence is currently loaded in Siril"))
                    raise SirilError(_("Interface error in _execute_command: {}").format(error_msg))
                raise SirilError(_("Error: unknown interface error"))

            return True

        except SirilError:
            raise # Simply raise any exceptions that are already subtypes of SirilError
        except Exception as e:
            raise SirilError(f"Error in _execute_command(): {e}") from e

    def _request_data(self, command: _Command, payload: Optional[bytes] = None, timeout: Optional[float] = DEFAULT_TIMEOUT) -> Optional[bytes]:
        """
        High-level method to request small-volume data from Siril. The
        payload limit is 63336 bytes. For commands expected to return
        larger volumes of data, SHM should be used.
        Private method, not for direct use in scripts.

        Args:
            command: The data request command
            payload: Optional request parameters
            timeout: Timeout for data request. None for indefinite timeout.

        Returns:
            Requested data or None if no data to return

        Raises:
            NoImageError: if no image is loaded when one should be,
            NoSequenceError: if no sequence is loaded when one should be,
            SirilConnectionError: for unknown errors to do with the connection.
        """
        if self.debug:
            timeout = None

        try:
            status, response = self._send_command(command, payload, timeout)

            if status is None:
                error_msg = response.decode('utf-8') if response else _("Error: no status returned")
                raise SirilConnectionError(error_msg)

            if status == _Status.NONE:
                # This indicates "allowed failure" - no data to return but not an error
                return None

            # Handle error responses
            if status == _Command.ERROR:
                if response:
                    error_msg = response.decode('utf-8', errors='replace').lower()
                    if "no image loaded" in error_msg:
                        raise NoImageError(_("No image is currently loaded in Siril"))
                    if "no seqeunce loaded" in error_msg:
                        raise NoSequenceError(_("No sequence is currently loaded in Siril"))
                    raise SirilConnectionError(_("Interface error in _request_data: {}").format(error_msg))
                raise SirilConnectionError(_("Error: unknown interface error"))

            return response

        except SirilError:
            raise # Simply raise any exceptions that are already subtypes of SirilError
        except Exception as e:
            raise SirilError(f"Error in _request_data(): {e}") from e # Wrap others in a SirilError

    def _request_shm(self, size: int) -> Optional[_SharedMemoryInfo]:
        """
        Request Siril to create a shared memory allocation and return an
        object containing the details (name, size). Private method, not for use in
        scripts.

        Args:
            size: int specifying the size of shm buffer to create.

        Returns:
            _SharedMemoryInfo: Details of the shared memory allocation, or None if allocation failed.

        Raises:
            SharedMemoryError: If an invalid response is received or no data is returned
            ValueError: If the shared memory information is invalid
            SirilError: For other errors
        """
        size_data = struct.pack('!Q', size)  # Pack the size as a uint64_t in network byte order.

        try:
            # Send the shared memory request command.
            status, response = self._send_command(_Command.REQUEST_SHM, size_data)

            # Check the status and handle accordingly.
            if status == _Status.ERROR:
                raise SirilConnectionError(_("Failed to initiate shared memory transfer: invalid response"))

            if status == _Status.NONE:
                # Return None if no shared memory was allocated.
                return None

            if not response:
                raise SirilConnectionError(_("Failed to initiate shared memory transfer: no data received"))

            try:
                # Attempt to parse the shared memory info from the response buffer.
                shm_info = _SharedMemoryInfo.from_buffer_copy(response)
                return shm_info

            except (AttributeError, BufferError, ValueError) as e:
                # Catch parsing errors and raise a descriptive exception.
                raise ValueError(_("Invalid shared memory information received: {}").format(e)) from e

        except SirilConnectionError as connerr:
            # Let specific SharedMemoryErrors propagate
            raise SharedMemoryError("Error creating shared memory allocation") from connerr

        except Exception as e:
            raise SirilError(
                _("Unexpected error during shared memory allocation request: {}").format(e)
            ) from e

    def _map_shared_memory(self, name: str, size: int) -> SharedMemoryWrapper:
        """
        Create or open a shared memory mapping using SharedMemoryWrapper.
        Private method, not for direct use in scripts.

        Args:
            name: Name of the shared memory segment,
            size: Size of the shared memory segment in bytes

        Returns:
            SharedMemoryWrapper: A wrapper object for the shared memory segment

        Raises:
            SharedMemoryError: If the shared memory mapping fails
        """
        try:
            return SharedMemoryWrapper(name=name, size=size)
        except Exception as e:
            raise SharedMemoryError(_("Failed to create shared memory mapping: {}").format(e)) from e

    def command_error_message(self, status_code: CommandStatus) -> str:
        """
        Provides a string describing the return status from a Siril command.

        Args:
            status_code: The status code returned by the Siril command handler, or by
                  the CommandError exception.
        Returns:
            str: A string providing a description of the error code returned by a Siril command, for use in exception handling.
        """
        error_messages = {
            CommandStatus.CMD_NOT_FOUND: _("Command not found"),
            CommandStatus.CMD_NO_WAIT: _("Command does not wait for completion"),
            CommandStatus.CMD_NO_CWD: _("Current working directory not set"),
            CommandStatus.CMD_NOT_SCRIPTABLE: _("Command not scriptable"),
            CommandStatus.CMD_WRONG_N_ARG: _("Wrong number of arguments"),
            CommandStatus.CMD_ARG_ERROR: _("Argument error"),
            CommandStatus.CMD_SELECTION_ERROR: _("Selection error"),
            CommandStatus.CMD_OK: _("Command succeeded"),
            CommandStatus.CMD_GENERIC_ERROR: _("Generic error"),
            CommandStatus.CMD_IMAGE_NOT_FOUND: _("Image not found"),
            CommandStatus.CMD_SEQUENCE_NOT_FOUND: _("Sequence not found"),
            CommandStatus.CMD_INVALID_IMAGE: _("Invalid image"),
            CommandStatus.CMD_LOAD_IMAGE_FIRST: _("Load image first"),
            CommandStatus.CMD_ONLY_SINGLE_IMAGE: _("Command requires a single image to be loaded"),
            CommandStatus.CMD_NOT_FOR_SINGLE: _("Command not for single images"),
            CommandStatus.CMD_NOT_FOR_MONO: _("Command not for monochrome images"),
            CommandStatus.CMD_NOT_FOR_RGB: _("Command not for RGB images"),
            CommandStatus.CMD_FOR_CFA_IMAGE: _("Command only for CFA images"),
            CommandStatus.CMD_FILE_NOT_FOUND: _("File not found"),
            CommandStatus.CMD_FOR_PLATE_SOLVED: _("Command requires plate-solved image"),
            CommandStatus.CMD_NEED_INIT_FIRST: _("Initialization required first"),
            CommandStatus.CMD_ALLOC_ERROR: _("Memory allocation error"),
            CommandStatus.CMD_THREAD_RUNNING: _("Command thread already running"),
            CommandStatus.CMD_DIR_NOT_FOUND: _("Directory not found")
        }
        return error_messages.get(status_code, f"Unknown error code: {status_code}")

    def _get_bundle_path(self) ->Optional[str]:
        """
        Request the bundle path directory. This is a private method used
        to ensure that the correct DLL paths are preconfigured on Windows:
        it is not for use by scriptwriters.

        **It is an error to call this method on non-Windows OSes**

        Returns:
            The Siril bundle path as a string.

        Raises:
            SirilError: if the bundle path is not received.
        """

        response = self._request_data(_Command.GET_BUNDLE_PATH)

        if response is None:
            raise SirilError("Error in _get_bundle_path(): Failed to get bundle path - received null response")

        # Let the decode raise UnicodeDecodeError if it fails
        # Let the rstrip operation pass through any string operation errors
        return response.decode('utf-8').rstrip('\x00')

    def log(self, my_string: str, color:LogColor=LogColor.DEFAULT) -> bool:
        """
        Send a log message to Siril. The maximum message length is
        1022 bytes: longer messages will be truncated.

        Args:
            my_string: The message to log
            color: Defines the text color, defaults to white. See the documentation
            for LogColor for an explanation of which colors should be used for which
            purposes.

        Raises:
            SirilError: if the command fails
        """

        try:
            # Append a newline character to the string
            truncated_string = truncate_utf8(my_string, 1021) + '\n'
            # Convert string to bytes using UTF-8 encoding
            message_bytes = truncated_string.encode('utf-8')
            # Prepend the color byte
            packed_message = bytes([color.value]) + message_bytes
            self._execute_command(_Command.LOG_MESSAGE, packed_message)

        except SirilError:
            raise
        except Exception as e:
            raise SirilError(f"Error sending log message {message_bytes}: {e}") from e

    def _claim_thread(self) -> None:
        """
        Claim the processing thread. This prevents other processes using the
        processing thread to operate on the current Siril image. The preferred
        method of thread control is to use the image_lock() context manager
        rather than using this function manually, therefore this is now a
        private method.

        This function **must** always be called before starting any processing
        that will end with ``SirilInterface.set_image_pixeldata()``. The
        sequence of operations should be:

        * Call ``SirilInterface.claim_thread()`` by entering image_lock() context
        * If an exception is raised, handle it appropriately:
        - ProcessingThreadBusyError: the thread is already in use
        - ImageDialogOpenError: an image processing dialog is open
        * If no exception is raised, you have the thread claimed.
        * Now you can call ``SirilInterface.get_image()`` or ``get_image_pixeldata()``
        * Carry out your image processing
        * Call ``SirilInterface.set_image_pixeldata()``
        * Call ``SirilInterface.release_thread()`` by exiting image_lock() context

        Note that the thread should only be claimed when the script itself is
        operating on the Siril image data. If the script is calling a Siril command
        to alter the Siril image then the thread **must not** be claimed or the
        Siril command will be unable to acquire it, and will fail.

        Raises:
            ProcessingThreadBusyError: If the thread is already in use
            ImageDialogOpenError: If an image processing dialog is open
            SirilError: If any other error occurs while claiming the thread
        """
        try:
            status, response = self._send_command(_Command.CLAIM_THREAD, None)
            if status in (_Status.NONE, _Status.ERROR):
                response_str = response.decode('utf-8')
                print(f"image_lock: {response_str}", file=sys.stderr)

                if "processing thread is locked" in response_str:
                    raise ProcessingThreadBusyError("Processing thread is already in use")
                if "processing dialog is open" in response_str:
                    raise ImageDialogOpenError("Image processing dialog is open")
                raise SirilError(f"Failed to claim processing thread: {response_str}")
        except (ProcessingThreadBusyError, ImageDialogOpenError):
            # Re-raise specific exceptions
            raise
        except Exception as e:
            raise SirilError(f"Error claiming processing thread: {e}") from e

    def _release_thread(self) -> None:
        """
        Release the processing thread. This permits other processes to use the
        processing thread to operate on the current Siril image. The preferred
        method of thread control is to use the image_lock() context manager
        rather than using this function manually, therefore this is now a
        private method.

        This function **must** always be called after completing any processing
        that has updated the image loaded in Siril. The sequence of operations
        should be:

        * Call ``SirilInterface.claim_thread()`` by entering image_lock() context
        * If an exception is raised, handle it appropriately
        * If no exception is raised, you have the thread claimed.
        * Now you can call ``SirilInterface.get_image()`` or ``get_image_pixeldata()``
        * Carry out your image processing
        * Call ``SirilInterface.set_image_pixeldata()``
        * Call ``SirilInterface.release_thread()`` by exiting image_lock() context

        Raises:
            SirilError: if an error occurred in releasing the thread
        """
        retval = self._execute_command(_Command.RELEASE_THREAD, None)
        if retval is not True:
            raise SirilError(_("Error trying to release the processing thread. "
                "It will be released when the script terminates."))

    def image_lock(self):
        """
        A context manager that handles claiming and releasing the processing thread.

        This method is designed to be used with a `with` statement to ensure that
        the thread is properly claimed before processing and released after processing,
        even if an exception occurs during processing. It is preferable to use this
        context manager rather than manually calling claim_thread() and
        release_thread() as the context manager will ensure correct cleanup if an
        exception occurs.

        Note that the image_lock() context should only be entered when the script itself
        is operating on the Siril image data. If the script is calling a Siril command
        to alter the Siril image then the context **must not** be entered or the Siril
        command will be unable to acquire the processing thread and will fail.

        Example:

        .. code-block:: python

            try:
                with siril.image_lock():
                    # Get image data
                    image_data = self.get_image_pixeldata()
                    # Process image data
                    processed_data = some_processing_function(image_data)
                    # Set the processed image data
                    siril.set_image_pixeldata(processed_data)
            except ProcessingThreadBusyError:
                # Handle busy thread case
                pass
            except ImageDialogOpenError:
                # Handle open dialog case
                pass

        Raises:
            ProcessingThreadBusyError: If the thread is already in use
            ImageDialogOpenError: If an image processing dialog is open
            SirilError: If the thread cannot be claimed or released for other reasons
        """
        class ImageLockContext:
            """ Class to manage the image_lock context """
            def __init__(self, outer_self: SirilInterface):
                self.outer_self = outer_self
                self.claimed = False

            def __enter__(self):
                # This will raise the appropriate exception if claim fails
                self.outer_self._claim_thread()
                self.claimed = True
                return self.outer_self

            def __exit__(self, exc_type, exc_val, exc_tb):
                if self.claimed:
                    # Always release the thread and raise any exceptions that occur
                    self.outer_self._release_thread()
                    self.claimed = False
                # Don't suppress exceptions
                return False

        return ImageLockContext(self)

    def confirm_messagebox(self, title: str, message: str, confirm_label: str) -> bool:
        """
        Create a modal confirmation message dialog in Siril and wait for the response.

        Args:
            title: The title to display in the message box (up to 256 characters)
            message: The message to display in the message box (up to 1021 characters)
            confirm_label: The label to display in the message box confirmation button (OK, Yes, Confirm etc.) (Up to 24 characters)

        Returns:
            bool: True if the message box confirmation button was clicked, False otherwise

        Raises:
            DataError: if no response was received,
            SirilError: if another error occurred.
        """
        try:
            # Truncate strings to byte-size limits minus null terminator
            encoded_title = truncate_utf8(title, 255).encode('utf-8') + b'\0'
            encoded_message = truncate_utf8(message, 1020).encode('utf-8') + b'\0'
            encoded_label = truncate_utf8(confirm_label, 23).encode('utf-8') + b'\0'

            # Final buffer
            message_bytes = encoded_title + encoded_message + encoded_label

            # Call the command with the encoded data
            response = self._request_data(_Command.CONFIRM_MESSAGEBOX, message_bytes, timeout=None)

            if response is None:
                raise DataError(_("Error in confirm_messagebox(): no response"))

            return bool(int.from_bytes(response, byteorder='little'))

        except Exception as e:
            raise SirilError(f"Error in confirm_messagebox(): {e}") from e

    def _messagebox(self, my_string: str, cmd_type: int, modal: Optional[bool] = False) -> bool:
        """
        Send a message to Siril for display in a messagebox.
        Private method: the public API is provided by error_messagebox(),
        warning_messagebox() and info_messagebox().

        Args:
            my_string: The message to display in the message box
            cmd_type: Sets whether to show an error, warning or info messagebox
            modal: Whether or not the message box is modal

        Returns:
            bool: True if the error was successfully displayed, False otherwise

        Raises:
            SirilError: if an error occurred.
         """

        # Append a newline character to the string
        message_bytes = (truncate_utf8(my_string, 1021) + '\n').encode('utf-8')
        if modal:
            return self._execute_command(cmd_type, message_bytes, timeout = None)
        return self._execute_command(cmd_type, message_bytes)

    def error_messagebox(self, my_string: str, modal: Optional[bool] = False) -> bool:
        """
        Send an error message to Siril. The maximum message length is
        1022 bytes: longer messages will be truncated (but this is more than
        enough for an error message box). Note that the error message box is
        not modal by default: this is intended for displaying an error message
        more prominently than using the Siril log prior to quitting the
        application.

        Args:
            my_string: The message to display in the error message box
            modal: Sets whether or not the message box should be modal and
                   wait for completion or non-modal and allow the script to
                   continue execution. Note that although a modal message box will
                   block execution of the script, if a TKinter main loop is
                   running events will continue to queue up, so if the message
                   box is triggered by clicking a button then the user may
                   click it while the message box is shown and trigger a second
                   message box which will display immediately the first one is
                   closed.

        Returns:
            bool: True if the error was successfully displayed, False otherwise

        Raises:
            SirilError: if an error occurred.
         """
        try:
            cmd_type = _Command.ERROR_MESSAGEBOX_MODAL if modal else _Command.ERROR_MESSAGEBOX
            return self._messagebox(my_string, cmd_type, modal)

        except Exception as e:
            raise SirilError(f"Error in error_messagebox(): {e}") from e

    def info_messagebox(self, my_string: str, modal: Optional[bool] = False) -> bool:
        """
        Send an information message to Siril. The maximum message length is
        1022 bytes: longer messages will be truncated. This is intended for
        displaying informational messages more prominently than using the Siril log.

        Args:
            my_string: The message to display in the info message box
            modal: Sets whether or not the message box should be modal and
                   wait for completion or non-modal and allow the script to
                   continue execution. Note that although a modal message box will
                   block execution of the script, if a TKinter main loop is
                   running events will continue to queue up, so if the message
                   box is triggered by clicking a button then the user may
                   click it while the message box is shown and trigger a second
                   message box which will display immediately the first one is
                   closed.

        Returns:
            bool: True if the info was successfully displayed, False otherwise

        Raises:
            SirilError: if an error occurred.
         """

        try:
            cmd_type = _Command.INFO_MESSAGEBOX_MODAL if modal else _Command.INFO_MESSAGEBOX
            return self._messagebox(my_string, cmd_type, modal)

        except Exception as e:
            raise SirilError(f"Error in info_messagebox(): {e}") from e

    def warning_messagebox(self, my_string: str, modal: Optional[bool] = False) -> bool:
        """
        Send a warning message to Siril. The maximum message length is
        1022 bytes: longer messages will be truncated. This is intended for
        displaying warning messages more prominently than using the Siril log.

        Args:
            my_string: The message to display in the warning message box
            modal: Sets whether or not the message box should be modal and
                   wait for completion or non-modal and allow the script to
                   continue execution. Note that although a modal message box will
                   block execution of the script, if a TKinter main loop is
                   running events will continue to queue up, so if the message
                   box is triggered by clicking a button then the user may
                   click it while the message box is shown and trigger a second
                   message box which will display immediately the first one is
                   closed.

        Returns:
            bool: True if the warning was successfully displayed, False otherwise

        Raises:
            SirilError: if an error occurred.
         """

        try:
            cmd_type = _Command.WARNING_MESSAGEBOX_MODAL if modal else _Command.WARNING_MESSAGEBOX
            return self._messagebox(my_string, cmd_type, modal)

        except Exception as e:
            raise SirilError(f"Error in warning_messagebox(): {e}") from e

    def undo_save_state(self, my_string: str) -> bool:
        """
        Saves an undo state. The maximum message length is 70 bytes: longer
        messages will be truncated. Requires a single image to be loaded.

        Args:
            my_string: The message to log in FITS HISTORY

        Returns:
            bool: True if the message was successfully logged, False otherwise

        Raises:
            SirilError: if an error occurred.

        """

        try:
            # Append a newline character to the string
            # Convert string to bytes using UTF-8 encoding
            message_bytes = my_string.encode('utf-8')
            return self._execute_command(_Command.UNDO_SAVE_STATE, message_bytes, timeout=60)

        except Exception as e:
            raise SirilError(f"Error in undo_save_state(): {e}") from e

    def update_progress(self, message: str, progress: float) -> bool:
        """
        Send a progress update to Siril with a message and completion percentage.

        Args:
            message: Status message to display,
            progress: Progress value in the range 0.0 to 1.0. The following special
                      values can be used: -1.0 will pulsate the progress bar, and
                      -2.0 will update the progress bar text but will not update
                      the progress shown in the bar.

        Raises:
            ValueError: If the progress argument is out of range,
            SirilError: For any other errors.
        """

        try:
            # Validate progress value
            if not (0.0 <= progress <= 1.0 or progress == -1.0 or progress == -2.0):
                raise ValueError(_("Progress value must be between 0.0 and 1.0"))

            # Convert string to UTF-8 bytes
            message_bytes = message.encode('utf-8')

            # Create payload: network-order float followed by string
            # '!f' for network byte order 32-bit float
            float_bytes = struct.pack('!f', progress)

            # Combine float and string bytes
            payload = float_bytes + message_bytes

            self._execute_command(_Command.UPDATE_PROGRESS, payload)

        except ValueError:
            raise
        except Exception as e:
            raise SirilError(f"Error in update_progress(): {e}") from e

    def reset_progress(self) -> bool:
        """
        Resets the Siril progress bar.

        Args:
            none

        Raises:
            SirilError: For any errors.
        """

        self.update_progress("", 0.0)

    def cmd(self, *args: str):
        """
        Send a command to Siril to be executed. The range of available commands can
        be found by checking the online documentation. The command and its arguments
        are provided as a list of strings.

        Args:
            *args: Variable number of string arguments to be combined into a command

        Raises:
            DataError: If no response (or an incorrect response) was received,
            CommandError: If the command returns an error status code,
            SirilError: If any other error occurs during execution.

        Example:
            .. code-block:: python

                siril.cmd("ght", "-D=0.5", "-b=2.0")
        """
        try:
            # Join arguments with spaces between them
            command_string = " ".join(str(arg) for arg in args)
            # Convert to bytes for transmission
            command_bytes = command_string.encode('utf-8')

            # Use _request_data instead of _execute_command
            response = self._request_data(_Command.SEND_COMMAND, command_bytes, timeout=None)

            if response is None:
                raise DataError(_(f"Error: _request_data({args}) failed."))

            # Convert response bytes to integer from network byte order
            if len(response) == 4:  # Valid response is int32_t ie 4 bytes
                status_code = int.from_bytes(response, byteorder='big')

                # If the status code is a no-error code, return
                if status_code in (CommandStatus.CMD_OK, CommandStatus.CMD_NO_WAIT):
                    return  # Command executed successfully
                # Else raise an exception with a helpful error message and the status code
                error_message = self.command_error_message(status_code)
                raise CommandError(_(f"Command '{args[0]}' failed: {error_message}"), status_code)

            # Handle case where response doesn't contain enough bytes for a status code
            raise DataError(_(f"Error: Response from {args[0]} incorrect size to contain a status code."))

        except CommandError:
            raise
        except Exception as e:
            raise SirilError(_(f"Error in cmd(): {e}")) from e

    def set_siril_selection(self,
                            x: Optional[int] = None,
                            y: Optional[int] = None,
                            w: Optional[int] = None,
                            h: Optional[int] = None,
                            selection: Optional[Tuple[int, int, int, int]] = None) -> bool:
        """
        Set the image selection in Siril using the provided coordinates and dimensions.

        Args:
            x: X-coordinate of the selection's top-left corner (must be provided with y, w, h)
            y: Y-coordinate of the selection's top-left corner (must be provided with x, w, h)
            w: Width of the selection (must be provided with x, y, h)
            h: Height of the selection (must be provided with x, y, w)
            selection: A tuple of (x, y, w, h) as returned by get_siril_selection()

        Raises:
            SirilError: if an error occurred.
            ValueError: if parameters are not properly provided.

        Returns:
            bool: True if the selection was set successfully
        """
        try:
            # Check if selection tuple is provided
            if selection is not None:
                # Make sure individual coordinates are not also provided
                if any(param is not None for param in (x, y, w, h)):
                    raise ValueError("Cannot provide both selection tuple and individual coordinates")

                if len(selection) != 4:
                    raise ValueError("Selection tuple must contain exactly 4 values (x, y, w, h)")

                x, y, w, h = selection
            # Otherwise check if all individual coordinates are provided
            elif all(param is not None for param in (x, y, w, h)):
                pass  # Just use the provided x, y, w, h values
            else:
                raise ValueError("Either provide a selection tuple or all four coordinate parameters (x, y, w, h)")

            # Pack the coordinates and dimensions into bytes using network byte order (!)
            payload = struct.pack('!IIII', x, y, w, h)
            self._execute_command(_Command.SET_SELECTION, payload)
            return True
        except Exception as e:
            raise SirilError(f"Error in set_siril_selection(): {e}") from e

    def get_siril_selection(self) -> Optional[Tuple[int, int, int, int]]:

        """
        Request the image selection from Siril.

        Returns:
            A tuple (x, y, height, width) representing the current selection, or
            None if no selection is made.

        Raises:
            SirilError: if an error occurred.
        """

        response = self._request_data(_Command.GET_SELECTION)
        if response is None:
            return None

        try:
            # Assuming the response is in the format: x (4 bytes), y (4 bytes), w (4 bytes), h (4 bytes)
            x, y, w, h = struct.unpack('!IIII', response)
            return x, y, w, h  # Returning as (x, y, w, h)
        except struct.error as e:
            raise SirilError(_("Error occurred in get_siril_selection(): {}").format(e)) from e

    def get_siril_active_vport(self) -> Optional[SirilVport]:

        """
        Request the active viewport from Siril.

        Returns:
            A SirilVport representing the active vport:

            - sirilpy.SirilVport.RED / sirilpy.SirilVport.MONO
            - sirilpy.SirilVport.GREEN,
            - sirilpy.SirilVport.BLUE,
            - sirilpy.SirilVport.RGB

            Note that RED and MONO share the
            same IntEnum value, so there is no difference between a test for
            one and the other; the two enum labels are provided solely to aid
            code legibility.

        Raises:
            DataError: if no response or an invalid response is received,
            SirilError: if an error occurred.
        """

        response = self._request_data(_Command.GET_ACTIVE_VPORT)
        if response is None:
            raise DataError(_("Error: no response or invalid response received in get_siril_active_vport()"))

        try:
            # Assuming the response is in the format: !I
            vport = struct.unpack('!I', response)[0]
            return vport
        except struct.error as e:
            raise SirilError(_("Error occurred in get_siril_active_vport()")) from e

    def get_image_shape(self) -> Optional[Tuple[int, int, int]]:

        """
        Request the shape of the image from Siril.

        Returns:
            A tuple (channels, height, width) representing the shape of the image,
            or None if no image shape is available to return.

        Raises: SirilError: if an error occurred.
        """

        try:
            response = self._request_data(_Command.GET_DIMENSIONS)

            if response is None:
                return None

            # Assuming the response is in the format: width (4 bytes), height (4 bytes), nb_channels (4 bytes)
            width, height, channels = struct.unpack('!III', response)
            return channels, height, width  # Returning as (channels, height, width)
        except Exception as e:
            raise SirilError(_("Error occurred in get_image_shape()")) from e

    def get_selection_star(self, shape: Optional[list[int]] = None, \
        channel: Optional[int] = None, assume_centred: bool = False) -> Optional[PSFStar]:
        """
        Retrieves a PSFStar star model from the current selection in Siril.
        Only a single PSFStar is returned: if there are more than one in the
        selection, the first one identified by Siril's internal star detection
        algorithm is returned.
        **Update**: from sirilpy 1.0.4 this method uses Siril's photometry functions to
        try to provide photometrically accurate values for PSFStar.mag, PSFStar.s_mag
        and PSFStar.SNR. If photometry succeeded and no saturated pixels were detected
        then PSFStar.phot_is_valid will be True, otherwise it will be False.
        Args:
            shape: Optional list of [x, y, w, h] specifying the selection to
                retrieve from. w x h must not exceed 300 px x 300 px.
                If provided, looks for a star in the specified selection
                If None, looks for a star in the selection already made in
                Siril, if one is made.
            channel: Optional int specifying the channel to retrieve from.
                    If provided 0 = Red / Mono, 1 = Green, 2 = Blue. If the
                    channel is omitted the current viewport will be used if
                    in GUI mode, or if not in GUI mode the method will fall back
                    to channel 0
            assume_centred: Optional bool specifying whether to assume the star
                        is already centred in the selection. Defaults to False.
        Returns:
            PSFStar: the PSFStar object representing the star model, or None if
                    no star is detected in the selection.
        Raises:
            ValueError: If an invalid shape is provided,
            NoImageError: If no image is loaded,
            SirilConnectionError: If a communication error occurs,
            SirilError: If any other error occurred during execution.
        """
        try:
            # Sentinel values for not-provided parameters
            SENTINEL_VALUE = 0xFFFFFFFF  # -1 as unsigned int

            # Default values
            x = y = w = h = SENTINEL_VALUE
            channel_val = SENTINEL_VALUE if channel is None else channel
            centred_flag = 1 if assume_centred else 0

            # Validate and set shape if provided
            if shape is not None:
                if len(shape) != 4:
                    raise ValueError(_("Shape must be a list of [x, y, w, h]"))
                if any(not isinstance(v, int) for v in shape):
                    raise ValueError(_("All shape values must be integers"))
                if any(v < 0 for v in shape):
                    raise ValueError(_("All shape values must be non-negative"))
                x, y, w, h = shape

            # Validate channel if provided
            if channel is not None and (channel < 0 or channel > 2):
                raise ValueError(_("Channel must be 0 (Red/Mono), 1 (Green), or 2 (Blue)"))

            # Always pack the same payload: x, y, w, h, channel, centred_flag (24 bytes)
            shape_data = struct.pack('!IIIIII', x, y, w, h, channel_val, centred_flag)

            status, response = self._send_command(_Command.GET_STAR_IN_SELECTION, shape_data)

            # Handle error responses
            if status == _Status.ERROR:
                if response:
                    error_msg = response.decode('utf-8', errors='replace')
                    if "no image loaded" in error_msg.lower():
                        raise NoImageError(_("No image is currently loaded in Siril"))
                    if "invalid selection" in error_msg.lower():
                        raise ValueError(_("No selection is currently made in Siril"))
                    raise SirilConnectionError(_("Server error: {}").format(error_msg))
                raise SirilConnectionError(_("Failed to transfer star data: Empty response"))

            if status == _Status.NONE:
                return None

            if not response:
                raise RuntimeError(_("Failed to transfer star data: No data received"))

            return PSFStar.deserialize(response)

        except (SirilError, ValueError):
            raise
        except Exception as e:
            raise SirilError(_("Error in get_selection_star(): {}").format(e)) from e

    def get_selection_stats(self, shape: Optional[list[int]] = None, \
        channel: Optional[int] = None) -> Optional[PSFStar]:
        """
        Retrieves statistics for the current selection in Siril. Requires a single
        image or a sequence to be loaded.

        Args:
            shape: Optional list of [x, y, w, h] specifying the selection to
                retrieve from.
                If provided, looks for a star in the specified selection
                If None, looks for a star in the selection already made in Siril,
                if one is made.
            channel: Optional int specifying the channel to retrieve from.
                If provided 0 = Red / Mono, 1 = Green, 2 = Blue. If the
                channel is omitted the current viewport will be used if
                in GUI mode, or if not in GUI mode the method will fall back
                to channel 0
        Returns:
            ImageStats: the ImageStats object representing the selection statistics.
        Raises:
            SirilError: If an error occurred during processing,
            ValueError: If an invalid shape is provided.
        """
        # Validate shape if provided
        shape_data = None
        if shape is not None:
            if len(shape) != 4:
                raise ValueError(_("Shape must be a list of [x, y, w, h]"))
            if any(not isinstance(v, int) for v in shape):
                raise ValueError(_("All shape values must be integers"))
            if any(v < 0 for v in shape):
                raise ValueError(_("All shape values must be non-negative"))
            # Pack shape data for the command
            if len(shape) == 4:
                if channel is None:
                    shape_data = struct.pack('!IIII', *shape)
                else:
                    shape_data = struct.pack('!IIIII', *shape, channel)
            else:
                raise ValueError(_("Incorrect shape data provided: must be (x,y,w,h)"))
        elif channel is not None:
            # Pack only channel data when no shape is provided but channel is specified
            shape_data = struct.pack('!I', channel)

        status, response = self._send_command(_Command.GET_STATS_FOR_SELECTION, shape_data)
        # Handle error responses
        if status == _Status.ERROR:
            if response:
                error_msg = response.decode('utf-8', errors='replace')
                if "no image loaded" in error_msg.lower():
                    raise NoImageError(_("No image is currently loaded in Siril"))
                if "invalid selection" in error_msg.lower():
                    raise ValueError(_("No selection is currently made in Siril"))
                raise SirilConnectionError(_("Server error: {}").format(error_msg))
            raise SirilError(_("Error in get_selection_stats(): Failed to transfer stats data: Empty response"))
        if status == _Status.NONE:
            return None
        if not response:
            raise SirilError(_("Error in get_selection_stats(): No data received"))
        try:
            return ImageStats.deserialize(response)
        except Exception as e:
            raise SirilError(_("Error in get_selection_stats(): {e}")) from e

    def pix2radec(self, x: float, y: float) -> Optional[Tuple[float, float]]:
        """
        Converts a pair of pixel coordinates into RA and dec coordinates using the
        WCS of the image loaded in Siril. This requires that an image is loaded in
        Siril and that it has been platesolved (i.e. it has a WCS solution).

        Args:
            x: float: provides the x coordinate to be converted
            y: float: provides the y coordinate to be converted

        Returns:
            Tuple[float, float]: (RA, Dec) as a Tuple of two floats.

        Raises:
            NoImageError: If no image or sequence is loaded,
            ValueError: If the image or loaded sequence frame is not plate solved,
            SirilError: For errors during pix2radec execution.
        """
        try:
            shape_data = struct.pack('!2d', x, y)
            status, response = self._send_command(_Command.PIX2WCS, shape_data)
            # Handle error responses
            if status == _Status.ERROR:
                if response:
                    error_msg = response.decode('utf-8', errors='replace')
                    if "no image loaded" in error_msg.lower():
                        raise NoImageError(_("No image is currently loaded in Siril"))
                    if "not plate solved" in error_msg.lower():
                        raise ValueError(_("Siril image is not plate solved"))
                    raise SirilConnectionError(_("Server error: {}").format(error_msg))
                raise SirilConnectionError(_("Failed to transfer coordinates: Empty response"))

            if status == _Status.NONE:
                return None

            if not response:
                raise SirilConnectionError(_("Failed to transfer coordinates: No data received"))
            # Define the format string for unpacking the C struct
            # '!' for network byte order (big-endian)
            # 'd' for double (all floating point values)
            format_string = '!2d'  # '!' ensures network byte order

            # Calculate expected size
            expected_size = struct.calcsize(format_string)

            # Verify we got the expected amount of data
            if len(response) != expected_size:
                raise SirilConnectionError(f"Received data size {len(response)} doesn't match expected size {expected_size}")

            # Unpack the binary data
            return struct.unpack(format_string, response)

        except (NoImageError, ValueError):
            raise
        except Exception as e:
            raise SirilError(_("Error in pix2radec(): {e}")) from e

    def radec2pix(self, ra: float, dec: float) -> Optional[Tuple[float, float]]:
        """
        Converts a pair of RA,dec coordinates into image pixel coordinates using the
        WCS of the image loaded in Siril. This requires that an image is loaded in
        Siril and that it has been platesolved (i.e. it has a WCS solution).

        Args:
            ra: float: provides the RA coordinate to be converted
            dec: float: provides the dec coordinate to be converted

        Returns:
            Tuple[float, float]: [x, y] as a Tuple of two floats.

        Raises:
            NoImageError: If no image or sequence is loaded,
            ValueError: If the image or loaded sequence frame is not plate solved,
            SirilError: For errors during radec2pix execution.
        """
        try:
            shape_data = struct.pack('!2d', ra, dec)
            status, response = self._send_command(_Command.WCS2PIX, shape_data)
            # Handle error responses
            if status == _Status.ERROR:
                if response:
                    error_msg = response.decode('utf-8', errors='replace')
                    if "no image loaded" in error_msg.lower():
                        raise NoImageError(_("No image is currently loaded in Siril"))
                    if "not plate solved" in error_msg.lower():
                        raise ValueError(_("Siril image is not plate solved"))
                    raise SirilConnectionError(_("Server error: {}").format(error_msg))
                raise SirilConnectionError(_("Failed to transfer coordinates: Empty response"))

            if status == _Status.NONE:
                return None

            if not response:
                raise SirilConnectionError(_("Failed to transfer coordinates: No data received"))
            try:
                # Define the format string for unpacking the C struct
                # '!' for network byte order (big-endian)
                # 'd' for double (all floating point values)
                format_string = '!2d'  # '!' ensures network byte order

                # Calculate expected size
                expected_size = struct.calcsize(format_string)

                # Verify we got the expected amount of data
                if len(response) != expected_size:
                    raise SirilConnectionError(f"Received data size {len(response)} doesn't match expected size {expected_size}")

                # Unpack the binary data
                values = struct.unpack(format_string, response)
            except struct.error as e:
                raise SirilError(_("Error unpacking data: {}").format(e)) from e
            return values

        except (ValueError, NoImageError):
            raise
        except Exception as e:
            raise SirilError(_("Error in radec2pix(): {e}")) from e

    def get_image_pixeldata(self, shape: Optional[list[int]] = None, preview: Optional[bool] = False, linked: Optional[bool] = False) -> Optional[np.ndarray]:

        """
        Retrieves the pixel data from the image currently loaded in Siril.

        Args:
            shape: Optional list of [x, y, w, h] specifying the region to retrieve.
                   If provided, gets pixeldata for just that region.
                   If None, gets pixeldata for the entire image.
            preview: optional bool specifying whether to get pixeldata as a preview
                     (i.e. 8-bit autostretched data) or as real image data. Defaults
                     to False (i.e. real image data)
            linked: optional bool specifying whether the autostretch preview should
                    be linked or unlinked. If preview == False then this option is
                    ignored.

        Returns:
            numpy.ndarray: The image data as a numpy array

        Raises:
            NoImageError: If no image is currently loaded,
            ValueError: If an invalid shape is provided,
            DataError: if the array cannot be reshaped to the correct dimensions,
            SirilError: For other errors during pixel data retrieval,
        """

        shm = None
        try:
            preview_data = struct.pack('!??', preview, linked)
            # Validate shape if provided
            if shape is not None:
                if len(shape) != 4:
                    raise ValueError(_("Shape must be a list of [x, y, w, h]"))
                if any(not isinstance(v, int) for v in shape):
                    raise ValueError(_("All shape values must be integers"))
                if any(v < 0 for v in shape):
                    raise ValueError(_("All shape values must be non-negative"))

                # Pack shape data for the command
                shape_data = struct.pack('!IIII', *shape)
                command = _Command.GET_PIXELDATA_REGION
            else:
                shape_data = None
                command = _Command.GET_PIXELDATA

            # Request shared memory setup
            payload = preview_data + shape_data if shape_data else preview_data
            status, response = self._send_command(command, payload)

            # Handle error responses
            if status == _Command.ERROR:
                if response:
                    error_msg = response.decode('utf-8', errors='replace')
                    if "no image loaded" in error_msg.lower():
                        raise NoImageError(_("Error in get_image_pixeldata(): no image is currently loaded in Siril"))
                    raise RuntimeError(_("Server error: {}").format(error_msg))
                raise SharedMemoryError(_("Failed to initiate shared memory transfer: Empty response"))

            if not response:
                raise SharedMemoryError(_("Failed to initiate shared memory transfer: No data received"))

            try:
                # Parse the shared memory information
                shm_info = _SharedMemoryInfo.from_buffer_copy(response)
            except (AttributeError, BufferError, ValueError) as e:
                raise SharedMemoryError(_("Invalid shared memory information received: {}").format(e)) from e

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
                raise SharedMemoryError(_("Failed to map shared memory: {}").format(e)) from e

            buffer = bytearray(shm.buf)[:shm_info.size]
            # Create numpy array from shared memory
            if preview:
                dtype = np.uint8
            else:
                dtype = np.float32 if shm_info.data_type == 1 else np.uint16
            try:
                arr = np.frombuffer(buffer, dtype=dtype)
            except (BufferError, ValueError, TypeError) as e:
                raise SharedMemoryError(_("Failed to create array from shared memory: {}").format(e)) from e

            # Validate array size matches expected dimensions
            expected_size = shm_info.width * shm_info.height * shm_info.channels
            if arr.size < expected_size:
                raise ValueError(
                    f"Data size mismatch: got {arr.size} elements, "
                    f"expected {expected_size} for dimensions "
                    f"{shm_info.width}x{shm_info.height}x{shm_info.channels}"
                )

            # Reshape the array according to the image dimensions
            try:
                if shm_info.channels > 1:
                    arr = arr.reshape((shm_info.channels, shm_info.height, shm_info.width))
                else:
                    arr = arr.reshape((shm_info.height, shm_info.width))
            except ValueError as e:
                raise DataError(_("Failed to reshape array to image dimensions: {}").format(e)) from e

            # Make a copy of the data since we'll be releasing the shared memory
            result = np.copy(arr)

            return result

        except (ValueError, SirilError):
            # Re-raise NoImageError without wrapping
            raise
        except Exception as e:
            # Wrap all other exceptions with context
            raise SirilError(_("Error in get_image_pixeldata(): {}").format(e)) from e
        finally:
            # Clean up shared memory using the wrapper's methods
            if shm is not None:
                try:
                    shm.close()  # First close the memory mapping as we have finished with it
                    # (We don't unlink it as C wll do that)

                    # Signal that Python is done with the shared memory and wait for C to finish
                    if not self._execute_command(_Command.RELEASE_SHM, shm_info):
                        raise SirilError(_("Failed to cleanup shared memory"))

                except Exception:
                    pass

    def get_seq_frame_pixeldata(self, frame: int, shape: Optional[List[int]] = None,
                                preview: Optional[bool] = False, linked: Optional[bool] = False) -> Optional[np.ndarray]:

        """
        Retrieves the pixel data from a frame in the sequence currently loaded in Siril.

        Args:
            frame: selects the frame to retrieve pixel data from. This
                uses a 0-based indexing scheme, i.e. the first frame is frame
                number 0, not frame numer 1.
            shape: Optional list of [x, y, w, h] specifying the region to retrieve.
                If provided, gets pixeldata for just that region.
                If None, gets pixeldata for the entire image.
            preview: optional bool specifying whether to get pixeldata as a preview
                (i.e. 8-bit autostretched data) or as real image data. Defaults
                to False (i.e. real image data).
            linked: optional bool specifying whether the autostretched preview should
                    be linked or unlinked. This option is ignored if preview is not True

        Returns:
            numpy.ndarray: The image data as a numpy array

        Raises:
            ValueError: If an invalid shape is provided,
            DataError: if the array cannot be reshaped to the correct dimensions,
            SirilError: For other errors during pixel data retrieval.
        """

        shm = None
        # Convert channel number to network byte order bytes

        try:
            preview_payload = struct.pack('!??', preview, linked)
            frame_payload = struct.pack('!I', frame)
            # Validate shape if provided
            if shape is not None:
                if len(shape) != 4:
                    raise ValueError(_("Shape must be a list of [x, y, w, h]"))
                if any(not isinstance(v, int) for v in shape):
                    raise ValueError(_("All shape values must be integers"))
                if any(v < 0 for v in shape):
                    raise ValueError(_("All shape values must be non-negative"))

                # Pack shape data for the command
                shape_data = struct.pack('!IIII', *shape)
            else:
                shape_data = None
            payload = preview_payload + frame_payload + shape_data if shape_data else preview_payload + frame_payload
            # Request shared memory setup
            status, response = self._send_command(_Command.GET_SEQ_PIXELDATA, payload)

            # Handle error responses
            if status == _Command.ERROR:
                if response:
                    error_msg = response.decode('utf-8', errors='replace')
                    if "no sequence loaded" in error_msg.lower():
                        raise NoSequenceError(_("No sequence is currently loaded in Siril"))
                    raise SirilConnectionError(_("Server error: {}").format(error_msg))
                raise SharedMemoryError(_("Failed to initiate shared memory transfer: Empty response"))

            if not response:
                raise SharedMemoryError(_("Failed to initiate shared memory transfer: No data received"))

            try:
                # Parse the shared memory information
                shm_info = _SharedMemoryInfo.from_buffer_copy(response)
            except (AttributeError, BufferError, ValueError) as e:
                raise SharedMemoryError(_("Invalid shared memory information received: {}").format(e)) from e

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
                raise SharedMemoryError(_("Failed to map shared memory: {}").format(e)) from e

            buffer = bytearray(shm.buf)[:shm_info.size]
            # Create numpy array from shared memory
            if preview:
                dtype = np.uint8
            else:
                dtype = np.float32 if shm_info.data_type == 1 else np.uint16
            try:
                arr = np.frombuffer(buffer, dtype=dtype)
            except (BufferError, ValueError, TypeError) as e:
                raise SharedMemoryError(_("Failed to create array from shared memory: {}").format(e)) from e

            # Validate array size matches expected dimensions
            expected_size = shm_info.width * shm_info.height * shm_info.channels
            if arr.size < expected_size:
                raise ValueError(
                    f"Data size mismatch: got {arr.size} elements, "
                    f"expected {expected_size} for dimensions "
                    f"{shm_info.width}x{shm_info.height}x{shm_info.channels}"
                )

            # Reshape the array according to the image dimensions
            try:
                if shm_info.channels > 1:
                    arr = arr.reshape((shm_info.channels, shm_info.height, shm_info.width))
                else:
                    arr = arr.reshape((shm_info.height, shm_info.width))
            except ValueError as e:
                raise SirilError(_("Failed to reshape array to image dimensions: {}").format(e)) from e

            # Make a copy of the data since we'll be releasing the shared memory
            result = np.copy(arr)

            return result

        except (ValueError, SirilError):
            raise
        except Exception as e:
            # Wrap all other exceptions with context
            raise SirilError(_("Error retrieving pixel data: {}").format(e)) from e
        finally:
            # Clean up shared memory using the wrapper's methods
            if shm is not None:
                try:
                    shm.close()  # First close the memory mapping as we have finished with it
                    # (We don't unlink it as C wll do that)

                    # Signal that Python is done with the shared memory and wait for C to finish
                    if not self._execute_command(_Command.RELEASE_SHM, shm_info):
                        raise SirilError(_("Failed to cleanup shared memory"))

                except Exception:
                    pass

    def xy_plot(self, plot_data: PlotData, display=True, save=False):
        """
        Serialize plot data and send via shared memory. See the sirilpy.plot submodule
        documentation for how to configure a PlotData object for use with SirilInterface.xy_plot()

        Args:
            plot_metadata: PlotMetadata object containing plot configuration
            display: bool indicating whether to display the plot on screen (defaults to True)
            save: bool indicating whether to save to the file specified in PlotData.savename
                  (defaults to False)

        Raises:
            DataError: if invalid xy_plot data is received via shared memory,
            SirilError: If an error occurs.
        """
        try:
            serialized_data, total_bytes = PlotData.serialize(plot_data)

            # Create shared memory using our wrapper
            try:
                # Request Siril to provide a big enough shared memory allocation
                shm_info = self._request_shm(total_bytes)

                # Map the shared memory
                try:
                    shm = self._map_shared_memory(
                        shm_info.shm_name.decode('utf-8'),
                        shm_info.size
                    )
                except (OSError, ValueError) as e:
                    raise SharedMemoryError(_("Failed to map shared memory: {}").format(e)) from e

            except Exception as e:
                raise SharedMemoryError(_("Failed to create shared memory: {}").format(e)) from e

            # Copy serialized data to shared memory
            try:
                buffer = memoryview(shm.buf).cast('B')
                buffer[:total_bytes] = serialized_data
                del buffer
            except Exception as e:
                raise DataError(_("Error in xy_plot(): {}").format(e)) from e

            # Pack the plot info structure
            info = struct.pack(
                '!IIIIQ256s',
                1 if save else 0,  # width (repurposed as a "save" flag)
                1 if display else 0,  # height (repurposed as a "display" flag)
                0,  # reserved/unused
                0,  # reserved/unused
                total_bytes,
                shm_info.shm_name
            )

            if not self._execute_command(_Command.SIRIL_PLOT, info):
                raise SirilError(_("Failed to send pixel data command"))

            return True

        except SirilError:
            raise
        except Exception as e:
            raise SirilError(f"Error in xy_plot(): {e}") from e
        finally:
            # Ensure shared memory is closed and unlinked
            if 'shm' in locals() and shm is not None:
                try:
                    shm.close()
                except Exception as e:
                    pass

    def clear_image_bgsamples(self):
        """
        Clears all background sample points from the image.
        """
        self._execute_command(_Command.CLEAR_BGSAMPLES, None)

    def set_image_bgsamples(self, points: Union[List[Tuple[float, float]], List[BGSample]], show_samples: bool = False, recalculate = True):
        """
        Serialize a set of background sample points and send via shared memory.
        Points can be provided either as:
        - List of (x,y) Tuples: BGSamples are created with these positions and Siril
        will automatically compute the statistics.
        - List of BGSample objects: The complete sample data is sent to Siril.
        By default Siril will recalculate statistics for the sample points
        on receipt, but this can be overridden with the argument recalculate=False

        Args:
            points: List of sample points, either as (x,y) tuples or BGSample objects
            show_samples: Whether to show the sample points in Siril
            recalculate: Whether to recalculate the sample points once set. This only
                         applies if the sample points are provided as a List of
                         BGSamples, in which case it defaults to True. If the sample
                         points are provided as a List of (x,y) Tuples then the
                         parameter has no effect. Setting recalculate=False is usually
                         a bad idea but the option is provided to support potential
                         advanced uses where the values are adjusted in python code to
                         manipulate the background fit.

        Returns: True if the command succeeded, otherwise False

        Raises:
            NoImageError: if no image is loaded in Siril,
            ValueError: if samples do not have valid positions,
            SirilError: if there was a Siril error in handling the command.
        """
        try:
            if not self.is_image_loaded():
                raise NoImageError(_("Error in set_image_bgsamples(): no image loaded"))

            # Convert tuples to BGSamples if needed
            recalc = True
            if points and not isinstance(points[0], BGSample):
                samples = [BGSample(position=pos) for pos in points]
            else:
                samples = points
                recalc = bool(recalculate)

            # Format string to match C struct layout:
            # 3 doubles for median array
            # 1 double for mean
            # 2 doubles for min,max
            # 1 unsigned long long for size_t
            # 2 doubles for position
            # 1 unsigned int for gboolean
            format_string = '3dd2dQ2dI' * len(samples)

            # Flatten all BGSample data
            flat_values = []
            for sample in samples:
                if sample.position is None:
                    raise ValueError("All BGSamples must have valid position values")
                flat_values.extend([
                    *sample.median,      # 3 doubles
                    sample.mean,         # 1 double
                    sample.min,          # 1 double
                    sample.max,          # 1 double
                    sample.size,         # unsigned long long
                    *sample.position,    # 2 doubles
                    1 if sample.valid else 0  # unsigned int (gboolean)
                ])

            # Pack all the values
            serialized_data = struct.pack(format_string, *flat_values)
            total_bytes = len(serialized_data)

            # Create shared memory using our wrapper
            try:
                # Request Siril to provide a big enough shared memory allocation
                shm_info = self._request_shm(total_bytes)
                # Map the shared memory
                try:
                    shm = self._map_shared_memory(
                        shm_info.shm_name.decode('utf-8'),
                        shm_info.size
                    )
                except (OSError, ValueError) as e:
                    raise SharedMemoryError(_("Failed to map shared memory: {}").format(e)) from e
            except Exception as e:
                raise SharedMemoryError(_("Failed to create shared memory: {}").format(e)) from e

            # Copy serialized data to shared memory
            try:
                buffer = memoryview(shm.buf).cast('B')
                buffer[:total_bytes] = serialized_data
                del buffer
            except Exception as e:
                raise SharedMemoryError(f"Failed to copy data to shared memory: {e}") from e

            # Pack the plot info structure
            info = struct.pack(
                '!IIIIQ256s',
                0,  # width (not used for plots)
                0,  # height (not used for plots)
                1 if recalc else 0,  # recalculate flag
                1 if show_samples else 0,  # show_samples flag
                total_bytes,
                shm_info.shm_name
            )

            if not self._execute_command(_Command.SET_BGSAMPLES, info):
                raise SirilError(_("Failed to send BG sample command"))
            return True

        except (ValueError, SirilError):
            raise
        except Exception as e:
            raise SirilError(f"Error in set_image_bgsamples(): {e}") from e
        finally:
            # Ensure shared memory is closed and unlinked
            if 'shm' in locals() and shm is not None:
                try:
                    shm.close()
                except Exception as e:
                    pass

    def set_image_pixeldata(self, image_data: np.ndarray) -> bool:
        """
        Send image data to Siril using shared memory.

        Args:
            image_data: numpy.ndarray containing the image data.
                        Must be 2D (single channel) or 3D (multi-channel) array
                        with dtype either np.float32 or np.uint16.

        Raises:
            NoImageError: if no image is loaded in Siril,
            ValueError: if the input array is invalid,
            SirilError: if there was an error in handling the command.
        """

        shm = None
        shm_info = None
        try:
            # Validate input array
            if not isinstance(image_data, np.ndarray):
                raise ValueError(_("Image data must be a numpy array"))

            if image_data.ndim not in (2, 3):
                raise ValueError(_("Image must be 2D or 3D array"))

            if image_data.dtype not in (np.float32, np.uint16):
                raise ValueError(_("Image data must be float32 or uint16"))

            # Check there is an image loaded in Siril
            if not self.is_image_loaded():
                raise NoImageError(_("Error in set_image_pixeldata(): no image loaded"))

            # Get dimensions
            if image_data.ndim == 2:
                height, width = image_data.shape
                channels = 1
            else:
                channels, height, width = image_data.shape

            if channels > 3:
                raise ValueError(_("Image cannot have more than 3 channels"))

            if any(dim <= 0 for dim in (width, height)):
                raise ValueError(_("Invalid image dimensions: {}x{}").format(width, height))

            # Calculate total size
            element_size = 4 if image_data.dtype == np.float32 else 2
            total_bytes = width * height * channels * element_size

            # Create shared memory using our wrapper
            try:
                # Request Siril to provide a big enough shared memory allocation
                shm_info = self._request_shm(total_bytes)

                # Map the shared memory
                try:
                    shm = self._map_shared_memory(
                        shm_info.shm_name.decode('utf-8'),
                        shm_info.size
                    )
                except (OSError, ValueError) as e:
                    raise SharedMemoryError(_("Failed to map shared memory: {}").format(e)) from e

            except Exception as e:
                raise SharedMemoryError(_("Failed to create shared memory: {}").format(e)) from e

            # Copy data to shared memory
            try:
                buffer = memoryview(shm.buf).cast('B')
                shared_array = np.frombuffer(buffer, dtype=image_data.dtype).reshape(image_data.shape)
                np.copyto(shared_array, image_data)
                # Delete transient objects used to structure copy
                del buffer
                del shared_array
            except Exception as e:
                raise SharedMemoryError(_("Failed to copy data to shared memory: {}").format(e)) from e

            # Pack the image info structure
            info = struct.pack(
                '!IIIIQ256s',
                width,
                height,
                channels,
                1 if image_data.dtype == np.float32 else 0,
                total_bytes,
                shm_info.shm_name
            )
            # Send command using the existing _execute_command method
            if not self._execute_command(_Command.SET_PIXELDATA, info):
                raise SirilError(_("Failed to send pixel data command"))

            return

        except (SirilError, ValueError):
            raise
        except Exception as e:
            raise SirilError(f"Error in set_image_pixeldata(): {e}") from e

        finally:
            if shm is not None:
                try:
                    shm.close()
                    self._execute_command(_Command.RELEASE_SHM, shm_info)
                except Exception as e:
                    pass

    def set_seq_frame_pixeldata(self, index: int, image_data: np.ndarray, prefix: str) -> bool:
        """
        Send sequence frame image data to Siril using shared memory. Note that this
        method only works with sequences of FITS images: it does **not** work with
        FITSEQ, SER or AVI single-file sequences. The image_lock() context manager
        is not required in order to use this method.

        Args:
            index: integer specifying which frame to set the pixeldata for. This
                uses a 0-based indexing scheme, i.e. the first frame is frame
                number 0, not frame numer 1.
            image_data: numpy.ndarray containing the image data.
                Must be 2D (single channel) or 3D (multi-channel) array
                with dtype either np.float32 or np.uint16.
            prefix: String prefix to use when saving the file to make a
                new sequence. Note that saving sequence frames with a new prefix
                does not by itself create a new sequence: once all the frames have
                been saved with the new sequence prefix,
                ``sirilpy.SirilInterface.create_new_seq()`` must be called to create
                the actual sequence file. Note that while it is permitted to pass
                prefix=None, this will overwrite the existing sequence and is not
                typically what is wanted, therefore the parameter is not optional
                and must be passed explicitly.

        Raises:
            NoSequenceError: if no sequence is loaded in Siril,
            ValueError: if the input array is invalid,
            SirilError: if there was an error in handling the command.
        """

        shm = None
        shm_info = None  # Initialize at the start so it's available in finally

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
                channels, height, width = image_data.shape

            if channels > 3:
                raise ValueError(_("Image cannot have more than 3 channels"))

            if any(dim <= 0 for dim in (width, height)):
                raise ValueError(_("Invalid image dimensions: {}x{}").format(width, height))

            if not self.is_sequence_loaded():
                raise NoSequenceError(_("Error in set_seq_frame_pixeldata(): no sequence loaded"))

            # Calculate total size
            element_size = 4 if image_data.dtype == np.float32 else 2
            total_bytes = width * height * channels * element_size

            # Create shared memory using our wrapper
            try:
                # Request Siril to provide a big enough shared memory allocation
                shm_info = self._request_shm(total_bytes)

                # Map the shared memory
                try:
                    shm = self._map_shared_memory(
                        shm_info.shm_name.decode('utf-8'),
                        shm_info.size
                    )
                except (OSError, ValueError) as e:
                    raise SharedMemoryError(_("Failed to map shared memory: {}").format(e)) from e

            except Exception as e:
                raise SharedMemoryError(_("Failed to create shared memory: {}").format(e)) from e

            # Copy data to shared memory
            try:
                buffer = memoryview(shm.buf).cast('B')
                shared_array = np.frombuffer(buffer, dtype=image_data.dtype).reshape(image_data.shape)
                np.copyto(shared_array, image_data)
                # Delete transient objects used to structure copy
                del buffer
                del shared_array
            except Exception as e:
                raise SharedMemoryError(_("Failed to copy data to shared memory: {}").format(e)) from e

            # Pack the image info structure
            info = struct.pack(
                '!IIIIQ256s',
                width,
                height,
                channels,
                1 if image_data.dtype == np.float32 else 0,
                total_bytes,
                shm_info.shm_name
            )

            # Prepare prefix bytes (256 bytes max, null-terminated)
            if prefix is None:
                prefix_bytes = b'\x00' * 256
            else:
                # Encode prefix and ensure it fits in 256 bytes (including null terminator)
                prefix_encoded = prefix.encode('utf-8')
                if len(prefix_encoded) > 255:
                    raise ValueError(_("Prefix too long (max 255 bytes)"))
                prefix_bytes = prefix_encoded + b'\x00' * (256 - len(prefix_encoded))

            # Create payload
            # We don't range check index here as it is done more efficiently in the C code
            index_bytes = struct.pack('!i', index)
            payload = index_bytes + info + prefix_bytes

            # Send command using the existing _execute_command method
            if not self._execute_command(_Command.SET_SEQ_FRAME_PIXELDATA, payload):
                raise SirilError(_("Failed to send set_seq_frame_pixeldata command"))

            return True

        except (ValueError, SirilError):
            raise
        except Exception as e:
            raise SirilError(f"Error in set_seq_frame_pixeldata(): {e}") from e

        finally:
            # Close the shared memory handle first
            if shm is not None:
                try:
                    shm.close()
                except Exception:
                    pass

            # Then send the release command to the C side
            if shm_info is not None:
                try:
                    shm_name_str = shm_info.shm_name.decode('utf-8').rstrip('\x00')
                    if not self._execute_command(_Command.RELEASE_SHM, shm_info):
                        print(f"Warning: Failed to release shared memory in set_seq_frame_pixeldata")
                except Exception as e:
                    print(f"Exception during shm cleanup in set_seq_frame_pixeldata: {e}")

    def get_image_iccprofile(self) -> Optional[bytes]:
        """
        Retrieve the ICC profile of the current Siril image using shared memory. Requires
        a single image to be loaded.

        Args:
        none.

        Returns:
            bytes: The image ICC profile as a byte array, or None if the current
            image has no ICC profile.

        Raises:
            NoImageError: If no image is currently loaded,
            SirilError: If any other error occurs.
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
                        raise NoImageError(_("Error in get_image_iccprofile(): no image is currently loaded in Siril"))
                    raise SirilConnectionError(_("Server error: {}").format(error_msg))
                raise SharedMemoryError(_("Failed to initiate shared memory transfer: Empty response"))

            if status == _Status.NONE:
                return None

            if not response:
                raise SharedMemoryError(_("Failed to initiate shared memory transfer: No data received"))

            if len(response) < 280: # Not a correct SharedMemoryInfo payload
                return None
            try:
                # Parse the shared memory information
                shm_info = _SharedMemoryInfo.from_buffer_copy(response)
            except (AttributeError, BufferError, ValueError) as e:
                raise SharedMemoryError(_("Invalid shared memory information received: {}").format(e)) from e

            # Map the shared memory
            try:
                shm = self._map_shared_memory(
                    shm_info.shm_name.decode('utf-8'),
                    shm_info.size
                )
            except (OSError, ValueError) as e:
                raise SharedMemoryError(_("Failed to map shared memory: {}").format(e)) from e

            try:
                buffer = bytearray(shm.buf)[:shm_info.size]
                result = bytes(buffer)
            except (BufferError, ValueError, TypeError) as e:
                raise SirilError(_("Failed to create bytes from shared memory: {}").format(e)) from e

            return result

        except SirilError:
            # Re-raise NoImageError without wrapping
            raise
        except Exception as e:
            # Wrap all other exceptions with context
            raise SirilError(_("Error in get_image_iccprofile(): {}").format(e)) from e

        finally:
            if shm is not None:
                try:
                    shm.close()  # First close the memory mapping as we have finished with it
                    # (We don't unlink it as C wll do that)

                    # Signal that Python is done with the shared memory and wait for C to finish
                    if not self._execute_command(_Command.RELEASE_SHM, shm_info):
                        raise SharedMemoryError(_("Failed to cleanup shared memory"))

                except Exception:
                    pass

    def get_image_fits_header(self, return_as = 'str') -> Union[str, dict, None]:
        """
        Retrieve the full FITS header of the current image loaded in Siril.
        Requires a single image to be loaded.

        Args:
            return_as: Optional string specifying the format of the returned header.
             Can be 'str' for a string or 'dict' for a dictionary.

        Returns:
            str: The image FITS header as a string, or None if there is no header.
            dict: The image FITS header as a dictionary, or None if there is no header.
            None: If the header is empty or not available.

        Raises:
            NoImageError: If no image is currently loaded,
            SirilError: For other errors during data retrieval,
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
                        raise NoImageError(_("Error in get_image_fits_header(): no image is currently loaded in Siril"))
                    raise SirilError(_("Server error: {}").format(error_msg))
                raise SharedMemoryError(_("Failed to initiate shared memory transfer: Empty response"))

            if not response:
                raise SharedMemoryError(_("Failed to initiate shared memory transfer: No data received"))

            if status == _Status.NONE:
                return None

            if len(response) < 25: # No payload
                return None

            try:
                # Parse the shared memory information
                shm_info = _SharedMemoryInfo.from_buffer_copy(response)
            except (AttributeError, BufferError, ValueError) as e:
                raise SharedMemoryError(_("Invalid shared memory information received: {}").format(e)) from e

            # Map the shared memory
            try:
                shm = self._map_shared_memory(
                    shm_info.shm_name.decode('utf-8'),
                    shm_info.size
                )
            except (OSError, ValueError) as e:
                raise SharedMemoryError(_("Failed to map shared memory: {}").format(e)) from e

            try:
                # Read entire buffer at once
                buffer = bytearray(shm.buf)[:shm_info.size]
                result = buffer.decode('utf-8', errors='ignore')
            except (BufferError, ValueError, TypeError) as e:
                raise SirilError(_("Failed to create string from shared memory: {}").format(e)) from e

            if return_as == 'dict':
                return parse_fits_header(result)
            if return_as == 'str':
                return result
            raise ValueError(_(f"Invalid return_as value: {return_as}"))

        except SirilError:
            # Re-raise NoImageError and other SirilErrors without wrapping
            raise
        except Exception as e:
            # Wrap all other exceptions with context
            raise SirilError(_("Error in get_image_fits_header(): {}").format(e)) from e
        finally:
            if shm is not None:
                try:
                    shm.close()  # First close the memory mapping as we have finished with it
                    # (We don't unlink it as C wll do that)

                    # Signal that Python is done with the shared memory and wait for C to finish
                    if not self._execute_command(_Command.RELEASE_SHM, shm_info):
                        raise SharedMemoryError(_("Failed to cleanup shared memory"))

                except Exception:
                    pass

    def get_image_unknown_keys(self) -> Optional[str]:
        """
        Retrieve the unknown key in a FITS header of the current loaded Siril
        image using shared memory. Requires a single image to be loaded.

        Args:
            none.

        Returns:
            bytes: The unknown keys as a string, or None if there are no unknown keys.

        Raises:
            NoImageError: If no image is currently loaded,
            SirilError: For other errors during data retrieval.
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
                        raise NoImageError(_("Error in get_image_unknown_keys(): no image is currently loaded in Siril"))
                    raise SirilError(_("Server error: {}").format(error_msg))
                raise SharedMemoryError(_("Failed to initiate shared memory transfer: Empty response"))

            if status == _Status.NONE:
                return None

            if status == _Status.NONE:
                return None

            if not response:
                raise SharedMemoryError(_("Failed to initiate shared memory transfer: No data received"))
            if len(response) < 25: # No payload
                return None

            try:
                # Parse the shared memory information
                shm_info = _SharedMemoryInfo.from_buffer_copy(response)
            except (AttributeError, BufferError, ValueError) as e:
                raise SharedMemoryError(_("Invalid shared memory information received: {}").format(e)) from e

            # Map the shared memory
            try:
                shm = self._map_shared_memory(
                    shm_info.shm_name.decode('utf-8'),
                    shm_info.size
                )
            except (OSError, ValueError) as e:
                raise SharedMemoryError(_("Failed to map shared memory: {}").format(e)) from e

            try:
                # Read entire buffer at once using memoryview
                buffer = bytearray(shm.buf)[:shm_info.size]
                result = buffer.decode('utf-8', errors='ignore')
            except (BufferError, ValueError, TypeError) as e:
                raise SirilError(_("Failed to create string from shared memory: {}").format(e)) from e

            return result

        except SirilError:
            # Re-raise NoImageError without wrapping
            raise
        except Exception as e:
            # Wrap all other exceptions with context
            raise SirilError(_("Error in get_image_unknown_keys(): {}").format(e)) from e
        finally:
            if shm is not None:
                try:
                    shm.close()  # First close the memory mapping as we have finished with it
                    # (We don't unlink it as C wll do that)

                    # Signal that Python is done with the shared memory and wait for C to finish
                    if not self._execute_command(_Command.RELEASE_SHM, shm_info):
                        raise SharedMemoryError(_("Failed to cleanup shared memory"))

                except Exception:
                    pass

    def get_image_history(self) -> Optional[list[str]]:
        """
        Retrieve history entries in the FITS header of the current loaded
        Siril image using shared memory. Requires a single image to be loaded.

        Args:
            none.

        Returns:
            list: The HISTORY entries in the FITS header as a list of strings, or
                  None if there are no HISTORY keywords.

        Raises:
            NoImageError: If no image is currently loaded,
            SirilError: For other errors during data retrieval.
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
                        raise NoImageError(_("Error in get_image_history(): no image is currently loaded in Siril"))
                    raise SirilError(_("Server error: {}").format(error_msg))
                raise SharedMemoryError(_("Failed to initiate shared memory transfer: Empty response"))

            if status == _Status.NONE:
                return None

            if not response:
                raise SharedMemoryError(_("Failed to initiate shared memory transfer: No data received"))

            if len(response) < 29:  # No payload
                return None

            try:
                # Parse the shared memory information
                shm_info = _SharedMemoryInfo.from_buffer_copy(response)
            except (AttributeError, BufferError, ValueError) as e:
                raise SharedMemoryError(_("Invalid shared memory information received: {}").format(e)) from e

            # Map the shared memory
            try:
                shm = self._map_shared_memory(
                    shm_info.shm_name.decode('utf-8'),
                    shm_info.size
                )
            except (OSError, ValueError) as e:
                raise SharedMemoryError(_("Failed to map shared memory: {}").format(e)) from e

            try:
                # Read entire buffer at once using memoryview
                buffer = bytearray(shm.buf)[:shm_info.size]
                string_data = buffer.decode('utf-8', errors='ignore')
                string_list = string_data.split('\x00')
            except (BufferError, ValueError, TypeError) as e:
                raise SirilError(_("Failed to create string from shared memory: {}").format(e)) from e

            return [s for s in string_list if s]

        except SirilError:
            # Raise SirilError without wrapping
            raise
        except Exception as e:
            # Wrap all other exceptions with context
            raise SirilError(_("Error in get_image_history(): {}").format(e)) from e
        finally:
            if shm is not None:
                try:
                    shm.close()  # First close the memory mapping as we have finished with it
                    # (We don't unlink it as C wll do that)

                    # Signal that Python is done with the shared memory and wait for C to finish
                    if not self._execute_command(_Command.RELEASE_SHM, shm_info):
                        raise SharedMemoryError(_("Failed to cleanup shared memory"))

                except Exception:
                    pass

    def get_siril_wd(self) -> str:
        """
        Request the current working directory from Siril.

        Returns:
            The current working directory as a string.

        Raises:
            DataError: if no response was obtained,
            SirilError: for all other errors.
        """

        response = self._request_data(_Command.GET_WORKING_DIRECTORY)

        if response is None:
            raise DataError(_("Error in get_siril_wd(): no response"))

        try:
            # Assuming the response is a null-terminated UTF-8 encoded string
            wd = response.decode('utf-8').rstrip('\x00')
            return wd
        except Exception as e:
            raise SirilError(_("Error in get_wd(): {}").format(e)) from e

    def get_siril_configdir(self) -> str:
        """
        Request the user config directory used by Siril.

        Returns:
            The user config directory as a string.

        Raises:
            DataError: if no response is received,
            SirilError: for all other errors.
        """

        response = self._request_data(_Command.GET_USERCONFIGDIR)

        if response is None:
            raise DataError(_("Error in get_siril_configdir(): no response"))

        try:
            # Assuming the response is a null-terminated UTF-8 encoded string
            wd = response.decode('utf-8').rstrip('\x00')
            return wd
        except Exception as e:
            raise SirilError(_("Error in get_siril_configdir(): {}").format(e)) from e

    def get_siril_userdatadir(self) -> str:
        """
        Request the user data directory used by Siril.

        Returns:
            The user data directory as a string.

        Raises:
            DataError: if no response is received,
            SirilError: for all other errors.
        """

        response = self._request_data(_Command.GET_USERDATADIR)

        if response is None:
            raise DataError(_("Error in get_siril_userdatadir(): no response"))

        try:
            # Assuming the response is a null-terminated UTF-8 encoded string
            path = response.decode('utf-8').rstrip('\x00')
            return path
        except Exception as e:
            raise SirilError(_("Error in get_siril_userdatadir(): {}").format(e)) from e

    def get_siril_systemdatadir(self) -> Optional[str]:
        """
        Request the system data directory used by Siril.

        Returns:
            The system data directory as a string.

        Raises:
            DataError: if no response is received,
            SirilError: for all other errors.
        """

        response = self._request_data(_Command.GET_SYSTEMDATADIR)

        if response is None:
            raise DataError(_("Error in get_siril_systemdatadir(): no response"))

        try:
            # Assuming the response is a null-terminated UTF-8 encoded string
            path = response.decode('utf-8').rstrip('\x00')
            return path
        except UnicodeDecodeError as e:
            raise SirilError(_("Error in get_siril_systemdatadir(): {}").format(e)) from e

    def is_image_loaded(self) -> bool:
        """
        Check if a single image is loaded in Siril.

        Returns:
            bool: True if a single image is loaded, False if a single image is
            not loaded

        Raises:
            DataError: if no response is received,
            SirilError: for all other errors.
        """
        response = self._request_data(_Command.GET_IS_IMAGE_LOADED)

        if response is None:
            raise DataError(_("Error in is_image_loaded(): no response received from Siril"))

        try:
            image_loaded = struct.unpack('!i', response)[0] != 0
            return bool(image_loaded)
        except Exception as e:
            raise SirilError("Error in is_image_loaded()") from e

    def is_sequence_loaded(self) -> bool:
        """
        Check if a sequence is loaded in Siril.

        Returns:
            bool: True if a sequence is loaded, False if a sequence is
            not loaded

        Raises:
            DataError: if no response is received,
            SirilError: for all other errors.
        """
        response = self._request_data(_Command.GET_IS_SEQUENCE_LOADED)

        if response is None:
            raise DataError(_("Error in is_sequence_loaded(): no response received from Siril"))

        try:
            sequence_loaded = struct.unpack('!i', response)[0] != 0
            return bool(sequence_loaded)
        except Exception as e:
            raise SirilError(_("Error in is_sequence_loaded()")) from e

    def get_image_filename(self) -> Optional[str]:
        """
        Request the filename of the loaded image from Siril. Requires a
        single image to be loaded.

        Returns:
            The filename as a string.

        Raises:
            NoImageError: if no image is loaded,
            SirilError: if a decoding error occurs.
        """

        response = self._request_data(_Command.GET_FILENAME)

        if response is None:
            return None

        try:
            # Assuming the response is a null-terminated UTF-8 encoded string
            wd = response.decode('utf-8').rstrip('\x00')
            return wd
        except Exception as e:
            raise SirilError(f"Error in get_image_filename(): {e}") from e

    def get_seq_frame_filename(self, frame: int) -> Optional[str]:
        """
        Request the filename of the specified frame of the loaded sequence from Siril.
        Requires a sequence to be loaded.

        Args:
            frame (int): Specifies the frame index. This
                uses a 0-based indexing scheme, i.e. the first frame is frame
                number 0, not frame numer 1.

        Returns:
            The filename as a string.

        Raises:
            NoSequenceError: if no sequence is loaded in Siril,
            SirilError: if a decoding error occurs.
        """

        if not self.is_sequence_loaded():
            raise NoSequenceError(_("Error in get_seq_frame_filename(): no sequence is loaded"))

        # Convert frame number to network byte order bytes
        frame_payload = struct.pack('!I', frame)  # '!I' for network byte order uint32_t

        response = self._request_data(_Command.GET_SEQ_FRAME_FILENAME, payload=frame_payload)

        if response is None:
            return None

        try:
            # Assuming the response is a null-terminated UTF-8 encoded string
            filename = response.decode('utf-8').rstrip('\x00')
            return filename
        except Exception as e:
            raise SirilError(f"Error in get_seq_frame_filename(): {e}") from e

    def get_image_stats(self, channel: int) -> Optional[ImageStats]:
        """
        Request image statistics from Siril for a specific channel. Requires a
        single image to be loaded.

        Args:
            channel: Integer specifying which channel to get statistics
                     for (typically 0, 1, or 2)

        Returns:
            ImageStats object containing the statistics, or None if no stats are available for the selected channel

        Raises:
            NoImageError: if no image is loaded,
            SirilError: if an error occurs.
        """

        # Convert channel number to network byte order bytes
        channel_payload = struct.pack('!I', channel)  # '!I' for network byte order uint32_t

        # Request data with the channel number as payload
        response = self._request_data(_Command.GET_IMAGE_STATS, payload=channel_payload)

        if response is None:
            return None

        try:
            return ImageStats.deserialize(response)

        except struct.error as e:
            raise SirilError(f"Error in get_image_stats(): {e}") from e

    def get_seq_regdata(self, frame: int, channel: int) -> Optional[RegData]:
        """
        Request sequence frame registration data from Siril.

        Args:
            frame: Integer specifying which frame in the sequence to get registration
                data for (between 0 and Sequence.number). This
                uses a 0-based indexing scheme, i.e. the first frame is frame
                number 0, not frame numer 1.
            channel: Integer specifying which channel to get registration data
                     for (typically 0, 1, or 2)

        Returns:
            RegData object containing the registration data, or None if
            no registration data is available for the specified frame
            and channel

        Raises:
            NoSequenceError: if no sequence is loaded in Siril,
            SirilError: if a decoding error occurs.
        """

        if not self.is_sequence_loaded():
            raise NoSequenceError(_("Error in get_seq_regdata(): no sequence is loaded"))

        data_payload = struct.pack('!II', frame, channel)  # '!I' for network byte order uint32_t

        # Request data with the channel number as payload
        response = self._request_data(_Command.GET_SEQ_REGDATA, payload=data_payload)

        if response is None:
            return None

        try:
            return RegData.deserialize(response)

        except struct.error as e:
            raise SirilError(f"Error in get_seq_regdata(): {e}") from e

    def get_seq_stats(self, frame: int, channel: int) -> Optional[ImageStats]:
        """
        Request sequence frame statistics from Siril.

        Args:
            frame: Integer specifying which frame in the sequence to get statistics
                data for (between 0 and Sequence.number). This
                uses a 0-based indexing scheme, i.e. the first frame is frame
                number 0, not frame numer 1.
            channel: Integer specifying which channel to get statistics
                     for (typically 0, 1, or 2)

        Returns:
            ImageStats object containing the statistics, or None if an error occurred

        Raises:
            NoSequenceError: if no sequence is loaded in Siril,
            SirilError: if a decoding error occurs.
        """

        if not self.is_sequence_loaded():
            raise NoSequenceError(_("Error in get_seq_stats(): no sequence is loaded"))

        data_payload = struct.pack('!II', frame, channel)  # '!I' for network byte order uint32_t

        # Request data with the channel number as payload
        response = self._request_data(_Command.GET_SEQ_STATS, payload=data_payload)

        if response is None:
            return None

        try:
            return ImageStats.deserialize(response)

        except struct.error as e:
            raise SirilError(f"Error in get_seq_stats(): {e}") from e

    def get_seq_imgdata(self, frame: int) -> Optional[ImgData]:
        """
        Request sequence frame metadata from Siril.

        Args:
            frame: Integer specifying which frame in the sequence to get image
                metadata for (between 0 and Sequence.number). This
                uses a 0-based indexing scheme, i.e. the first frame is frame
                number 0, not frame numer 1.

        Returns:
            ImgData object containing the frame metadata, or None if an error occurred

        Raises:
            NoSequenceError: if no sequence is loaded in Siril,
            SirilError: if a decoding error occurs.
        """

        if not self.is_sequence_loaded():
            raise NoSequenceError(_("Error in get_seq_imgdata(): no sequence is loaded"))

        data_payload = struct.pack('!I', frame)  # '!I' for network byte order uint32_t

        # Request data with the channel number as payload
        response = self._request_data(_Command.GET_SEQ_IMGDATA, payload=data_payload)
        if response is None:
            return None

        try:
            return ImgData.deserialize(response)
        except struct.error as e:
            raise SirilError(f"Error in get_seq_imgdata(): {e}") from e

    def get_seq_distodata(self, channel: int) -> Optional[DistoData]:
        """
        Request sequence distortion data from Siril

        channel: Integer specifying which channel to get registration data
        for (typically 0, 1, or 2)

        Returns:
            DistoData object containing the channel distortion parameters, or None if an error occurred

        Raises:
            NoSequenceError: if no sequence is loaded in Siril,
            SirilError: if a decoding error occurs.
        """

        if not self.is_sequence_loaded():
            raise NoSequenceError(_("Error in get_seq_distodata(): no sequence is loaded"))

        # Request data with the channel number as payload
        data_payload = struct.pack('!I', channel)
        response = self._request_data(_Command.GET_SEQ_DISTODATA, payload=data_payload)
        if response is None:
            return None
        try:
            format_string = '!q2d'

            fixed_length = struct.calcsize(format_string)
            values = struct.unpack(format_string, response[:fixed_length])
            # Extract remaining bytes for the null-terminated string
            if len(response) > fixed_length:
                remaining_data = response[fixed_length:]
                distofilename_string = remaining_data.decode('utf-8').rstrip('\x00')
            else:
                distofilename_string = ''

            return DistoData (
                index = values[0],
                velocity = (values[1], values[2]),
                filename = distofilename_string
            )
        except struct.error as e:
            raise SirilError(f"Error in get_seq_distodata(): {e}") from e

    def get_seq(self) -> Optional[Sequence]:
        """
        Request metadata for the current sequence loaded in Siril.

        Returns:
            Sequence object containing the current sequence metadata, or None
            if an error occurred

        Raises:
            NoSequenceError: if no sequence is loaded in Siril,
            SirilError: if a decoding error occurs.
        """

        if not self.is_sequence_loaded():
            raise NoSequenceError(_("Error in get_seq(): no sequence is loaded"))

        # Request data with the channel number as payload
        response = self._request_data(_Command.GET_SEQ)

        if response is None:
            return None

        try:
            format_string = '!4q3Q4qdQqQq'
            fixed_length = struct.calcsize(format_string)

            values = struct.unpack(format_string, response[:fixed_length])
            # Extract remaining bytes for the null-terminated string
            remaining_data = response[fixed_length:]
            seqname_string = remaining_data.decode('utf-8').rstrip('\x00')

            number = values[0]
            nb_layers = values[3]

            imgparam_list = [self.get_seq_imgdata(frame)
                             for frame in range(number)]

            regdata_list = [[self.get_seq_regdata(frame, channel)
                            for frame in range(number)]
                            for channel in range(nb_layers)]

            stats_list = [[self.get_seq_stats(frame, channel)
                          for frame in range(number)]
                          for channel in range(nb_layers)]

            disto_list = [self.get_seq_distodata(channel)
                             for channel in range(nb_layers)]

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
                distoparam = disto_list,
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
            raise SirilError(f"Error in get_seq(): {e}") from e

    def get_image_keywords(self) -> Optional[FKeywords]:
        """
        Request FITS keywords data from Siril as a FKeywords object. Requires
        a single image to be loaded.

        Returns:
            FKeywords object containing the FITS keywords, or None if an error occurred

        Raises:
            SirilError: if a decoding error occurs.
        """

        response = self._request_data(_Command.GET_KEYWORDS)
        if response is None:
            return None

        try:
            return FKeywords.deserialize(response)

        except struct.error as e:
            raise SirilError(f"Error in get_image_keywords(): {e}") from e

    def get_image(self, with_pixels: Optional[bool] = True, preview: Optional[bool] = False) -> Optional[FFit]:
        """
        Request a copy of the current image open in Siril. Requires a single
        image to be loaded.

        Args:
            with_pixels: optional bool specifying whether to get pixel data as a
                         NumPy array, or only the image metadata. Defaults to True
            preview: optional bool specifying whether to get pixeldata as a preview
                     (i.e. 8-bit autostretched data) or as real image data. Defaults
                     to False (i.e. real image data)

        Returns:
            FFit object containing the image metadata and (optionally)
            pixel data, or None if an error occurred

        Raises:
            NoImageError: if no image is loaded in Siril
            SirilError: if a decoding error occurs
        """

        # Request data with the channel number as payload
        response = self._request_data(_Command.GET_IMAGE)
        if response is None:
            raise NoImageError(_("Error in get_image(): no image currently loaded in Siril"))

        try:
            # Build format string for struct unpacking
            # Network byte order for all values
            format_parts = [
                'q',  # rx padded to 64bit
                'q',  # ry padded to 64bit
                'q',  # naxes[2] padded to 64bit
                'q',  # bitpix padded to 64bit
                'q',  # orig_bitpix padded to 64bit
                'Q',  # gboolean checksum padded to 64bit
                'd',  # mini
                'd',  # maxi
                'd',  # neg_ratio padded to 64bit
                'Q',  # gboolean top_down (padded to uint64_t)
                'Q',  # gboolean focalkey (padded to uint64_t)
                'Q',  # gboolean pixelkey (padded to uint64_t)
                'Q',  # gboolean color_managed (padded to uint64_t)
            ]

            format_string = '!' + ''.join(format_parts)

            # Verify data size
            expected_size = struct.calcsize(format_string)
            if len(response) != expected_size:
                raise ValueError(f"Received image data size {len(response)} doesn't match expected size {expected_size}")

            # Unpack the binary data
            values = struct.unpack(format_string, response)

            # Create FFit object
            try:
                img_header = self.get_image_fits_header()
            except Exception:
                img_header = ""

            try:
                img_unknown_keys = self.get_image_unknown_keys()
            except Exception:
                img_unknown_keys = ""

            try:
                img_history = self.get_image_history()
            except Exception:
                img_history = []

            try:
                img_icc_profile = self.get_image_iccprofile()
            except Exception as e:
                print(f"Failed to retrieve ICC profile: {e}")
                img_icc_profile = None

            try:
                img_keywords = self.get_image_keywords()
            except Exception as e:
                print(f"Error: failed to retrieve FITS keywords: {e}")
                img_keywords = None

            img_pixeldata = None
            if with_pixels:
                try:
                    img_pixeldata = self.get_image_pixeldata(preview = preview)
                except Exception as e:
                    print(f"Error: failed to retrieve pixel data: {e}")

            return FFit(
                _naxes = (values[0], values[1], values[2]),
                naxis = 2 if values[2] == 1 else 3,
                bitpix=values[3],
                orig_bitpix=values[4],
                checksum=bool(values[5]),
                mini=values[6],
                maxi=values[7],
                neg_ratio=values[8],
                top_down=bool(values[9]),
                _focalkey=bool(values[10]),
                _pixelkey=bool(values[11]),
                color_managed=bool(values[12]),
                _data = img_pixeldata,
                stats=[
                    self.get_image_stats(0),
                    self.get_image_stats(1) if values[2] > 1 else None,
                    self.get_image_stats(2) if values[2] > 1 else None,
                ],
                keywords = img_keywords,
                _icc_profile = img_icc_profile,
                header = img_header,
                unknown_keys = img_unknown_keys,
                history = img_history
            )

        except Exception as e:
            raise SirilError(f"Error in get_image(): {e}") from e

    def get_seq_frame(self, frame: int, with_pixels: Optional[bool] = True, preview: Optional[bool] = False, linked: Optional[bool] = False) -> Optional[FFit]:
        """
        Request sequence frame as a FFit from Siril. The keywords, statistics, header and
        other metadata are always returned: if an ICC profile is present, this will also
        be populated in the resulting FFit, and optionally the pixel data can also be
        returned.

        Args:
            frame: Integer specifying which frame in the sequence to retrieve data for
                (between 0 and Sequence.number - 1). This uses a 0-based indexing scheme,
                i.e. the first frame is frame number 0, not frame numer 1.
            with_pixels: bool specifying whether or not to return the pixel data for the
                frame (default is True).
            preview: bool specifying whether or not to return the real pixel data or an
                autostretched uint8_t preview version. Only has an effect in
                conjunction with with_pixels = True

        Returns:
            FFit object containing the frame data

        Raises:
            NoSequenceError: if no sequence is loaded in Siril,
            DataError: on receipt of incorrect data,
            SirilError: if an error occurs.
        """
        if not self.is_sequence_loaded():
            raise NoSequenceError(_("Error in get_seq_frame(): no sequence is loaded"))

        shm_pixels = None
        shm_header = None
        shm_icc_profile = None

        # Initialize info variables at the start so they're available in finally
        shm_pixels_info = None
        shm_header_info = None
        shm_icc_info = None

        data_payload = struct.pack('!I???', frame, with_pixels, preview, linked)
        response = self._request_data(_Command.GET_SEQ_IMAGE, data_payload, timeout = None)
        if response is None:
            return None

        try:
            # --- Unpack the core FFit (the fields before the keyword block) ---
            core_format_parts = [
                'q',  # rx
                'q',  # ry
                'q',  # naxes[2]
                'q',  # bitpix
                'q',  # orig_bitpix
                'Q',  # checksum (padded boolean)
                'd',  # mini
                'd',  # maxi
                'd',  # neg_ratio
                'Q',  # top_down
                'Q',  # focalkey
                'Q',  # pixelkey
                'Q',  # color_managed
            ]
            core_format = '!' + ''.join(core_format_parts)
            core_size = struct.calcsize(core_format)

            # Calculate expected size with shared memory info structs
            # Header and ICC profile are always included, pixels only if requested
            shminfo_size = 0
            shminfo_format = '!Qiiii256s'  # One shm_info_t struct
            single_shminfo_size = struct.calcsize(shminfo_format)

            # Always include header and ICC profile shared memory info structs
            shminfo_size = single_shminfo_size * 2  # header + icc_profile
            if with_pixels:
                shminfo_size += single_shminfo_size  # + pixels

            # Make sure response contains core + keywords + shared memory info
            expected_total = core_size + FKeywords.KEYWORDS_SIZE + shminfo_size
            if len(response) != expected_total:
                raise DataError(
                    f"Received image data size {len(response)} doesn't match expected size {expected_total}"
                )

            # Unpack core part
            core_vals = struct.unpack_from(core_format, response, 0)

            # Extract keyword block bytes and deserialize them
            kw_offset = core_size
            kw_bytes = response[kw_offset: kw_offset + FKeywords.KEYWORDS_SIZE]
            fits_keywords = FKeywords.deserialize(kw_bytes)

            # After keywords comes shared-memory info blocks
            shminfo_offset = kw_offset + FKeywords.KEYWORDS_SIZE

            pixeldata = None
            header_data = None
            icc_profile_data = None

            # Determine the order and offsets for shared memory info structs
            current_offset = shminfo_offset

            # Extract pixel data shm_info (only if with_pixels=True)
            if with_pixels:
                (size_be, data_type_be, width_be, height_be, channels_be, shm_name_bytes) = struct.unpack_from(shminfo_format, response, current_offset)

                shm_pixels_info = _SharedMemoryInfo(
                    size = size_be,
                    data_type = data_type_be,
                    width = width_be,
                    height = height_be,
                    channels = channels_be,
                    shm_name = shm_name_bytes
                )

                # Check if pixel data is available (non-zero size and valid shm_name)
                shm_name_str = shm_pixels_info.shm_name.decode('utf-8').rstrip('\x00')
                if shm_pixels_info.size > 0 and shm_name_str:
                    # Validate dims
                    if any(dim <= 0 for dim in (shm_pixels_info.width, shm_pixels_info.height)):
                        raise DataError(_("Invalid image dimensions: {}x{}").format(shm_pixels_info.width, shm_pixels_info.height))

                    if shm_pixels_info.channels <= 0 or shm_pixels_info.channels > 3:
                        raise DataError(_("Invalid number of channels: {}").format(shm_pixels_info.channels))

                    # Map shared memory and build numpy array
                    try:
                        shm_pixels = self._map_shared_memory(shm_name_str, shm_pixels_info.size)
                    except (OSError, ValueError) as e:
                        raise SharedMemoryError(_("Failed to map pixel shared memory: {}").format(e)) from e

                    buffer = bytearray(shm_pixels.buf)[:shm_pixels_info.size]

                    if preview:
                        dtype = np.uint8
                    else:
                        dtype = np.float32 if shm_pixels_info.data_type == 1 else np.uint16

                    try:
                        arr = np.frombuffer(buffer, dtype=dtype)
                    except (BufferError, ValueError, TypeError) as e:
                        raise SharedMemoryError(_("Failed to create array from shared memory: {}").format(e)) from e

                    expected_elems = shm_pixels_info.width * shm_pixels_info.height * shm_pixels_info.channels
                    if arr.size < expected_elems:
                        raise DataError(
                            f"Data size mismatch: got {arr.size} elements, expected {expected_elems} for dimensions {shm_pixels_info.width}x{shm_pixels_info.height}x{shm_pixels_info.channels}"
                        )

                    try:
                        if shm_pixels_info.channels > 1:
                            arr = arr.reshape((shm_pixels_info.channels, shm_pixels_info.height, shm_pixels_info.width))
                        else:
                            arr = arr.reshape((shm_pixels_info.height, shm_pixels_info.width))
                    except ValueError as e:
                        raise DataError(_("Error in get_seq_frame(): Failed to reshape array to image dimensions: {}").format(e)) from e

                    pixeldata = np.copy(arr)

                current_offset += single_shminfo_size

            # Extract header data shm_info (always present)
            (size_be, data_type_be, width_be, height_be, channels_be, shm_name_bytes) = struct.unpack_from(shminfo_format, response, current_offset)

            shm_header_info = _SharedMemoryInfo(
                size = size_be,
                data_type = data_type_be,
                width = width_be,
                height = height_be,
                channels = channels_be,
                shm_name = shm_name_bytes
            )

            # Check if header data is available
            shm_name_str = shm_header_info.shm_name.decode('utf-8').rstrip('\x00')
            if shm_header_info.size > 0 and shm_name_str:
                try:
                    shm_header = self._map_shared_memory(shm_name_str, shm_header_info.size)
                except (OSError, ValueError) as e:
                    raise SharedMemoryError(_("Failed to map header shared memory: {}").format(e)) from e

                # CRITICAL: Must copy data out of shared memory before it's released
                buffer = bytearray(shm_header.buf[:shm_header_info.size])
                try:
                    header_data = buffer.decode('utf-8').rstrip('\x00')
                except UnicodeDecodeError as e:
                    raise DataError(_("Failed to decode header data: {}").format(e)) from e

            current_offset += single_shminfo_size

            # Extract ICC profile data shm_info (always present)
            (size_be, data_type_be, width_be, height_be, channels_be, shm_name_bytes) = struct.unpack_from(shminfo_format, response, current_offset)

            shm_icc_info = _SharedMemoryInfo(
                size = size_be,
                data_type = data_type_be,
                width = width_be,
                height = height_be,
                channels = channels_be,
                shm_name = shm_name_bytes
            )

            # Check if ICC profile data is available
            shm_name_str = shm_icc_info.shm_name.decode('utf-8').rstrip('\x00')
            if shm_icc_info.size > 0 and shm_name_str:
                try:
                    shm_icc_profile = self._map_shared_memory(shm_name_str, shm_icc_info.size)
                except (OSError, ValueError) as e:
                    raise SharedMemoryError(_("Failed to map ICC profile shared memory: {}").format(e)) from e

                # Copy data out of shared memory before it's released
                buffer = bytearray(shm_icc_profile.buf[:shm_icc_info.size])
                icc_profile_data = bytes(buffer)  # Convert bytearray copy to immutable bytes

            # Build FFit object using core_vals and fits_keywords
            fit = FFit(
                _naxes = (core_vals[0], core_vals[1], core_vals[2]),
                naxis = 2 if core_vals[2] == 1 else 3,
                bitpix = core_vals[3],
                orig_bitpix = core_vals[4],
                checksum = bool(core_vals[5]),
                mini = core_vals[6],
                maxi = core_vals[7],
                neg_ratio = core_vals[8],
                top_down = bool(core_vals[9]),
                _focalkey = bool(core_vals[10]),
                _pixelkey = bool(core_vals[11]),
                color_managed = bool(core_vals[12]),
                _data = pixeldata,
                stats = [
                    self.get_seq_stats(frame, 0),
                    self.get_seq_stats(frame, 1) if core_vals[2] > 1 else None,
                    self.get_seq_stats(frame, 2) if core_vals[2] > 1 else None,
                ],
                keywords = fits_keywords,
                _icc_profile = icc_profile_data,
                header = header_data,
                unknown_keys = None,
                history = None
            )
            return fit

        except Exception as e:
            raise SirilError(f"Error in get_seq_frame(): {e}") from e

        finally:
            # Clean up all shared memory mappings
            # Must release ALL shared memory segments that the C code sent us,
            # even if we didn't successfully open them. The C code tracks ALL allocations
            # and expects release commands for each one.

            # Close any successfully opened shared memory handles
            if shm_pixels is not None:
                try:
                    shm_pixels.close()
                except Exception:
                    pass

            if shm_header is not None:
                try:
                    shm_header.close()
                except Exception:
                    pass

            if shm_icc_profile is not None:
                try:
                    shm_icc_profile.close()
                except Exception:
                    pass

            # Now send RELEASE_SHM commands for ALL segments that were reported in the response,
            # regardless of whether we successfully opened them
            shm_infos_to_release = []

            if shm_pixels_info is not None:
                shm_infos_to_release.append(('pixels', shm_pixels_info))

            if shm_header_info is not None:
                shm_infos_to_release.append(('header', shm_header_info))

            if shm_icc_info is not None:
                shm_infos_to_release.append(('icc_profile', shm_icc_info))

            # Signal that Python is done with all shared memory segments
            for shm_type, shm_info in shm_infos_to_release:
                try:
                    shm_name_str = shm_info.shm_name.decode('utf-8').rstrip('\x00')
                    # Only send release command if there was actually a shared memory allocation
                    # (indicated by non-empty name and non-zero size)
                    if shm_name_str and shm_info.size > 0:
                        if not self._execute_command(_Command.RELEASE_SHM, shm_info):
                            print(f"Warning: Failed to cleanup {shm_type} shared memory")
                except Exception as e:
                    print(f"Exception during {shm_type} shm cleanup: {e}")
                    pass

    def get_seq_frame_header(self, frame: int, return_as = 'str') -> Union[str, dict, None]:
        """
        Retrieve the full FITS header of an image from the sequence loaded in Siril.

        Args:
            frame: Integer specifying which frame in the sequence to retrieve data for
                (between 0 and Sequence.number - 1). This uses a 0-based indexing scheme, i.e.
                the first frame is frame number 0, not frame numer 1.
            return_as: Optional string specifying the format of the returned header.
                Can be 'str' for a string or 'dict' for a dictionary.

        Returns:
            str: The image FITS header as a string, or None if there is no header.
            dict: The image FITS header as a dictionary, or None if there is no header.
            None: If the header is empty or not available.

        Raises:
            NoSequenceError: If no sequence is currently loaded,
            SirilError: For other errors during data retrieval,
        """

        if not self.is_sequence_loaded():
            raise NoSequenceError(_("Error in get_seq_frame_header(): no sequence is loaded"))

        try:
            payload = struct.pack('!I', frame)
            # Request shared memory setup
            status, response = self._send_command(_Command.GET_SEQ_FRAME_HEADER, payload)

            # Handle error responses
            if status == _Status.ERROR:
                if response:
                    error_msg = response.decode('utf-8', errors='replace')
                    if "no sequence loaded" in error_msg.lower():
                        raise NoSequenceError(_("Error in get_seq_frame_header(): no sequence is currently loaded in Siril"))
                    raise SirilError(_("Server error: {}").format(error_msg))
                raise SharedMemoryError(_("Failed to initiate shared memory transfer: Empty response"))

            if not response:
                raise SharedMemoryError(_("Failed to initiate shared memory transfer: No data received"))

            if status == _Status.NONE:
                return None

            if len(response) < 25: # No payload
                return None

            try:
                # Parse the shared memory information
                shm_info = _SharedMemoryInfo.from_buffer_copy(response)
            except (AttributeError, BufferError, ValueError) as e:
                raise SharedMemoryError(_("Invalid shared memory information received: {}").format(e)) from e

            # Map the shared memory
            try:
                shm = self._map_shared_memory(
                    shm_info.shm_name.decode('utf-8'),
                    shm_info.size
                )
            except (OSError, ValueError) as e:
                raise SharedMemoryError(_("Failed to map shared memory: {}").format(e)) from e

            try:
                # Read entire buffer at once
                buffer = bytearray(shm.buf)[:shm_info.size]
                result = buffer.decode('utf-8', errors='ignore')
            except (BufferError, ValueError, TypeError) as e:
                raise SirilError(_("Failed to create string from shared memory: {}").format(e)) from e

            if return_as == 'dict':
                return parse_fits_header(result)
            if return_as == 'str':
                return result
            raise ValueError(_(f"Invalid return_as value: {return_as}"))

        except SirilError:
            # Re-raise NoSequenceError and other SirilErrors without wrapping
            raise
        except Exception as e:
            # Wrap all other exceptions with context
            raise SirilError(_("Error in get_seq_frame_header(): {}").format(e)) from e
        finally:
            if shm is not None:
                try:
                    shm.close()  # First close the memory mapping as we have finished with it
                    # (We don't unlink it as C will do that)

                    # Signal that Python is done with the shared memory and wait for C to finish
                    if not self._execute_command(_Command.RELEASE_SHM, shm_info):
                        raise SharedMemoryError(_("Failed to cleanup shared memory"))

                except Exception:
                    pass

    def get_image_stars(self, channel: Optional[int] = None) -> List[PSFStar]:
        """
        Request star model PSF data from Siril.

        Args:
            channel: Optional int specifying the channel to retrieve from.
                    If provided 0 = Red / Mono, 1 = Green, 2 = Blue. If the
                    channel is omitted the default behavior will be used:
                    channel 0 for mono images, channel 1 (green) for color images.
                    **channel requires sirilpy v1.0.8 or higher.**

        Returns:
            List of PSFStar objects containing the star data, or None if
            no stars can be found. If stars have already been detected using
            the `findstar` command then this list will be returned, otherwise
            automatic star detection will be attempted with the current
            star finder settings.

        Raises:
            NoImageError: If no image is currently loaded,
            ValueError: If an invalid channel is provided,
            SirilError: For other errors during data retrieval,
        """

        stars = []
        shm = None
        shm_info = None

        try:
            # Sentinel value for not-provided channel
            SENTINEL_VALUE = 0xFFFFFFFF  # -1 as unsigned int

            # Validate channel if provided
            if channel is not None and (channel < 0 or channel > 2):
                raise ValueError(_("Channel must be 0 (Red/Mono), 1 (Green), or 2 (Blue)"))

            channel_val = SENTINEL_VALUE if channel is None else channel

            # Pack channel data for the command
            channel_data = struct.pack('!I', channel_val)

            # Request shared memory setup
            status, response = self._send_command(_Command.GET_PSFSTARS, channel_data)

            # Handle error responses
            if status == _Status.ERROR:
                if response:
                    error_msg = response.decode('utf-8', errors='replace')
                    if "no image loaded" in error_msg.lower():
                        raise NoImageError(_("No image is currently loaded in Siril"))
                    raise SirilError(_("Server error: {}").format(error_msg))
                raise SharedMemoryError(_("Failed to initiate shared memory transfer: Empty response"))

            if status == _Status.NONE:
                return None

            if not response:
                raise SharedMemoryError(_("Failed to initiate shared memory transfer: No data received"))

            try:
                # Parse the shared memory information
                shm_info = _SharedMemoryInfo.from_buffer_copy(response)
            except (AttributeError, BufferError, ValueError) as e:
                raise SharedMemoryError(_("Invalid shared memory information received: {}").format(e)) from e

            # Map the shared memory
            try:
                shm = self._map_shared_memory(
                    shm_info.shm_name.decode('utf-8'),
                    shm_info.size
                )
            except (OSError, ValueError) as e:
                raise SharedMemoryError(_("Failed to map shared memory: {}").format(e)) from e

            format_string = '!13d2qdq7d q d8d q 2d'  # Define the format string based on PSFStar structure
            fixed_size = struct.calcsize(format_string)

            # Read entire buffer at once using memoryview
            buffer = bytearray(shm.buf)[:shm_info.size]

            # Validate buffer size
            if shm_info.size % fixed_size != 0:
                raise ValueError(_("Buffer size {} is not a multiple "
                                "of struct size {}").format(len(buffer), fixed_size))

            num_stars = len(buffer) // fixed_size

            # Sanity check for number of stars
            if num_stars <= 0:
                raise ValueError(_("Invalid number of stars: {}").format(num_stars))

            if num_stars > 200000:  # to match the #define MAX_STARS
                num_stars = 200000
                self.log(_("Limiting stars to max 200000"))

            for i in range(num_stars):
                # Calculate start and end positions for each struct
                start = i * fixed_size
                end = start + fixed_size

                try:
                    star = PSFStar.deserialize(buffer[start:end])
                    stars.append(star)

                except struct.error as e:
                    print(f"Error unpacking star data for index {i}: {e}", file=sys.stderr)
                    break

            return stars

        except SirilError:
            raise
        except Exception as e:
            raise SirilError(_("Error in get_image_stars(): {}").format(e)) from e

        finally:
            if shm is not None:
                try:
                    shm.close()  # First close the memory mapping as we have finished with it
                    # (We don't unlink it as C will do that)

                    # Signal that Python is done with the shared memory and wait for C to finish
                    if shm_info is not None and not self._execute_command(_Command.RELEASE_SHM, shm_info):
                        raise SharedMemoryError(_("Failed to cleanup shared memory"))

                except Exception:
                    pass

    def get_siril_config(self, group: str, key: str) -> Optional[Union[bool, int, float, str, List[str]]]:
        """
        Request a configuration value from Siril.

        Args:
            group: Configuration group name,
            key: Configuration key name within the group
                 (Available values for group and key can be determined using
                 the "get -A" command)

        Returns:
            The configuration value with appropriate Python type, or None if an
            error occurred.

        Raises:
            DataError: if an unknown config type is encountered,
            SirilError: if an error occurred getting the requested config value
        """

        try:
            # Encode both group and key with null terminator
            group_bytes = group.encode('utf-8') + b'\0'
            key_bytes = key.encode('utf-8') + b'\0'
            payload = group_bytes + key_bytes

            # Request the config value
            response = self._request_data(_Command.GET_CONFIG, payload=payload)
            if response is None:
                return None

            # First byte should be the type
            config_type = _ConfigType(response[0])
            value_data = response[1:]

            # Parse based on type
            if config_type == _ConfigType.BOOL:
                return bool(struct.unpack('!I', value_data)[0])
            if config_type == _ConfigType.INT:
                return struct.unpack('!i', value_data)[0]
            if config_type == _ConfigType.DOUBLE:
                return struct.unpack('!d', value_data)[0]
            if config_type in (_ConfigType.STR, _ConfigType.STRDIR):
                # Assume null-terminated string
                string_value = value_data.split(b'\0')[0].decode('utf-8')
                return string_value
            if config_type == _ConfigType.STRLIST:
                # Split on null bytes and decode each string
                strings = value_data.split(b'\0')
                # Remove empty strings at the end
                while strings and not strings[-1]:
                    strings.pop()
                return [s.decode('utf-8') for s in strings]
            raise DataError(_("Error: unknown config type"))

        except SirilError:
            raise
        except Exception as e:
            raise SirilError(_("Error in get_siril_config(): {}").format(e)) from e

    def set_seq_frame_incl(self, index: Union[int, List[int]], incl: bool):
        """
        Set whether given frame(s) are included in the currently loaded sequence
        in Siril. This method is intended for use in creating custom sequence
        filters.

        Args:
            index: integer or list of integers specifying which frame(s) to set the 
                inclusion status for. This uses a 0-based indexing scheme, i.e. the 
                first frame is frame number 0, not frame number 1.
                Passing a list is available since sirilpy 1.0.16
            incl: bool specifying whether the frame(s) are included or not.

        Raises:
            NoSequenceError: if no sequence is loaded in Siril,
            TypeError: if index is not an int or list of ints,
            SirilError: on failure.
        """
        try:
            if not self.is_sequence_loaded():
                raise NoSequenceError(_("Error in set_seq_frame_incl(): no sequence loaded"))
            
            # Validate and normalize index parameter
            if isinstance(index, int):
                indices = [index]
            elif isinstance(index, list):
                if not index:
                    raise ValueError(_("Index list cannot be empty"))
                if not all(isinstance(i, int) for i in index):
                    raise TypeError(_("All elements in index list must be integers"))
                indices = index
            else:
                raise TypeError(_("Index must be an integer or list of integers"))
            
            # Build payload: count, indices, incl
            # Format: count (I) + indices (I * count) + incl (I)
            count = len(indices)
            format_string = f'!I{count}II'  # count + indices + incl
            payload = struct.pack(format_string, count, *indices, incl)
            
            self._execute_command(_Command.SET_SEQ_FRAME_INCL, payload)
            return
        except SirilError:
            raise
        except Exception as e:
            raise SirilError(f"Error in set_seq_frame_incl(): {e}") from e

    def get_image_bgsamples(self) -> Optional[List[BGSample]]:
        """
        Request background samples data from Siril.

        Returns:
            List of BGSamples background samples, with each set of coordinates
            expressed as a tuple[float, float], or None if no background
            samples have been set.

        Raises:
            NoImageError: If no image is currently loaded,
            DataError: on receipt of bad data,
            SirilError: For other errors during  data retrieval,
        """

        samples = []
        shm = None

        try:
            # Request shared memory setup
            status, response = self._send_command(_Command.GET_BGSAMPLES)

            # Handle error responses
            if status == _Status.ERROR:
                if response:
                    error_msg = response.decode('utf-8', errors='replace')
                    if "no image loaded" in error_msg.lower():
                        raise NoImageError(_("No image is currently loaded in Siril"))
                    raise SirilConnectionError(_("Server error: {}").format(error_msg))
                raise SharedMemoryError(_("Failed to initiate shared memory transfer: Empty response"))

            if status == _Status.NONE:
                return None # no samples have been set

            if not response:
                raise SharedMemoryError(_("Failed to initiate shared memory transfer: No data received"))

            try:
                # Parse the shared memory information
                shm_info = _SharedMemoryInfo.from_buffer_copy(response)
            except (AttributeError, BufferError, ValueError) as e:
                raise SharedMemoryError(_("Invalid shared memory information received: {}").format(e)) from e

            # Map the shared memory
            try:
                shm = self._map_shared_memory(
                    shm_info.shm_name.decode('utf-8'),
                    shm_info.size
                )
            except (OSError, ValueError) as e:
                raise SharedMemoryError(_("Failed to map shared memory: {}").format(e)) from e

            format_string = '!6dQ2dQ' # Define the format string based on background_sample structure
            fixed_size = struct.calcsize(format_string)

            # Read entire buffer at once using memoryview
            buffer = bytearray(shm.buf)[:shm_info.size]

            # Validate buffer size
            if shm_info.size % fixed_size != 0:
                raise ValueError(_("Buffer size {} is not a multiple "
                                   "of struct size {}").format(len(buffer), fixed_size))

            num_samples = len(buffer) // fixed_size

            # Sanity check for number of stars
            if num_samples < 1:
                raise ValueError(_("Invalid number of samples: {}").format(num_samples))

            if num_samples > 10000:   # to match the max samples per line settable in the UI
                                    # and assuming equal dimensions. More than enough!
                num_samples = 10000
                self.log(_("Limiting stars to max 200000"))

            for i in range(num_samples):
                # Calculate start and end positions for each struct
                start = i * fixed_size
                end = start + fixed_size

                try:
                    sample = BGSample.deserialize(buffer[start:end])
                    samples.append(sample)

                except struct.error as e:
                    raise DataError(_("Error in get_image_bgsamples(): {}").format(e)) from e

            return samples

        except SirilError:
            raise
        except Exception as e:
            raise SirilError(_("Error in get_image_bgsamples(): {}").format(e)) from e

        finally:
            if shm is not None:
                try:
                    shm.close()  # First close the memory mapping as we have finished with it
                    # (We don't unlink it as C wll do that)

                    # Signal that Python is done with the shared memory and wait for C to finish
                    if not self._execute_command(_Command.RELEASE_SHM, shm_info):
                        raise SirilError(_("Failed to cleanup shared memory"))

                except Exception:
                    pass

    def set_image_metadata_from_header_string(self, header: str) -> bool:
        """
        Send image metadata to Siril from a FITS header. The header can be
        obtained from a sirilpy FFit.header or alternatively from a FITS file
        obened from disk using astropy.fits.

        Example:

        .. code-block:: python

            hdul = fits.open('your_fits_file.fits')
            # Get the header from the primary HDU (or any other HDU you want)
            header = hdul[0].header
            # Convert the header to string
            header_string = header.tostring(sep='\\\\n')
            # Send the metadata to Siril
            siril.set_image_metadata_from_header_string(header_string)

        Args:
            header: string containing the FITS header data

        Returns:
            bool: True if successful, False otherwise

        Raises:
            TypeError: invalid parameter provided,
            NoImageError: if no image is loaded in Siril,
            SirilError: if an error occurs.

        """

        shm = None
        shm_info = None
        try:
            # Validate input array
            if not isinstance(header, str):
                raise TypeError(_("Header data must be a string"))

            if not self.is_image_loaded():
                raise NoImageError(_("Error in set_image_metadata_from_header_string(): no image loaded"))

            header_bytes = header.encode('utf-8')
            # Calculate total size
            total_bytes = len(header_bytes) + 1

            # Create shared memory using our wrapper
            try:
                # Request Siril to provide a big enough shared memory allocation
                shm_info = self._request_shm(total_bytes)

                # Map the shared memory
                try:
                    shm = self._map_shared_memory(
                        shm_info.shm_name.decode('utf-8'),
                        shm_info.size
                    )
                except (OSError, ValueError) as e:
                    raise SharedMemoryError(_("Failed to map shared memory: {}").format(e)) from e

            except Exception as e:
                raise SharedMemoryError(_("Failed to create shared memory: {}").format(e)) from e

            # Copy data to shared memory
            try:
                buffer = memoryview(shm.buf).cast('B')
                buffer[:len(header_bytes)] = header_bytes

                # Delete transient objects used to structure copy
                del buffer
            except Exception as e:
                raise SirilError(_("Failed to copy data to shared memory: {}").format(e)) from e

            # Pack the image info structure
            info = struct.pack(
                '!IIIIQ256s',
                0,
                0,
                0,
                0,
                total_bytes,
                shm_info.shm_name
            )
            # Send command using the existing _execute_command method
            if not self._execute_command(_Command.SET_IMAGE_HEADER, info):
                raise SirilConnectionError(_("_execute_command failed"))

            return True

        except (TypeError, SirilError):
            raise
        except Exception as e:
            raise SirilError(f"Error in set_image_metadata_from_header_string(): {e}") from e

        finally:
            if shm is not None:
                try:
                    shm.close()
                    self._execute_command(_Command.RELEASE_SHM, shm_info)
                except Exception as e:
                    pass

    def overlay_add_polygon(self, polygon: Polygon) -> Polygon:
        """
        Adds a user polygon to the Siril display overlay
        Args:
            polygon: Polygon defining the polygon to be added
        Returns:
            Polygon: the input updated with the ID assigned by Siril
        Raises:
            TypeError: invalid data provided,
            SirilConnectionError: on a connection failure,
            DataError: on receipt of invalid data,
            SirilError: on any failure
        """
        shm = None
        shm_info = None
        try:
            # Serialize the provided polygon
            polygon_bytes = polygon.serialize()

            # Calculate total size needed
            total_bytes = len(polygon_bytes)

            # Create shared memory using our wrapper
            try:
                # Request Siril to provide a big enough shared memory allocation
                shm_info = self._request_shm(total_bytes)
                # Map the shared memory
                try:
                    shm = self._map_shared_memory(
                        shm_info.shm_name.decode('utf-8'),
                        shm_info.size
                    )
                except (OSError, ValueError) as e:
                    raise SharedMemoryError(_("Failed to map shared memory: {}").format(e)) from e

            except Exception as e:
                raise SharedMemoryError(_("Failed to create shared memory: {}").format(e)) from e

            # Copy polygon data to shared memory
            try:
                buffer = memoryview(shm.buf).cast('B')
                buffer[:len(polygon_bytes)] = polygon_bytes

                # Delete transient objects used to structure copy
                del buffer
            except Exception as e:
                raise SirilError(_("Failed to copy data to shared memory: {}").format(e)) from e

            # Pack the polygon info structure
            info = struct.pack(
                '!IIIIQ256s',
                0,  # unused field 1
                0,  # unused field 2
                0,  # unused field 3
                0,  # unused field 4
                total_bytes,  # size of polygon data
                shm_info.shm_name  # shared memory name
            )

            # Send command using the existing _execute_command method
            response = self._request_data(_Command.ADD_USER_POLYGON, info)
            if response is None:
                raise SirilConnectionError(("Failed to get a response from the SirilInterface"))
            try:
                polygon_id = struct.unpack('!i', response[:4])[0]
                polygon.polygon_id = polygon_id
                return polygon

            except struct.error as e:
                raise DataError(_("Error unpacking polygon ID")) from e

        except (TypeError, SirilError, DataError, SirilConnectionError, SharedMemoryError):
            raise
        except Exception as e:
            raise SirilError(f"Error in overlay_add_polygon(): {e}") from e
        finally:
            if shm is not None:
                try:
                    shm.close()
                    if shm_info is not None:
                        self._execute_command(_Command.RELEASE_SHM, shm_info)
                except Exception:
                    pass

    def overlay_delete_polygon(self, polygon_id: int):
        """
        Deletes a single user polygon from the Siril overlay, specified by ID

        Args:
            id: int specifying the polygon ID to be deleted

        Raises:
            SirilError: on failure
        """
        try:
            # Create payload: network-order int followed by string
            # '!I' for network byte order 32-bit int
            payload = struct.pack('!i', polygon_id)

            self._execute_command(_Command.DELETE_USER_POLYGON, payload)
            return

        except Exception as e:
            raise SirilError(f"Error in overlay_delete_polygon(): {e}") from e

    def overlay_clear_polygons(self) -> bool:
        """
        Clears all user polygons from the Siril overlay

        Returns:
            bool: True if the command succeeded, False otherwise

        Raises:
            SirilError: if an error occurred
            SharedMemoryError: if a shared memory error occurred
        """

        try:
            self._execute_command(_Command.CLEAR_USER_POLYGONS, None)
            return

        except Exception as e:
            raise SirilError(f"Error in overlay_clear_polygons(): {e}") from e

    def overlay_draw_polygon(self, color=0x00FF0040, fill=False):
        """
        Enters a mode where the user can draw a Polygon in the Siril window
        by clicking the main mouse button and dragging. Releasing the mouse
        button finalises and closes the polygon.

        Args:
            color: uint32 specifying packed RGBA values. Default: 0x00FF0040, 75% transparent green)
            fill: bool specifying whether or not to fill the polygon (default: False)
        """
        payload = struct.pack('!I?', color, fill)
        status, response = self._send_command(_Command.DRAW_POLYGON, payload)
        if status == _Status.ERROR:
            if response:
                error_msg = response.decode('utf-8', errors='replace')
                if "no image loaded" in error_msg.lower():
                    raise NoImageError(_("No image is currently loaded in Siril"))
                raise SirilConnectionError(_("Server error: {}").format(error_msg))
            raise SirilConnectionError(_("Failed to confirm command: Empty response"))

        if status == _Status.NONE:
            # Not a serious error but raising this allows it to be handled
            raise MouseModeError(_("Cannot draw a polygon at present"))

    def overlay_get_polygon(self, polygon_id: int) -> 'Polygon':
        """
        Gets a single user polygon from the Siril overlay, specified by ID

        Args:
            id: int specifying the polygon ID to be retrieved. The special ID -1 will
        retrieve the most recently added polygon.

        Returns:
            Polygon: the specified Polygon if it exists, None otherwise

        Raises:
            SirilError: on failure
        """
        shm = None
        try:
            payload = struct.pack('!i', polygon_id)
            # Send it using _request_data
            status, response = self._send_command(_Command.GET_USER_POLYGON, payload)
            # Handle error responses
            if status == _Status.ERROR:
                if response:
                    error_msg = response.decode('utf-8', errors='replace')
                    if "no image loaded" in error_msg.lower():
                        raise NoImageError(_("No image is currently loaded in Siril"))
                    raise SirilConnectionError(_("Server error: {}").format(error_msg))
                raise SirilConnectionError(_("Failed to confirm command: Empty response"))

            if status == _Status.NONE:
                return None # May be correct if there are no user polygons defined

            if not response:
                raise SharedMemoryError(_("Failed to initiate shared memory transfer: No data received"))

            try:
                # Parse the shared memory information
                shm_info = _SharedMemoryInfo.from_buffer_copy(response)
            except (AttributeError, BufferError, ValueError) as e:
                raise ValueError(_("Invalid shared memory information received: {}").format(e)) from e

            # Map the shared memory
            try:
                shm = self._map_shared_memory(
                    shm_info.shm_name.decode('utf-8'),
                    shm_info.size
                )
            except (OSError, ValueError) as e:
                raise SharedMemoryError(_("Failed to map shared memory: {}").format(e)) from e

            # Read entire buffer at once using memoryview
            buffer = bytearray(shm.buf)[:shm_info.size]

            polygon = Polygon.deserialize_polygon(buffer)[0]

            return polygon

        except Exception as e:
            raise SirilError(_("Error in overlay_get_polygon(): {}").format(e)) from e

        finally:
            if shm is not None:
                try:
                    shm.close()  # First close the memory mapping as we have finished with it
                    # (We don't unlink it as C wll do that)

                    # Signal that Python is done with the shared memory and wait for C to finish
                    if not self._execute_command(_Command.RELEASE_SHM, shm_info):
                        raise SirilError(_("Failed to cleanup shared memory"))

                except Exception:
                    pass

    def overlay_get_polygons_list(self) -> List['Polygon']:
        """
        Gets a List of all user polygons from the Siril overlay

        Returns:
            List[Polygon]: the list of Polygon if some exist, None otherwise

        Raises:
            SirilError: if an error occurred
            SharedMemoryError: if a shared memory error occurred
        """
        shm = None
        try:
            status, response = self._send_command(_Command.GET_USER_POLYGON_LIST)
            # Handle error responses
            if status == _Status.ERROR:
                if response:
                    error_msg = response.decode('utf-8', errors='replace')
                    if "no image loaded" in error_msg.lower():
                        raise NoImageError(_("No image is currently loaded in Siril"))
                    raise SirilConnectionError(_("Server error: {}").format(error_msg))
                raise SharedMemoryError(_("Failed to initiate shared memory transfer: Empty response"))

            if status == _Status.NONE:
                return None # May be correct if there are no user polygons defined

            if not response:
                raise SharedMemoryError(_("Failed to initiate shared memory transfer: No data received"))

            try:
                # Parse the shared memory information
                shm_info = _SharedMemoryInfo.from_buffer_copy(response)
            except (AttributeError, BufferError, ValueError) as e:
                raise ValueError(_("Invalid shared memory information received: {}").format(e)) from e

            # Map the shared memory
            try:
                shm = self._map_shared_memory(
                    shm_info.shm_name.decode('utf-8'),
                    shm_info.size
                )
            except (OSError, ValueError) as e:
                raise SharedMemoryError(_("Failed to map shared memory: {}").format(e)) from e

            # Read entire buffer at once using memoryview
            buffer = bytearray(shm.buf)[:shm_info.size]

            polygon_list = Polygon.deserialize_polygon_list(buffer)

            return polygon_list

        except Exception as e:
            raise SirilError(_("Error in overlay_get_polygons_list(): {}").format(e)) from e

        finally:
            if shm is not None:
                try:
                    shm.close()  # First close the memory mapping as we have finished with it
                    # (We don't unlink it as C wll do that)

                    # Signal that Python is done with the shared memory and wait for C to finish
                    if not self._execute_command(_Command.RELEASE_SHM, shm_info):
                        raise SirilError(_("Failed to cleanup shared memory"))

                except Exception:
                    pass

    def create_new_seq(self, seq_root: str) -> bool:
        """
        Creates a new .seq file with all images named seq_rootXXXXX.ext located in
        the current home folder. If a sequence with the same name is already loaded
        in Siril, it will not be recreated. This only works for FITS files, not FITSEQ nor SER.
        The newly created sequence is not loaded in Siril.

        Args:
            seq_root: The root name of the sequence to be created.

        Returns:
            bool: True if the sequence was successfully created, False otherwise.

        Raises:
            SirilError: if an error occurred.

        """

        try:
            if self.is_sequence_loaded():
                seq = self.get_seq()
                seq1_stripped = seq.seqname.rstrip('_')
                seq2_stripped = seq_root.rstrip('_')
                if seq1_stripped == seq2_stripped:
                    self.log(_('A sequence with the same name is already loaded in Siril, aborting'), LogColor.RED)
                    return False
            home_folder = self.get_siril_wd()
            all_files = os.listdir(home_folder)
            ext = self.get_siril_config('core','extension')
            pattern = fr'^{re.escape(seq_root)}\d{{5}}{re.escape(ext)}$'
            regex = re.compile(pattern)
            exact_matches = [os.path.join(home_folder, f) for f in all_files if regex.match(f)]
            if len(exact_matches) == 0:
                self.log(_(f'No files matching {seq_root} pattern found in the Home folder'), LogColor.RED)
                return False
            if len(exact_matches) == 1:
                self.log(_(f'Only one file matching {seq_root} found in the Home folder, cannot create sequence'), LogColor.RED)
            seq_root_dummy_ext = seq_root + ".ext" # We have to add this to work around remove_ext_from_filename in create_one_regular_seq
            message_bytes = seq_root_dummy_ext.encode('utf-8')
            return self._execute_command(_Command.CREATE_NEW_SEQ, message_bytes)

        except Exception as e:
            raise SirilError(f"Error in create_new_seq(): {e}") from e

    def is_cli(self) -> bool:
        """
        Check if the current instance is running in CLI mode. This method is useful
        to detect how the script was invoked and whether to show or not a GUI.
        This is False when the script is called by clicking in the Script menu,
        True otherwise.

        Returns:
            bool: True if running in CLI mode, False otherwise.
        """
        return self._is_cli

    def load_image_from_file(self, filepath: str, with_pixels: Optional[bool] = True,
                            preview: Optional[bool] = False,
                            linked: Optional[bool] = False) -> Optional[FFit]:
        """
        Request Siril to load an image from a file and transfer it to sirilpy. This
        method does not change the image currently loaded in Siril. Any image format
        supported by Siril is supported. This may be used as an alternative to loading
        an image using astropy.io.fits, however perhaps the main benefit to using it
        is that it supports the preview option which can be used to obtain an 8-bit
        autostretched rendering of the image more quickly than is possible using astropy
        and applying an autostretch using numpy.

        Args:
            filepath: String specifying the path to the image file to load.
            with_pixels: bool specifying whether or not to return the pixel data for the
                image (default is True).
            preview: bool specifying whether or not to return the real pixel data or an
                autostretched uint8_t preview version. Only has an effect in
                conjunction with with_pixels = True
            linked: bool specifying whether the autostretch preview should be linked or
                    not. Has no effect unless preview is True.

        Returns:
            FFit object containing the image data

        Raises:
            FileNotFoundError: if the specified file does not exist,
            DataError: on receipt of incorrect data,
            SirilError: if an error occurs.
        """

        if not filepath or not isinstance(filepath, str):
            raise ValueError("filepath must be a non-empty string")

        # Convert filepath to bytes for transmission
        filepath_bytes = filepath.encode('utf-8')
        if len(filepath_bytes) > 4000:  # Reasonable limit for filepath length
            raise ValueError("filepath is too long")

        shm_pixels = None
        shm_header = None
        shm_icc_profile = None

        # Pack boolean flags + filepath
        data_payload = struct.pack('!???', with_pixels, preview, linked) + filepath_bytes
        response = self._request_data(_Command.GET_IMAGE_FILE, data_payload, timeout=None)
        if response is None:
            return None

        try:
            # Build format string for struct unpacking
            # Network byte order for all values
            FLEN_VALUE = 71  # Standard FITS keyword length
            format_parts = [
                # Core FFfit (start at index 0)
                'q',  # rx padded to 64bit
                'q',  # ry padded to 64bit
                'q',  # naxes[2] padded to 64bit
                'q',  # bitpix padded to 64bit
                'q',  # orig_bitpix padded to 64bit
                'Q',  # gboolean checksum padded to 64bit
                'd',  # mini
                'd',  # maxi
                'd',  # neg_ratio padded to 64bit
                'Q',  # gboolean top_down (padded to uint64_t)
                'Q',  # gboolean focalkey (padded to uint64_t)
                'Q',  # gboolean pixelkey (padded to uint64_t)
                'Q']  # gboolean color_managed (padded to uint64_t)
            # Add keywords format (skip the '!' prefix)
            format_parts.extend(FKeywords._KEYWORD_FORMAT_PARTS)
            keywords_count = len(FKeywords._KEYWORD_FORMAT_PARTS)

            # Add stats for 3 channels (14 doubles each)
            for i in range(3):
                format_parts.extend(['d'] * 14)  # 14 doubles per channel

            # Add shared memory info structs
            # Pixel data shm_info (always present, may be zeroed)
            format_parts.extend([
                'Q',  # size (size_t)
                'i',  # data_type
                'i',  # width
                'i',  # height
                'i',  # channels
                '256s'  # shm_name (char[256])
            ])

            # Header shm_info (always present, may be zeroed)
            format_parts.extend([
                'Q',  # size (size_t)
                'i',  # data_type
                'i',  # width
                'i',  # height
                'i',  # channels
                '256s'  # shm_name (char[256])
            ])

            # ICC Profile shm_info (always present, may be zeroed)
            format_parts.extend([
                'Q',  # size (size_t)
                'i',  # data_type
                'i',  # width
                'i',  # height
                'i',  # channels
                '256s'  # shm_name (char[256])
            ])

            format_string = '!' + ''.join(format_parts)

            # Verify data size
            expected_size = struct.calcsize(format_string)
            if len(response) != expected_size:
                raise DataError(f"Received image data size {len(response)} doesn't match expected size {expected_size}")

            # Unpack the binary data
            values = struct.unpack(format_string, response)

            # Helper function to decode and strip null-terminated strings
            def decode_string(s: bytes) -> str:
                return s.decode('utf-8').rstrip('\x00')

            # Helper function to convert timestamp to datetime
            def timestamp_to_datetime(timestamp: int) -> Optional[datetime]:
                return datetime.fromtimestamp(timestamp) if timestamp != 0 else None

            # Extract stats data (starts after core fields + keywords)
            stats = [None, None, None]
            stats_start_idx = 13 + keywords_count
            for channel in range(3):
                start_idx = stats_start_idx + (channel * 14)
                # Check if this channel has valid stats (non-zero values)
                channel_stats = values[start_idx:start_idx + 14]
                if any(stat != 0.0 for stat in channel_stats):
                    # Create ImageStats object directly from values
                    stats[channel] = ImageStats(
                        total=int(channel_stats[0]),
                        ngoodpix=int(channel_stats[1]),
                        mean=channel_stats[2],
                        median=channel_stats[3],
                        sigma=channel_stats[4],
                        avgDev=channel_stats[5],
                        mad=channel_stats[6],
                        sqrtbwmv=channel_stats[7],
                        location=channel_stats[8],
                        scale=channel_stats[9],
                        min=channel_stats[10],
                        max=channel_stats[11],
                        normValue=channel_stats[12],
                        bgnoise=channel_stats[13]
                    )

            # Extract shared memory info structs
            # Calculate the correct starting indices based on the actual format
            keywords_count = len(FKeywords._KEYWORD_FORMAT_PARTS)
            pixel_start_idx = 13 + keywords_count + 42  # 13 core fields + keywords + 42 stats fields

            # Header shm_info starts after pixel shm_info (6 fields)
            header_start_idx = pixel_start_idx + 6

            # ICC Profile shm_info starts after header shm_info (6 fields)
            icc_start_idx = header_start_idx + 6

            # Get pixeldata if requested and available
            pixeldata = None
            if with_pixels:
                try:
                    shm_pixels_info = _SharedMemoryInfo(
                        size=values[pixel_start_idx],
                        data_type=values[pixel_start_idx + 1],
                        width=values[pixel_start_idx + 2],
                        height=values[pixel_start_idx + 3],
                        channels=values[pixel_start_idx + 4],
                        shm_name=values[pixel_start_idx + 5]
                    )

                    # Check if pixel data is available (non-zero size and valid shm_name)
                    shm_name_str = shm_pixels_info.shm_name.decode('utf-8').rstrip('\x00')
                    if shm_pixels_info.size > 0 and shm_name_str:
                        # Validate dimensions
                        if any(dim <= 0 for dim in (shm_pixels_info.width, shm_pixels_info.height)):
                            raise DataError(_("Invalid image dimensions: {}x{}").format(shm_pixels_info.width, shm_pixels_info.height))

                        if shm_pixels_info.channels <= 0 or shm_pixels_info.channels > 3:
                            raise DataError(_("Invalid number of channels: {}").format(shm_pixels_info.channels))

                        # Map the shared memory
                        try:
                            shm_pixels = self._map_shared_memory(shm_name_str, shm_pixels_info.size)
                        except (OSError, ValueError) as e:
                            raise SharedMemoryError(_("Failed to map pixel shared memory: {}").format(e)) from e

                        buffer = bytearray(shm_pixels.buf)[:shm_pixels_info.size]
                        # Create numpy array from shared memory
                        if preview:
                            dtype = np.uint8
                        else:
                            dtype = np.float32 if shm_pixels_info.data_type == 1 else np.uint16
                        try:
                            arr = np.frombuffer(buffer, dtype=dtype)
                        except (BufferError, ValueError, TypeError) as e:
                            raise SharedMemoryError(_("Failed to create array from shared memory: {}").format(e)) from e

                        # Validate array size matches expected dimensions
                        expected_size = shm_pixels_info.width * shm_pixels_info.height * shm_pixels_info.channels
                        if arr.size < expected_size:
                            raise DataError(
                                f"Data size mismatch: got {arr.size} elements, "
                                f"expected {expected_size} for dimensions "
                                f"{shm_pixels_info.width}x{shm_pixels_info.height}x{shm_pixels_info.channels}"
                            )

                        # Reshape the array according to the image dimensions
                        try:
                            if shm_pixels_info.channels > 1:
                                arr = arr.reshape((shm_pixels_info.channels, shm_pixels_info.height, shm_pixels_info.width))
                            else:
                                arr = arr.reshape((shm_pixels_info.height, shm_pixels_info.width))
                        except ValueError as e:
                            raise DataError(_("Error in get_image_file(): Failed to reshape array to image dimensions: {}").format(e)) from e

                        # Make a copy of the data since we'll be releasing the shared memory
                        pixeldata = np.copy(arr)

                except SirilError:
                    raise
                except Exception as e:
                    raise SirilError(f"Error obtaining pixeldata: {e}") from e

            # Get header data if available
            header_data = None
            try:
                shm_header_info = _SharedMemoryInfo(
                    size=values[header_start_idx],
                    data_type=values[header_start_idx + 1],
                    width=values[header_start_idx + 2],
                    height=values[header_start_idx + 3],
                    channels=values[header_start_idx + 4],
                    shm_name=values[header_start_idx + 5]
                )

                # Check if header data is available (non-zero size and valid shm_name)
                shm_name_str = shm_header_info.shm_name.decode('utf-8').rstrip('\x00')
                if shm_header_info.size > 0 and shm_name_str:
                    try:
                        shm_header = self._map_shared_memory(shm_name_str, shm_header_info.size)
                    except (OSError, ValueError) as e:
                        raise SharedMemoryError(_("Failed to map header shared memory: {}").format(e)) from e

                    buffer = bytes(shm_header.buf)[:shm_header_info.size]
                    # Header data should be text, decode as UTF-8
                    try:
                        header_data = buffer.decode('utf-8').rstrip('\x00')
                    except UnicodeDecodeError as e:
                        raise DataError(_("Failed to decode header data: {}").format(e)) from e

            except SirilError:
                raise
            except Exception as e:
                raise SirilError(f"Error obtaining header data: {e}") from e

            # Get ICC profile data if available
            icc_profile_data = None
            try:
                shm_icc_info = _SharedMemoryInfo(
                    size=values[icc_start_idx],
                    data_type=values[icc_start_idx + 1],
                    width=values[icc_start_idx + 2],
                    height=values[icc_start_idx + 3],
                    channels=values[icc_start_idx + 4],
                    shm_name=values[icc_start_idx + 5]
                )

                # Check if ICC profile data is available (non-zero size and valid shm_name)
                shm_name_str = shm_icc_info.shm_name.decode('utf-8').rstrip('\x00')
                if shm_icc_info.size > 0 and shm_name_str:
                    try:
                        shm_icc_profile = self._map_shared_memory(shm_name_str, shm_icc_info.size)
                    except (OSError, ValueError) as e:
                        raise SharedMemoryError(_("Failed to map ICC profile shared memory: {}").format(e)) from e

                    buffer = bytes(shm_icc_profile.buf)[:shm_icc_info.size]
                    # Make a copy of the ICC profile data
                    icc_profile_data = bytes(buffer)

            except SirilError:
                raise
            except Exception as e:
                raise SirilError(f"Error obtaining ICC profile data: {e}") from e

            keywords_start_idx = 13
            keywords_end_idx = 13 + keywords_count
            keys = values[keywords_start_idx:keywords_end_idx]

            # Pack the keywords back into bytes for deserialization
            keywords_format = '!' + ''.join(FKeywords._KEYWORD_FORMAT_PARTS)
            keywords_bytes = struct.pack(keywords_format, *keys)
            fits_keywords = FKeywords.deserialize(keywords_bytes)

            fit = FFit(
                _naxes=(values[0], values[1], values[2]),
                naxis=2 if values[2] == 1 else 3,
                bitpix=values[3],
                orig_bitpix=values[4],
                checksum=bool(values[5]),
                mini=values[6],
                maxi=values[7],
                neg_ratio=values[8],
                top_down=bool(values[9]),
                _focalkey=bool(values[10]),
                _pixelkey=bool(values[11]),
                color_managed=bool(values[12]),
                _data=pixeldata,
                stats=stats,  # Now populated with ImageStats objects
                keywords=fits_keywords,
                _icc_profile=icc_profile_data,
                header=header_data,  # Now populated from shared memory
                unknown_keys=None,
                history=None
            )
            return fit

        except Exception as e:
            raise SirilError(f"Error in get_image_file(): {e}") from e

        finally:
            # Clean up all shared memory mappings
            shm_infos_to_cleanup = []

            if shm_pixels is not None:
                try:
                    shm_pixels.close()
                    shm_infos_to_cleanup.append(('pixels', shm_pixels_info))
                except Exception:
                    pass

            if shm_header is not None:
                try:
                    shm_header.close()
                    shm_infos_to_cleanup.append(('header', shm_header_info))
                except Exception:
                    pass

            if shm_icc_profile is not None:
                try:
                    shm_icc_profile.close()
                    shm_infos_to_cleanup.append(('icc_profile', shm_icc_profile))
                except Exception:
                    pass

            # Signal that Python is done with all shared memory segments
            for shm_type, shm_info in shm_infos_to_cleanup:
                try:
                    if not self._execute_command(_Command.RELEASE_SHM, shm_info):
                        raise SirilError(_("Failed to cleanup {} shared memory").format(shm_type))
                except Exception:
                    pass

    def analyse_image_from_file(self, filepath: str) -> ImageAnalysis:
        """
        Request Siril to load an image from a file and analyse it. This
        method does not change the image currently loaded in Siril. Any image format
        supported by Siril is supported. An ImageAnalysis object is returned, containing
        parameters that may be used to assess the quality of an image for use in culling.

        Args:
            filepath: String specifying the path to the image file to load.

        Returns:
            ImageAnalysis object containing quality metrics for the analysed image.

        Raises:
            FileNotFoundError: if the specified file does not exist,
            DataError: on receipt of incorrect data,
            SirilError: if an error occurs.
        """

        if not filepath or not isinstance(filepath, str):
            raise ValueError("filepath must be a non-empty string")

        if not os.path.isfile(filepath):
            raise FileNotFoundError(f"No such file: {filepath}")

        # Convert filepath to bytes for transmission
        filepath_bytes = filepath.encode('utf-8')
        if len(filepath_bytes) > 4000:  # Reasonable limit for filepath length
            raise ValueError("filepath is too long")

        response = self._request_data(_Command.ANALYSE_IMAGE_FILE, filepath_bytes, timeout=None)
        if response is None:
            raise SirilError(f"Error: no response received in SirilInterface.analyse_image_from_file()")

        try:
            analysis = ImageAnalysis.deserialize(response)
            return analysis
        except Exception as e:
            raise SirilError(f"Error unpacking data in SirilInterface.analyse_image_from_file(): {e}") from e

    def undo(self):
        """
        Undoes the last operation, if there is an undo history available.
        """

        try:
            if self._execute_command(_Command.UNDO, None):
                return True
            else:
                raise SirilError(f"Error in undo(): {e}") from e

        except Exception as e:
            raise SirilError(f"Error in undo(): {e}") from e

    def redo(self):
        """
        Redoes the last undone operation, if there is an undo history and
        an undone operation available to be redone.
        """

        try:
            if self._execute_command(_Command.REDO, None):
                return True
            else:
                raise SirilError(f"Error in redo(): {e}") from e

        except Exception as e:
            raise SirilError(f"Error in redo(): {e}") from e

    def clear_undo_history(self):
        """
        Clears the image undo history, if there is undo history to be cleared
        """

        try:
            if self._execute_command(_Command.CLEAR_UNDO_HISTORY, None):
                return True
            else:
                raise SirilError(f"Error in clear_undo_history(): {e}") from e

        except Exception as e:
            raise SirilError(f"Error in clear_undo_history(): {e}") from e

    def set_image_iccprofile(self, iccprofile: bytes):
        """
        Set the image ICC profile in Siril

        Args:
            iccprofile (bytes): The ICC profile to send to Siril. This will
            replace an existing ICC profile, if one is set. If None,
            any existing image ICC profile will be removed.

        Returns: True if the command succeeded, otherwise False

        Raises:
            NoImageError: if no image is loaded in Siril,
            ValueError: if iccprofile is not valid type,
            SirilError: if there was a Siril error in handling the command.
        """
        try:
            if not self.is_image_loaded():
                raise NoImageError(_("Error in set_image_bgsamples(): no image loaded"))

            if iccprofile is None: # Remove the ICC profile
                return self.cmd("icc_remove")

            total_bytes = len(iccprofile)

            # Create shared memory using our wrapper
            try:
                # Request Siril to provide a big enough shared memory allocation
                shm_info = self._request_shm(total_bytes)
                # Map the shared memory
                try:
                    shm = self._map_shared_memory(
                        shm_info.shm_name.decode('utf-8'),
                        shm_info.size
                    )
                except (OSError, ValueError) as e:
                    raise SharedMemoryError(_("Failed to map shared memory: {}").format(e)) from e
            except Exception as e:
                raise SharedMemoryError(_("Failed to create shared memory: {}").format(e)) from e

            # Copy serialized data to shared memory
            try:
                buffer = memoryview(shm.buf).cast('B')
                buffer[:total_bytes] = iccprofile
                del buffer
            except Exception as e:
                raise SharedMemoryError(f"Failed to copy data to shared memory: {e}") from e

            # Pack the plot info structure
            info = struct.pack(
                '!IIIIQ256s',
                0,  # width (not used for ICC profile)
                0,  # height (not used for ICC profile)
                0, # not used
                0, # not used
                total_bytes,
                shm_info.shm_name
            )

            if not self._execute_command(_Command.SET_IMAGE_ICCPROFILE, info):
                raise SirilError(_("Failed to send set_image_iccprofile() command"))
            return True

        except (ValueError, SirilError):
            raise
        except Exception as e:
            raise SirilError(f"Error in set_image_iccprofile(): {e}") from e
        finally:
            # Ensure shared memory is closed and unlinked
            if 'shm' in locals() and shm is not None:
                try:
                    shm.close()
                except Exception as e:
                    pass

    def get_siril_slider_state(self, float_range: Optional[bool] = False) -> Tuple[int, int, 'SlidersMode']:

        """
        Request the display slider state from Siril.

        Returns:
            A tuple (min, max, slider mode) representing the current slider state.

        Raises:
            SirilError: if an error occurred.
        """

        response = self._request_data(_Command.GET_SLIDER_STATE)
        if response is None:
            return None

        try:
            lo, hi, mode = struct.unpack('!HHI', response)
            if float_range:
                lo = lo / 65535.0
                hi = hi / 65535.0
            return lo, hi, mode  # Returning as (lo, hi, mode)
        except struct.error as e:
            raise SirilError(_("Error occurred in get_siril_slider_state(): {}").format(e)) from e

    def set_siril_slider_mode(self, mode: 'SlidersMode') -> bool:
        """
        Set the slider state in Siril using the provided lo, hi and mode values.
        Args:
            mode: SlidersMode enum value

        Raises:
            SirilError: if an error occurred.
            ValueError: if parameters are not properly provided.

        Returns:
            bool: True if the slider state was set successfully
        """
        try:
            # Check for valid argument combinations

            # Validate that mode is a SlidersMode enum if provided
            if not isinstance(mode, SlidersMode):
                raise ValueError("Mode must be a SlidersMode enum value")

            # Pack the values into bytes using network byte order (!)
            # Format depends on which arguments are provided
            payload = struct.pack('!I', mode.value)
            self._execute_command(_Command.SET_SLIDER_STATE, payload)
            return True
        except Exception as e:
            raise SirilError(f"Error in set_siril_slider_state(): {e}") from e

    def set_siril_slider_lohi(self,
                            lo: Union[float, int] = None,
                            hi: Union[float, int] = None) -> bool:
        """
        Set the slider values in Siril using the provided lo and hi values.
        If the sliders mode is not already set to USER, it is set to that mode
        as setting slider values is only relevant in that mode.
        Args:
            lo: Low value for the slider (float [0,1] or uint16)
            hi: High value for the slider (float [0,1] or uint16)

        Raises:
            SirilError: if an error occurred.
            ValueError: if parameters are not properly provided.

        Returns:
            bool: True if the slider state was set successfully
        """
        try:

            # Convert float values to uint16 if necessary (but pack as uint32)
            def convert_value(value):
                if value is None:
                    return None
                if isinstance(value, float):
                    if not 0.0 <= value <= 1.0:
                        raise ValueError("Float values must be in range [0,1]")
                    return int(value * 65535)
                elif isinstance(value, int):
                    if not 0 <= value <= 65535:
                        raise ValueError("Integer values must be uint16 (0-65535)")
                    return value
                else:
                    raise ValueError("Values must be float or int")

            converted_lo = convert_value(lo)
            converted_hi = convert_value(hi)

            # Pack the values into bytes using network byte order (!)
            # All values packed as 32-bit integers for consistency
            payload = struct.pack('!II', converted_lo, converted_hi)
            self._execute_command(_Command.SET_SLIDER_LOHI, payload)
            return True
        except Exception as e:
            raise SirilError(f"Error in set_siril_slider_state(): {e}") from e

    def get_siril_stf(self) -> 'STFType':

        """
        Request the Screen Transfer Function in use in Siril.

        Returns:
            A STFType representing the current STF state.

        Raises:
            SirilError: if an error occurred.
        """

        response = self._request_data(_Command.GET_STFMODE)
        if response is None:
            return None

        try:
            mode = struct.unpack('!I', response)
            return mode[0]
        except struct.error as e:
            raise SirilError(_("Error occurred in get_siril_stf(): {}").format(e)) from e

    def set_siril_stf(self, mode: 'SlidersMode') -> bool:
        """
        Set the screen transfer function in Siril using the provided enum value.
        Args:
            mode: STFType enum value

        Raises:
            SirilError: if an error occurred.
            ValueError: if parameters are not properly provided.

        Returns:
            bool: True if the slider state was set successfully
        """
        try:
            # Check for valid argument combinations

            # Validate that mode is a STFType enum
            if not isinstance(mode, STFType):
                raise ValueError("Mode must be a STFType enum value")

            # Pack the values into bytes using network byte order (!)
            # Format depends on which arguments are provided
            payload = struct.pack('!I', mode.value)
            self._execute_command(_Command.SET_STFMODE, payload)
            return True
        except Exception as e:
            raise SirilError(f"Error in set_siril_stftype(): {e}") from e

    def get_siril_panzoom(self) -> Tuple[float, float, float]:

        """
        Request the pan and zoom state from Siril.

        Returns:
            A tuple (x offset, y offset, zoomlevel) representing the current image pan and zoom state.

        Raises:
            SirilError: if an error occurred.
        """

        response = self._request_data(_Command.GET_PANZOOM)
        if response is None:
            return None

        try:
            x_off, y_off, zoomlevel = struct.unpack('!ddd', response)
            return x_off, y_off, zoomlevel
        except struct.error as e:
            raise SirilError(_("Error occurred in get_siril_slider_state(): {}").format(e)) from e

    def set_siril_pan(self,
                            xoff: float,
                            yoff: float) -> bool:
        """
        Set the display offset (pan) in Siril to (xoff, yoff).
        Args:
            xoff (float): x display offset
            yoff (float): y display offset

        Raises:
            SirilError: if an error occurred.
            ValueError: if parameters are not properly provided.

        Returns:
            bool: True if the display offset was set successfully
        """
        try:

            # Pack the values into bytes using network byte order (!)
            # All values packed as 32-bit integers for consistency
            payload = struct.pack('!dd', xoff, yoff)
            self._execute_command(_Command.SET_PAN, payload)
            return True
        except Exception as e:
            raise SirilError(f"Error in set_siril_pan(): {e}") from e

    def set_siril_zoom(self, zoom: float) -> bool:
        """
        Set the display zoom. Passing any negative value will set 'zoom to fit'.
        Args:
            zoom (float): zoom level

        Raises:
            SirilError: if an error occurred.
            ValueError: if parameters are not properly provided.

        Returns:
            bool: True if the display offset was set successfully
        """
        try:

            # Pack the values into bytes using network byte order (!)
            # All values packed as 32-bit integers for consistency
            payload = struct.pack('!d', zoom)
            self._execute_command(_Command.SET_ZOOM, payload)
            return True
        except Exception as e:
            raise SirilError(f"Error in set_siril_zoom(): {e}") from e

    def get_siril_display_iccprofile(self) -> Optional[bytes]:
        """
        Retrieve the Siril display ICC profile.

        Args:
        none.

        Returns:
            bytes: The display ICC profile as a byte array, or None
            if Siril is running in headless mode.

        Raises:
            NoImageError: If no image is currently loaded,
            SirilError: If any other error occurs.
        """

        shm = None
        try:
            # Request shared memory setup
            status, response = self._send_command(_Command.GET_DISPLAY_ICC_PROFILE)

            # Handle error responses
            if status == _Status.ERROR:
                if response:
                    error_msg = response.decode('utf-8', errors='replace')
                    if "no image loaded" in error_msg.lower():
                        raise NoImageError(_("Error in get_siril_display_iccprofile(): no image is currently loaded in Siril"))
                    raise SirilConnectionError(_("Server error: {}").format(error_msg))
                raise SharedMemoryError(_("Failed to initiate shared memory transfer: Empty response"))

            if status == _Status.NONE:
                return None

            if not response:
                raise SharedMemoryError(_("Failed to initiate shared memory transfer: No data received"))

            if len(response) < 280: # Not a correct SharedMemoryInfo payload
                raise SharedMemoryError(_("Invalid shared memory information received (size < 280 bytes)"))
            try:
                # Parse the shared memory information
                shm_info = _SharedMemoryInfo.from_buffer_copy(response)
            except (AttributeError, BufferError, ValueError) as e:
                raise SharedMemoryError(_("Invalid shared memory information received: {}").format(e)) from e

            # Map the shared memory
            try:
                shm = self._map_shared_memory(
                    shm_info.shm_name.decode('utf-8'),
                    shm_info.size
                )
            except (OSError, ValueError) as e:
                raise SharedMemoryError(_("Failed to map shared memory: {}").format(e)) from e

            try:
                buffer = bytearray(shm.buf)[:shm_info.size]
                result = bytes(buffer)
            except (BufferError, ValueError, TypeError) as e:
                raise SirilError(_("Failed to create bytes from shared memory: {}").format(e)) from e

            return result

        except SirilError:
            # Re-raise NoImageError without wrapping
            raise
        except Exception as e:
            # Wrap all other exceptions with context
            raise SirilError(_("Error in get_siril_display_iccprofile(): {}").format(e)) from e

        finally:
            if shm is not None:
                try:
                    shm.close()  # First close the memory mapping as we have finished with it
                    # (We don't unlink it as C wll do that)

                    # Signal that Python is done with the shared memory and wait for C to finish
                    if not self._execute_command(_Command.RELEASE_SHM, shm_info):
                        raise SharedMemoryError(_("Failed to cleanup shared memory"))

                except Exception:
                    pass

    def get_siril_stf_linked(self) -> bool:

        """
        Determine whether the AutoStretch STF is configured to be linked or not

        Returns:
            A bool representing the current STF channel-linked state.

        Raises:
            SirilError: if an error occurred.
        """

        response = self._request_data(_Command.GET_STF_LINKED)
        if response is None:
            return None

        try:
            mode = struct.unpack('!I', response)
            return bool(mode[0])
        except struct.error as e:
            raise SirilError(_("Error occurred in get_siril_stf(): {}").format(e)) from e

    def set_siril_stf_linked(self, state: bool) -> bool:
        """
        Set the screen transfer function in Siril using the provided enum value.
        Args:
            state: bool

        Raises:
            SirilError: if an error occurred.
            ValueError: if parameters are not properly provided.

        Returns:
            bool: True if the slider state was set successfully
        """
        try:
            # Check for valid argument combinations

            # Validate that mode is a STFType enum
            if not isinstance(state, bool):
                raise ValueError("Mode must be True or False")

            # Pack the value into a byte
            # Format depends on which arguments are provided
            payload = struct.pack('?', state)
            self._execute_command(_Command.SET_STF_LINKED, payload)
            return True
        except Exception as e:
            raise SirilError(f"Error in set_siril_stf_linked(): {e}") from e

    def set_image_filename(self,
                        filename: str) -> bool:
        """
        Set the image filename in Siril.
        Args:
            filename (str): the image filename to set
        Raises:
            SirilError: if an error occurred.
            ValueError: if parameters are not properly provided.
        Returns:
            bool: True if the filename was set successfully
        """
        try:
            # Encode the filename as UTF-8 and pack it as bytes
            filename_bytes = filename.encode('utf-8')
            payload = filename_bytes
            self._execute_command(_Command.SET_IMAGE_FILENAME, payload)
            return True
        except Exception as e:
            raise SirilError(f"Error in set_image_filename(): {e}") from e

    def get_siril_log(self) -> str:
        """
        Retrieve the full Siril log as a text string.

        Returns:
            str: The text of the Siril log.

        Raises:
            SirilError: For errors during data retrieval,
        """

        shm = None
        try:
            # Request shared memory setup
            status, response = self._send_command(_Command.GET_SIRIL_LOG)

            # Handle error responses
            if status == _Status.ERROR:
                if response:
                    error_msg = response.decode('utf-8', errors='replace')
                    raise SirilError(_("Server error: {}").format(error_msg))
                raise SharedMemoryError(_("Failed to initiate shared memory transfer: Empty response"))

            if not response:
                raise SharedMemoryError(_("Failed to initiate shared memory transfer: No data received"))

            if status == _Status.NONE:
                return None

            if len(response) < 25: # No payload
                return None

            try:
                # Parse the shared memory information
                shm_info = _SharedMemoryInfo.from_buffer_copy(response)
            except (AttributeError, BufferError, ValueError) as e:
                raise SharedMemoryError(_("Invalid shared memory information received: {}").format(e)) from e

            # Map the shared memory
            try:
                shm = self._map_shared_memory(
                    shm_info.shm_name.decode('utf-8'),
                    shm_info.size
                )
            except (OSError, ValueError) as e:
                raise SharedMemoryError(_("Failed to map shared memory: {}").format(e)) from e

            try:
                # Read entire buffer at once
                buffer = bytearray(shm.buf)[:shm_info.size]
                result = buffer.decode('utf-8', errors='ignore')
            except (BufferError, ValueError, TypeError) as e:
                raise SirilError(_("Failed to create string from shared memory: {}").format(e)) from e

            return result

        except SirilError:
            raise

        except Exception as e:
            # Wrap all other exceptions with context
            raise SirilError(_("Error in get_siril_log(): {}").format(e)) from e
        finally:
            if shm is not None:
                try:
                    shm.close()  # First close the memory mapping as we have finished with it
                    # (We don't unlink it as C wll do that)

                    # Signal that Python is done with the shared memory and wait for C to finish
                    if not self._execute_command(_Command.RELEASE_SHM, shm_info):
                        raise SharedMemoryError(_("Failed to cleanup shared memory"))

                except Exception:
                    pass

    def save_image_file(self, data, header=None, filename=None) -> bool:
        """
        Save image pixeldata and metadata to a FITS file. This uses Siril to
        save the image and can therefore be used to avoid a script dependency on
        astropy if it is only required for saving an image. This allows images to
        be saved directly from python without having to use set_image_pixeldata()
        and set_image_metadata_from_header_string() and saving the Siril image. It
        support processing multiple images without affecting the image currently
        loaded in Siril.

        Args:
            data: Either a numpy.ndarray containing the image data (must be 2D or 3D
                array with dtype float32 or uint16), OR a FFit object containing
                both data and header.
            header: string containing the FITS header data. Required if data is a
                    numpy array, ignored if data is a FFit object.
            filename: string containing the path where the file should be saved.
                    Required if data is a numpy array. If data is a FFit object
                    and filename is None, will use fit.filename.

        Returns:
            bool: True if successful, False otherwise

        Raises:
            ValueError: if the input array or header is invalid,
            TypeError: if invalid parameter types are provided,
            SirilError: if there was an error in handling the command.

        Examples:
            # Using numpy array and header string
            siril.save_image_file(data_array, header_string, "output.fit")

            # Using FFit object
            siril.save_image_file(fit, filename="output.fit")

            # Using FFit object with its own filename
            siril.save_image_file(fit)
        """

        shm_data = None
        shm_header = None
        shm_data_info = None
        shm_header_info = None

        try:
            # Check if data is a FFit object
            is_ffit = hasattr(data, 'data') and hasattr(data, 'header')

            if is_ffit:
                # Extract data and header from FFit object
                image_data = data.data
                header_str = data.header
                # Use FFit's filename if none provided
                if filename is None:
                    if hasattr(data, 'filename') and data.filename:
                        filename = data.filename
                    else:
                        raise ValueError(_("No filename provided and FFit object has no filename"))
            else:
                # data is a numpy array
                image_data = data
                header_str = header

                if header is None:
                    raise ValueError(_("Header must be provided when data is a numpy array"))
                if filename is None:
                    raise ValueError(_("Filename must be provided when data is a numpy array"))

            # Validate input array
            if not isinstance(image_data, np.ndarray):
                raise ValueError(_("Image data must be a numpy array"))

            if image_data.ndim not in (2, 3):
                raise ValueError(_("Image must be 2D or 3D array"))

            if image_data.dtype not in (np.float32, np.uint16):
                raise ValueError(_("Image data must be float32 or uint16"))

            if not isinstance(header_str, str):
                raise TypeError(_("Header data must be a string"))

            if not isinstance(filename, str):
                raise TypeError(_("Filename must be a string"))

            # Get image dimensions
            if image_data.ndim == 2:
                height, width = image_data.shape
                channels = 1
            else:
                channels, height, width = image_data.shape

            if channels > 3:
                raise ValueError(_("Image cannot have more than 3 channels"))

            if any(dim <= 0 for dim in (width, height)):
                raise ValueError(_("Invalid image dimensions: {}x{}").format(width, height))

            # Calculate image data size
            element_size = 4 if image_data.dtype == np.float32 else 2
            image_bytes = width * height * channels * element_size

            # Prepare header data
            header_bytes = header_str.encode('utf-8')
            header_size = len(header_bytes) + 1

            # Calculate total payload size
            total_bytes = image_bytes + header_size

            # Create shared memory for image data
            try:
                shm_data_info = self._request_shm(image_bytes)
                shm_data = self._map_shared_memory(
                    shm_data_info.shm_name.decode('utf-8'),
                    shm_data_info.size
                )
            except (OSError, ValueError) as e:
                raise SharedMemoryError(_("Failed to map shared memory for image data: {}").format(e)) from e

            # Create shared memory for header
            try:
                shm_header_info = self._request_shm(header_size)
                shm_header = self._map_shared_memory(
                    shm_header_info.shm_name.decode('utf-8'),
                    shm_header_info.size
                )
            except (OSError, ValueError) as e:
                raise SharedMemoryError(_("Failed to map shared memory for header: {}").format(e)) from e

            # Copy image data to shared memory
            try:
                buffer = memoryview(shm_data.buf).cast('B')
                shared_array = np.frombuffer(buffer, dtype=image_data.dtype).reshape(image_data.shape)
                np.copyto(shared_array, image_data)
                del buffer
                del shared_array
            except Exception as e:
                raise SharedMemoryError(_("Failed to copy image data to shared memory: {}").format(e)) from e

            # Copy header to shared memory
            try:
                buffer = memoryview(shm_header.buf).cast('B')
                buffer[:len(header_bytes)] = header_bytes
                del buffer
            except Exception as e:
                raise SharedMemoryError(_("Failed to copy header to shared memory: {}").format(e)) from e

            # Encode filename
            filename_bytes = filename.encode('utf-8')
            if len(filename_bytes) > 255:
                raise ValueError(_("Filename too long (max 255 bytes)"))

            # Pack the save image info structure
            # Format: width, height, channels, data_type, image_size, image_shm_name,
            #         header_size, header_shm_name, filename
            info = struct.pack(
                '!IIIIQ256sQ256s256s',
                width,
                height,
                channels,
                1 if image_data.dtype == np.float32 else 0,
                image_bytes,
                shm_data_info.shm_name,
                header_size,
                shm_header_info.shm_name,
                filename_bytes + b'\0' * (256 - len(filename_bytes))
            )

            # Send command using the existing _execute_command method
            if not self._execute_command(_Command.SAVE_IMAGE_FILE, info):
                raise SirilError(_("Failed to send save image file command"))

            return True

        except (SirilError, ValueError, TypeError):
            raise
        except Exception as e:
            raise SirilError(f"Error in save_image_file(): {e}") from e

        finally:
            # Cleanup shared memory for image data
            if shm_data is not None:
                try:
                    shm_data.close()
                    self._execute_command(_Command.RELEASE_SHM, shm_data_info)
                except Exception:
                    pass

            # Cleanup shared memory for header
            if shm_header is not None:
                try:
                    shm_header.close()
                    self._execute_command(_Command.RELEASE_SHM, shm_header_info)
                except Exception:
                    pass
