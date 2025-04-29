# Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
# Copyright (C) 2012-2025 team free-astro (see more in AUTHORS file)
# Reference site is https://siril.org
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Connection module for Siril, providing the ability to connect to a running
Siril instance and communicate with it. Includes an extensive range
of methods that can be used to get and set data from / to Siril.
"""

import os
import sys
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
from .exceptions import SirilError, DataError, SirilConnectionError, CommandError, NoImageError, NoSequenceError, SharedMemoryError
from .models import ImageStats, FKeywords, FFit, PSFStar, BGSample, RegData, ImgData, DistoData, Sequence, SequenceType, Polygon
from .enums import _Command, _Status, CommandStatus, _ConfigType, LogColor, SirilVport
from .utility import truncate_utf8

DEFAULT_TIMEOUT = 5.

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

    def connect(self) -> bool:
        """
        Establish a connection to Siril based on the pipe or socket path.

        Returns:
            True on success

        Raises:
            SirilConnectionError: if a connection error occurred
        """

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
        try:
            self.reset_progress()
        except Exception:
            print("Warning: unable to reset progress bar in disconnect()", file=sys.stderr)

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
                SirilInterface._connected = False
                return
            raise SirilError(_("No pipe connection to close"))
        if hasattr(self, 'sock'):
            self.sock.close()
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
                    raise SirilError(_("Interface error: {}").format(error_msg))
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
                    raise SirilConnectionError(_("Interface error: {}").format(error_msg))
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

        except SirilConnectionError as re:
            # Let specific SharedMemoryErrors propagate
            raise SharedMemoryError("Error creating shared memory allocation") from re

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

    def _claim_thread(self) -> bool:
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
        * If the result is False, alert the user and await further input: the
          thread is already in use, or an image processing dialog is open.
        * If the result is True, you have the thread claimed.
        * Now you can call ``SirilInterface.get_image()`` or ``get_image_pixeldata()``
        * Carry out your image processing
        * Call ``SirilInterface.set_image_pixeldata()``
        * Call ``SirilInterface.release_thread()`` by exiting image_lock() context

        Note that the thread should only be claimed when the script itself is
        operating on the Siril image data. If the script is calling a Siril command
        to alter the Siril image then the thread **must not** be claimed or the
        Siril command will be unable to acquire it, and will fail.

        Returns:
            bool: True if the processing thread is claimed successfully, or
                  False if the thread is in use or if an error occured. In either case
                  processing cannot continue, though the script can wait and allow the
                  user to try again once the thread is free.
        """
        try:
            retval = self._execute_command(_Command.CLAIM_THREAD, None)
            if retval is True:
                return True
            if retval is None:
                print(_("The processing thread is locked. Wait for the current "
                    "processing task to finish."))
                return False
            print(_("Error trying to claim the processing thread. Thread is "
                "in use or an image processing dialog is open."))
            return False

        except Exception as e:
            print(f"Error claiming processing thread: {e}", file=sys.stderr)
            return False

    def _release_thread(self) -> bool:
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
        * If the result is False, alert the user and await further input: the
          thread is already in use, or an image processing dialog is open.
        * If the result is True, you have the thread claimed.
        * Now you can call ``SirilInterface.get_image()`` or ``get_image_pixeldata()``
        * Carry out your image processing
        * Call ``SirilInterface.set_image_pixeldata()``
        * Call ``SirilInterface.release_thread()`` by exiting image_lock() context

        Returns:
            True if the thread was successfully released

        Raises:
            SirilError: if an error occurred in releasing the thread
        """
        try:
            retval = self._execute_command(_Command.RELEASE_THREAD, None)
            if retval is True:
                return True
            raise SirilError(_("Error trying to release the processing thread. "
                "It will be released when the script terminates."))

        except Exception as e:
            print(f"Error releasing the processing thread: {e}", file=sys.stderr)
            return False

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

            with siril.image_lock():
                # Get image data
                image_data = self.get_image_pixeldata()
                # Process image data
                processed_data = some_processing_function(image_data)
                # Set the processed image data
                siril.set_image_pixeldata(processed_data)

        Raises:
            SirilError: If the thread cannot be claimed.
        """
        class ImageLockContext:
            """ Class to manage the image_lock context """
            def __init__(self, outer_self: SirilInterface):
                self.outer_self = outer_self
                self.claimed = False

            def __enter__(self):
                if not self.outer_self._claim_thread():
                    raise SirilError("Error in image_lock(). Thread may be in use or an image processing dialog is open.")
                self.claimed = True
                return self.outer_self

            def __exit__(self, exc_type, exc_val, exc_tb):
                if self.claimed:
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
        messages will be truncated.

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
            return self._execute_command(_Command.UNDO_SAVE_STATE, message_bytes)

        except Exception as e:
            raise SirilError(f"Error in undo_save_state(): {e}") from e

    def update_progress(self, message: str, progress: float) -> bool:
        """
        Send a progress update to Siril with a message and completion percentage.

        Args:
            message: Status message to display,
            progress: Progress value in the range 0.0 to 1.0

        Raises:
            ValueError: If the progress argument is out of range,
            SirilError: For any other errors.
        """

        try:
            # Validate progress value
            if not 0.0 <= progress <= 1.0:
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
            raise SirilError(_("Error in cmd(): {e}")) from e

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
        channel: Optional[int] = None)-> Optional[PSFStar]:

        """
        Retrieves a PSFStar star model from the current selection in Siril.
        Only a single PSFStar is returned: if there are more than one in the
        selection, the first one identified by Siril's internal star detection
        algorithm is returned.

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
                if channel is None:
                    shape_data = struct.pack('!IIII', *shape)
                else:
                    shape_data = struct.pack('!IIIII', *shape, channel)

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

            format_string = '!13d2qdq16dqdd'  # Define the format string based on PSFStar structure

            # Extract the bytes for this struct and unpack
            values = struct.unpack(format_string, response)

            return PSFStar.deserialize(response)

        except (SirilError, ValueError):
            raise
        except Exception as e:
            raise SirilError(_("Error in get_selection_star(): {}").format(e)) from e

    def get_selection_stats(self, shape: Optional[list[int]] = None, \
        channel: Optional[int] = None) -> Optional[PSFStar]:

        """
        Retrieves statistics for the current selection in Siril.

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

    def get_image_pixeldata(self, shape: Optional[list[int]] = None, preview: Optional[bool] = False) -> Optional[np.ndarray]:

        """
        Retrieves the pixel data from the image currently loaded in Siril.

        Args:
            shape: Optional list of [x, y, w, h] specifying the region to retrieve.
                   If provided, gets pixeldata for just that region.
                   If None, gets pixeldata for the entire image.
            preview: optional bool specifying whether to get pixeldata as a preview
                     (i.e. 8-bit autostretched data) or as real image data. Defaults
                     to False (i.e. real image data)

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
            preview_data = struct.pack('!?', preview)
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
                    finish_info = struct.pack('256s', shm_info.shm_name)
                    if not self._execute_command(_Command.RELEASE_SHM, finish_info):
                        raise SirilError(_("Failed to cleanup shared memory"))

                except Exception:
                    pass

    def get_seq_frame_pixeldata(self, frame: int, shape: Optional[List[int]] = None, preview: Optional[bool] = False) -> Optional[np.ndarray]:

        """
        Retrieves the pixel data from a frame in the sequence currently loaded in Siril.

        Args:
            frame: selects the frame to retrieve pixel data from
            shape: Optional list of [x, y, w, h] specifying the region to retrieve.
                   If provided, gets pixeldata for just that region.
                   If None, gets pixeldata for the entire image.
            preview: optional bool specifying whether to get pixeldata as a preview
                     (i.e. 8-bit autostretched data) or as real image data. Defaults
                     to False (i.e. real image data).

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
            preview_payload = struct.pack('!?', preview)
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
                    finish_info = struct.pack('256s', shm_info.shm_name)
                    if not self._execute_command(_Command.RELEASE_SHM, finish_info):
                        raise SirilError(_("Failed to cleanup shared memory"))

                except Exception:
                    pass

    def xy_plot(self, plot_data: PlotData):
        """
        Serialize plot data and send via shared memory. See the sirilpy.plot submodule
        documentation for how to configure a PlotData object for use with SirilInterface.xy_plot()

        Args:
            plot_metadata: PlotMetadata object containing plot configuration

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
                0,  # width (not used for plots)
                0,  # height (not used for plots)
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

    def set_seq_frame_pixeldata(self, index: int, image_data: np.ndarray) -> bool:
        """
        Send image data to Siril using shared memory.

        Args:
            index: integer specifying which frame to set the pixeldata for.
            image_data: numpy.ndarray containing the image data.
                        Must be 2D (single channel) or 3D (multi-channel) array
                        with dtype either np.float32 or np.uint16.

        Raises:
            NoSequenceError: if no sequence is loaded in Siril,
            ValueError: if the input array is invalid,
            SirilError: if there was an error in handling the command.
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

            # Create payload
            index_bytes = struct.pack('!i', index)
            payload = index_bytes + info

            # Send command using the existing _execute_command method
            if not self._execute_command(_Command.SET_SEQ_FRAME_PIXELDATA, payload):
                raise SirilError(_("Failed to send set_seq_frame_pixeldata command"))

            return

        except (ValueError, SirilError):
            raise
        except Exception as e:
            raise SirilError("Error in set_seq_frame_pixeldata(): {e}") from e

        finally:
            if shm is not None:
                try:
                    shm.close()
                except Exception as e:
                    pass

    def get_image_iccprofile(self) -> Optional[bytes]:
        """
        Retrieve the ICC profile of the current Siril image using shared memory.

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
                    finish_info = struct.pack('256s', shm_info.shm_name)
                    if not self._execute_command(_Command.RELEASE_SHM, finish_info):
                        raise SharedMemoryError(_("Failed to cleanup shared memory"))

                except Exception:
                    pass

    def get_image_fits_header(self) -> Optional[str]:
        """
        Retrieve the full FITS header of the current image loaded in Siril.

        Args:
            none.

        Returns:
            bytes: The image FITS header as a string, or None if there is no header.

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

            return result

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
                    finish_info = struct.pack('256s', shm_info.shm_name)
                    if not self._execute_command(_Command.RELEASE_SHM, finish_info):
                        raise SharedMemoryError(_("Failed to cleanup shared memory"))

                except Exception:
                    pass

    def get_image_unknown_keys(self) -> Optional[str]:
        """
        Retrieve the unknown key in a FITS header of the current loaded Siril
        image using shared memory.

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
                    finish_info = struct.pack('256s', shm_info.shm_name)
                    if not self._execute_command(_Command.RELEASE_SHM, finish_info):
                        raise SharedMemoryError(_("Failed to cleanup shared memory"))

                except Exception:
                    pass

    def get_image_history(self) -> Optional[list[str]]:
        """
        Retrieve history entries in the FITS header of the current loaded
        Siril image using shared memory.

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
                    finish_info = struct.pack('256s', shm_info.shm_name)
                    if not self._execute_command(_Command.RELEASE_SHM, finish_info):
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
        Request the filename of the loaded image from Siril.

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
        Request image statistics from Siril for a specific channel.

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
                   data for (between 0 and Sequence.number),
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
                   data for (between 0 and Sequence.number)
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
                   metadata for (between 0 and Sequence.number)

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
        Request FITS keywords data from Siril as a FKeywords object.

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
        Request a copy of the current image open in Siril.

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

    def get_seq_frame(self, frame: int, with_pixels: Optional[bool] = True, preview: Optional[bool] = False) -> Optional[FFit]:
        """
        Request sequence frame as a FFit from Siril.

        Args:
            frame: Integer specifying which frame in the sequence to retrieve data for
                   (between 0 and Sequence.number)
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

        shm = None
        data_payload = struct.pack('!I??', frame, with_pixels, preview)
        response = self._request_data(_Command.GET_SEQ_IMAGE, data_payload, timeout = None)
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
                'Q',  # gboolean color_managed (padded to uint64_t)
                # Keywords (start at index 14)
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
            if with_pixels:
                # Starts at index 65
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

            # Create FFit object:
            # Get pixeldata if requested
            pixeldata = None
            if with_pixels:
                try:
                    shm_info = _SharedMemoryInfo(
                        size = values[65],
                        data_type = values[66],
                        width = values[67],
                        height = values[68],
                        channels = values[69],
                        shm_name = values[70]
                    )
                    # Validate dimensions
                    if any(dim <= 0 for dim in (shm_info.width, shm_info.height)):
                        raise DataError(_("Invalid image dimensions: {}x{}").format(shm_info.width, shm_info.height))

                    if shm_info.channels <= 0 or shm_info.channels > 3:
                        raise DataError(_("Invalid number of channels: {}").format(shm_info.channels))

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
                        raise DataError(
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
                        raise DataError(_("Error in get_seq_frame(): Failed to reshape array to image dimensions: {}").format(e)) from e

                    # Make a copy of the data since we'll be releasing the shared memory
                    pixeldata = np.copy(arr)

                except SirilError:
                    raise
                except Exception as e:
                    raise SirilError(f"Error obtaining pixeldata: {e}") from e

            fits_keywords = FKeywords(
                program=decode_string(values[13]),
                filename=decode_string(values[14]),
                row_order=decode_string(values[15]),
                filter=decode_string(values[16]),
                image_type=decode_string(values[17]),
                object=decode_string(values[18]),
                instrume=decode_string(values[19]),
                telescop=decode_string(values[20]),
                observer=decode_string(values[21]),
                sitelat_str=decode_string(values[22]),
                sitelong_str=decode_string(values[23]),
                bayer_pattern=decode_string(values[24]),
                focname=decode_string(values[25]),
                bscale=values[26],
                bzero=values[27],
                lo=values[28],
                hi=values[29],
                flo=values[30],
                fhi=values[31],
                data_max=values[32],
                data_min=values[33],
                pixel_size_x=values[34],
                pixel_size_y=values[35],
                binning_x=values[36],
                binning_y=values[37],
                expstart=values[38],
                expend=values[39],
                centalt=values[40],
                centaz=values[41],
                sitelat=values[42],
                sitelong=values[43],
                siteelev=values[44],
                bayer_xoffset=values[45],
                bayer_yoffset=values[46],
                airmass=values[47],
                focal_length=values[48],
                flength=values[49],
                iso_speed=values[50],
                exposure=values[51],
                aperture=values[52],
                ccd_temp=values[53],
                set_temp=values[54],
                livetime=values[55],
                stackcnt=values[56],
                cvf=values[57],
                gain=values[58],
                offset=values[59],
                focuspos=values[60],
                focussz=values[61],
                foctemp=values[62],
                date=timestamp_to_datetime(values[63]),
                date_obs=timestamp_to_datetime(values[64])
            )
            fit = FFit(
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
                _data = pixeldata,
                stats=[
                    self.get_seq_stats(frame, 0),
                    self.get_seq_stats(frame, 1) if values[2] > 1 else None,
                    self.get_seq_stats(frame, 2) if values[2] > 1 else None,
                ],
                keywords = fits_keywords,
                _icc_profile = None,
                header = None,
                unknown_keys = None,
                history = None
            )
            return fit

        except Exception as e:
            raise SirilError(f"Error in get_seq_frame(): {e}") from e

        finally:
            if shm is not None:
                try:
                    shm.close()  # First close the memory mapping as we have finished with it
                    # (We don't unlink it as C wll do that)

                    # Signal that Python is done with the shared memory and wait for C to finish
                    finish_info = struct.pack('256s', shm_info.shm_name)
                    if not self._execute_command(_Command.RELEASE_SHM, finish_info):
                        raise SirilError(_("Failed to cleanup shared memory"))

                except Exception:
                    pass

    def get_seq_frame_header(self, frame: int) -> Optional[str]:
        """
        Retrieve the full FITS header of an image from the sequence loaded in Siril.

        Args:
            frame: Integer specifying which frame in the sequence to retrieve data for
            (between 0 and Sequence.number)

        Returns:
            str: The image FITS header as a string, or None if there is no header.

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

            return result

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
                    finish_info = struct.pack('256s', shm_info.shm_name)
                    if not self._execute_command(_Command.RELEASE_SHM, finish_info):
                        raise SharedMemoryError(_("Failed to cleanup shared memory"))

                except Exception:
                    pass

    def get_image_stars(self) -> List[PSFStar]:
        """
        Request star model PSF data from Siril.

        Returns:
            List of PSFStar objects containing the star data, or None if
            no stars can be found. If stars have already been detected using
            the `findstar` command then this list will be returned, otherwise
            automatic star detection will be attempted with the current
            star finder settings.

        Raises:
            NoImageError: If no image is currently loaded,
            SirilError: For other errors during  data retrieval,
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
                        raise NoImageError(_("No image is currently loaded in Siril"))
                    raise SirilError(_("Server error: {}").format(error_msg))
                raise SharedMemoryError(_("Failed to initiate shared memory transfer: Empty response"))

            if status == _Status.NONE:
                return None

            if not response:
                raise SharedMemoryError(_("Failed to initiate shared memory transfer: No data received"))

            if status == _Status.NONE:
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

            format_string = '!13d2qdq16dqdd'  # Define the format string based on PSFStar structure
            fixed_size = struct.calcsize(format_string)

            # Read entire buffer at once using memoryview
            buffer = bytearray(shm.buf)[:shm_info.size]

            # Validate buffer size
            if shm_info.size % fixed_size != 0:
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
                    # (We don't unlink it as C wll do that)

                    # Signal that Python is done with the shared memory and wait for C to finish
                    finish_info = struct.pack('256s', shm_info.shm_name)
                    if not self._execute_command(_Command.RELEASE_SHM, finish_info):
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

    def set_seq_frame_incl(self, index: int, incl: bool):
        """
        Set whether a given frame is included in the currently loaded sequence
        in Siril. This method is intended for use in creating custom sequence
        filters.

        Args:
            index: integer specifying which frame to set the pixeldata for.
            incl: bool specifying whether the frame is included or not.

        Raises:
            NoSequenceError: if no sequence is loaded in Siril,
            SirilError: on failure.
        """
        try:
            if not self.is_image_loaded():
                raise NoSequenceError(_("Error in set_seq_frame_incl(): no sequence loaded"))
            # Pack the index and incl into bytes using network byte order (!)
            payload = struct.pack('!II', index, incl)
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
                    finish_info = struct.pack('256s', shm_info.shm_name)
                    if not self._execute_command(_Command.RELEASE_SHM, finish_info):
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

    def overlay_add_polygon(self, polygon: Polygon)-> Polygon:
        """
        Adds a user polygon to the Siril display overlay

        Args:
            polyon: Polygon defining the polygon to be added

        Returns:
            Polygon: the input updated with the ID assigned by Siril

        Raises:
            SirilConnectionError: on a connection failure,
            DataError: on receipt of invalid data,
            SirilError: on any failure
        """

        try:
            # Serialize the provided polygon
            polygon_bytes = polygon.serialize()
            # Send it using _execute_command
            response = self._request_data(_Command.ADD_USER_POLYGON, polygon_bytes)
            if response is None:
                raise SirilConnectionError(_("Failed to get a response from the SirilInterface"))

            try:
                # Assuming the response is in the format: !i (ID) (4 bytes)
                polygon_id = struct.unpack('!i', response)[0]
                polygon.polygon_id = polygon_id
                return polygon
            except struct.error as e:
                raise DataError(_("Error unpacking polygon ID")) from e

        except Exception as e:
            raise SirilError(f"Error in overlay_add_polygon(): {e}") from e

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
            raise(f"Error in overlay_clear_polygons(): {e}") from e

    def overlay_get_polygon(self, polygon_id: int) -> 'Polygon':
        """
        Gets a single user polygon from the Siril overlay, specified by ID

        Args:
            id: int specifying the polygon ID to be deleted

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
                    finish_info = struct.pack('256s', shm_info.shm_name)
                    if not self._execute_command(_Command.RELEASE_SHM, finish_info):
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
                    finish_info = struct.pack('256s', shm_info.shm_name)
                    if not self._execute_command(_Command.RELEASE_SHM, finish_info):
                        raise SirilError(_("Failed to cleanup shared memory"))

                except Exception:
                    pass
