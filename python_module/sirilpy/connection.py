# Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
# Copyright (C) 2012-2025 team free-astro (see more in AUTHORS file)
# Reference site is https://siril.org
# SPDX-License-Identifier: GPL-3.0-or-later

import os
import sys
import ctypes
import socket
import struct
import threading
from enum import IntEnum, unique
from datetime import datetime
from typing import Tuple, Optional, List, Union
import numpy as np
from .translations import _
from .shm import SharedMemoryWrapper
from .plot import PlotData, _PlotSerializer
from .exceptions import SirilError, SirilConnectionError, CommandError, NoImageError
from .models import ImageStats, FKeywords, FFit, Homography, PSFStar, BGSample, RegData, ImgData, DistoData, Sequence, SequenceType, SirilPoint, UserPolygon

DEFAULT_TIMEOUT = 5.

if os.name == 'nt':
    import win32file
    import win32event
    import pywintypes
    import winerror

@unique
class _Status(IntEnum):
    """
    Returns the status of a command. NONE is for commands that
    may legitimately fail to return data but which should not be
    regarded as an error, instead this triggers the command processor
    to return the special python value None
    Internal class: this is not intended for use in scripts.
    """

    OK = 0
    NONE = 1
    ERROR = 0xFF

@unique
class _Command(IntEnum):
    """
    Enumerates the commands. This enum MUST match the one in
    siril_pythonmodule.h. Internal class: this is not intended for
    use in scripts.
    """
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
    GET_USERCONFIGDIR = 26
    GET_IS_IMAGE_LOADED = 27
    GET_IS_SEQUENCE_LOADED = 28
    GET_SELECTION = 29
    SET_SELECTION = 30
    GET_ACTIVE_VPORT = 31
    GET_STAR_IN_SELECTION = 32
    GET_STATS_FOR_SELECTION = 33
    PIX2WCS = 34
    WCS2PIX = 35
    UNDO_SAVE_STATE = 36
    GET_BUNDLE_PATH = 37
    ERROR_MESSAGEBOX = 38
    ERROR_MESSAGEBOX_MODAL = 39
    SIRIL_PLOT = 40
    CLAIM_THREAD = 41
    RELEASE_THREAD = 42
    SET_SEQ_FRAME_PIXELDATA = 43
    REQUEST_SHM = 44
    SET_SEQ_FRAME_INCL = 45
    GET_USERDATADIR = 46
    GET_SYSTEMDATADIR = 47
    GET_BGSAMPLES = 48
    SET_BGSAMPLES = 49
    GET_SEQ_FRAME_FILENAME = 50
    INFO_MESSAGEBOX = 51
    INFO_MESSAGEBOX_MODAL = 52
    WARNING_MESSAGEBOX = 53
    WARNING_MESSAGEBOX_MODAL = 54
    GET_SEQ_DISTODATA = 55
    SET_IMAGE_HEADER = 56
    ADD_USER_POLYGON = 57
    DELETE_USER_POLYGON = 58
    CLEAR_USER_POLYGONS = 59
    GET_USER_POLYGON = 60
    GET_USER_POLYGON_LIST = 61
    CONFIRM_MESSAGEBOX = 62
    ERROR = 0xFF

@unique
class _CommandStatus(IntEnum):
    CMD_NOT_FOUND = 1
    CMD_NO_WAIT = 1 << 1
    CMD_NO_CWD = 1 << 2
    CMD_NOT_SCRIPTABLE = 1 << 3
    CMD_WRONG_N_ARG = 1 << 4
    CMD_ARG_ERROR = 1 << 5
    CMD_SELECTION_ERROR = 1 << 6
    CMD_OK = 0
    CMD_GENERIC_ERROR = 1 << 7
    CMD_IMAGE_NOT_FOUND = 1 << 8
    CMD_SEQUENCE_NOT_FOUND = 1 << 9
    CMD_INVALID_IMAGE = 1 << 10
    CMD_LOAD_IMAGE_FIRST = 1 << 11
    CMD_ONLY_SINGLE_IMAGE = 1 << 12
    CMD_NOT_FOR_SINGLE = 1 << 13
    CMD_NOT_FOR_MONO = 1 << 14
    CMD_NOT_FOR_RGB = 1 << 15
    CMD_FOR_CFA_IMAGE = 1 << 16
    CMD_FILE_NOT_FOUND = 1 << 17
    CMD_FOR_PLATE_SOLVED = 1 << 18
    CMD_NEED_INIT_FIRST = 1 << 19
    CMD_ALLOC_ERROR = 1 << 20
    CMD_THREAD_RUNNING = 1 << 21
    CMD_DIR_NOT_FOUND = 1 << 22

@unique
class LogColor (IntEnum):
    """
    Defines colors available for use with ``SirilInterface.log()``
    For consistency ``LogColor.Default`` should be used for normal messages,
    ``LogColor.Red`` should be used for error messages, ``LogColor.Salmon``
    should be used for warning messages, LogColor.Green should  be used
    for completion notifications, and ``LogColor.Blue`` should be used for
    technical messages such as equations, coefficients etc.
    """
    DEFAULT = 0
    RED = 1
    SALMON = 2
    GREEN = 3
    BLUE = 4

class _Defaults:
    """
    Contains default values for different datatypes, matching Siril
    """
    DEFAULT_DOUBLE_VALUE = -999.0
    DEFAULT_FLOAT_VALUE = -999.0
    DEFAULT_INT_VALUE = -2147483647
    DEFAULT_UINT_VALUE = 2147483647
    VALUES = {DEFAULT_DOUBLE_VALUE, DEFAULT_FLOAT_VALUE, DEFAULT_INT_VALUE, DEFAULT_UINT_VALUE}

@unique
class _ConfigType(IntEnum):
    """
    Enumerates config variable types for use with the
    ``get_siril_config()`` method. Internal class: this is not intended
    for use in scripts.
    """
    BOOL = 0
    INT = 1
    DOUBLE = 2
    STR = 3
    STRDIR = 4
    STRLIST = 5

class _SharedMemoryInfo(ctypes.Structure):
    """
    Structure matching the C-side shared memory info. Internal class:
    this is not intended for use in scripts.
    """
    _fields_ = [
        ("size", ctypes.c_size_t),
        ("data_type", ctypes.c_int),  # 0 for WORD, 1 for float
        ("width", ctypes.c_int),
        ("height", ctypes.c_int),
        ("channels", ctypes.c_int),
        ("shm_name", ctypes.c_char * 256)
    ]

class SirilInterface:
    """
    SirilInterface is the main class providing an interface to a running
    Siril instance and access to methods to interact with it through
    Siril's inbuilt command system and accessing image and sequence data.
    """

    def __init__(self):
        """
        Initialize the SirilInterface, automatically determining the
        correct pipe or socket path based on the environment variable and
        operating system. Internal method.
        """
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

    def connect(self) -> Optional[bool]:
        """
        Establish a connection to Siril based on the pipe or socket path.

        Returns:
            True if the connection is successful, otherwise False.

        Raises:
            SirilConnectionError: if a connection error occurred
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
                    if self.debug:
                        current_pid = os.getpid()
                        print(f'Current ProcessID is {current_pid}')
                        self.info_messagebox(f'Current ProcessID is {current_pid}', False)
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
                return True

        except Exception as e:
            raise SirilConnectionError(_("Failed to connect: {}").format(e)) from e

    def disconnect(self):
        """
        Closes the established socket or pipe connection.

        Returns:
            True if the connection is closed successfully.

        Raises:
            SirilConnectionError: if the connection cannot be closed because the
                             pipe / socket cannot be found.
        """

        if os.name == 'nt':
            if hasattr(self, 'pipe_handle'):
                # Close the pipe handle
                win32file.CloseHandle(self.pipe_handle)
                # Close the event handles
                if hasattr(self, 'overlap_read'):
                    win32file.CloseHandle(self.overlap_read.hEvent)
                if hasattr(self, 'overlap_write'):
                    win32file.CloseHandle(self.overlap_write.hEvent)
                return True
            raise SirilError(_("No pipe connection to close"))
        if hasattr(self, 'sock'):
            self.sock.close()
            return True
        raise SirilError(_("No socket connection to close"))

    def _recv_exact(self, n: int, timeout: Optional[float] = DEFAULT_TIMEOUT) -> Optional[bytes]:
        """
        Helper method to receive exactly n bytes from the socket or pipe.
        Internal method, not for direct use in scripts.

        Args:
            n: Number of bytes to receive
            timeout: Timeout in seconds. None for indefinite timeout.
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
        Send a command and receive response with optional timeout.

        Args:
            command: Command to send
            data: Optional data payload
            timeout: Timeout for receive operations. None for indefinite timeout.
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

        except Exception as e:
            raise CommandError(_("Error sending command: {}").format(e)) from e

    def _execute_command(self, command: _Command, payload: Optional[bytes] = None, timeout: Optional[float] = DEFAULT_TIMEOUT) -> bool:
        """
        High-level method to execute a command and handle the response.
        Internal method, not for end-user use.

        Args:
            command: The command to execute
            payload: Optional command payload
            timeout: Timeout for command execution. None for indefinite timeout.

        Returns:
            True if command was successful, False otherwise
        """
        if self.debug:
            timeout = None

        status, response = self._send_command(command, payload, timeout)

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

    def _request_data(self, command: _Command, payload: Optional[bytes] = None, timeout: Optional[float] = DEFAULT_TIMEOUT) -> Optional[bytes]:
        """
        High-level method to request small-volume data from Siril. The
        payload limit is 63336 bytes. For commands expected to return
        larger volumes of data, SHM should be used.
        Internal method, not for direct use in scripts.

        Args:
            command: The data request command
            payload: Optional request parameters
            timeout: Timeout for data request. None for indefinite timeout.

        Returns:
            Requested data or None if error
        """
        if self.debug:
            timeout = None

        status, response = self._send_command(command, payload, timeout)

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

    def _request_shm(self, size: int) -> Optional[_SharedMemoryInfo]:
        """
        Request Siril to create a shared memory allocation and return an
        object containing the details (name, size).

        Args:
            size: int specifying the size of shm buffer to create.

        Returns:
            _SharedMemoryInfo: Details of the shared memory allocation, or None if allocation failed.

        Raises:
            RuntimeError: If an invalid response is received or no data is returned.
            ValueError: If the shared memory information is invalid.
        """
        size_data = struct.pack('!Q', size)  # Pack the size as a uint64_t in network byte order.

        try:
            # Send the shared memory request command.
            status, response = self._send_command(_Command.REQUEST_SHM, size_data)

            # Check the status and handle accordingly.
            if status == _Status.ERROR:
                raise RuntimeError(_("Failed to initiate shared memory transfer: invalid response"))

            if status == _Status.NONE:
                # Return None if no shared memory was allocated.
                return None

            if not response:
                raise RuntimeError(_("Failed to initiate shared memory transfer: no data received"))

            try:
                # Attempt to parse the shared memory info from the response buffer.
                shm_info = _SharedMemoryInfo.from_buffer_copy(response)
                return shm_info

            except (AttributeError, BufferError, ValueError) as e:
                # Catch parsing errors and raise a descriptive exception.
                raise ValueError(_("Invalid shared memory information received: {}").format(e)) from e

        except RuntimeError as re:
            # Let specific RuntimeErrors propagate with contextual information.
            raise RuntimeError(
                _("Runtime error during shared memory allocation request: {}").format(re)
            ) from re

        except Exception as e:
            # Catch any unexpected errors to improve debuggability.
            raise RuntimeError(
                _("Unexpected error during shared memory allocation request: {}").format(e)
            ) from e

    def _map_shared_memory(self, name: str, size: int) -> SharedMemoryWrapper:
        """
        Create or open a shared memory mapping using SharedMemoryWrapper.
        Internal method, not for direct use in scripts.

        Args:
            name: Name of the shared memory segment,
            size: Size of the shared memory segment in bytes

        Returns:
            SharedMemoryWrapper: A wrapper object for the shared memory segment

        Raises:
            RuntimeError: If the shared memory mapping fails
        """
        try:
            return SharedMemoryWrapper(name=name, size=size)
        except Exception as e:
            raise RuntimeError(_("Failed to create shared memory mapping: {}").format(e)) from e

    def _get_bundle_path(self) ->Optional[str]:
        """
        Request the bundle path directory. This is an internal method used
        to ensure that the correct DLL paths are preconfigured on Windows:
        it is not for use by scriptwriters.

        **It is an error to call this method on non-Windows OSes**

        Returns:
            The Siril bundle path as a string.

        Raises:
            All exceptions are raised to the caller.
        """

        response = self._request_data(_Command.GET_BUNDLE_PATH)

        if response is None:
            raise RuntimeError("Failed to get bundle path - received null response")

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

        Returns:
            bool: True if the message was successfully logged, False otherwise
        """

        try:
            # Append a newline character to the string
            truncated_string = my_string[:1021] + '\n'
            # Convert string to bytes using UTF-8 encoding
            message_bytes = truncated_string.encode('utf-8')
            # Prepend the color byte
            packed_message = bytes([color.value]) + message_bytes
            return self._execute_command(_Command.LOG_MESSAGE, packed_message)

        except Exception as e:
            print(f"Error sending log message: {e}", file=sys.stderr)
            return False

    def claim_thread(self) -> bool:
        """
        Claim the processing thread. This prevents other processes using the
        processing thread to operate on the current Siril image. The preferred
        method of thread control is to use the image_lock() context manager
        rather than using this function manually.

        This function **must** always be called before starting any processing
        that will end with ``SirilInterface.set_image_pixeldata()``. The
        sequence of operations should be:

        * Call ``SirilInterface.claim_thread()``
        * If the result is False, alert the user and await further input: the
          thread is already in use, or an image processing dialog is open.
        * If the result is True, you have the thread claimed.
        * Now you can call ``SirilInterface.get_image()`` or ``get_image_pixeldata()``
        * Carry out your image processing
        * Call ``SirilInterface.set_image_pixeldata()``
        * Call ``SirilInterface.release_thread()``

        As a precaution, the thread will be released automatically if it is still
        held at the point the script process terminates, but that should not be
        seen as an excuse for failing to call ``SirilInterface.release_thread()``

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

    def release_thread(self) -> bool:
        """
        Release the processing thread. This permits other processes to use the
        processing thread to operate on the current Siril image. The preferred
        method of thread control is to use the image_lock() context manager
        rather than using this function manually.

        This function **must** always be called after completing any processing
        that has updated the image loaded in Siril. The sequence of operations
        should be:

        * Call ``SirilInterface.claim_thread()``
        * If the result is False, alert the user and await further input: the
          thread is already in use, or an image processing dialog is open.
        * If the result is True, you have the thread claimed.
        * Now you can call ``SirilInterface.get_image()`` or ``get_image_pixeldata()``
        * Carry out your image processing
        * Call ``SirilInterface.set_image_pixeldata()``
        * Call ``SirilInterface.release_thread()``

        As a precaution, the thread will be released automatically if it is still
        held at the point the script process terminates, but that should not be
        seen as an excuse for failing to call this method.

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

        Example usage:

        .. code-block:: python

            with self.image_lock():
                # Get image data
                image_data = self.get_image_pixeldata()
                # Process image data
                processed_data = some_processing_function(image_data)
                # Set the processed image data
                self.set_image_pixeldata(processed_data)

        Raises:
            RuntimeError: If the thread cannot be claimed.
        """
        class ImageLockContext:
            def __init__(self, outer_self):
                self.outer_self = outer_self
                self.claimed = False

            def __enter__(self):
                if not self.outer_self.claim_thread():
                    raise RuntimeError("Failed to claim processing thread. Thread may be in use or an image processing dialog is open.")
                self.claimed = True
                return self.outer_self

            def __exit__(self, exc_type, exc_val, exc_tb):
                if self.claimed:
                    self.outer_self.release_thread()
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
            RuntimeError: if an error occurred.
        """
        try:
            # Truncate strings to allowed lengths
            truncated_title = title[:256]
            truncated_message = message[:1021]
            truncated_label = confirm_label[:24]

            # Encode strings and add null terminators
            encoded_title = truncated_title.encode('utf-8') + b'\0'
            encoded_message = truncated_message.encode('utf-8') + b'\0'
            encoded_label = truncated_label.encode('utf-8') + b'\0'

            # Concatenate into one buffer
            message_bytes = encoded_title + encoded_message + encoded_label

            # Call the command with the encoded data
            response = self._request_data(_Command.CONFIRM_MESSAGEBOX, message_bytes, timeout=None)

            if response is None:
                raise SirilError(_("Error sending confirm_messagebox command"))

            return bool(int.from_bytes(response, byteorder='little'))

        except Exception as e:
            print(f"Error sending confirmation message: {e}", file=sys.stderr)
            return False

    def _messagebox(self, my_string: str, cmd_type: int, modal: Optional[bool] = False) -> bool:
        """
        Send a message to Siril for display in a messagebox.
        Helper method for error_messagebox, warning_messagebox, info_messagebox.

        Args:
            my_string: The message to display in the message box
            type: Sets whether to show an error, warning or info messagebox
            modal: Whether or not the message box is modal

        Returns:
            bool: True if the error was successfully displayed, False otherwise
        """

        try:
            # Append a newline character to the string
            truncated_string = my_string[:1021] + '\n'
            # Convert string to bytes using UTF-8 encoding
            message_bytes = truncated_string.encode('utf-8')
            if modal:
                return self._execute_command(cmd_type, message_bytes, timeout = None)
            return self._execute_command(cmd_type, message_bytes)

        except Exception as e:
            print(f"Error sending log message: {e}", file=sys.stderr)
            return False

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
        """
        cmd_type = _Command.ERROR_MESSAGEBOX_MODAL if modal else _Command.ERROR_MESSAGEBOX
        return self._messagebox(my_string, cmd_type, modal)

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
        """

        cmd_type = _Command.INFO_MESSAGEBOX_MODAL if modal else _Command.INFO_MESSAGEBOX
        return self._messagebox(my_string, cmd_type, modal)

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
        """

        cmd_type = _Command.WARNING_MESSAGEBOX_MODAL if modal else _Command.WARNING_MESSAGEBOX
        return self._messagebox(my_string, cmd_type, modal)


    def undo_save_state(self, my_string: str) -> bool:
        """
        Saves an undo state. The maximum message length is 70 bytes: longer
        messages will be truncated.

        Args:
            my_string: The message to log in FITS HISTORY

        Returns:
            bool: True if the message was successfully logged, False otherwise
        """

        try:
            # Append a newline character to the string
            # Convert string to bytes using UTF-8 encoding
            message_bytes = my_string.encode('utf-8')
            return self._execute_command(_Command.UNDO_SAVE_STATE, message_bytes)

        except Exception as e:
            print(f"Error saving undo state: {e}", file=sys.stderr)
            return False

    def update_progress(self, message: str, progress: float) -> bool:
        """
        Send a progress update to Siril with a message and completion percentage.

        Args:
            message: Status message to display,
            progress: Progress value in the range 0.0 to 1.0

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

            return self._execute_command(_Command.UPDATE_PROGRESS, payload)

        except Exception as e:
            print(f"Error sending progress update: {e}", file=sys.stderr)
            return False

    def reset_progress(self) -> bool:
        """
        Resets the Siril progress bar.

        Args:
            none

        Returns:
            bool: True if the progress update was successfully sent, False otherwise
        """

        return self.update_progress("", 0.0)

    def cmd(self, *args: str):
        """
        Send a command to Siril to be executed. The range of available commands can
        be found by checking the online documentation. The command and its arguments
        are provided as a list of strings.

        Args:
            *args: Variable number of string arguments to be combined into a command

        Raises:
            CommandError: If the command fails with a specific error code.
            SirilError: If another error occurs during execution.

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
                raise SirilError(_(f"Error: _request_data({args}) failed."))

            # Convert response bytes to integer from network byte order
            if len(response) == 4:  # Valid response is int32_t ie 4 bytes
                status_code = int.from_bytes(response, byteorder='big')

                # Check against _CommandStatus enum
                if status_code == _CommandStatus.CMD_OK or status_code == _CommandStatus.CMD_NO_WAIT:
                    return  # Command executed successfully
                # ERROR HANDLING
                # Map status code to error message
                error_messages = {
                    _CommandStatus.CMD_NOT_FOUND: "Command not found",
                    _CommandStatus.CMD_NO_WAIT: "Command does not wait for completion",
                    _CommandStatus.CMD_NO_CWD: "Current working directory not set",
                    _CommandStatus.CMD_NOT_SCRIPTABLE: "Command not scriptable",
                    _CommandStatus.CMD_WRONG_N_ARG: "Wrong number of arguments",
                    _CommandStatus.CMD_ARG_ERROR: "Argument error",
                    _CommandStatus.CMD_SELECTION_ERROR: "Selection error",
                    _CommandStatus.CMD_GENERIC_ERROR: "Generic error",
                    _CommandStatus.CMD_IMAGE_NOT_FOUND: "Image not found",
                    _CommandStatus.CMD_SEQUENCE_NOT_FOUND: "Sequence not found",
                    _CommandStatus.CMD_INVALID_IMAGE: "Invalid image",
                    _CommandStatus.CMD_LOAD_IMAGE_FIRST: "Load image first",
                    _CommandStatus.CMD_ONLY_SINGLE_IMAGE: "Command requires a single image to be loaded",
                    _CommandStatus.CMD_NOT_FOR_SINGLE: "Command not for single images",
                    _CommandStatus.CMD_NOT_FOR_MONO: "Command not for monochrome images",
                    _CommandStatus.CMD_NOT_FOR_RGB: "Command not for RGB images",
                    _CommandStatus.CMD_FOR_CFA_IMAGE: "Command only for CFA images",
                    _CommandStatus.CMD_FILE_NOT_FOUND: "File not found",
                    _CommandStatus.CMD_FOR_PLATE_SOLVED: "Command requires plate-solved image",
                    _CommandStatus.CMD_NEED_INIT_FIRST: "Initialization required first",
                    _CommandStatus.CMD_ALLOC_ERROR: "Memory allocation error",
                    _CommandStatus.CMD_THREAD_RUNNING: "Command thread already running",
                    _CommandStatus.CMD_DIR_NOT_FOUND: "Directory not found"
                }
                print(f"Status code: {status_code}")
                error_message = error_messages.get(status_code, f"Unknown error code: {status_code}")
                raise CommandError(_(f"Command '{args[0]}' failed: {error_message}"), status_code)
            else:
                # Handle case where response doesn't contain enough bytes for a status code
                raise SirilError(_(f"Error: Response from {args[0]} incorrect size to contain a status code."))

        except Exception as e:
            raise  # Re-raise without wrapping

    def set_siril_selection(self, x: int, y: int, w: int, h: int) -> bool:
        """
        Set the image selection in Siril using the provided coordinates and dimensions.

        Args:
            x: X-coordinate of the selection's top-left corner
            y: Y-coordinate of the selection's top-left corner
            w: Width of the selection
            h: Height of the selection

        Returns:
            bool: True if the selection was successfully set, False otherwise
        """
        try:
            # Pack the coordinates and dimensions into bytes using network byte order (!)
            payload = struct.pack('!IIII', x, y, w, h)
            return self._execute_command(_Command.SET_SELECTION, payload)
        except Exception as e:
            print(f"Error setting selection: {e}", file=sys.stderr)
            return False

    def get_siril_selection(self) -> Optional[Tuple[int, int, int, int]]:

        """
        Request the image selection from Siril.

        Returns:
            A tuple (height, width, channels) representing the current selection,
            or None if an error occurred.
        """

        response = self._request_data(_Command.GET_SELECTION)

        if response is None:
            return None

        try:
            # Assuming the response is in the format: x (4 bytes), y (4 bytes), w (4 bytes), h (4 bytes)
            x, y, w, h = struct.unpack('!IIII', response)
            return x, y, w, h  # Returning as (x, y, w, h)
        except struct.error as e:
            print(f"Error unpacking image selection: {e}", file=sys.stderr)
            return None

    def get_siril_active_vport(self) -> Optional[int]:

        """
        Request the active viewport from Siril.

        Returns:
            An int representing the active vport:
            0 = Red (or Mono),
            1 = Green,
            2 = Blue,
            3 = RGB,
            or None if an error occurred.
        """

        response = self._request_data(_Command.GET_ACTIVE_VPORT)

        if response is None:
            return None

        try:
            # Assuming the response is in the format: !I
            vport = struct.unpack('!I', response)[0]
            return vport
        except struct.error as e:
            print(f"Error unpacking data: {e}", file=sys.stderr)
            return None

    def get_image_shape(self) -> Optional[Tuple[int, int, int]]:

        """
        Request the shape of the image from Siril.

        Returns:
            A tuple (channels, height, width) representing the shape of the image,
            or None if an error occurred.
        """

        response = self._request_data(_Command.GET_DIMENSIONS)

        if response is None:
            return None

        try:
            # Assuming the response is in the format: width (4 bytes), height (4 bytes), nb_channels (4 bytes)
            width, height, channels = struct.unpack('!III', response)
            return channels, height, width  # Returning as (channels, height, width)
        except struct.error as e:
            print(f"Error unpacking image dimensions: {e}", file=sys.stderr)
            return None

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
            PSFStar: the PSFStar object representing the star model.

        Raises:
            NoImageError: If no image is currently loaded,
            RuntimeError: For other errors during data retrieval,
            ValueError: If the received data format is invalid or no selection can
                        be determined
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
                if len(shape) == 4:
                    if channel is None:
                        shape_data = struct.pack('!IIII', *shape)
                    else:
                        shape_data = struct.pack('!IIIII', *shape, channel)
                shape_data = None

            status, response = self._send_command(_Command.GET_STAR_IN_SELECTION, shape_data)

            # Handle error responses
            if status == _Status.ERROR:
                if response:
                    error_msg = response.decode('utf-8', errors='replace')
                    if "no image loaded" in error_msg.lower():
                        raise NoImageError(_("No image is currently loaded in Siril"))
                    if "invalid selection" in error_msg.lower():
                        raise ValueError(_("No selection is currently made in Siril"))
                    raise RuntimeError(_("Server error: {}").format(error_msg))
                raise RuntimeError(_("Failed to transfer star data: Empty response"))

            if status == _Status.NONE:
                return None

            if not response:
                raise RuntimeError(_("Failed to transfer star data: No data received"))

            format_string = '!13d2qdq16dqdd'  # Define the format string based on PSFStar structure

            # Extract the bytes for this struct and unpack
            values = struct.unpack(format_string, response)

            try:
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
            except struct.error as e:
                print(f"Error unpacking star data: {e}", file=sys.stderr)
                raise RuntimeError(_("Error processing star data: {}").format(e)) from e

            return star

        except Exception as e:
            raise RuntimeError(_("Error processing star data: {}").format(e)) from e

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
            NoImageError: If no image is currently loaded,
            RuntimeError: For other errors during data retrieval,
            ValueError: If the received data format is invalid or no selection can
                        be determined
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
                if len(shape) == 4:
                    if channel is None:
                        shape_data = struct.pack('!IIII', *shape)
                    else:
                        shape_data = struct.pack('!IIIII', *shape, channel)
                else:
                    raise RuntimeError(_("Incorrect shape data provided: must be (x,y,w,h)"))

            status, response = self._send_command(_Command.GET_STATS_FOR_SELECTION, shape_data)

            # Handle error responses
            if status == _Status.ERROR:
                if response:
                    error_msg = response.decode('utf-8', errors='replace')
                    if "no image loaded" in error_msg.lower():
                        raise NoImageError(_("No image is currently loaded in Siril"))
                    if "invalid selection" in error_msg.lower():
                        raise ValueError(_("No selection is currently made in Siril"))
                    raise RuntimeError(_("Server error: {}").format(error_msg))
                raise RuntimeError(_("Failed to transfer stats data: Empty response"))

            if status == _Status.NONE:
                return None

            if not response:
                raise RuntimeError(_("Failed to transfer stats data: No data received"))
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

            except struct.error as e:
                print(f"Error unpacking star data: {e}", file=sys.stderr)
                raise RuntimeError(_("Error processing star data: {}").format(e)) from e
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
        except Exception as e:
            raise RuntimeError(_("Failed to transfer stats data: error occurred")) from e

    def pix2radec(self, x: float, y: float) -> Optional[Tuple[float, float]]:
        """
        Converts a pair of pixel coordinates into RA and dec coordinates using the
        WCS of the image loaded in Siril. This requires that an image is loaded in
        Siril and that it has been platesolved (i.e. it has a WCS solution).

        Args:
            x: float: provides the x coordinate to be converted
            y: float: provides the y coordinate to be converted

        Returns:
            Tuple[float, float]: [RA, dec] as a Tuple of two floats.

        Raises:
            NoImageError: If no image is currently loaded,
            RuntimeError: For other errors during pixel data retrieval,
            ValueError: If the received data format is invalid or no WCS is found
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
                    raise RuntimeError(_("Server error: {}").format(error_msg))
                raise RuntimeError(_("Failed to transfer coordinates: Empty response"))

            if status == _Status.NONE:
                return None

            if not response:
                raise RuntimeError(_("Failed to transfer coordinates: No data received"))
            try:
                # Define the format string for unpacking the C struct
                # '!' for network byte order (big-endian)
                # 'd' for double (all floating point values)
                format_string = '!2d'  # '!' ensures network byte order

                # Calculate expected size
                expected_size = struct.calcsize(format_string)

                # Verify we got the expected amount of data
                if len(response) != expected_size:
                    print(f"Received data size {len(response)} doesn't match expected size {expected_size}",
                        file=sys.stderr)
                    return None

                # Unpack the binary data
                values = struct.unpack(format_string, response)
            except struct.error as e:
                print(f"Error unpacking data: {e}", file=sys.stderr)
                raise RuntimeError(_("Error processing data: {}").format(e)) from e
            return values

        except Exception as e:
            raise RuntimeError(_("Failed to transfer stats data: error occurred")) from e

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
            NoImageError: If no image is currently loaded,
            RuntimeError: For other errors during pixel data retrieval,
            ValueError: If the received data format is invalid or no WCS is found
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
                    raise RuntimeError(_("Server error: {}").format(error_msg))
                raise RuntimeError(_("Failed to transfer coordinates: Empty response"))

            if status == _Status.NONE:
                return None

            if not response:
                raise RuntimeError(_("Failed to transfer coordinates: No data received"))
            try:
                # Define the format string for unpacking the C struct
                # '!' for network byte order (big-endian)
                # 'd' for double (all floating point values)
                format_string = '!2d'  # '!' ensures network byte order

                # Calculate expected size
                expected_size = struct.calcsize(format_string)

                # Verify we got the expected amount of data
                if len(response) != expected_size:
                    print(f"Received data size {len(response)} doesn't match expected size {expected_size}",
                        file=sys.stderr)
                    return None

                # Unpack the binary data
                values = struct.unpack(format_string, response)
            except struct.error as e:
                print(f"Error unpacking data: {e}", file=sys.stderr)
                raise RuntimeError(_("Error processing data: {}").format(e)) from e
            return values
        except Exception as e:
            raise RuntimeError(_("Failed to transfer stats data: error occurred")) from e

    def get_image_pixeldata(self, shape: Optional[list[int]] = None) -> Optional[np.ndarray]:

        """
        Retrieves the pixel data from the image currently loaded in Siril.

        Args:
            shape: Optional list of [x, y, w, h] specifying the region to retrieve.
                   If provided, gets pixeldata for just that region.
                   If None, gets pixeldata for the entire image.

        Returns:
            numpy.ndarray: The image data as a numpy array

        Raises:
            NoImageError: If no image is currently loaded,
            RuntimeError: For other errors during pixel data retrieval,
            ValueError: If the received data format is invalid or shape is invalid
        """

        shm = None
        try:
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
            status, response = self._send_command(command, shape_data)

            # Handle error responses
            if status == _Command.ERROR:
                if response:
                    error_msg = response.decode('utf-8', errors='replace')
                    if "no image loaded" in error_msg.lower():
                        raise NoImageError(_("No image is currently loaded in Siril"))
                    raise RuntimeError(_("Server error: {}").format(error_msg))
                raise RuntimeError(_("Failed to initiate shared memory transfer: Empty response"))

            if not response:
                raise RuntimeError(_("Failed to initiate shared memory transfer: No data received"))

            try:
                # Parse the shared memory information
                shm_info = _SharedMemoryInfo.from_buffer_copy(response)
            except (AttributeError, BufferError, ValueError) as e:
                raise ValueError(_("Invalid shared memory information received: {}").format(e)) from e

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
                raise RuntimeError(_("Failed to map shared memory: {}").format(e)) from e

            buffer = bytearray(shm.buf)[:shm_info.size]
            # Create numpy array from shared memory
            dtype = np.float32 if shm_info.data_type == 1 else np.uint16
            try:
                arr = np.frombuffer(buffer, dtype=dtype)
            except (BufferError, ValueError, TypeError) as e:
                raise RuntimeError(_("Failed to create array from shared memory: {}").format(e)) from e

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
                raise ValueError(_("Failed to reshape array to image dimensions: {}").format(e)) from e

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
                    shm.close()  # First close the memory mapping as we have finished with it
                    # (We don't unlink it as C wll do that)

                    # Signal that Python is done with the shared memory and wait for C to finish
                    finish_info = struct.pack('256s', shm_info.shm_name)
                    if not self._execute_command(_Command.RELEASE_SHM, finish_info):
                        raise RuntimeError(_("Failed to cleanup shared memory"))

                except Exception:
                    pass

    def get_seq_frame_pixeldata(self, frame: int, shape: Optional[List[int]] = None) -> Optional[np.ndarray]:

        """
        Retrieves the pixel data from the image currently loaded in Siril.

        Args:
            frame: selects the frame to retrieve pixel data from

        Returns:
            numpy.ndarray: The image data as a numpy array

        Raises:
            NoImageError: If no image is currently loaded,
            RuntimeError: For other errors during pixel data retrieval,
            ValueError: If the received data format is invalid or shape is invalid
        """

        shm = None
        # Convert channel number to network byte order bytes

        try:
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
            payload = frame_payload + shape_data if shape_data else frame_payload
            # Request shared memory setup
            status, response = self._send_command(_Command.GET_SEQ_PIXELDATA, payload)

            # Handle error responses
            if status == _Command.ERROR:
                if response:
                    error_msg = response.decode('utf-8', errors='replace')
                    if "no image loaded" in error_msg.lower():
                        raise NoImageError(_("No image is currently loaded in Siril"))
                    raise RuntimeError(_("Server error: {}").format(error_msg))
                raise RuntimeError(_("Failed to initiate shared memory transfer: Empty response"))

            if not response:
                raise RuntimeError(_("Failed to initiate shared memory transfer: No data received"))

            try:
                # Parse the shared memory information
                shm_info = _SharedMemoryInfo.from_buffer_copy(response)
            except (AttributeError, BufferError, ValueError) as e:
                raise ValueError(_("Invalid shared memory information received: {}").format(e)) from e

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
                raise RuntimeError(_("Failed to map shared memory: {}").format(e)) from e

            buffer = bytearray(shm.buf)[:shm_info.size]
            # Create numpy array from shared memory
            dtype = np.float32 if shm_info.data_type == 1 else np.uint16
            try:
                arr = np.frombuffer(buffer, dtype=dtype)
            except (BufferError, ValueError, TypeError) as e:
                raise RuntimeError(_("Failed to create array from shared memory: {}").format(e)) from e

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
                raise ValueError(_("Failed to reshape array to image dimensions: {}").format(e)) from e

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
                    shm.close()  # First close the memory mapping as we have finished with it
                    # (We don't unlink it as C wll do that)

                    # Signal that Python is done with the shared memory and wait for C to finish
                    finish_info = struct.pack('256s', shm_info.shm_name)
                    if not self._execute_command(_Command.RELEASE_SHM, finish_info):
                        raise RuntimeError(_("Failed to cleanup shared memory"))

                except Exception:
                    pass

    def xy_plot(self, plot_data: PlotData):
        """
        Serialize plot data and send via shared memory. See the sirilpy.plot submodule
        documentation for how to configure a PlotData object for use with SirilInterface.xy_plot()

        Args:
            plot_metadata: PlotMetadata object containing plot configuration
        """
        try:
            serialized_data, total_bytes = _PlotSerializer._serialize_plot_data(plot_data)

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
                    raise RuntimeError(_("Failed to map shared memory: {}").format(e)) from e

            except Exception as e:
                raise RuntimeError(_("Failed to create shared memory: {}").format(e)) from e

            # Copy serialized data to shared memory
            try:
                buffer = memoryview(shm.buf).cast('B')
                buffer[:total_bytes] = serialized_data
                del buffer
            except Exception as e:
                print(f"Failed to copy data to shared memory: {e}", file=sys.stderr)
                return False

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
                raise RuntimeError(_("Failed to send pixel data command"))

            return True

        except Exception as e:
            print(f"Error sending plot data: {e}", file=sys.stderr)
            return False
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
        """
        try:
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
            format_string = '3d d 2d Q 2d I' * len(samples)

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
                    raise RuntimeError(_("Failed to map shared memory: {}").format(e)) from e
            except Exception as e:
                raise RuntimeError(_("Failed to create shared memory: {}").format(e)) from e

            # Copy serialized data to shared memory
            try:
                buffer = memoryview(shm.buf).cast('B')
                buffer[:total_bytes] = serialized_data
                del buffer
            except Exception as e:
                print(f"Failed to copy data to shared memory: {e}", file=sys.stderr)
                return False

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
                raise RuntimeError(_("Failed to send BG sample command"))
            return True

        except Exception as e:
            print(f"Error sending background samples: {e}", file=sys.stderr)
            return False
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

        Returns:
            bool: True if successful, False otherwise
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
                    raise RuntimeError(_("Failed to map shared memory: {}").format(e)) from e

            except Exception as e:
                raise RuntimeError(_("Failed to create shared memory: {}").format(e)) from e

            # Copy data to shared memory
            try:
                buffer = memoryview(shm.buf).cast('B')
                shared_array = np.frombuffer(buffer, dtype=image_data.dtype).reshape(image_data.shape)
                np.copyto(shared_array, image_data)
                # Delete transient objects used to structure copy
                del buffer
                del shared_array
            except Exception as e:
                raise RuntimeError(_("Failed to copy data to shared memory: {}").format(e)) from e

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
                raise RuntimeError(_("Failed to send pixel data command"))

            return True

        except Exception as e:
            print(f"Error sending pixel data: {e}", file=sys.stderr)
            return False

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
                    raise RuntimeError(_("Failed to map shared memory: {}").format(e)) from e

            except Exception as e:
                raise RuntimeError(_("Failed to create shared memory: {}").format(e)) from e

            # Copy data to shared memory
            try:
                buffer = memoryview(shm.buf).cast('B')
                shared_array = np.frombuffer(buffer, dtype=image_data.dtype).reshape(image_data.shape)
                np.copyto(shared_array, image_data)
                # Delete transient objects used to structure copy
                del buffer
                del shared_array
            except Exception as e:
                raise RuntimeError(_("Failed to copy data to shared memory: {}").format(e)) from e

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
                raise RuntimeError(_("Failed to send set_seq_frame_pixeldata command"))

            return True

        except Exception as e:
            print("Error sending pixel data: {e}", file=sys.stderr)
            return False

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
            RuntimeError: For other errors during  data retrieval,
            ValueError: If the shared memory data is invalid
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
                    raise RuntimeError(_("Server error: {}").format(error_msg))
                raise RuntimeError(_("Failed to initiate shared memory transfer: Empty response"))

            if status == _Status.NONE:
                return None

            if not response:
                raise RuntimeError(_("Failed to initiate shared memory transfer: No data received"))

            if len(response) < 280: # Not a correct SharedMemoryInfo payload
                return None
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
                raise RuntimeError(_("Failed to map shared memory: {}").format(e)) from e

            try:
                buffer = bytearray(shm.buf)[:shm_info.size]
                result = bytes(buffer)
            except (BufferError, ValueError, TypeError) as e:
                raise RuntimeError(_("Failed to create bytes from shared memory: {}").format(e)) from e

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
                    shm.close()  # First close the memory mapping as we have finished with it
                    # (We don't unlink it as C wll do that)

                    # Signal that Python is done with the shared memory and wait for C to finish
                    finish_info = struct.pack('256s', shm_info.shm_name)
                    if not self._execute_command(_Command.RELEASE_SHM, finish_info):
                        raise RuntimeError(_("Failed to cleanup shared memory"))

                except Exception:
                    pass

    def get_image_fits_header(self) -> Optional[str]:
        """
        Retrieve the full FITS header of the current image loaded in Siril.

        Args:
            none.

        Returns:
            bytes: The image FITS header as a string.

        Raises:
            NoImageError: If no image is currently loaded,
            RuntimeError: For other errors during  data retrieval,
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
                        raise NoImageError(_("No image is currently loaded in Siril"))
                    raise RuntimeError(_("Server error: {}").format(error_msg))
                raise RuntimeError(_("Failed to initiate shared memory transfer: Empty response"))

            if not response:
                raise RuntimeError(_("Failed to initiate shared memory transfer: No data received"))

            if status == _Status.NONE:
                return None

            if len(response) < 25: # No payload
                return None

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
                raise RuntimeError(_("Failed to map shared memory: {}").format(e)) from e

            try:
                # Read entire buffer at once
                buffer = bytearray(shm.buf)[:shm_info.size]
                result = buffer.decode('utf-8', errors='ignore')
            except (BufferError, ValueError, TypeError) as e:
                raise RuntimeError(_("Failed to create string from shared memory: {}").format(e)) from e

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
                    shm.close()  # First close the memory mapping as we have finished with it
                    # (We don't unlink it as C wll do that)

                    # Signal that Python is done with the shared memory and wait for C to finish
                    finish_info = struct.pack('256s', shm_info.shm_name)
                    if not self._execute_command(_Command.RELEASE_SHM, finish_info):
                        raise RuntimeError(_("Failed to cleanup shared memory"))

                except Exception:
                    pass

    def get_image_unknown_keys(self) -> Optional[str]:
        """
        Retrieve the unknown key in a FITS header of the current loaded Siril
        image using shared memory.

        Args:
            none.

        Returns:
            bytes: The unknown keys as a string.

        Raises:
            NoImageError: If no image is currently loaded,
            RuntimeError: For other errors during  data retrieval,
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
                        raise NoImageError(_("No image is currently loaded in Siril"))
                    raise RuntimeError(_("Server error: {}").format(error_msg))
                raise RuntimeError(_("Failed to initiate shared memory transfer: Empty response"))

            if status == _Status.NONE:
                return None

            if status == _Status.NONE:
                return None

            if not response:
                raise RuntimeError(_("Failed to initiate shared memory transfer: No data received"))
            if len(response) < 25: # No payload
                return None

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
                raise RuntimeError(_("Failed to map shared memory: {}").format(e)) from e

            try:
                # Read entire buffer at once using memoryview
                buffer = bytearray(shm.buf)[:shm_info.size]
                result = buffer.decode('utf-8', errors='ignore')
            except (BufferError, ValueError, TypeError) as e:
                raise RuntimeError(_("Failed to create string from shared memory: {}").format(e)) from e

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
                    shm.close()  # First close the memory mapping as we have finished with it
                    # (We don't unlink it as C wll do that)

                    # Signal that Python is done with the shared memory and wait for C to finish
                    finish_info = struct.pack('256s', shm_info.shm_name)
                    if not self._execute_command(_Command.RELEASE_SHM, finish_info):
                        raise RuntimeError(_("Failed to cleanup shared memory"))

                except Exception:
                    pass

    def get_image_history(self) -> Optional[list[str]]:
        """
        Retrieve history entries in the FITS header of the current loaded
        Siril image using shared memory.

        Args:
            none.

        Returns:
            list: The history entries in the FITS header as a list of strings.

        Raises:
            NoImageError: If no image is currently loaded,
            RuntimeError: For other errors during data retrieval,
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
                        raise NoImageError(_("No image is currently loaded in Siril"))
                    raise RuntimeError(_("Server error: {}").format(error_msg))
                raise RuntimeError(_("Failed to initiate shared memory transfer: Empty response"))

            if status == _Status.NONE:
                return None

            if not response:
                raise RuntimeError(_("Failed to initiate shared memory transfer: No data received"))

            if len(response) < 29:  # No payload
                return None

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
                raise RuntimeError(_("Failed to map shared memory: {}").format(e)) from e

            try:
                # Read entire buffer at once using memoryview
                buffer = bytearray(shm.buf)[:shm_info.size]
                string_data = buffer.decode('utf-8', errors='ignore')
                string_list = string_data.split('\x00')
            except (BufferError, ValueError, TypeError) as e:
                raise RuntimeError(_("Failed to create string from shared memory: {}").format(e)) from e

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
                    shm.close()  # First close the memory mapping as we have finished with it
                    # (We don't unlink it as C wll do that)

                    # Signal that Python is done with the shared memory and wait for C to finish
                    finish_info = struct.pack('256s', shm_info.shm_name)
                    if not self._execute_command(_Command.RELEASE_SHM, finish_info):
                        raise RuntimeError(_("Failed to cleanup shared memory"))

                except Exception:
                    pass

    def get_siril_wd(self) -> Optional[str]:
        """
        Request the current working directory from Siril.

        Returns:
            The current working directory as a string, or None if an error occurred.
        """

        response = self._request_data(_Command.GET_WORKING_DIRECTORY)

        if response is None:
            return None

        try:
            # Assuming the response is a null-terminated UTF-8 encoded string
            wd = response.decode('utf-8').rstrip('\x00')
            return wd
        except UnicodeDecodeError as e:
            print(f"Error decoding working directory: {e}", file=sys.stderr)
            return None

    def get_siril_configdir(self) -> Optional[str]:
        """
        Request the user config directory used by Siril.

        Returns:
            The user config directory as a string, or None if an error occurred.
        """

        response = self._request_data(_Command.GET_USERCONFIGDIR)

        if response is None:
            return None

        try:
            # Assuming the response is a null-terminated UTF-8 encoded string
            wd = response.decode('utf-8').rstrip('\x00')
            return wd
        except UnicodeDecodeError as e:
            print(f"Error decoding user config directory: {e}", file=sys.stderr)
            return None

    def get_siril_userdatadir(self) -> Optional[str]:
        """
        Request the user data directory used by Siril.

        Returns:
            The user data directory as a string, or None if an error occurred.
        """

        response = self._request_data(_Command.GET_USERDATADIR)

        if response is None:
            return None

        try:
            # Assuming the response is a null-terminated UTF-8 encoded string
            path = response.decode('utf-8').rstrip('\x00')
            return path
        except UnicodeDecodeError as e:
            print(f"Error decoding user data directory: {e}", file=sys.stderr)
            return None

    def get_siril_systemdatadir(self) -> Optional[str]:
        """
        Request the system data directory used by Siril.

        Returns:
            The system data directory as a string, or None if an error occurred.
        """

        response = self._request_data(_Command.GET_SYSTEMDATADIR)

        if response is None:
            return None

        try:
            # Assuming the response is a null-terminated UTF-8 encoded string
            path = response.decode('utf-8').rstrip('\x00')
            return path
        except UnicodeDecodeError as e:
            print(f"Error decoding user data directory: {e}", file=sys.stderr)
            return None

    def is_image_loaded(self) -> Optional[bool]:
        """
        Check if a single image is loaded in Siril.

        Returns:
            bool: True if a single image is loaded, False if a single image is
            not loaded, or None if an error occurred.
        """
        response = self._request_data(_Command.GET_IS_IMAGE_LOADED)

        if response is None:
            return None

        image_loaded = struct.unpack('!i', response)[0] != 0
        return image_loaded

    def is_sequence_loaded(self) -> Optional[bool]:
        """
        Check if a sequence is loaded in Siril.

        Returns:
            bool: True if a sequence is loaded, False if a sequence is not loaded,
            or None if an error occurred.
        """
        response = self._request_data(_Command.GET_IS_SEQUENCE_LOADED)

        if response is None:
            return None

        sequence_loaded = struct.unpack('!i', response)[0] != 0
        return sequence_loaded

    def get_image_filename(self) -> Optional[str]:
        """
        Request the filename of the loaded image from Siril.

        Returns:
            The filename as a string, or None if an error occurred.
        """

        response = self._request_data(_Command.GET_FILENAME)

        if response is None:
            return None

        try:
            # Assuming the response is a null-terminated UTF-8 encoded string
            wd = response.decode('utf-8').rstrip('\x00')
            return wd
        except UnicodeDecodeError as e:
            print(f"Error decoding loaded image filename: {e}", file=sys.stderr)
            return None

    def get_seq_frame_filename(self, frame: int) -> Optional[str]:
        """
        Request the filename of the specified frame of the loaded sequence from Siril.

        Returns:
            The filename as a string, or None if an error occurred.
        """

        # Convert frame number to network byte order bytes
        frame_payload = struct.pack('!I', frame)  # '!I' for network byte order uint32_t

        response = self._request_data(_Command.GET_SEQ_FRAME_FILENAME, payload=frame_payload)

        if response is None:
            return None

        try:
            # Assuming the response is a null-terminated UTF-8 encoded string
            filename = response.decode('utf-8').rstrip('\x00')
            return filename
        except UnicodeDecodeError as e:
            print(f"Error decoding loaded image filename: {e}", file=sys.stderr)
            return None

    def get_image_stats(self, channel: int) -> Optional[ImageStats]:
        """
        Request image statistics from Siril for a specific channel.

        Args:
            channel: Integer specifying which channel to get statistics
                     for (typically 0, 1, or 2)

        Returns:
            ImageStats object containing the statistics, or None if an error occurred
        """

        # Convert channel number to network byte order bytes
        channel_payload = struct.pack('!I', channel)  # '!I' for network byte order uint32_t

        # Request data with the channel number as payload
        response = self._request_data(_Command.GET_IMAGE_STATS, payload=channel_payload)

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

        Args:
            frame: Integer specifying which frame in the sequence to get registration
                   data for (between 0 and Sequence.number),
            channel: Integer specifying which channel to get registration data
                     for (typically 0, 1, or 2)

        Returns:
            RegData object containing the registration data, or None if an error occurred
        """

        data_payload = struct.pack('!II', frame, channel)  # '!I' for network byte order uint32_t

        # Request data with the channel number as payload
        response = self._request_data(_Command.GET_SEQ_REGDATA, payload=data_payload)

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
        """

        data_payload = struct.pack('!II', frame, channel)  # '!I' for network byte order uint32_t

        # Request data with the channel number as payload
        response = self._request_data(_Command.GET_SEQ_STATS, payload=data_payload)

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
        Request sequence frame metadata from Siril.

        Args:
            frame: Integer specifying which frame in the sequence to get image
                   metadata for (between 0 and Sequence.number)

        Returns:
            ImgData object containing the frame metadata, or None if an error occurred
        """

        data_payload = struct.pack('!I', frame)  # '!I' for network byte order uint32_t

        # Request data with the channel number as payload
        response = self._request_data(_Command.GET_SEQ_IMGDATA, payload=data_payload)
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

    def get_seq_distodata(self, channel: int) -> Optional[DistoData]:
        """
        Request sequence distortion data from Siril

        channel: Integer specifying which channel to get registration data
        for (typically 0, 1, or 2)

        Returns:
            DistoData object containing the channel distortion parameters, or None if an error occurred
        """

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
                distofilename_string = remaining_data.decode('utf-8')
            else:
                distofilename_string = ''

            return DistoData (
                index = values[0],
                velocity = (values[1], values[2]),
                filename = distofilename_string
            )
        except struct.error as e:
            print(f"Error unpacking distortion data: {e}", file=sys.stderr)
            return None

    def get_seq(self) -> Optional[Sequence]:
        """
        Request metadata for the current sequence loaded in Siril.

        Returns:
            Sequence object containing the current sequence metadata, or None
            if an error occurred
        """

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
            seqname_string = remaining_data.decode('utf-8')

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
            print(f"Error unpacking sequence data: {e}", file=sys.stderr)
            return None

    def get_image_keywords(self) -> Optional[FKeywords]:
        """
        Request FITS keywords data from Siril as a FKeywords object.

        Returns:
            FKeywords object containing the FITS keywords, or None if an error occurred
        """

        response = self._request_data(_Command.GET_KEYWORDS)
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

            # Replace default values and unphysical values
            def decimal_to_dms(decimal, is_latitude=True):
                """Convert decimal degrees to degrees, minutes, seconds string."""
                # Get the absolute value and direction
                absolute = abs(decimal)
                if is_latitude:
                    direction = 'N' if decimal >= 0 else 'S'
                else:
                    direction = 'E' if decimal >= 0 else 'W'

                # Calculate degrees, minutes, seconds
                degrees = int(absolute)
                minutes_decimal = (absolute - degrees) * 60
                minutes = int(minutes_decimal)
                seconds = round((minutes_decimal - minutes) * 60, 2)

                # Format as string
                return f"{degrees}°{minutes}'{seconds}\"{direction}"

            values = [None if val in _Defaults.VALUES else val for val in values]
            if values[9] == "" and values[29]: # sitelat_str
                values[9] = decimal_to_dms(values[29])
            if values[10] == "" and values[30]: # sitelong_str
                values[10] = decimal_to_dms(values[30])

            # Helper function to decode and strip null-terminated strings
            def decode_string(s: bytes) -> str:
                return s.decode('utf-8').rstrip('\x00')

            # Helper function to convert timestamp to datetime
            def timestamp_to_datetime(timestamp: int) -> Optional[datetime]:
                return datetime.fromtimestamp(timestamp) if timestamp != 0 else None

            # Create FKeywords object
            kw = FKeywords(
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
                # if fhi is 0.0, set both fhi and flo to None
                flo=values[17] if values[18] != 0.0 else None,
                fhi=values[18] if values[18] != 0.0 else None,
                data_max=values[19],
                data_min=values[20],
                pixel_size_x=values[21] if values[21] and values[21] > 0.0 else None,
                pixel_size_y=values[22] if values[22] and values[21] > 0.0 else None,
                binning_x=values[23] if values[23] and values[24] > 1 else 1,
                binning_y=values[24] if values[24] and values[24] > 1 else 1,
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
                focal_length=values[35] if values[35] and values[35] > 0.0 else None,
                flength=values[36] if values[36] and values[36] > 0.0 else None,
                iso_speed=values[37],
                exposure=values[38],
                aperture=values[39],
                ccd_temp=values[40],
                set_temp=values[41],
                livetime=values[42],
                stackcnt=values[43],
                cvf=values[44],
                gain=values[45],
                offset=values[46],
                focuspos=values[47],
                focussz=values[48],
                foctemp=values[49],
                date=timestamp_to_datetime(values[50]),
                date_obs=timestamp_to_datetime(values[51])
            )

            return kw

        except struct.error as e:
            print(f"Error unpacking FITS keywords data: {e}", file=sys.stderr)
            return None
        except Exception as e:
            print(f"Error processing FITS keywords data: {e}", file=sys.stderr)
            return None

    def get_image(self, with_pixels: Optional[bool] = True) -> Optional[FFit]:
        """
        Request a copy of the current image open in Siril.

        Args:
            with_pixels: optional bool specifying whether to get pixel data as a
                         NumPy array, or only the image metadata. Defaults to True

        Returns:
            FFit object containing the image metadata and (optionally)
            pixel data, or None if an error occurred
        """

        # Request data with the channel number as payload
        response = self._request_data(_Command.GET_IMAGE)
        if response is None:
            return None

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
                    img_pixeldata = self.get_image_pixeldata()
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
                top_down=bool(values[10]),
                _focalkey=bool(values[11]),
                _pixelkey=bool(values[12]),
                color_managed=bool(values[13]),
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

        except struct.error as e:
            print(f"Error unpacking FITS metadata: {e}", file=sys.stderr)
            return None
        except Exception as e:
            print(f"Error processing FITS metadata: {e}", file=sys.stderr)
            return None

    def get_seq_frame(self, frame: int, with_pixels: Optional[bool] = True) -> Optional[FFit]:
        """
        Request sequence frame as a FFit from Siril.

        Args:
            frame: Integer specifying which frame in the sequence to retrieve data for
            (between 0 and Sequence.number)
            with_pixels: bool specifying whether or not to return the pixel data for the
            frame (default is True).

        Returns:
            FFit object containing the frame data, or None if an error occurred
        """

        shm = None
        data_payload = struct.pack('!I?', frame, with_pixels)
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
                'Q',  # data_type (padded to uint64_t)
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
                print(f"Received image data size {len(response)} doesn't match expected size {expected_size}",
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

            # Create FFit object:
            # Get pixeldata if requested
            pixeldata = None
            if with_pixels:
                try:
                    shm_info = _SharedMemoryInfo(
                        size = values[66],
                        data_type = values[67],
                        width = values[68],
                        height = values[69],
                        channels = values[70],
                        shm_name = values[71]
                    )
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
                        raise RuntimeError(_("Failed to map shared memory: {}").format(e)) from e

                    buffer = bytearray(shm.buf)[:shm_info.size]
                    # Create numpy array from shared memory
                    dtype = np.float32 if shm_info.data_type == 1 else np.uint16
                    try:
                        arr = np.frombuffer(buffer, dtype=dtype)
                    except (BufferError, ValueError, TypeError) as e:
                        raise RuntimeError(_("Failed to create array from shared memory: {}").format(e)) from e

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
                        raise ValueError(_("Failed to reshape array to image dimensions: {}").format(e)) from e

                    # Make a copy of the data since we'll be releasing the shared memory
                    pixeldata = np.copy(arr)

                except Exception as e:
                    print(f"Error obtaining pixeldata: {e}")

            fits_keywords = FKeywords(
                program=decode_string(values[14]),
                filename=decode_string(values[15]),
                row_order=decode_string(values[16]),
                filter=decode_string(values[17]),
                image_type=decode_string(values[18]),
                object=decode_string(values[19]),
                instrume=decode_string(values[20]),
                telescop=decode_string(values[21]),
                observer=decode_string(values[22]),
                sitelat_str=decode_string(values[23]),
                sitelong_str=decode_string(values[24]),
                bayer_pattern=decode_string(values[25]),
                focname=decode_string(values[26]),
                bscale=values[27],
                bzero=values[28],
                lo=values[29],
                hi=values[30],
                flo=values[31],
                fhi=values[32],
                data_max=values[33],
                data_min=values[34],
                pixel_size_x=values[35],
                pixel_size_y=values[36],
                binning_x=values[37],
                binning_y=values[38],
                expstart=values[39],
                expend=values[40],
                centalt=values[41],
                centaz=values[42],
                sitelat=values[43],
                sitelong=values[44],
                siteelev=values[45],
                bayer_xoffset=values[46],
                bayer_yoffset=values[47],
                airmass=values[48],
                focal_length=values[49],
                flength=values[50],
                iso_speed=values[51],
                exposure=values[52],
                aperture=values[53],
                ccd_temp=values[54],
                set_temp=values[55],
                livetime=values[56],
                stackcnt=values[57],
                cvf=values[58],
                gain=values[59],
                offset=values[60],
                focuspos=values[61],
                focussz=values[62],
                foctemp=values[63],
                date=timestamp_to_datetime(values[64]),
                date_obs=timestamp_to_datetime(values[65])
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
                top_down=bool(values[10]),
                _focalkey=bool(values[11]),
                _pixelkey=bool(values[12]),
                color_managed=bool(values[13]),
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

        except struct.error as e:
            print(f"Error unpacking FITS metadata: {e}", file=sys.stderr)
            return None
        except Exception as e:
            print(f"Error processing FITS metadata: {e}", file=sys.stderr)
            return None
        finally:
            if shm is not None:
                try:
                    shm.close()  # First close the memory mapping as we have finished with it
                    # (We don't unlink it as C wll do that)

                    # Signal that Python is done with the shared memory and wait for C to finish
                    finish_info = struct.pack('256s', shm_info.shm_name)
                    if not self._execute_command(_Command.RELEASE_SHM, finish_info):
                        raise RuntimeError(_("Failed to cleanup shared memory"))

                except Exception:
                    pass

    def get_image_stars(self) -> List[PSFStar]:
        """
        Request star model PSF data from Siril.

        Returns:
            List of PSFStar objects containing the star data, or None if
            no stars have been detected. (The "findstar" command should
            be run first to detect stars in the image.)

        Raises:
            NoImageError: If no image is currently loaded,
            RuntimeError: For other errors during  data retrieval,
            ValueError: If the received data format is invalid
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
                    raise RuntimeError(_("Server error: {}").format(error_msg))
                raise RuntimeError(_("Failed to initiate shared memory transfer: Empty response"))

            if status == _Status.NONE:
                return None

            if not response:
                raise RuntimeError(_("Failed to initiate shared memory transfer: No data received"))

            if status == _Status.NONE:
                return None

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
                raise RuntimeError(_("Failed to map shared memory: {}").format(e)) from e

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
                    # Extract the bytes for this struct and unpack
                    values = struct.unpack(format_string, buffer[start:end])

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
            raise RuntimeError(_("Error processing star data: {}").format(e)) from e

        finally:
            if shm is not None:
                try:
                    shm.close()  # First close the memory mapping as we have finished with it
                    # (We don't unlink it as C wll do that)

                    # Signal that Python is done with the shared memory and wait for C to finish
                    finish_info = struct.pack('256s', shm_info.shm_name)
                    if not self._execute_command(_Command.RELEASE_SHM, finish_info):
                        raise RuntimeError(_("Failed to cleanup shared memory"))

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
            error occurred

        Raises:
            RuntimeError: if an error occurred getting the requested config value
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
            return None

        except Exception as e:
            raise RuntimeError(_("Error getting config value: {}").format(e)) from e

    def set_seq_frame_incl(self, index: int, incl: bool) -> bool:
        """
        Set whether a given frame is included in the currently loaded sequence
        in Siril. This method is intended for use in creating custom sequence
        filters.

        Args:
            index: integer specifying which frame to set the pixeldata for.
            incl: bool specifying whether the frame is included or not.

        Returns:
            bool: True if successful, False otherwise
        """
        try:
            # Pack the index and incl into bytes using network byte order (!)
            payload = struct.pack('!II', index, incl)
            return self._execute_command(_Command.SET_SEQ_FRAME_INCL, payload)
        except Exception as e:
            print(f"Error setting selection: {e}", file=sys.stderr)
            return False

    def get_image_bgsamples(self) -> Optional[List[BGSample]]:
        """
        Request background samples data from Siril.

        Returns:
            List of BGSamples background samples, with each set of coordinates
            expressed as a tuple[float, float], or None if no background
            samples have been set.

        Raises:
            NoImageError: If no image is currently loaded,
            RuntimeError: For other errors during  data retrieval,
            ValueError: If the received data format is invalid
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
                    raise RuntimeError(_("Server error: {}").format(error_msg))
                raise RuntimeError(_("Failed to initiate shared memory transfer: Empty response"))

            if status == _Status.NONE:
                return None

            if not response:
                raise RuntimeError(_("Failed to initiate shared memory transfer: No data received"))

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
                raise RuntimeError(_("Failed to map shared memory: {}").format(e)) from e

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
                    # Extract the bytes for this struct and unpack
                    values = struct.unpack(format_string, buffer[start:end])

                    sample = BGSample(
                        median = (values[0], values[1], values[2]),
                        mean = values[3],
                        min = values[4],
                        max = values[5],
                        size = values[6],
                        position = (values[7], values[8]),
                        valid = bool(values[9])
                    )
                    samples.append(sample)

                except struct.error as e:
                    print(f"Error unpacking sample data for index {i}: {e}", file=sys.stderr)
                    break

            return samples

        except Exception as e:
            raise RuntimeError(_("Error processing sample data: {}").format(e)) from e

        finally:
            if shm is not None:
                try:
                    shm.close()  # First close the memory mapping as we have finished with it
                    # (We don't unlink it as C wll do that)

                    # Signal that Python is done with the shared memory and wait for C to finish
                    finish_info = struct.pack('256s', shm_info.shm_name)
                    if not self._execute_command(_Command.RELEASE_SHM, finish_info):
                        raise RuntimeError(_("Failed to cleanup shared memory"))

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
            header_string = header.tostring(sep='\n')
            # Send the metadata to Siril
            siril.set_image_metadata_from_header_string(header_string)

        Args:
            header: string containing the FITS header data

        Returns:
            bool: True if successful, False otherwise
        """

        shm = None
        shm_info = None
        try:
            # Validate input array
            if not isinstance(header, str):
                raise ValueError(_("Header data must be a string"))

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
                    raise RuntimeError(_("Failed to map shared memory: {}").format(e)) from e

            except Exception as e:
                raise RuntimeError(_("Failed to create shared memory: {}").format(e)) from e

            # Copy data to shared memory
            try:
                buffer = memoryview(shm.buf).cast('B')
                buffer[:len(header_bytes)] = header_bytes

                # Delete transient objects used to structure copy
                del buffer
            except Exception as e:
                raise RuntimeError(_("Failed to copy data to shared memory: {}").format(e)) from e

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
                raise RuntimeError(_("_execute_command failed"))

            return True

        except Exception as e:
            print(f"Error sending FITS header metadata: {e}", file=sys.stderr)
            return False

        finally:
            if shm is not None:
                try:
                    shm.close()
                    self._execute_command(_Command.RELEASE_SHM, shm_info)
                except Exception as e:
                    pass

    def add_user_polygon(self, polygon: UserPolygon)-> UserPolygon:
        """
        Adds a user polygon to the Siril display overlay

        Args:
            polyon: UserPolygon defining the polygon to be added

        Returns:
            UserPolygon: the input updated with the ID assigned by Siril
        """

        try:
            # Serialize the provided polygon
            polygon_bytes = polygon.serialize()
            # Send it using _execute_command
            response = self._request_data(_Command.ADD_USER_POLYGON, polygon_bytes)
            if response is None:
                print("Error sending user polygon", file=sys.stderr)
                return None

            try:
                # Assuming the response is in the format: !i (ID) (4 bytes)
                id = struct.unpack('!i', response)[0]
                polygon.polygon_id = id
                return polygon
            except struct.error as e:
                print(f"Error retrieving assigned polygon ID: {e}", file=sys.stderr)
                return None

        except Exception as e:
            print(f"Error sending user polygon: {e}", file=sys.stderr)
            return None

    def delete_user_polygon(self, polygon_id: int) -> bool:
        """
        Deletes a single user polygon from the Siril overlay, specified by ID

        Args:
            id: int specifying the polygon ID to be deleted
        Returns:
            bool: True if the command succeeded, False otherwise
        """
        try:
            # Create payload: network-order int followed by string
            # '!I' for network byte order 32-bit int
            payload = struct.pack('!i', polygon_id)

            return self._execute_command(_Command.DELETE_USER_POLYGON, payload)

        except Exception as e:
            print(f"Error sending progress update: {e}", file=sys.stderr)
            return False

    def clear_user_polygons(self) -> bool:
        """
        Clears all user polygons from the Siril overlay

        Returns:
            bool: True if the command succeeded, False otherwise
        """

        try:
            success = self._execute_command(_Command.CLEAR_USER_POLYGONS, None)
            return success

        except Exception as e:
            print(f"Error clearing user polygons: {e}", file=sys.stderr)
            return False

    def get_user_polygon(self, polygon_id: int) -> 'UserPolygon':
        """
        Gets a single user polygon from the Siril overlay, specified by ID

        Args:
            id: int specifying the polygon ID to be deleted

        Returns:
            UserPolygon: the specified UserPolygon if it exists, None otherwise

        Raises:
            RuntimeError: if an error occurred processing the command
        """
        try:
            payload = struct.pack('!i', polygon_id)
            # Send it using _request_data
            response = self._request_data(_Command.GET_USER_POLYGON, payload)
            if response is None:
                return None

            # Catch the polygon and disregard leftover bytes
            polygon, _ = UserPolygon.deserialize_polygon(response)
            return polygon
        except Exception as e:
            raise RuntimeError(_("Failed to get user polygon: {}").format(e)) from e

    def get_user_polygon_list(self) -> List['UserPolygon']:
        """
        Gets a List of all user polygons from the Siril overlay

        Returns:
            List[UserPolygon]: the list of UserPolygon if some exist, None otherwise

        Raises:
            RuntimeError: if an error occurred processing the command
        """
        try:
            status, response = self._send_command(_Command.GET_USER_POLYGON_LIST)
            # Handle error responses
            if status == _Status.ERROR:
                if response:
                    error_msg = response.decode('utf-8', errors='replace')
                    if "no image loaded" in error_msg.lower():
                        raise NoImageError(_("No image is currently loaded in Siril"))
                    raise RuntimeError(_("Server error: {}").format(error_msg))
                raise RuntimeError(_("Failed to initiate shared memory transfer: Empty response"))

            if status == _Status.NONE:
                return None

            if not response:
                raise RuntimeError(_("Failed to initiate shared memory transfer: No data received"))

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
                raise RuntimeError(_("Failed to map shared memory: {}").format(e)) from e

            # Read entire buffer at once using memoryview
            buffer = bytearray(shm.buf)[:shm_info.size]

            polygon_list = UserPolygon.deserialize_polygon_list(buffer)

            return polygon_list

        except Exception as e:
            raise RuntimeError(_("Error processing polygon data: {}").format(e)) from e

        finally:
            if shm is not None:
                try:
                    shm.close()  # First close the memory mapping as we have finished with it
                    # (We don't unlink it as C wll do that)

                    # Signal that Python is done with the shared memory and wait for C to finish
                    finish_info = struct.pack('256s', shm_info.shm_name)
                    if not self._execute_command(_Command.RELEASE_SHM, finish_info):
                        raise RuntimeError(_("Failed to cleanup shared memory"))

                except Exception:
                    pass
