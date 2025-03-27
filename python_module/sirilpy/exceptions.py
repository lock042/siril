# Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
# Copyright (C) 2012-2025 team free-astro (see more in AUTHORS file)
# Reference site is https://siril.org
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Exceptions submodule for Siril, providing exception classes for use
in exception raising within the sirilpy module.
"""

from .translations import _
from .enums import CommandStatus

class SirilError(Exception):
    """
    Base exception class for all Siril-related errors.

    All other Siril exceptions inherit from this class, making it easy
    to catch any Siril-related error with a single except clause.
    """
    def __init__(self, message: str = _("An error occurred")):
        self.message = message
        super().__init__(self.message)

class SirilConnectionError(SirilError):
    """Raised when there are problems connecting to or
    communicating with Siril.
    This includes cases like:

    - Siril not running

    - Socket connection failures

    - Communication protocol errors

    - Unexpected disconnections"""

    def __init__(self, message: str = _("Failed to connect to Siril")):
        super().__init__(message)

class SharedMemoryError(SirilError):
    """
    Raised when there are problems connecting to or
    communicating with Siril using shared memory.
    """

    def __init__(self, message: str = _("Siril shared memory error")):
        super().__init__(message)

class CommandError(SirilError):
    """
    Raised when a command sent to Siril fails to execute properly.
    (Note: 'command' in this case refers to internal commands sent
    from the python module to the Siril python handler, not Siril
    commands of the type that might be entered in the Siril command
    entry.) This includes cases like:

    - Invalid command parameters

    - Command execution failures

    - Unexpected command responses

    - Command timeout
    """

    def __init__(self, message: str = _("Command execution failed"),
                 status_code=CommandStatus.CMD_GENERIC_ERROR):
        super().__init__(message)
        self.status_code = status_code

class DataError(SirilError):
    """
    Raised when there are problems with data handling.
    This includes cases like:

    - Invalid image data

    - Data conversion errors

    - Memory allocation failures

    - Buffer overflows
    """
    def __init__(self, message: str = _("Error handling data")):
        super().__init__(message)

class NoImageError(SirilError):
    """
    Raised when a method requires an image to be loaded
    but no image is loaded.
    """
    def __init__(self, message: str = _("No Siril image loaded")):
        super().__init__(message)

class NoSequenceError(SirilError):
    """
    Raised when a method requires a sequence to be loaded
    but no sequence is loaded.
    """
    def __init__(self, message: str = _("No Siril sequence loaded")):
        super().__init__(message)
