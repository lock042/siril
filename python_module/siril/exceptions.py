from .translations import _

class SirilError(Exception):
    """
    Base exception class for all Siril-related errors.

    All other Siril exceptions inherit from this class, making it easy
    to catch any Siril-related error with a single except clause.
    """
    def __init__(self, message: str = _("An error occurred")):
        self.message = message
        super().__init__(self.message)

class ConnectionError(SirilError):
    """Raised when there are problems connecting to or
    communicating with Siril.
    This includes cases like:

    - Siril not running

    - Socket connection failures

    - Communication protocol errors

    - Unexpected disconnections"""

    def __init__(self, message: str = _("Failed to connect to Siril")):
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
    def __init__(self, message: str = _("Command execution failed")):
        super().__init__(message)

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