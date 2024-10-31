from .translations import _, N_

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
    __doc__ = N_("Raised when there are problems connecting to or "
                 "communicating with Siril.\n"
                 "This includes cases like:\n"
                 "- Siril not running\n"
                 "- Socket connection failures\n"
                 "- Communication protocol errors\n"
                 "- Unexpected disconnections")

    def __init__(self, message: str = _("Failed to connect to Siril")):
        super().__init__(message)

class CommandError(SirilError):
    __doc__ = N_("Raised when a Siril command fails to execute properly.\n"
                 "This includes cases like:\n"
                 "- Invalid command parameters\n"
                 "- Command execution failures\n"
                 "- Unexpected command responses\n"
                 "- Command timeout\n")
    def __init__(self, message: str = _("Command execution failed")):
        super().__init__(message)

class DataError(SirilError):
    __doc__ = N_("Raised when there are problems with data handling.\n"
                 "This includes cases like:\n"
                 "- Invalid image data\n"
                 "- Data conversion errors\n"
                 "- Memory allocation failures\n"
                 "- Buffer overflows")
    def __init__(self, message: str = _("Error handling data")):
        super().__init__(message)

class NoImageError(SirilError):
    __doc__ = N_("Raised when a method requires an image to be loaded "
                 " but no image is loaded.\n")
    def __init__(self, message: str = _("No Siril image loaded")):
        super().__init__(message)
