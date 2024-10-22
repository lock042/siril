# siril/__init__.py - Package initialization

from .connection import SirilInterface

from .commands import (
    get_ffit,
)

from .models import (
    DataType,
    ImageStats,
    FKeywords,
    FFit,
)

from .exceptions import SirilError, ConnectionError, CommandError, DataError

__version__ = "0.1.0"
__all__ = [
#    'get_image_info',
#    'get_region',
#    'find_stars',
    'SirilInterface',
    'DataType',
    'ImageStats',
    'FKeywords',
    'FFit',
    'SirilError',
    'ConnectionError',
    'CommandError',
    'DataError',
]
