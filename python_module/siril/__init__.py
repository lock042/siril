# siril/__init__.py - Package initialization
from typing import TYPE_CHECKING

# Import translation functions first
from .translations import _

# TYPE_CHECKING is False at runtime but True during type checking.
if TYPE_CHECKING:
    from .connection import SirilInterface, import_or_install
    from .models import (
        DataType,
        ImageStats,
        FKeywords,
        FFit,
        Homography,
        StarProfile,
        SequenceType,
        PSFStar,
        RegData,
        ImgData,
        Sequence
    )
    from .shm import SharedMemoryWrapper
    from .exceptions import (
        SirilError,
        ConnectionError,
        CommandError,
        DataError,
        NoImageError,
        NoSequenceError
    )

# Runtime imports
from .connection import SirilInterface, import_or_install, ensure_installed
from .models import (
    DataType,
    ImageStats,
    FKeywords,
    FFit,
    Homography,
    StarProfile,
    SequenceType,
    PSFStar,
    RegData,
    ImgData,
    Sequence
)
from .shm import SharedMemoryWrapper
from .exceptions import (
    SirilError,
    ConnectionError,
    CommandError,
    DataError,
    NoImageError,
    NoSequenceError
)

# Package metadata
__version__ = "0.1.0"
__author__ = "Team free-astro"
__license__ = "GPLv3+"
__copyright__ = "(c) Team free-astro 2024"

# Define public API
__all__ = [
    'import_or_install',
    'ensure_installed',
    'SirilInterface',
    'DataType',
    'ImageStats',
    'FKeywords',
    'FFit',
    'Homography',
    'StarProfile',
    'SequenceType',
    'PSFStar',
    'RegData',
    'ImgData',
    'Sequence',
    'SirilError',
    'ConnectionError',
    'CommandError',
    'DataError',
    'NoImageError',
    'NoSequenceError',
    'SharedMemoryWrapper',
    '_'
]
