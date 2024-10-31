# siril/__init__.py - Package initialization
from typing import TYPE_CHECKING

# Import translation functions first
from .translations import _, N_

# TYPE_CHECKING is False at runtime but True during type checking.
if TYPE_CHECKING:
    from .connection import SirilInterface
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
        NoImageError
    )

# Runtime imports
from .connection import SirilInterface
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
    NoImageError
)

# Package metadata
__version__ = "0.1.0"
__author__ = "Team free-astro"
__license__ = "GPLv3+"
__copyright__ = "(c) Team free-astro 2024"

# Define public API
__all__ = [
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
    'SharedMemoryWrapper',
    '_',
    'N_'
]