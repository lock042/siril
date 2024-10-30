# siril/__init__.py - Package initialization

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

from .exceptions import SirilError, ConnectionError, CommandError, DataError, NoImageError

__version__ = "0.1.0"
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
    'SharedMemoryWrapper'
]
