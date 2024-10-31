# siril/__init__.py - Package initialization

from typing import TYPE_CHECKING

# TYPE_CHECKING is False at runtime but True during type checking.
# This allows us to have circular imports in type hints without causing
# runtime import cycles, while still maintaining type safety.
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
__license__ = "GPLv3+"  # Add if desired
__copyright__ = "(c) Team free-astro 2024"  # Add if desired

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
    'SharedMemoryWrapper'
]

# Set up i18n for the package
from pathlib import Path
import gettext
import locale
import os

def setup_i18n(domain: str = 'siril', localedir: str | None = None) -> None:
    """
    Initialize internationalization for the package.

    Args:
        domain: Translation domain name
        localedir: Directory containing translation files. If None,
                  defaults to the 'locale' directory in the package root.
    """
    if localedir is None:
        localedir = Path(__file__).parent.parent / 'locale'

    try:
        # Try to use system locale
        locale_str = locale.getdefaultlocale()[0]
        if locale_str:
            locale.setlocale(locale.LC_ALL, f"{locale_str}.UTF-8")
    except (locale.Error, TypeError):
        # Fall back to environment variable or default
        locale.setlocale(locale.LC_ALL, '')

    translation = gettext.translation(domain, localedir=localedir, fallback=True)
    translation.install()

    # Make translation functions available package-wide
    global _
    global N_
    _ = translation.gettext
    N_ = gettext.gettext

# Initialize translations when package is imported
setup_i18n()

# Type hints for the translation functions
from typing import Callable
_: Callable[[str], str]
N_: Callable[[str], str]
