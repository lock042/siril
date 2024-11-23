# Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
# Copyright (C) 2012-2024 team free-astro (see more in AUTHORS file)
# Reference site is https://siril.org
# SPDX-License-Identifier: GPL-3.0-or-later

from typing import TYPE_CHECKING
import importlib

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

try: # import from the packaging specification
    from importlib.metadata import metadata
    meta = metadata("siril")
    __version__ = meta.get("version", "unknown")
    __author__ = meta.get("author", "unknown")
    __license__ = meta.get("license", "unknown")
except Exception:
    pass

__copyright__ = " (c) Team free-astro 2024" # not a standard metadata

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

if sys.platform == "win32":
    s = SirilInterface()
    s.connect()
    bundlepath = s._get_bundle_path()
    os.add_dll_directory(os.path.join(bundlepath, 'bin'))
    s.disconnect()
