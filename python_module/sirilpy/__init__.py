# Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
# Copyright (C) 2012-2025 team free-astro (see more in AUTHORS file)
# Reference site is https://siril.org
# SPDX-License-Identifier: GPL-3.0-or-later

from typing import TYPE_CHECKING

# Import translation functions first
from .translations import _

# TYPE_CHECKING is False at runtime but True during type checking.
if TYPE_CHECKING:
    from .connection import LogColor, SirilInterface
    from .models import (
        DataType,
        ImageStats,
        FKeywords,
        FFit,
        Homography,
        StarProfile,
        SequenceType,
        PSFStar,
        BGSample,
        RegData,
        ImgData,
        DistoData,
        Sequence
    )
    from .plot import PlotType, SeriesData, PlotData, _PlotSerializer
    from .shm import SharedMemoryWrapper
    from .utility import (
        human_readable_size,
        download_with_progress,
        ensure_installed,
        check_module_version,
        SuppressedStdout,
        SuppressedStderr
    )
    from .exceptions import (
        SirilError,
        ConnectionError,
        CommandError,
        DataError,
        NoImageError,
        NoSequenceError
    )

# Runtime imports
from .models import (
    DataType,
    ImageStats,
    FKeywords,
    FFit,
    Homography,
    StarProfile,
    SequenceType,
    PSFStar,
    BGSample,
    RegData,
    ImgData,
    DistoData,
    Sequence
)
from .plot import PlotType, SeriesData, PlotData, _PlotSerializer
from .shm import SharedMemoryWrapper
from .utility import (
    human_readable_size,
    download_with_progress,
    ensure_installed,
    check_module_version,
    SuppressedStdout,
    SuppressedStderr
)
from .exceptions import (
    SirilError,
    ConnectionError,
    CommandError,
    DataError,
    NoImageError,
    NoSequenceError
)
from .connection import LogColor, SirilInterface

try: # import from the packaging specification
    from importlib.metadata import metadata
    meta = metadata("sirilpy")
    __version__ = meta.get("version", "unknown")
    __author__ = meta.get("author", "unknown")
    __license__ = meta.get("license", "unknown")
except Exception:
    pass

__copyright__ = " (c) Team free-astro 2024-2025" # not a standard metadata

# Define public API
__all__ = [
    'ensure_installed',
    'check_module_version',
    'SirilInterface',
    'LogColor',
    'DataType',
    'ImageStats',
    'FKeywords',
    'FFit',
    'Homography',
    'StarProfile',
    'SequenceType',
    'PSFStar',
    'BGSample',
    'RegData',
    'ImgData',
    'DistoData',
    'Sequence',
    'PlotType',
    'SeriesData',
    'PlotData',
    '_PlotSerializer',
    'SirilError',
    'ConnectionError',
    'CommandError',
    'DataError',
    'NoImageError',
    'NoSequenceError',
    'SharedMemoryWrapper',
    'SuppressedStdout',
    'SuppressedStderr',
    'human_readable_size',
    'download_with_progress',
    '_'
]
