# Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
# Copyright (C) 2012-2025 team free-astro (see more in AUTHORS file)
# Reference site is https://siril.org
# SPDX-License-Identifier: GPL-3.0-or-later
"""
SirilPy - Python interface for Siril astronomical image processing software.

This module provides bindings and utilities for interacting with Siril
from Python, enabling advanced astronomical image processing workflows.
"""

# Import translation functions first
# TODO: this is currently unused (there are no actual translations yet)
from .translations import _

# Regular imports - all modules needed at runtime
from .enums import (
    LogColor,
    CommandStatus,
    DataType,
    StarProfile,
    SequenceType,
    DistoType,
    PlotType,
    SirilVport
)
from .models import (
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
    Sequence,
    SirilPoint,
    UserPolygon
)
from .plot import PlotType, SeriesData, PlotData
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
    SirilConnectionError,
    SharedMemoryError,
    CommandError,
    DataError,
    NoImageError,
    NoSequenceError
)
from .connection import SirilInterface

try:  # import from the packaging specification
    from importlib.metadata import metadata, PackageNotFoundError
    meta = metadata("sirilpy")
    __version__ = meta.get("version", "unknown")
    __author__ = meta.get("author", "unknown")
    __license__ = meta.get("license", "unknown")
except (ImportError, PackageNotFoundError):
    # Specific exceptions rather than general Exception
    __version__ = "unknown"
    __author__ = "unknown"
    __license__ = "unknown"

__copyright__ = " (c) Team free-astro 2024-2025"  # not a standard metadata

# Define public API
__all__ = [
    'ensure_installed',
    'check_module_version',
    'SirilInterface',
    'DataType',
    'ImageStats',
    'FKeywords',
    'FFit',
    'Homography',
    'PSFStar',
    'BGSample',
    'RegData',
    'ImgData',
    'DistoData',
    'Sequence',
    'SirilPoint',
    'UserPolygon',
    'PlotType',
    'SeriesData',
    'PlotData',
    'SirilError',
    'SirilConnectionError',
    'SharedMemoryError',
    'CommandError',
    'DataError',
    'NoImageError',
    'NoSequenceError',
    'SharedMemoryWrapper',
    'SuppressedStdout',
    'SuppressedStderr',
    'human_readable_size',
    'download_with_progress',
    'LogColor',
    'CommandStatus',
    'DataType',
    'StarProfile',
    'SequenceType',
    'DistoType',
    'PlotType',
    'SirilVport',
    '_'
]
