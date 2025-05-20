# Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
# Copyright (C) 2012-2025 team free-astro (see more in AUTHORS file)
# Reference site is https://siril.org
# SPDX-License-Identifier: GPL-3.0-or-later
"""
SirilPy - Python interface for Siril astronomical image processing software.

This module provides bindings and utilities for interacting with Siril
from Python, enabling advanced astronomical image processing workflows.
"""

# Core imports required to configure stdout and stderr to accept utf-8
import io
import sys

# Import translation functions first
from .translations import _

# Regular imports - all modules needed at runtime
from .enums import (
    LogColor,
    CommandStatus,
    BitpixType,
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
    PSFStar,
    BGSample,
    RegData,
    ImgData,
    DistoData,
    Sequence,
    FPoint,
    Polygon
)
from .plot import SeriesData, PlotData
from .shm import SharedMemoryWrapper
from .utility import (
    truncate_utf8,
    human_readable_size,
    download_with_progress,
    ensure_installed,
    check_module_version,
    SuppressedStdout,
    SuppressedStderr,
    ONNXHelper
)
from .exceptions import (
    SirilError,
    SirilConnectionError,
    SharedMemoryError,
    CommandError,
    DataError,
    NoImageError,
    NoSequenceError,
    ProcessingThreadBusyError,
    ImageDialogOpenError
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
    'FPoint',
    'Polygon',
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
    'ProcessingThreadBusyError',
    'ImageDialogOpenError',
    'SharedMemoryWrapper',
    'truncate_utf8',
    'SuppressedStdout',
    'SuppressedStderr',
    'ONNXHelper',
    'human_readable_size',
    'download_with_progress',
    'LogColor',
    'CommandStatus',
    'BitpixType',
    'StarProfile',
    'SequenceType',
    'DistoType',
    'PlotType',
    'SirilVport',
    '_'
]

def safe_reconfigure_stream(stream_name):
    stream = getattr(sys, stream_name)
    if hasattr(stream, "reconfigure"):
        try:
            stream.reconfigure(encoding="utf-8", errors="strict")
            return
        except Exception:
            pass  # Fallback below
    # Replace stream with a UTF-8 TextIOWrapper that uses errors='replace'
    buffer = getattr(sys, f"{stream_name}.buffer", None)
    if buffer is not None:
        wrapper = io.TextIOWrapper(buffer, encoding="utf-8", errors="replace")
        setattr(sys, stream_name, wrapper)

safe_reconfigure_stream("stdout")
safe_reconfigure_stream("stderr")
