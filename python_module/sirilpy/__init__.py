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
import os
import sys

import locale
def _fix_locale():
    """
    Aggressive locale fix - clear all locale environment variables and set to en_US.UTF-8 or C
    """
    # First, clear ALL locale-related environment variables
    locale_vars = ['LC_ALL', 'LC_CTYPE', 'LC_NUMERIC', 'LC_TIME', 'LC_COLLATE',
                   'LC_MONETARY', 'LC_MESSAGES', 'LC_PAPER', 'LC_NAME',
                   'LC_ADDRESS', 'LC_TELEPHONE', 'LC_MEASUREMENT',
                   'LC_IDENTIFICATION', 'LANG', 'LANGUAGE']

    for var in locale_vars:
        if var in os.environ:
            del os.environ[var]

    # Try to set to a proper locale that returns valid language info
    fallback_locales = ['en_US.UTF-8', 'en_US', 'C.UTF-8', 'C']

    for loc in fallback_locales:
        try:
            os.environ['LC_ALL'] = loc
            os.environ['LANG'] = loc
            locale.setlocale(locale.LC_ALL, loc)

            # Test that it works and returns valid language info
            lang, encoding = locale.getdefaultlocale()
            if lang is not None:  # Make sure we get a valid language
                break

        except Exception:
            continue
    else:
        # If all locales fail, set basic fallback
        os.environ['LC_ALL'] = 'en_US.UTF-8'
        os.environ['LANG'] = 'en_US.UTF-8'
        try:
            locale.setlocale(locale.LC_ALL, 'en_US.UTF-8')
        except Exception as e:
            print(f"Warning: Could not fix locale: {e}", file=sys.stderr)
_fix_locale()

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
    SirilVport,
    ImageType,
    STFType,
    SlidersMode
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
    Polygon,
    ImageAnalysis
)
from .plot import SeriesData, PlotData
from .shm import SharedMemoryWrapper
from .utility import (
    siril_header_to_dict,
    truncate_utf8,
    human_readable_size,
    download_with_progress,
    ensure_installed,
    check_module_version,
    needs_module_version,
    SuppressedStdout,
    SuppressedStderr,
)
from .gpuhelper import (
    ONNXHelper,
    TorchHelper,
    JaxHelper,
)
from .exceptions import (
    SirilError,
    SirilConnectionError,
    SharedMemoryError,
    CommandError,
    DataError,
    NoImageError,
    MouseModeError,
    NoSequenceError,
    ProcessingThreadBusyError,
    ImageDialogOpenError
)
from .connection import SirilInterface

from .version import __version__, __author__, __license__, __copyright__

# Define public API
__all__ = [
    'ensure_installed',
    'check_module_version',
    'needs_module_version',
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
    'MouseModeError',
    'ProcessingThreadBusyError',
    'ImageDialogOpenError',
    'SharedMemoryWrapper',
    'truncate_utf8',
    'SuppressedStdout',
    'SuppressedStderr',
    'ONNXHelper',
    'TorchHelper',
    'JaxHelper',
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
    'ImageAnalysis',
    'ImageType',
    'STFType',
    'SlidersMode',
    'siril_header_to_dict',
    '_'
]

def _safe_reconfigure_stream(stream_name):
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

_safe_reconfigure_stream("stdout")
_safe_reconfigure_stream("stderr")
