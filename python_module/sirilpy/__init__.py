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
    Validate and fix the process locale:
    - Keep the current locale if it's supported by Python.
    - Otherwise, fall back to a known good UTF-8 locale (en_US.UTF-8).
    - Never fall back to 'C' or 'C.UTF-8'.
    """

    def is_locale_supported(loc_name: str) -> bool:
        """Return True if loc_name works in both Python and Babel."""
        try:
            # Try to activate in Python's locale system
            locale.setlocale(locale.LC_ALL, loc_name)
            lang, encoding = locale.getlocale()
            if not lang or not encoding:
                return False
            return True
        except Exception:
            return False

    # 1: Detect current environment locale
    current_locale = os.environ.get('LC_ALL') or os.environ.get('LANG') or ''
    current_locale = current_locale.strip() or None

    # 2: Try the current locale first
    if current_locale and is_locale_supported(current_locale):
        # Ensure environment variables are consistent
        os.environ['LC_ALL'] = current_locale
        os.environ['LANG'] = current_locale
        locale.setlocale(locale.LC_ALL, current_locale)
        return

    # 3: Fallbacks to known good locales (all UTF-8, all English)
    fallback_locales = [
        'en_US.UTF-8',
        'en_GB.UTF-8',
        'en_AU.UTF-8',
        'en_CA.UTF-8',
        'en.UTF-8',
    ]

    for loc in fallback_locales:
        if is_locale_supported(loc):
            os.environ['LC_ALL'] = loc
            os.environ['LANG'] = loc
            locale.setlocale(locale.LC_ALL, loc)
            return

    # 4: Final safety fallback if all else fails
    fallback = 'en_US.UTF-8'
    os.environ['LC_ALL'] = fallback
    os.environ['LANG'] = fallback
    try:
        locale.setlocale(locale.LC_ALL, fallback)
    except Exception as e:
        print(f"Warning: Could not set locale to {fallback}: {e}", file=sys.stderr)

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
