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

def _fix_locale():
    """Fix problematic locales by testing compatibility and falling back to en_US."""
    try:
        current_locale = locale.getlocale()

        # Test if current locale might cause issues by trying to use it
        try:
            # Test the locale by trying to set it again (this will fail for problematic locales)
            if current_locale[0]:
                locale.setlocale(locale.LC_ALL, current_locale)

            # Additional test: try to format a number
            test_val = locale.format_string("%.2f", 3.14)

        except (locale.Error, ValueError, TypeError):
            # Current locale is problematic, switch to en_US.UTF-8 (typical known-good locale)
            try:
                locale.setlocale(locale.LC_ALL, 'en_US.UTF-8')
                os.environ['LC_ALL'] = 'en_US.UTF-8'
                os.environ['LANG'] = 'en_US.UTF-8'
            except locale.Error:
                # If en_US.UTF-8 isn't available, try en_US
                try:
                    locale.setlocale(locale.LC_ALL, 'en_GB')
                    os.environ['LC_ALL'] = 'en_US'
                    os.environ['LANG'] = 'en_US'
                except locale.Error:
                    # Final fallback to C locale
                    try:
                        locale.setlocale(locale.LC_ALL, 'C')
                        os.environ['LC_ALL'] = 'C'
                        os.environ['LANG'] = 'C'
                    except locale.Error:
                        print("Warning: locale appears invalid and unable to set "
                                "safe fallback", file=sys.stderr)
                        pass  # If even C locale fails, leave as is, but with a warning

    except (locale.Error, AttributeError, TypeError):
        # If we can't detect or test locale, don't change anything
        pass

_fix_locale()
