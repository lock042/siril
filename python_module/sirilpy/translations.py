# Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
# Copyright (C) 2012-2025 team free-astro (see more in AUTHORS file)
# Reference site is https://siril.org
# SPDX-License-Identifier: GPL-3.0-or-later

import os
import gettext
import locale
from functools import wraps
from typing import Callable, TypeVar

T = TypeVar('T')

def setup_translations(domain: str = 'sirilpy', localedir: str = None) -> Callable:
    """
    Set up translations for the module based on the system locale.
    """
    # If no localedir is specified, find the module's parent directory
    if localedir is None:
        try:
            current_module_dir = os.path.dirname(os.path.abspath(__file__))
            localedir = os.path.join(os.path.dirname(current_module_dir), 'locale')
        except Exception:
            # Fallback to module root if locale directory can't be found
            localedir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    # Bind the text domain to the locale directory
    gettext.bindtextdomain(domain, localedir)
    gettext.textdomain(domain)

    # Try to set locale, with graceful error handling
    try:
        locale.setlocale(locale.LC_ALL, '')
    except locale.Error:
        pass

    # Get default locale
    lang, encoding = locale.getdefaultlocale()

    # Try to set locale with progressive fallback
    try:
        locale.setlocale(locale.LC_ALL, lang)
        return gettext.gettext
    except locale.Error:
        try:
            locale.setlocale(locale.LC_ALL, lang.split('_')[0])
            return gettext.gettext
        except locale.Error:
            return gettext.gettext

# Initialize the translation function
_ = setup_translations()

# Export these symbols
__all__ = ['_', 'setup_translations']
