import os
import gettext
import locale
from functools import wraps
from typing import Callable, TypeVar

T = TypeVar('T')

def setup_translations(domain: str = 'siril', localedir: str = 'locale') -> Callable:
    """
    Set up translations for the module based on the system locale.
    """
    if not os.path.exists(localedir):
        os.makedirs(localedir)

    gettext.bindtextdomain(domain, localedir)
    gettext.textdomain(domain)

    try:
        locale.setlocale(locale.LC_ALL, '')
    except locale.Error:
        pass

    lang, encoding = locale.getdefaultlocale()

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
