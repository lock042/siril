# siril/translations.py
from typing import Callable
import gettext
import locale
from pathlib import Path

# Initialize with default implementations
_: Callable[[str], str] = gettext.gettext
N_: Callable[[str], str] = gettext.gettext

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

    # Update the global translation functions
    global _
    global N_
    _ = translation.gettext
    N_ = gettext.gettext

# Initialize translations when module is imported
setup_i18n()
