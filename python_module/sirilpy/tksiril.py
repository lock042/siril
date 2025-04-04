# Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
# Copyright (C) 2012-2025 team free-astro (see more in AUTHORS file)
# Reference site is https://siril.org
# SPDX-License-Identifier: GPL-3.0-or-later

"""
TKsiril submodule for Siril, providing utility methods to achieve consistent
script GUI appearance using the TKinter toolkit.
"""

import tkinter as tk
from tkinter import ttk
from .connection import SirilInterface
from .exceptions import SirilError

def create_tooltip(widget, text, wrap_length=250):
    """
    Create a tooltip for a given Tkinter widget.

    Args:
        widget (tk.Widget): The widget to attach the tooltip to
        text (str): The tooltip text to display
        max_width (int, optional): Maximum width of the tooltip. Defaults to 300.
        wrap_length (int, optional): Length at which text wraps. Defaults to 250.

    Raises:
        TypeError: If text is not a string or the provided widget is not a
                   valid Tkinter widget
    """
    # Validate widget argument
    if not isinstance(widget, (tk.Widget, tk.Tk, tk.Toplevel)):
        raise TypeError(f"Invalid widget type. Expected a Tkinter widget, got {type(widget)}")

    # Validate text argument
    if not isinstance(text, str):
        raise TypeError(f"Tooltip text must be a string, got {type(text)}")

    def show_tooltip(event):
        try:
            # Ensure the widget still exists
            widget.winfo_exists()
        except tk.TclError:
            return  # Widget has been destroyed

        tooltip = tk.Toplevel(widget)
        tooltip.wm_overrideredirect(True)
        tooltip.wm_geometry(f"+{event.x_root+10}+{event.y_root+10}")

        # Configure tooltip style
        tooltip.configure(bg='lightyellow')

        # Create label with text wrapping
        label = ttk.Label(tooltip, text=text,
                        justify=tk.LEFT,
                        relief=tk.SOLID,
                        borderwidth=1,
                        wraplength=wrap_length,
                        background='lightyellow')
        label.pack(ipadx=5, ipady=5)

        def hide_tooltip():
            try:
                tooltip.destroy()
            except tk.TclError:
                pass  # Ignore if tooltip is already destroyed

        widget.tooltip = tooltip
        widget.bind('<Leave>', lambda e: hide_tooltip())
        tooltip.bind('<Leave>', lambda e: hide_tooltip())

    widget.bind('<Enter>', show_tooltip)
    
def set_as_transient_child(window, parent_id):
    """
    Configure a Tkinter window as transient child of a parent window
    identified by its ID.
    
    Args:
        window: Tkinter window
        parent_id: Parent window ID (provided by SIRIL_PARENT_WINDOW)
        
    Returns:
        bool: True if successful, False otherwise
    """
    import platform
    
    try:
        system = platform.system()
        
        # For all platforms, set dialog-like properties
        try:
            # Use the Tk command interface directly which is more reliable
            window.tk.call('wm', 'attributes', window._w, '-topmost', '1')
            window.tk.call('focus', '-force', window._w)
            
            # Try to set transient using the lower-level Tk command
            if parent_id:
                try:
                    window.tk.call('wm', 'transient', window._w, parent_id)
                except Exception as e:
                    print(f"Transient setting failed: {e}")
                    
                    # Try platform-specific approaches
                    if system == "Linux":
                        try:
                            # Set window type to dialog
                            window.tk.call('wm', 'attributes', window._w, '-type', 'dialog')
                        except Exception:
                            pass
                    elif system == "Darwin":
                        try:
                            # For macOS, try to set utility style
                            window.tk.call('::tk::unsupported::MacWindowStyle', 'style', window._w, 'utility')
                        except Exception:
                            pass
            
            return True
        except Exception as e:
            print(f"Window property setting failed: {e}")
            return False
            
    except Exception as e:
        print(f"Failed to set window properties: {e}")
        return False
    
def match_theme_to_siril(themed_tk, s, on_top=True):
    """
    Match the Tkinter theme to the Siril configuration and set the script dialog
    to have topmost status, meaning that it will remain in front of other
    non-topmost windows.

    Args:
        s (SirilInterface): sirilpy.SirilInterface class to provide the
                            Siril theme (light or dark) to match
        themed_tk (ThemedTk): ThemedTk instance to apply the theme to
        on_top: whether the script window should be always on top of other windows

    Raises:
        TypeError: If input arguments are of incorrect type
        ValueError: If the theme configuration is not 0 or 1
        AttributeError: If required methods are not available
        RuntimeError: If there are errors installing or setting the theme
    """
    # Strict type checking for s
    if not isinstance(s, SirilInterface):
        raise TypeError(f"First argument must be a SirilInterface. Got {type(s)}")

    # Check if s is an instance of expected SirilInterface class
    try:
        s.__class__.__name__  # Ensure the object is instantiated
    except Exception as e:
        raise TypeError(f"Invalid SirilInterface object: {e}") from e

    # Check if themed_tk has the required method
    if not (hasattr(themed_tk, 'set_theme') or
            (hasattr(themed_tk, 'configure') and hasattr(themed_tk, 'winfo_class'))):
        raise TypeError(f"Second argument must be a Tkinter-like object with theme-setting capabilities. Got {type(themed_tk)}")

    # Get parent window ID from environment variable
    import os
    parent_window_id = os.environ.get("SIRIL_PARENT_WINDOW")
    
    # Try to set the window as transient if parent ID is available
    if parent_window_id:
        set_as_transient_child(themed_tk, parent_window_id)
    
    #if on_top is True:
        # Settings to keep the script window above others
        #themed_tk.focus_force()
        #themed_tk.attributes('-topmost', True)

    # Get theme configuration from Siril
    #try:
    theme_value = s.get_siril_config("gui", "theme")
    #except Exception as e:
        #raise AttributeError(f"Unable to retrieve theme configuration: {e}") from e

    # Map theme values to theme names
    theme_map = {
        0: "equilux",
        1: "arc"
    }

    # Check if theme value is valid
    if theme_value not in theme_map:
        raise ValueError(f"Invalid theme value: {theme_value}. Expected 0 or 1.")

    # Attempt to set the theme
    try:
        # Try set_theme method first (for ttkbootstrap/ttkthemes)
        if hasattr(themed_tk, 'set_theme'):
            themed_tk.set_theme(theme_map[theme_value])
        # Fallback for other theming methods
        elif hasattr(themed_tk, 'configure'):
            themed_tk.configure(theme=theme_map[theme_value])
        else:
            raise RuntimeError("No valid theme-setting method found")
    except Exception as e:
        raise RuntimeError(f"Failed to set theme: {e}") from e

def standard_style():
    """
    Provide a standardised ttk style to allow consistent visual appearance
    between different Siril python scripts.

    Args:
        none

    Raises:
        SirilError: If the style creation or configuration fails
    """
    try:
        style = ttk.Style()
        # Configure style
        style.configure("TFrame", padding=5)
        style.configure("TButton", padding=5)
        style.configure("TCheckbutton", padding=2)
        style.configure("TLabel", padding=2)
        style.configure("Header.TLabel", font=("Helvetica", 12, "bold"))
        style.configure("Value.TLabel", font=("Helvetica", 9))
        return style

    except Exception as e:
        raise SirilError(f"Failed to configure style: {e}") from e
