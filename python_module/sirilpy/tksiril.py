# Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
# Copyright (C) 2012-2025 team free-astro (see more in AUTHORS file)
# Reference site is https://siril.org
# SPDX-License-Identifier: GPL-3.0-or-later

"""
TKsiril submodule for Siril, providing utility methods to achieve consistent
script GUI appearance using the TKinter toolkit.
"""

import tkinter as tk
import sys
from tkinter import ttk
from .connection import SirilInterface
from .exceptions import SirilError

def create_tooltip(widget, text, wrap_length=250, parent_window=None):
    """
    Create a tooltip for a given Tkinter widget.

    Args:
        widget (tk.Widget): The widget to attach the tooltip to
        text (str): The tooltip text to display
        wrap_length (int, optional): Length at which text wraps. Defaults to 250.
        parent_window (tk.Tk or tk.Toplevel, optional): The parent window that might have topmost attribute

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

        # If parent window has topmost attribute, temporarily disable it
        is_topmost = False
        if parent_window and hasattr(parent_window, 'attributes'):
            try:
                is_topmost = parent_window.attributes('-topmost')
                if is_topmost:
                    parent_window.attributes('-topmost', False)
            except tk.TclError:
                pass  # Ignore if we can't get/set attributes

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
                # Restore topmost attribute if it was enabled
                if is_topmost and parent_window and hasattr(parent_window, 'attributes'):
                    parent_window.attributes('-topmost', True)
            except tk.TclError:
                pass  # Ignore if tooltip is already destroyed

        widget.tooltip = tooltip
        widget.bind('<Leave>', lambda e: hide_tooltip())
        tooltip.bind('<Leave>', lambda e: hide_tooltip())

    widget.bind('<Enter>', show_tooltip)

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

    # Get theme configuration from Siril
    try:
        theme_value = s.get_siril_config("gui", "theme")
    except Exception as e:
        raise AttributeError(f"Unable to retrieve theme configuration: {e}") from e

    # Map theme values to theme names
    theme_map = {
        0: "equilux",
        1: "arc"
    }

    if on_top is True:
        # Settings to keep the script window above others
        themed_tk.focus_force()
        themed_tk.attributes('-topmost', True)

    # This allows temporarily disabling topmost when opening dialogs
    original_filedialog_functions = {}
    original_messagebox_functions = {}
      
    # Function to wrap dialog functions
    def wrap_dialog(original_func):
        def wrapper(*args, **kwargs):
            if on_top is True:
                 # Temporarily disable topmost for the main window
                 themed_tk.attributes('-topmost', False)
            # Call the original dialog function
            result = original_func(*args, **kwargs)
            if on_top is True:
                 # Re-enable topmost after the dialog closes
                 themed_tk.attributes('-topmost', True)
            return result
        return wrapper
    
    # Replace standard dialog functions if tkinter.filedialog is used
    if 'tkinter.filedialog' in sys.modules:
        import tkinter.filedialog as fd
        for func_name in ['askopenfilename', 'askopenfilenames', 'asksaveasfilename', 
                         'askdirectory', 'askopenfile', 'askopenfiles', 'asksaveasfile']:
            if hasattr(fd, func_name):
                original_filedialog_functions[func_name] = getattr(fd, func_name)
                setattr(fd, func_name, wrap_dialog(getattr(fd, func_name)))
                
    # Replace standard messagebox functions if tkinter.messagebox is used
    if 'tkinter.messagebox' in sys.modules:
        import tkinter.messagebox as mb
        for func_name in ['showinfo', 'showwarning', 'showerror', 'askquestion', 
                          'askokcancel', 'askyesno', 'askretrycancel', 'askyesnocancel']:
            if hasattr(mb, func_name):
                original_messagebox_functions[func_name] = getattr(mb, func_name)
                setattr(mb, func_name, wrap_dialog(getattr(mb, func_name)))

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
