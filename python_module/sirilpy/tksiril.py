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

    # Main window reference for transient relationship
    main_window = themed_tk

    # This allows making dialogs transient of the main window
    original_filedialog_functions = {}
    original_messagebox_functions = {}

    # Function to wrap dialog functions to make them transient
    def wrap_dialog_transient(original_func):
        def wrapper(*args, **kwargs):
            # Call the original dialog function
            result = original_func(*args, **kwargs)

            # For dialogs that create a Toplevel window and return it
            if isinstance(result, tk.Toplevel):
                try:
                    # Make it transient of the main window
                    result.transient(main_window)
                    # Bring it to front
                    result.lift()
                except tk.TclError:
                    pass  # Handle edge cases where window is already destroyed

            return result
        return wrapper

    # Hook into Toplevel creation for dialogs that don't return their window
    original_toplevel_init = tk.Toplevel.__init__
    def patched_toplevel_init(self, master=None, **kw):
        # Call the original constructor
        original_toplevel_init(self, master, **kw)

        # Skip transient if this is a special window (like our main window)
        if kw.get('class_') == 'Toplevel' and not kw.get('transient'):
            return

        # For standard dialogs, make them transient of the main window
        if on_top is True and main_window.winfo_exists():
            try:
                # Make this Toplevel transient of the main window
                self.transient(main_window)

                # Ensure it appears on top
                self.after(10, self.lift)
            except tk.TclError:
                pass  # Window might be gone

    # Set the patched initializer
    tk.Toplevel.__init__ = patched_toplevel_init

    # Replace standard dialog functions if tkinter.filedialog is used
    if 'tkinter.filedialog' in sys.modules:
        import tkinter.filedialog as fd
        for func_name in ['askopenfilename', 'askopenfilenames', 'asksaveasfilename', 
                         'askdirectory', 'askopenfile', 'askopenfiles', 'asksaveasfile']:
            if hasattr(fd, func_name):
                # We don't wrap these since they don't return Toplevel objects
                # The fix applies through our Toplevel.__init__ patch
                pass

    # Replace standard messagebox functions if tkinter.messagebox is used
    if 'tkinter.messagebox' in sys.modules:
        import tkinter.messagebox as mb
        for func_name in ['showinfo', 'showwarning', 'showerror', 'askquestion', 
                          'askokcancel', 'askyesno', 'askretrycancel', 'askyesnocancel']:
            if hasattr(mb, func_name):
                # We don't wrap these since they don't return Toplevel objects
                # The fix applies through our Toplevel.__init__ patch
                pass

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
