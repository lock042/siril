# Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
# Copyright (C) 2012-2026 team free-astro (see more in AUTHORS file)
# Reference site is https://siril.org
# SPDX-License-Identifier: GPL-3.0-or-later

"""
TKsiril submodule for Siril, providing utility methods to achieve consistent
script GUI appearance using the TKinter toolkit.
"""

import platform
from time import sleep
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

def match_theme_to_siril(themed_tk, s, on_top=False):
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

def elevate(root):
    """
    Raises the Tk root window to the top of the window stack. Useful after
    calls to sirilpy methods that present child windows of the main Siril
    window such as info_messagebox().

    NOTE: For this to work on KDE desktops, focus-stealing prevention must
    be disabled.
    """
    root.lift()
    root.focus_force()

class ScrollableFrame(ttk.Frame):
    """
    A scrollable frame widget that can contain other widgets.

    This class creates a frame with vertical scrolling capability using a Canvas
    widget and Scrollbar. It supports both scrollbar and mouse wheel scrolling
    across all platforms (Windows, Mac, Linux).

    Usage:
        scrollable = ScrollableFrame(parent)
        scrollable.pack(fill="both", expand=True)

        # Add widgets to scrollable.scrollable_frame
        label = ttk.Label(scrollable.scrollable_frame, text="Hello")
        label.pack()

        # Optionally bind mouse wheel to child widgets
        scrollable.add_mousewheel_binding(label)
    """

    def __init__(self, container, *args, **kwargs):
        """
        Initialize the ScrollableFrame.

        Args:
            container: The parent widget
            *args: Additional arguments passed to ttk.Frame
            **kwargs: Additional keyword arguments passed to ttk.Frame
        """
        super().__init__(container, *args, **kwargs)

        # Create canvas and scrollbar
        self.canvas = tk.Canvas(self, highlightthickness=0)
        self.scrollbar = ttk.Scrollbar(self, orient="vertical", command=self.canvas.yview)
        self.scrollable_frame = ttk.Frame(self.canvas)

        # Configure canvas to work with scrollbar
        self.canvas.configure(yscrollcommand=self.scrollbar.set)

        # Pack scrollbar and canvas
        self.scrollbar.pack(side="right", fill="y")
        self.canvas.pack(side="left", fill="both", expand=True)

        # Create window in canvas for the scrollable frame
        self.canvas_window = self.canvas.create_window((0, 0), window=self.scrollable_frame, anchor="nw")

        # Bind events
        self.scrollable_frame.bind("<Configure>", self._on_frame_configure)
        self.canvas.bind("<Configure>", self._on_canvas_configure)

        # Mouse wheel binding - this approach works more reliably
        self._setup_mousewheel()

    def _on_frame_configure(self, event):
        """
        Handle frame configuration changes.

        This method is called when the scrollable frame's size changes.
        It updates the canvas scroll region to encompass all the content.

        Args:
            event: The tkinter event object
        """
        self.canvas.configure(scrollregion=self.canvas.bbox("all"))

    def _on_canvas_configure(self, event):
        """
        Handle canvas configuration changes.

        This method is called when the canvas size changes (e.g., window resize).
        It adjusts the width of the scrollable frame to match the canvas width,
        preventing unwanted horizontal scrolling.

        Args:
            event: The tkinter event object containing the new canvas dimensions
        """
        canvas_width = event.width
        self.canvas.itemconfig(self.canvas_window, width=canvas_width)

    def _setup_mousewheel(self):
        """
        Setup cross-platform mouse wheel scrolling.

        This method configures mouse wheel event bindings that work across
        Windows, Mac, and Unix-like platforms. It uses platform detection to
        apply the appropriate event bindings:
        - Windows/Mac: <MouseWheel> with event.delta
        - Linux/BSD/Unix: <Button-4> (scroll up) and <Button-5> (scroll down)

        Platform-specific binding prevents conflicts where Button-4/5 might
        represent different mouse buttons on Windows/Mac systems.

        The bindings are applied when the mouse enters the widget area
        and removed when it leaves to avoid conflicts with other scrollable widgets.
        """
        system = platform.system()

        def _on_mousewheel_windows_mac(event):
            """Handle mouse wheel events on Windows and Mac"""
            self.canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")

        def _on_mousewheel_linux(event):
            """Handle mouse wheel events on Linux, BSD, and other Unix-like systems"""
            if event.num == 4:
                self.canvas.yview_scroll(-1, "units")
            elif event.num == 5:
                self.canvas.yview_scroll(1, "units")

        def _bind_mousewheel(event):
            """Bind platform-appropriate mouse wheel events"""
            if system in ("Windows", "Darwin"):  # Darwin is macOS
                self.canvas.bind_all("<MouseWheel>", _on_mousewheel_windows_mac)
            else:  # Linux, BSD, and other Unix-like systems
                self.canvas.bind_all("<Button-4>", _on_mousewheel_linux)
                self.canvas.bind_all("<Button-5>", _on_mousewheel_linux)

        def _unbind_mousewheel(event):
            """Unbind platform-appropriate mouse wheel events"""
            if system in ("Windows", "Darwin"):
                self.canvas.unbind_all("<MouseWheel>")
            else:  # Linux, BSD, and other Unix-like systems
                self.canvas.unbind_all("<Button-4>")
                self.canvas.unbind_all("<Button-5>")

        # Bind to canvas enter/leave events
        self.canvas.bind('<Enter>', _bind_mousewheel)
        self.canvas.bind('<Leave>', _unbind_mousewheel)

        # Also bind to the main frame
        self.bind('<Enter>', _bind_mousewheel)
        self.bind('<Leave>', _unbind_mousewheel)

    def add_mousewheel_binding(self, widget=None):
        """
        Add mouse wheel scrolling support to a widget and its children.

        This method recursively binds mouse wheel events to the specified widget
        and all its child widgets. It uses platform detection to apply the
        appropriate event bindings for each operating system.

        Args:
            widget: The tkinter widget to bind mouse wheel events to.
                   The binding will be applied recursively to all children.
                   If no widget is specified it will default to the ScrollableFrame
                   itself.

        Example:
            # Add a complex widget to the scrollable frame
            frame = ttk.Frame(scrollable.scrollable_frame)
            label = ttk.Label(frame, text="Hello")
            button = ttk.Button(frame, text="Click me")

            # Bind mouse wheel to the entire widget hierarchy
            scrollable.add_mousewheel_binding(frame)
        """
        if widget is None:
            widget = self

        system = platform.system()

        def _on_mousewheel_windows_mac(event):
            """Handle mouse wheel events on Windows and Mac"""
            self.canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")

        def _on_mousewheel_linux(event):
            """Handle mouse wheel events on Linux, BSD, and other Unix-like systems"""
            if event.num == 4:
                self.canvas.yview_scroll(-1, "units")
            elif event.num == 5:
                self.canvas.yview_scroll(1, "units")

        # Bind platform-appropriate events
        if system in ("Windows", "Darwin"):  # Darwin is macOS
            widget.bind("<MouseWheel>", _on_mousewheel_windows_mac)
        else:  # Linux, BSD, and other Unix-like systems
            widget.bind("<Button-4>", _on_mousewheel_linux)
            widget.bind("<Button-5>", _on_mousewheel_linux)

        # Recursively bind to all children
        for child in widget.winfo_children():
            self.add_mousewheel_binding(child)
