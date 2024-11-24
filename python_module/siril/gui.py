# Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
# Copyright (C) 2012-2024 team free-astro (see more in AUTHORS file)
# Reference site is https://siril.org
# SPDX-License-Identifier: GPL-3.0-or-later

import sys
if sys.platform == "darwin"
    import objc
    from Cocoa import NSAppplication, NSWindow

def check_gi_requirements(version: float = None):
   """
   Check if gi (GObject Introspection) and optionally specific GTK/GDK versions are available.

   Args:
       version (float, optional): GTK/GDK version to check for, e.g. 3.0 or 4.0.
           If provided, checks for specific version availability.
           If None, only checks basic gi import capability.

   Returns:
       bool: True if required gi components are available, False otherwise.

   Example:
       >>> check_gi_requirements()  # Basic gi check
       True
       >>> check_gi_requirements(3.0)  # Check GTK3/GDK3 availability
       True
   """
   try:
       import gi
       if version is not None:
           gi.require_version('Gtk', f'{version}')
           gi.require_version('Gdk', f'{version}')
       from gi import repository
       return True
   except (ImportError, ValueError):
       return False

def set_parent_window(self):
    """
    Sets the parent window for the dialog using the window ID from environment variables.
    Must be called after the dialog is initialized.
    This method retrieves the parent window ID from the SIRIL_PARENT_WINDOW environment
    variable and sets up the proper parent-child window relationship. It handles both
    Windows and Unix/Linux platforms differently:
    - Windows: Converts HWND from hex string and uses GdkWin32
    - Unix/Linux: Converts XID from decimal string and uses GdkX11
    Required environment variables:
        SIRIL_PARENT_WINDOW: Window ID of the parent window (HWND on Windows, XID on Unix)
    Required GTK/GDK version: 3.0
    Raises:
        SirilException: If GTK3/GDK3 components are not available
    Returns:
        None: The method returns silently if no parent window ID is found
    """
    if not check_gi_requirements(3.0):
        raise SirilException("Required GTK3/GDK3 components not available")
    import os
    from gi import repository as Gir
    parent_window_id = os.environ.get("SIRIL_PARENT_WINDOW")
    if not parent_window_id:
        return

    # Get the GDK window for our dialog
    dialog_window = self.get_window()
    if dialog_window is None:
        # Realize the window if it hasn't been yet
        self.realize()
        dialog_window = self.get_window()

    display = Gir.Gdk.Display.get_default()
    if os.name == 'nt':  # Windows
        # Convert hex string to integer
        hwnd = int(parent_window_id, 16)
        # Get GDK window from HWND
        parent = Gir.GdkWin32.Window.foreign_new_for_display(display, hwnd)
    else if sys.platform == 'darwin':
        parent_nswindow = objc.objc_object(long(parent_window_id))
        parent_gtk_window = Gtk.Window.get_toplevel(parent_nswindow)
        dialog.set_transient_for(parent_gtk_window)

    else if display.get_name == "x11":
        # Unix/Linux
        # Convert string to integer
        xid = int(parent_window_id)
        # Get GDK window from XID
        parent = Gir.GdkX11.X11Window.foreign_new_for_display(display, xid)
    else if display.get_name == "wayland":
        import GdkWayland
        # Retrieve wl_surface handle from environment variable
        wl_surface_handle = os.getenv("WL_SURFACE_HANDLE") if wl_surface_handle: # Convert string to pointer
        wl_surface_pointer = int(wl_surface_handle, 16)
        parent = Gir.GdkWayland.WaylandWindow.foreign_new_for_display(GdkWayland.Display.get_default(), wl_surface_pointer)
    if parent:
        # Set the parent window directly
        dialog_window.set_transient_for(parent)
