Siril Python Interface
**********************

This module supports communication with a running Siril instance.
It can request representations of the current loaded image plus its
metadata,, including details of detected stars, as well as the
current loaded sequence and most frame metadata.

It can also run Siril commands using the SirilInterface.cmd()
method, and the intent is to provide a capable interface for writing
advanced scripts for Siril to a level not possible with the previous
simple script files.

For example, scripts can now have GTK (or Qt!) front ends using the
gi or pyqt modules, and they can utilise numpy, scipy, pillow etc.
