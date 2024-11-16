The siril python module supports communication with a running Siril
instance. It can request representations of the current loaded image
plus its metadata, including details of detected stars, as well as the
current loaded sequence and most frame metadata.

It can also run Siril commands using the SirilInterface.cmd()
method, and the intent is to provide a capable interface for writing
advanced scripts for Siril to a level not possible with the previous
simple script files.

For example, scripts can now have GTK front ends using the Gtk module
from the gi package, and they can utilise the entire ecosystem of python
modules including numpy, scipy, pillow and many more.

In the initial module release, most methods relating to the image or
sequence loaded in Siril are read-only. The intent is that the parameters
of the loaded image or sequence can be obtained and used as inputs to
scripts, for example for calculations or input to conditionals, but the
existing mature Siril command set should in most cases be used to act on
the loaded image. Thus header keywords can be set using
``cmd("update_key", "key", "value")``, and most built-in image operations
can be carried out using the appropriate Siril command. The main exception
to the rule of python methods providing read-only access is the
``set_pixeldata()`` method, which allows for setting the pixel data in the
loaded image from a numpy array. This means that new pixel processing
algorithms can be added using python, initially getting the pixel data from
the loaded image using ``get_pixeldata()`` and setting it on completion
using ``set_pixeldata()``.
