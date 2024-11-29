# Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
# Copyright (C) 2012-2024 team free-astro (see more in AUTHORS file)
# Reference site is https://siril.org
# SPDX-License-Identifier: GPL-3.0-or-later

from dataclasses import dataclass, field
from datetime import datetime
from typing import Optional, Tuple, List
import numpy as np
from enum import IntEnum
from .translations import _

class DataType(IntEnum):
    """
    Mimics the Siril data_type enum. Note that although Siril can
    handle opening FITS files of any data type, internally it processes
    images only as USHORT_IMG (uint16) or FLOAT_IMG (float32).
    """
    BYTE_IMG = 8
    SHORT_IMG = 16
    USHORT_IMG = 16
    LONG_IMG = 32
    FLOAT_IMG = 32
    DOUBLE_IMG = 64

@dataclass
class ImageStats:
    """
    Python equivalent of Siril imstats structure. Contains statistics for
    a particular channel of a Siril image.
    """

    total: int = 0  #: total number of pixels
    ngoodpix: int = 0   #: number of non-zero pixels
    mean: float = 0.0   #: mean value of pixels
    median: float = 0.0 #: median value of pixels
    sigma: float = 0.0  #: standard deviation of pixels
    avgDev: float = 0.0 #: average deviation of pixels
    mad: float = 0.0    #: mean average deviation of pixels
    sqrtbwmv: float = 0.0   #: square root of the biweight midvariance of pixel values
    location: float = 0.0   #: location of pixel values
    scale: float = 0.0  #: scale value of the pixels
    min: float = 0.0    #: minimum pixel value
    max: float = 0.0    #: maximum pixel value
    normValue: float = 0.0  #: norm value of the pixels
    bgnoise: float = 0.0    #: RMS background noise

from dataclasses import dataclass
from typing import Optional
from datetime import datetime
from gettext import gettext as _

@dataclass
class FKeywords:
    """
    Python equivalent of Siril fkeywords structure. Contains the FITS
    header keyword values converted to suitable datatypes.
    """

    # FITS file data
    bscale: float = 1.0
    bzero: float = 0.0
    lo: int = 0 #: MIPS-LO key in FITS file, "Lower visualization cutoff"
    hi: int = 0 #: MIPS-HI key in FITS file, "Upper visualization cutoff"
    flo: float = 0.0 #: MIPS-LO key in FITS file, "Lower visualization cutoff (float)"
    fhi: float = 0.0  #: MIPS-Hi key in FITS file, "Upper visualization cutoff (float)"

    # Private string attributes with corresponding properties
    _program: str = ""
    _filename: str = ""
    _row_order: str = ""
    _filter: str = ""
    _image_type: str = ""
    _object: str = ""
    _instrume: str = ""
    _telescop: str = ""
    _observer: str = ""
    _bayer_pattern: str = ""
    _sitelat_str: str = ""
    _sitelong_str: str = ""
    _focname: str = ""

    # Datetime attributes
    _date: Optional[datetime] = None
    _date_obs: Optional[datetime] = None

    # Remaining attributes
    data_max: float = 0.0 #: used to check if 32b float is in the [0, 1] range
    data_min: float = 0.0  #: used to check if 32b float is in the [0, 1] range
    pixel_size_x: float = 0.0 #: XPIXSZ FITS header card as a float
    pixel_size_y: float = 0.0 #: YPIXSZ FITS header card as a float
    binning_x: int = 1 #: XBINNING FITS header card as an int
    binning_y: int = 1 #: YBINNING FITS header card as an int

    expstart: float = 0.0 #: Exposure start as a Julian date
    expend: float = 0.0 #: Exposure end as a Julian date

    # Camera settings
    bayer_xoffset: int = 0 #: X offset of the Bayer pattern
    bayer_yoffset: int = 0 #: Y offset of the Bayer pattern
    _airmass: float = 1.0
    focal_length: float = 0.0 #: focal length
    flength: float = 0.0
    iso_speed: float = 0.0 #: ISO speed value as a float
    exposure: float = 0.0 #: Exposure time as a float (s)
    aperture: float = 0.0 #: Aperture value as a float
    ccd_temp: float = 0.0 #: CCD temperature as a float
    set_temp: float = 0.0 #: CCD set temperature as a float
    livetime: float = 0.0 #: Sum of exposure times (s)
    stackcnt: int = 0 #: Number of stacked frames
    cvf: float = 0.0 #: Conversion factor (e- / ADU)
    key_gain: int = 0 #: Gain factor read in camera
    key_offset: int = 0 #: Offset value read in camera

    # Focuser data
    focuspos: int = 0 #: Focuser position
    focussz: int = 0
    foctemp: float = 0.0 #: Focuser temperature

    # Position data
    _centalt: float = 0.0
    _centaz: float = 0.0
    _sitelat: float = 0.0
    _sitelong: float = 0.0
    siteelev: float = 0.0

    # Properties with setters

    @property
    def program(self) -> str:
        """
        Gets the PROGRAM keyword.
        """
        return self._program

    @program.setter
    def program(self, value: str) -> None:
        """
        Sets the PROGRAM keyword.
        """
        self._program = value[:70]

    @property
    def filename(self) -> str:
        """
        Gets the FILENAME keyword.
        """
        return self._filename

    @filename.setter
    def filename(self, value: str) -> None:
        """
        Sets the FILENAME keyword.
        """
        self._filename = value[:70]

    @property
    def row_order(self) -> str:
        """
        Gets the ROWORDER keyword. This sets the row order in
        which the sensor indicates rows and should be either 'TOP-DOWN'
        or 'BOTTOM-UP'
        """
        return self._row_order

    @row_order.setter
    def row_order(self, value: str) -> None:
        """
        Sets the ROWORDER keyword. This sets the row order in
        which the sensor indicates rows and should be either 'TOP-DOWN'
        or 'BOTTOM-UP'
        """
        self._row_order = value[:70]

    @property
    def date(self) -> Optional[datetime]:
        """
        Gets the FITS DATE keyword, which represents the date on
        which the HDU was created.
        """
        return self._date

    @date.setter
    def date(self, value: Optional[datetime]) -> None:
        """
        Sets the FITS DATE keyword, which represents the date on
        which the HDU was created.
        """
        self._date = value

    @property
    def date_obs(self) -> Optional[datetime]:
        """
        Gets the FITS DATE-OBS keyword, which represents the date
        on which the observation was made.
        """
        return self._date_obs

    @date_obs.setter
    def date_obs(self, value: Optional[datetime]) -> None:
        """
        Sets the FITS DATE-OBS keyword, which represents the date
        on which the observation was made.
        """
        self._date_obs = value

    @property
    def filter(self) -> str:
        """
        Gets the FITS FILTER keyword. Will be truncated
        to 70 characters according to the FITS standard.
        """
        return self._filter

    @filter.setter
    def filter(self, value: str) -> None:
        """
        Sets the FITS FILTER keyword. Will be truncated
        to 70 characters according to the FITS standard.
        """
        self._filter = value[:70]

    @property
    def image_type(self) -> str:
        """
        Gets the FITS IMAGETYP keyword. Will be truncated
        to 70 characters according to the FITS standard.
        """
        return self._image_type

    @image_type.setter
    def image_type(self, value: str) -> None:
        """
        Sets the FITS IMAGETYP keyword. Will be truncated
        to 70 characters according to the FITS standard.
        """
        self._image_type = value[:70]

    @property
    def object(self) -> str:
        """
        Gets the FITS OBJECT keyword. Will be truncated
        to 70 characters according to the FITS standard.
        """
        return self._object

    @object.setter
    def object(self, value: str) -> None:
        """
        Sets the FITS OBJECT keyword. Will be truncated
        to 70 characters according to the FITS standard.
        """
        self._object = value[:70]

    @property
    def instrume(self) -> str:
        """
        Gets the FITS INSTRUME keyword. Will be truncated "
        to 70 characters according to the FITS standard.
        """
        return self._instrume

    @instrume.setter
    def instrume(self, value: str) -> None:
        """
        Sets the FITS INSTRUME keyword. Will be truncated "
        to 70 characters according to the FITS standard.
        """
        self._instrume = value[:70]

    @property
    def telescop(self) -> str:
        """
        Gets the FITS TELESCOP keyword. Will be truncated
        to 70 characters according to the FITS standard.
        """
        return self._telescop

    @telescop.setter
    def telescop(self, value: str) -> None:
        """
        Sets the FITS TELESCOP keyword. Will be truncated
        to 70 characters according to the FITS standard.
        """
        self._telescop = value[:70]

    @property
    def observer(self) -> str:
        """
        Gets the FITS OBSERVER keyword. Will be truncated
        to 70 characters according to the FITS standard.
        """
        return self._observer

    @observer.setter
    def observer(self, value: str) -> None:
        """
        Sets the FITS OBSERVER keyword. Will be truncated
        to 70 characters according to the FITS standard.
        """
        self._observer = value[:70]

    @property
    def centalt(self) -> float:
        """
        Gets the centre altitude of the image as a float.
        """
        return self._centalt

    def set_centalt(self, value: float) -> None:
        """
        Sets the centre altitude of the image as a float.
        """
        if value <= 90:
            self._centalt = value
        else:
            raise ValueError(_("centalt must be less than or equal to 90 degrees"))

    @property
    def centaz(self) -> float:
        """
        Gets the centre azimuth of the image as a float.
        """
        return self._centaz

    def set_centaz(self, value: float) -> None:
        """
        Sets the centre azimuth of the image as a float.
        """
        if 0 <= value < 360:
            self._centaz = value
        else:
            raise ValueError(_("centaz must be between 0 and 360 degrees (exclusive)"))

    @property
    def sitelat(self) -> float:
        """
        Gets the site latitude keyword value as a float.
        """
        return self._sitelat

    def set_sitelat(self, value: float) -> None:
        """
        Sets the site latitude keyword value as a float.
        """
        if -90 <= value <= 90:
            self._sitelat = value
        else:
            raise ValueError(_("sitelat must be between -90 and 90 degrees"))

    @property
    def sitelong(self) -> float:
        """
        Gets the site longitude keyword value as a float.
        """
        return self._sitelong

    def set_sitelong(self, value: float) -> None:
        """
        Sets the site longitude keyword value as a float.
        """
        if 0 <= value < 360:
            self._sitelong = value
        else:
            raise ValueError(_("sitelong must be between 0 and 360 degrees (exclusive)"))

    @property
    def sitelat_str(self) -> str:
        """
        Gets the site latitude keyword value of the image as a string.
        """
        return self._sitelat_str

    @sitelat_str.setter
    def sitelat_str(self, value: str) -> None:
        """
        Sets the site latitude keyword value of the image as a string.
        """
        self._sitelat_str = value[:70]

    @property
    def sitelong_str(self) -> str:
        """Gets the site longitude keyword value of the image as a string."""
        return self._sitelong_str

    @sitelong_str.setter
    def sitelong_str(self, value: str) -> None:
        """Sets the site longitude keyword value of the image as a string."""
        self._sitelong_str = value[:70]

    @property
    def bayer_pattern(self) -> str:
        """
        Gets the image CFA pattern. This may be a Bayer pattern
        or X-TRANS pattern. Normally this is set by the capture software.
        """
        return self._bayer_pattern

    @bayer_pattern.setter
    def bayer_pattern(self, value: str) -> None:
        """
        Sets the image CFA pattern. This may be a Bayer pattern
        or X-TRANS pattern. Normally this is set by the capture software.
        """
        self._bayer_pattern = value[:70]

    @property
    def airmass(self) -> float:
        """Gets the airmass keyword value as a floating point value."""
        return self._airmass

    @airmass.setter
    def set_airmass(self, value: float) -> None:
        """Sets the airmass keyword value as a floating point value."""
        if value >= 1:
            self._airmass = value
        else:
            raise ValueError(_("airmass must be greater than or equal to 1"))

    @property
    def focname(self) -> str:
        """
        Gets the focuser name. Will be truncated
        to 70 characters according to the FITS standard.
        """
        return self._focname

    @focname.setter
    def focname(self, value: str) -> None:
        """
        Sets the focuser name. Will be truncated
        to 70 characters according to the FITS standard.
        """
        self._focname = value[:70]

@dataclass
class FFit:
    """
    Python equivalent of Siril ffit (FITS) structure, holding image
    pixel data and metadata.
    """
    bitpix: int = 0
    orig_bitpix: int = 0
    naxis: int = 0
    _naxes: Tuple[int, int, int] = (0, 0, 0)

    keywords: FKeywords = field(default_factory=FKeywords)
    checksum: bool = False

    header: Optional[str] = None
    unknown_keys: Optional[str] = None

    stats: Tuple[Optional[ImageStats], Optional[ImageStats], Optional[ImageStats]] = \
        (None, None, None)
    mini: float = 0.0
    maxi: float = 0.0
    neg_ratio: float = 0.0

    type: DataType = DataType.FLOAT_IMG
    _data: Optional[np.ndarray] = None

    top_down: bool = False
    focalkey: bool = False
    pixelkey: bool = False

    history: list[str] = field(default_factory=list)

    color_managed: bool = False
    _icc_profile: Optional[bytes] = None

    def __post_init__(self):
        """Initialize after creation"""
        if self.history is None:
            self.history = []
        if all(stat is None for stat in self.stats):
            self.stats = tuple(ImageStats() for _ in range(3))

    @property
    def data(self) -> Optional[np.ndarray]:
        """
        Get the pixel data of the current image loaded in Siril
        as a NumPy array
        """
        return self._data

    @data.setter
    def data(self, value: Optional[np.ndarray]):
        """
        Set the pixel data of the FFit to the provided NumPy array. Note: this
        does not update the image loaded in Siril - ``SirilInterface.set_pixeldata()``
        must be used to achieve this.

        Args:
            numpy.ndarray: NumPy array representing the pixel data. This must
                           be in either uint16 or float32 planar format

        Raises:
            ValueError: if the array shape is incompatible
        """

        if value is not None:
            shape = value.shape
            if len(shape) == 2:
                self._naxes = (shape[1], shape[0], 1)  # width, height, channels
            elif len(shape) == 3:
                if shape[2] not in (1, 3):
                    raise ValueError(_("Third dimension must be 1 or 3"))
                self._naxes = (shape[1], shape[0], shape[2])  # width, height, channels
            else:
                raise ValueError(_("Data must be 2D or 3D"))
            self._update_naxis()
        self._data = value

    def _update_naxis(self):
        """Update naxis based on naxes value"""
        self.naxis = 3 if self._naxes[2] == 3 else 2

    @property
    def naxes(self) -> Tuple[int, int, int]:
        """Get the naxes tuple"""
        return self._naxes

    @property
    def width(self) -> int:
        """Get image width"""
        return self._naxes[0]

    @property
    def height(self) -> int:
        """Get image height"""
        return self._naxes[1]

    @property
    def channels(self) -> int:
        """Get image channels"""
        return self._naxes[2]

    @property
    def icc_profile(self) -> Optional[bytes]:
        """
        Get the ICC profile as raw bytes data. This may be converted
        for use by modules such as pillow which can handle ICC profiles.
        """
        return self._icc_profile

    @icc_profile.setter
    def icc_profile(self, value: Optional[bytes]):
        """
        Set ICC profile and update color_managed flag. Note this only
        updates the python FFit structure: the API does not currently
        support setting the ICC profile of the currently loaded FITS image.
        """
        self._icc_profile = value
        self.color_managed = value is not None

    @property
    def dtype(self) -> np.dtype:
        """Get the NumPy dtype based on the current type"""
        return np.uint16 if self.type == DataType.USHORT_IMG else np.float32

    def allocate_data(self):
        """
        Allocate memory for image data with appropriate type. self.rx, self.ry,
        self.naxis, self.naxes and self.dtype must be set before calling this
        method.
        """
        shape = (self.ry, self.rx) if self.naxis == 2 else (self.ry, self.rx, self.naxes[2])
        self.data = np.zeros(shape, dtype=self.dtype)

    def ensure_data_type(self, target_type=None):
        """
        Ensure data is in the correct type with proper scaling

        Args:
            target_type: Optional type to convert to. Can be either DataType or np.dtype.
                         If None, uses self.type

        Raises:
            ValueError: if the conversion is between data types that are not
                        internally used by Siril for calculation
        """
        if self.data is None:
            return

        # Handle input type and determine target dtype
        if target_type is None:
            type_to_use = self.type
            dtype_to_use = np.float32 if type_to_use is DataType.FLOAT_IMG else np.uint16
        elif target_type in (np.float32, np.uint16):
            dtype_to_use = target_type
            type_to_use = DataType.FLOAT_IMG if dtype_to_use == np.float32 else DataType.USHORT_IMG
        elif isinstance(target_type, np.dtype):
            raise ValueError(f"Unsupported type conversion from {self.data.dtype} to {dtype_to_use}")
        else:  # Assume DataType
            type_to_use = target_type
            dtype_to_use = np.float32 if type_to_use is DataType.FLOAT_IMG else np.uint16

        if self.data.dtype == dtype_to_use:
            return

        if dtype_to_use == np.float32 and self.data.dtype == np.uint16:
            # Convert from USHORT (0-65535) to FLOAT (0.0-1.0)
            self._data = self.data.astype(np.float32) / 65535.0
            self.type = DataType.FLOAT_IMG
        elif dtype_to_use == np.uint16 and self.data.dtype == np.float32:
            # Convert from FLOAT (0.0-1.0) to USHORT (0-65535)
            self._data = (self.data * 65535.0).clip(0, 65535).astype(np.uint16)
            self.type = DataType.USHORT_IMG
        else:
            raise ValueError(f"Unsupported type conversion from {self.data.dtype} to {dtype_to_use}")

    def get_channel(self, channel: int) -> np.ndarray:
        """
        Get a specific channel of the pixel data. Note that this does
        not pull pixel data directly from the image loaded in Siril: that must
        previously have been obtained using get_pixel_data() or get_image()
        """
        if self.data is None:
            raise ValueError(_("No data allocated"))
        if self.naxis == 2:
            if channel != 0:
                raise ValueError(_("Cannot get channel > 0 for 2D data"))
            return self.data
        return self.data[:, :, channel]

    def update_stats(self):
        """
        Update image statistics for all channels. Note that this only
        updates the statistics based on the NumPy array representing pixel data
        in the python FFit object, it does not update the statistics of the
        image in Siril.
        """
        if self.data is None:
            raise ValueError(_("No data allocated"))

        for i in range(self.naxes[2]):
            channel_data = self.get_channel(i)
            stats = self.stats[i] or ImageStats()

            # Update basic statistics
            stats.total = channel_data.size
            stats.ngoodpix = np.count_nonzero(channel_data)
            stats.mean = float(np.mean(channel_data))
            stats.median = float(np.median(channel_data))
            stats.sigma = float(np.std(channel_data))
            stats.min = float(np.min(channel_data))
            stats.max = float(np.max(channel_data))

            # More complex statistics
            deviations = np.abs(channel_data - stats.median)
            stats.mad = float(np.median(deviations))
            stats.avgDev = float(np.mean(deviations))

            self.stats[i] = stats

@dataclass
class Homography:
    """
    Python equivalent of the Siril Homography structure. Contains coefficients
    for the Homography matrix that maps a sequence frame onto the reference
    frame.
    """
    h00: float = 0.0
    h01: float = 0.0
    h02: float = 0.0
    h10: float = 0.0
    h11: float = 0.0
    h12: float = 0.0
    h20: float = 0.0
    h21: float = 0.0
    h22: float = 0.0
    pair_matched: int = 0
    Inliers: int = 0

class StarProfile(IntEnum):
    """
    Python equivalent of the Siril starprofile enum. Used to identify the type
    of fit used to model a star in the image. Note that MOFFAT_FIXED is currently
    not used in Siril, but is reserved for future use for Moffat stars modelled
    with a fixed beta parameter
    """
    GAUSSIAN = 0
    MOFFAT = 1
    MOFFAT_FIXED = 2

class SequenceType(IntEnum):
    """Python equivalent of the Siril sequence_type enum"""
    SEQ_REGULAR = 0
    SEQ_SER = 1
    SEQ_FITSEQ = 2
    SEQ_AVI = 3
    SEQ_INTERNAL = 4

@dataclass
class PSFStar:
    """
    Python equivalent of the Siril fwhm_struct structure. Contains
    data on a modelled fit to a star identified in the image.
    """
    star_name: Optional[str] = field(
        default=None,
        metadata={"doc": "Name or identifier of the star"}
    )
    B: float = 0.0              #: average sky background value
    A: float = 0.0              #: amplitude
    x0: float = 0.0            #: coordinates of the peak
    y0: float = 0.0
    sx: float = 0.0            #: Size of the fitted function on the x and y axis in PSF coordinates
    sy: float = 0.0
    fwhmx: float = 0.0         #: FWHM in x and y axis
    fwhmy: float = 0.0
    fwhmx_arcsec: float = 0.0  #: FWHM in x and y axis in arc second
    fwhmy_arcsec: float = 0.0
    angle: float = 0.0         #: angle of the axis x,y with respect to the image's
    rmse: float = 0.0          #: RMSE of the minimization
    sat: float = 0.0           #: Level above which pixels have satured
    R: int = 0                 #: Optimized box size to enclose sufficient pixels in the background
    has_saturated: bool = False

    # Moffat parameters
    beta: float = 0.0          #: Moffat equation beta parameter
    profile: StarProfile = StarProfile.GAUSSIAN  # Whether profile is Gaussian or Moffat

    xpos: float = 0.0          #: position of the star in the image
    ypos: float = 0.0

    # photometry data
    mag: float = 0.0           #: (V)magnitude, approximate or accurate
    Bmag: float = 0.0          #: B magnitude
    s_mag: float = 999.99      #: error on the (V)magnitude
    s_Bmag: float = 999.99     #: error on the B magnitude
    SNR: float = 0.0           #: SNR of the star
    BV: float = 0.0            #: only used to pass data in photometric color calibration

    # uncertainties
    B_err: float = 0.0 #: error in B
    A_err: float = 0.0 #: error in A
    x_err: float = 0.0 #: error in x
    y_err: float = 0.0 #: error in y
    sx_err: float = 0.0 #: error in sx
    sy_err: float = 0.0 #: error in sy
    ang_err: float = 0.0 #: error in angle
    beta_err: float = 0.0 #: error in beta

    layer: int = 0  #: image channel on which the star modelling was carried out
    units: Optional[str] = None
    ra: float = 0.0            #: Right Ascension
    dec: float = 0.0           #: Declination

@dataclass
class RegData:
    """Python equivalent of Siril regdata structure"""
    fwhm: float = 0.0                    #: copy of fwhm->fwhmx, used as quality indicator
    weighted_fwhm: float = 0.0           #: used to exclude spurious images
    roundness: float = 0.0               #: fwhm->fwhmy / fwhm->fwhmx, 0 when uninit, ]0, 1] when set
    quality: float = 0.0                 #: measure of image quality
    background_lvl: float = 0.0          #: background level
    number_of_stars: int = 0             #: number of stars detected in the image
    H: Homography = field(default_factory=Homography)

@dataclass
class ImgData:
    """Python equivalent of Siril imgdata structure"""
    filenum: int = 0              #: real file index in the sequence
    incl: bool = False           #: selected in the sequence
    _date_obs: Optional[datetime] = None  #: date of the observation
    _airmass: float = 0.0         #: airmass of the image
    rx: int = 0                 #: width
    ry: int = 0                 #: height

    @property
    def airmass(self) -> float:
        """Get the airmass value of the Siril imgdata structure"""
        return self._airmass

    @airmass.setter
    def set_airmass(self, value: float) -> None:
        """
        Get the airmass value of the Siril imgdata structure
        The value of airmass must be >= 1.0
        """
        if value >= 1:
            self._airmass = value
        else:
            raise ValueError(_("airmass must be greater than or equal to 1"))

    @property
    def date_obs(self) -> Optional[datetime]:
        """
        Gets the FITS DATE-OBS keyword, which represents the date
        on which the observation was made.
        """
        return self._date_obs

    @date_obs.setter
    def date_obs(self, value: Optional[datetime]) -> None:
        """
        Sets the FITS DATE-OBS keyword, which represents the date
        on which the observation was made.
        """
        self._date_obs = value

@dataclass
class Sequence:
    """Python equivalent of Siril sequ structure"""
    seqname: str = ""                    #: name of the sequence
    number: int = 0                      #: number of images in the sequence
    selnum: int = 0                      #: number of selected images
    fixed: int = 0                       #: fixed length of image index in filename
    nb_layers: int = -1                  #: number of layers embedded in each image file
    rx: int = 0                          #: first image width
    ry: int = 0                          #: first image height
    is_variable: bool = False            #: sequence has images of different sizes
    bitpix: int = 0                      #: image pixel format, from fits
    reference_image: int = 0             #: reference image for registration
    imgparam: List[ImgData] = None       #: a structure for each image of the sequence
    regparam: List[List[RegData]] = None #: registration parameters for each layer
    stats: List[List[ImageStats]] = None #: statistics of the images for each layer
    beg: int = 0                         #: imgparam[0]->filenum
    end: int = 0                         #: imgparam[number-1]->filenum
    exposure: float = 0.0                #: exposure of frames
    fz: bool = False
    type: SequenceType = None
    cfa_opened_monochrome: bool = False  #: CFA SER opened in monochrome mode
    current: int = 0                     #: file number currently loaded
    # The following fields are not currently implemented:
    # photometry: List[List[PSFStar]] = None  # psf for multiple stars
    # reference_star: int = 0              # reference star for apparent magnitude
    # reference_mag: float = 0.0           # reference magnitude for reference star
    # photometry_colors: List[List[float]] = None  # colors for each photometry curve

    def __post_init__(self):
        """Initialize lists that were set to None by default"""
        if self.imgparam is None:
            self.imgparam = []
        if self.regparam is None:
            self.regparam = []
        if self.stats is None:
            self.stats = []
