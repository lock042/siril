# Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
# Copyright (C) 2012-2024 team free-astro (see more in AUTHORS file)
# Reference site is https://siril.org
# SPDX-License-Identifier: GPL-3.0-or-later

from dataclasses import dataclass, field
from datetime import datetime
from typing import Optional, Tuple, Union, List
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

    total: int = 0
    ngoodpix: int = 0
    mean: float = 0.0
    median: float = 0.0
    sigma: float = 0.0
    avgDev: float = 0.0
    mad: float = 0.0
    sqrtbwmv: float = 0.0
    location: float = 0.0
    scale: float = 0.0
    min: float = 0.0
    max: float = 0.0
    normValue: float = 0.0
    bgnoise: float = 0.0

@dataclass
class FKeywords:
    """
    Python equivalent of Siril fkeywords structure. Contains the FITS
    header keyword values converted to suitable datatypes.
    """

    # FITS file data
    bscale: float = 1.0
    bzero: float = 0.0
    lo: int = 0
    hi: int = 0
    flo: float = 0.0
    fhi: float = 0.0
    program: str = ""
    filename: str = ""
    data_max: float = 0.0
    data_min: float = 0.0
    pixel_size_x: float = 0.0
    pixel_size_y: float = 0.0
    binning_x: int = 1
    binning_y: int = 1
    row_order: str = ""

    # Date / time fields
    date: Optional[datetime] = None
    date_obs: Optional[datetime] = None
    expstart: float = 0.0
    expend: float = 0.0

    # Image metadata
    filter: str = ""
    image_type: str = ""
    object: str = ""
    instrume: str = ""
    telescop: str = ""
    observer: str = ""

    # Position data
    centalt: float = 0.0
    centaz: float = 0.0
    sitelat: float = 0.0
    sitelong: float = 0.0
    sitelat_str: str = ""
    sitelong_str: str = ""
    siteelev: float = 0.0

    # Camera settings
    bayer_pattern: str = ""
    bayer_xoffset: int = 0
    bayer_yoffset: int = 0
    airmass: float = 1.0
    focal_length: float = 0.0
    flength: float = 0.0
    iso_speed: float = 0.0
    exposure: float = 0.0
    aperture: float = 0.0
    ccd_temp: float = 0.0
    set_temp: float = 0.0
    livetime: float = 0.0
    stackcnt: int = 0
    cvf: float = 0.0
    key_gain: int = 0
    key_offset: int = 0

    # Focuser data
    focname: str = ""
    focuspos: int = 0
    focussz: int = 0
    foctemp: float = 0.0

    # Setters
    @bscale.setter
    def bscale(self, value: float) -> None:
        """
        Set FITS BSCALE value. This keyword shall be used,
        along with the BZERO keyword, when the array pixel values are
        not the true physical values, to transform the primary data
        array values to the true physical values they represent, using
        the equation: physical_value = BZERO + BSCALE * array_value.
        The value field shall contain a floating point number
        representing the coefficient of the linear term in the scaling
        equation, the ratio of physical value to array value at zero
        offset. The default value for this keyword is 1.0.
        """
        self.bscale = value

    @bzero.setter
    def bzero(self, value: float) -> None:
        """
        Sets FITS BZERO value. This keyword shall be used,
        along with the BSCALE keyword, when the array pixel values are
        not the true physical values, to transform the primary data
        array values to the true values using the equation:
        physical_value = BZERO + BSCALE * array_value. The value field
        shall contain a floating point number representing the physical
        value corresponding to an array value of zero.  The default value
        for this keyword is 0.0.
        """
        self.bzero = value

    @lo.setter
    def lo(self, value: int) -> None:
        """Sets lo value. This is the lower visualization cutoff."""
        self.lo = value

    @hi.setter
    def hi(self, value: int) -> None:
        """
        Sets hi value. This is the upper visualization cutoff.
        """
        self.hi = value

    @flo.setter
    def flo(self, value: float) -> None:
        """
        Sets floating point lo value. This is the lower visualization cutoff.
        """
        self.flo = value

    @fhi.setter
    def fhi(self, value: float) -> None:
        """
        Sets floating point hi value. This is the lower visualization cutoff.
        """
        self.fhi = value

    @program.setter
    def program(self, value: str) -> None:
        """Sets the PROGRAM keyword."""
        self.program = value[:70]

    @filename.setter
    def filename(self, value: str) -> None:
        """Sets the FILENAME keyword."""
        self.filename = value[:70]

    @data_max.setter
    def data_max(self, value: float) -> None:
        """
        Sets data_max. This is used to check if 32b float is
        in the [0, 1] range; it should not normally need to be adjusted
        from python scripts, as it is an indicator of the value range of
        the data in a saved FITS.
        """
        self.data_max = value

    @data_min.setter
    def data_min(self, value: float) -> None:
        """
        Sets data_min. This is used to check if 32b float is
        in the [0, 1] range; it should not normally need to be adjusted
        from python scripts, as it is an indicator of the value range of
        the data in a saved FITS.
        """
        self.data_min = value

    @pixel_size_x.setter
    def pixel_size_x(self, value: float) -> None:
        """
        Sets the pixel size in the x dimension. This is normally
        set by the capture software.
        """
        if value > 0:
            self.pixel_size_x = value
        else:
            raise ValueError(_("pixel_size_x must be greater than 0"))

    @pixel_size_y.setter
    def pixel_size_y(self, value: float) -> None:
        """
        Sets the pixel size in the y dimension. This is normally
        set by the capture software.
        """
        if value > 0:
            self.pixel_size_y = value
        else:
            raise ValueError(_("pixel_size_y must be greater than 0"))

    @binning_x.setter
    def binning_x(self, value: int) -> None:
        """
        Sets the binning in the x dimension. This is normally
        set by the capture software.
        """
        if value >= 1:
            self.binning_x = value
        else:
            raise ValueError(_("binning_x must be greater than or equal to 1"))

    @binning_y.setter
    def binning_y(self, value: int) -> None:
        """
        Sets the binning in the y dimension. This is normally
        set by the capture software.
        """
        if value >= 1:
            self.binning_y = value
        else:
            raise ValueError(_("binning_y must be greater than or equal to 1"))

    @row_order.setter
    def row_order(self, value: str) -> None:
        """
        Sets the ROWORDER keyword. This sets the row order in
        which the sensor indicates rows and should be either 'TOP-DOWN'
        or 'BOTTOM-UP'
        """
        self.row_order = value[:70]

    @date.setter
    def date(self, value: Optional[datetime]) -> None:
        """
        Sets the FITS DATE keyword, which represents the date on
        which the HDU was created.
        """
        self.date = value

    @date_obs.setter
    def date_obs(self, value: Optional[datetime]) -> None:
        """
        Sets the FITS DATE-OBS keyword, which represents the date
        on which the observation was made.
        """
        self.date_obs = value

    @expstart.setter
    def expstart(self, value: float) -> None:
        """Sets the start of the exposure as a Julian date."""
        self.expstart = value

    @expend.setter
    def expend(self, value: float) -> None:
        """Sets the end of the exposure as a Julian date."""
        self.expend = value

    @filter.setter
    def filter(self, value: str) -> None:
        """
        Sets the FITS FILTER keyword. Will be truncated
        to 70 characters according to the FITS standard.
        """
        self.filter = value[:70]

    @image_type.setter
    def image_type(self, value: str) -> None:
        """
        Sets the FITS IMAGETYP keyword. Will be truncated
        to 70 characters according to the FITS standard.
        """
        self.image_type = value[:70]

    @object.setter
    def object(self, value: str) -> None:
        """
        Sets the FITS OBJECT keyword. Will be truncated
        to 70 characters according to the FITS standard.
        """
        self.object = value[:70]

    @instrume.setter
    def set_instrume(self, value: str) -> None:
        """
        Sets the FITS INSTRUME keyword. Will be truncated "
        to 70 characters according to the FITS standard.
        """
        self.instrume = value[:70]

    @telescop.setter
    def telescop(self, value: str) -> None:
        """
        Sets the FITS TELESCOP keyword. Will be truncated
        to 70 characters according to the FITS standard.
        """
        self.telescop = value[:70]

    @observer.setter
    def observer(self, value: str) -> None:
        """
        Sets the FITS OBSERVER keyword. Will be truncated
        to 70 characters according to the FITS standard.
        """
        self.observer = value[:70]

    @centalt.setter
    def set_centalt(self, value: float) -> None:
        """Sets the centre altitude of the image."""
        if value <= 90:
            self.centalt = value
        else:
            raise ValueError(_("centalt must be less than or equal to 90 degrees"))

    @centaz.setter
    def centaz(self, value: float) -> None:
        """Sets the centre azimuth of the image."""
        if 0 <= value < 360:
            self.centaz = value
        else:
            raise ValueError(_("centaz must be between 0 and 360 degrees (exclusive)"))

    @sitelat.setter
    def sitelat(self, value: float) -> None:
        """Sets the site latitude."""
        if -90 <= value <= 90:
            self.sitelat = value
        else:
            raise ValueError(_("sitelat must be between -90 and 90 degrees"))

    @sitelong.setter
    def sitelong(self, value: float) -> None:
        """Sets the site longitude of the image."""
        if 0 <= value < 360:
            self.sitelong = value
        else:
            raise ValueError(_("sitelong must be between 0 and 360 degrees (exclusive)"))

    @sitelat_str.setter
    def sitelat_str(self, value: str) -> None:
        """Sets the site latitude of the image as a string."""
        self.sitelat_str = value[:70]

    @sitelong_str.setter
    def sitelong_str(self, value: str) -> None:
        """Sets the site longitude of the image as a string."""
        self.sitelong_str = value[:70]

    @siteelev.setter
    def siteelev(self, value: float) -> None:
        """Sets the site elevation of the image."""
        self.siteelev = value

    @bayer_pattern.setter
    def bayer_pattern(self, value: str) -> None:
        """
        Sets the image CFA pattern. This may be a Bayer pattern
        or X-TRANS pattern. Normally this is set by the capture software.
        """
        self.bayer_pattern = value[:70]

    @bayer_xoffset.setter
    def bayer_xoffset(self, value: int) -> None:
        """Sets the image Bayer x offset."""
        self.bayer_xoffset = value

    @bayer_yoffset.setter
    def bayer_yoffset(self, value: int) -> None:
        """Sets the image Bayer y offset."""
        self.bayer_yoffset = value

    @airmass.setter
    def airmass(self, value: float) -> None:
        """Sets the airmass as a floating point value."""
        if value >= 1:
            self.airmass = value
        else:
            raise ValueError(_("airmass must be greater than or equal to 1"))

    @focal_length.setter
    def focal_length(self, value: float) -> None:
        """Sets the focal length as a floating point value."""
        self.focal_length = value

    @flength.setter
    def flength(self, value: float) -> None:
        """Sets the flength as a floating point value."""
        self.flength = value

    @iso_speed.setter
    def iso_speed(self, value: float) -> None:
        """Sets the ISO speed as a floating point value."""
        self.iso_speed = value

    @exposure.setter
    def exposure(self, value: float) -> None:
        """Sets the exposure time as a floating point value."""
        self.exposure = value

    @aperture.setter
    def aperture(self, value: float) -> None:
        """Sets the aperture as a floating point value."""
        self.aperture = value

    @ccd_temp.setter
    def ccd_temp(self, value: float) -> None:
        """Sets the CCD temperature as a floating point value."""
        self.ccd_temp = value

    @set_temp.setter
    def set_temp(self, value: float) -> None:
        """Sets the CCD set temperature as a floating point value."""
        self.set_temp = value

    @livetime.setter
    def livetime(self, value: float) -> None:
        """
        Sets the sum of the exposure times (s) as a floating
        point value. Only relevant to stacked images.
        """
        self.livetime = value

    @stackcnt.setter
    def stackcnt(self, value: int) -> None:
        """
        Sets the count of images in a stack. Only relevant to
        stacked images.
        """
        self.stackcnt = value

    @cvf.setter
    def cvf(self, value: float) -> None:
        """Sets the conversion factor (e-/ADU)."""
        self.cvf = value

    @key_gain.setter
    def key_gain(self, value: int) -> None:
        """
        Sets the gain value read in camera headers. Normally
        this is set by capture software.
        """
        self.key_gain = value

    @key_offset.setter
    def key_offset(self, value: int) -> None:
        """
        Sets the offset value read in camera headers. Normally
        this is set by capture software.
        """
        self.key_offset = value

    @focname.setter
    def focname(self, value: str) -> None:
        """
        Sets the focuser name. Will be truncated
        to 70 characters according to the FITS standard.
        """
        self.focname = value[:70]

    @focuspos.setter
    def focuspos(self, value: int) -> None:
        """Sets the focus position. Normally this is set by
        capture software."""
        self.focuspos = value

    @focussz.setter
    def focussz(self, value: int) -> None:
        """Sets the focuser step size. Normally this is set by
        capture software."""
        self.focussz = value

    @foctemp.setter
    def foctemp(self, value: float) -> None:
        """Sets the focuser temperature. Normally this is set by
        capture software."""
        self.foctemp = value

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
    B: float = 0.0              # average sky background value
    A: float = 0.0              # amplitude
    x0: float = 0.0            # coordinates of the peak
    y0: float = 0.0
    sx: float = 0.0            # Size of the fitted function on the x and y axis in PSF coordinates
    sy: float = 0.0
    fwhmx: float = 0.0         # FWHM in x and y axis
    fwhmy: float = 0.0
    fwhmx_arcsec: float = 0.0  # FWHM in x and y axis in arc second
    fwhmy_arcsec: float = 0.0
    angle: float = 0.0         # angle of the axis x,y with respect to the image's
    rmse: float = 0.0          # RMSE of the minimization
    sat: float = 0.0           # Level above which pixels have satured
    R: int = 0                 # Optimized box size to enclose sufficient pixels in the background
    has_saturated: bool = False

    # Moffat parameters
    beta: float = 0.0          # Moffat equation beta parameter
    profile: StarProfile = StarProfile.GAUSSIAN  # Whether profile is Gaussian or Moffat

    xpos: float = 0.0          # position of the star in the image
    ypos: float = 0.0

    # photometry data
    mag: float = 0.0           # (V)magnitude, approximate or accurate
    Bmag: float = 0.0          # B magnitude
    s_mag: float = 999.99      # error on the (V)magnitude
    s_Bmag: float = 999.99     # error on the B magnitude
    SNR: float = 0.0           # SNR of the star
    BV: float = 0.0            # only used to pass data in photometric color calibration

    # uncertainties
    B_err: float = 0.0
    A_err: float = 0.0
    x_err: float = 0.0
    y_err: float = 0.0
    sx_err: float = 0.0
    sy_err: float = 0.0
    ang_err: float = 0.0
    beta_err: float = 0.0

    layer: int = 0
    units: Optional[str] = None
    ra: float = 0.0            # Right Ascension
    dec: float = 0.0           # Declination

@dataclass
class RegData:
    """Python equivalent of Siril regdata structure"""
    fwhm: float = 0.0                    # copy of fwhm->fwhmx, used as quality indicator
    weighted_fwhm: float = 0.0           # used to exclude spurious images
    roundness: float = 0.0               # fwhm->fwhmy / fwhm->fwhmx, 0 when uninit, ]0, 1] when set
    quality: float = 0.0
    background_lvl: float = 0.0
    number_of_stars: int = 0
    H: Homography = field(default_factory=Homography)

@dataclass
class ImgData:
    """Python equivalent of Siril imgdata structure"""
    filenum: int = 0              # real file index in the sequence
    incl: bool = False           # selected in the sequence
    date_obs: Optional[datetime] = None  # date of the observation
    airmass: float = 0.0         # airmass of the image
    rx: int = 0
    ry: int = 0

@dataclass
class Sequence:
    """Python equivalent of Siril sequ structure"""
    seqname: str = ""                    # name of the sequence
    number: int = 0                      # number of images in the sequence
    selnum: int = 0                      # number of selected images
    fixed: int = 0                       # fixed length of image index in filename
    nb_layers: int = -1                  # number of layers embedded in each image file
    rx: int = 0                          # first image width
    ry: int = 0                          # first image height
    is_variable: bool = False            # sequence has images of different sizes
    bitpix: int = 0                      # image pixel format, from fits
    reference_image: int = 0             # reference image for registration
    imgparam: List[ImgData] = None       # a structure for each image of the sequence
    regparam: List[List[RegData]] = None # registration parameters for each layer
    stats: List[List[ImageStats]] = None  # statistics of the images for each layer
    beg: int = 0                         # imgparam[0]->filenum
    end: int = 0                         # imgparam[number-1]->filenum
    exposure: float = 0.0                # exposure of frames
    fz: bool = False
    type: SequenceType = None
    cfa_opened_monochrome: bool = False  # CFA SER opened in monochrome mode
    current: int = 0                     # file number currently loaded
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
