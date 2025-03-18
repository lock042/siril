# Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
# Copyright (C) 2012-2025 team free-astro (see more in AUTHORS file)
# Reference site is https://siril.org
# SPDX-License-Identifier: GPL-3.0-or-later

from dataclasses import dataclass, field
from datetime import datetime
from typing import Optional, Tuple, List
from enum import IntEnum
import struct
import logging
import numpy as np
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

@dataclass
class FKeywords:
    """
    Python equivalent of Siril fkeywords structure. Contains the FITS
    header keyword values converted to suitable datatypes.
    """

    # FITS file data
    bscale: float = 1.0 #: Offset data range to that of unsigned short
    bzero: float = 0.0 #: Default scaling factor
    lo: int = 0 #: MIPS-LO key in FITS file, "Lower visualization cutoff"
    hi: int = 0 #: MIPS-HI key in FITS file, "Upper visualization cutoff"
    flo: np.float32 = 0.0 #: MIPS-LO key in FITS file, "Lower visualization cutoff (float)"
    fhi: np.float32 = 0.0  #: MIPS-Hi key in FITS file, "Upper visualization cutoff (float)"

    # string attributes with corresponding properties
    program: str = "" #: Software that created this HDU
    filename: str = "" #: Original Filename
    row_order: str = "" #: Order of the rows in image array
    filter: str = "" #: Active filter name
    image_type: str = "" #: Type of image
    object: str = "" #: Name of the object of interest
    instrume: str = "" #: Instrument name
    telescop: str = "" #: Telescope used to acquire this image
    observer: str = "" #: Observer name
    bayer_pattern: str = "" #: Bayer color pattern
    sitelat_str: str = "" # [deg] Observation site latitude
    sitelong_str: str = "" # [deg] Observation site longitude
    focname: str = "" #: Focusing equipment name

    # Datetime attributes
    date: Optional[datetime] = None #: UTC date that FITS file was created
    date_obs: Optional[datetime] = None #: YYYY-MM-DDThh:mm:ss observation start, UT

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
    airmass: float = 1.0 #: Airmass at frame center (Gueymard 1993)
    focal_length: float = 0.0 #: [mm] Focal length
    flength: float = 0.0 #: [mm] Focal length
    iso_speed: float = 0.0 #: ISO speed value as a float
    exposure: float = 0.0 #: Exposure time as a float (s)
    aperture: float = 0.0 #: Aperture value as a float
    ccd_temp: float = 0.0 #: CCD temperature as a float
    set_temp: float = 0.0 #: CCD set temperature as a float
    livetime: float = 0.0 #: Sum of exposure times (s)
    stackcnt: int = 0 #: Number of stacked frames
    cvf: float = 0.0 #: Conversion factor (e- / ADU)
    gain: int = 0 #: Gain factor read in camera
    offset: int = 0 #: Offset value read in camera

    # Focuser data
    focuspos: int = 0 #: Focuser position
    focussz: int = 0 #: [um] Focuser step size
    foctemp: float = 0.0 #: Focuser temperature

    # Position data
    centalt: float = 0.0 #: [deg] Altitude of telescope
    centaz: float = 0.0 #: [deg] Azimuth of telescope
    sitelat: float = 0.0 # [deg] Observation site latitude
    sitelong: float = 0.0 #: [deg] Observation site longitude
    siteelev: float = 0.0 #: [m] Observation site elevation



@dataclass
class FFit:
    """
    Python equivalent of Siril ffit (FITS) structure, holding image
    pixel data and metadata.
    """
    bitpix: int = 0 #: Bits per pixel
    orig_bitpix: int = 0 #: Original bits per pixel (accounts for changes from original file).
    naxis: int = 0 #: The number of axes (2 for a mono image, 3 for a RGB image). Corresponds to the FITS kwyword NAXIS.
    _naxes: Tuple[int, int, int] = (0, 0, 0) #: A tuple holding the image dimensions.

    keywords: FKeywords = field(default_factory=FKeywords) #: A FKeywords object containing FITS header keywords.
    checksum: bool = False #: Whether Siril will write FITS data checksums for this file.

    header: Optional[str] = None #: The FITS header as a string.
    unknown_keys: Optional[str] = None #: All unknown FITS header keys as a string. This gives access to header cards that Siril does not use internally.

    stats: List[Optional[ImageStats]] = field(default_factory=lambda: [ImageStats() for _ in range(3)]) #: A list of ImageStats objects, one for each channel.
    mini: float = 0.0 #: The minimum value across all image channels.
    maxi: float = 0.0 #: The maximum value across all image channels.
    neg_ratio: np.float32 = 0.0 #: The ratio of negative pixels to the total pixels.

    type: DataType = DataType.FLOAT_IMG #: Specifies the image data type.
    _data: Optional[np.ndarray] = None #: Holds the image data as a numpy array.

    top_down: bool = False #: Specifies the ROWORDER for this image. The FITS specification directs that FITS should be stored bottom-up, but many CMOS sensors are natively TOP_DOWN and capture software tends to save FITS images captured by these sensors as TOP_DOWN.
    _focalkey: bool = False 
    _pixelkey: bool = False

    history: list[str] = field(default_factory=list) #: Contains a list of strings holding the HISTORY entries for this image.

    color_managed: bool = False #: Specifies whether the image is color managed or not.
    _icc_profile: Optional[bytes] = None #: Holds the ICC profile for the image as Bytes. This can be used by some modules including pillow.

    def __post_init__(self):
        """Initialize after creation"""
        if self.history is None:
            self.history = []

    @property
    def data(self) -> Optional[np.ndarray]:
        """
        The pixel data of the current image loaded in Siril, stored as a NumPy array
        """
        return self._data

    @data.setter
    def data(self, value: Optional[np.ndarray]):
        """
        Set the pixel data of the FFit to the provided NumPy array. Note: this
        does not update the image loaded in Siril - ``SirilInterface.set_image_pixeldata()``
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
                self._naxes = (shape[2], shape[1], 1)  # width, height, 1 channel
            elif len(shape) == 3:
                if shape[2] not in (1, 3):
                    raise ValueError(_("Third dimension must be 1 or 3"))
                self._naxes = (shape[2], shape[1], shape[0])  # width, height, channels
            else:
                raise ValueError(_("Data must be 2D or 3D"))
            self._update_naxis()
        self._data = value

    def _update_naxis(self):
        """Update naxis based on naxes value"""
        self.naxis = 3 if self._naxes[2] == 3 else 2

    @property
    def naxes(self) -> Tuple[int, int, int]:
        """
        The naxes tuple holds the image dimensions as width x height x channels. Note that the axis ordering differs between Siril representation as held in naxes and numpy representation as held in _data.shape (which is channels x height x width)
        """
        return self._naxes

    @property
    def width(self) -> int:
        """Image width"""
        return self._naxes[0]

    @property
    def height(self) -> int:
        """Image height"""
        return self._naxes[1]

    @property
    def channels(self) -> int:
        """Image channels"""
        return self._naxes[2]

    @property
    def icc_profile(self) -> Optional[bytes]:
        """
        The ICC profile as raw bytes data. This may be converted
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
        """The NumPy dtype based on the current type"""
        return np.uint16 if self.type == DataType.USHORT_IMG else np.float32

    def allocate_data(self):
        """
        Allocate memory for image data with appropriate type. self.width, self.height,
        self.naxis, self.naxes and self.dtype must be set before calling this
        method.
        """
        shape = (self.height, self.width) if self.naxis == 2 else (self.channels, self.height, self.width)
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
        previously have been obtained using get_image_pixeldata() or get_image()
        """
        if self.data is None:
            raise ValueError(_("No data allocated"))
        if self.naxis == 2:
            if channel != 0:
                raise ValueError(_("Cannot get channel > 0 for 2D data"))
            return self.data
        return self.data[channel, ...]

    def _estimate_noise(self, array: np.ndarray, nullcheck: Optional[bool] = True, nullvalue: Optional[float] = 0.0) -> float:
        """
        Estimate the background noise in the input image using the sigma of first-order differences.

        noise = 1.0 / sqrt(2) * RMS of (flux[i] - flux[i-1])

        Parameters:
            array (np.ndarray): 2D array of image pixels (np.uint16 or np.float32).
            nullcheck (bool): If True, check for null values.
            nullvalue: The value of null pixels (only used if nullcheck is True).

        Returns:
            float: Estimated noise value.
        """
        farray = array.astype(np.float32)
        if farray.ndim != 2:
            raise ValueError("Input array must be a 2D array.")

        ny, nx = farray.shape
        if nx < 3:
            return 0.0  # Not enough pixels in a row to compute differences

        diffs = []
        scale_factor = 1.0 / np.sqrt(2.0)

        for row in farray:
            # Handle null values if required
            if nullcheck:
                row = row[row != nullvalue]

            # Skip row if it has less than 2 valid pixels
            if len(row) < 2:
                continue

            # Compute first-order differences
            differences = np.diff(row)

            # Compute standard deviation with iterative sigma clipping
            for _ in range(3):  # NITER = 3
                mean = np.mean(differences)
                stdev = np.std(differences)
                if stdev == 0:
                    break
                differences = differences[np.abs(differences - mean) < 3 * stdev]  # SIGMA_CLIP = 3

            if len(differences) > 0:
                diffs.append(np.std(differences))

        if not diffs:
            return 0.0

        # Compute median of standard deviations
        median_stdev = np.median(diffs)
        return scale_factor * median_stdev

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
            try:
                channel_data = self.get_channel(i)
                if not isinstance(channel_data, np.ndarray):
                    raise TypeError(f"Channel data must be a NumPy array, got {type(channel_data)}")

                stats = self.stats[i] or ImageStats()

                # Handle case where channel_data is all NaN
                if np.all(np.isnan(channel_data)):
                    raise ValueError(f"Channel {i} contains only NaN values")

                # Update basic statistics - these should work even with all zeros
                stats.total = channel_data.size
                stats.ngoodpix = np.count_nonzero(~np.isnan(channel_data))

                # Remove both zeros and NaN values for remaining calculations
                nonzero = channel_data[~np.isnan(channel_data) & (channel_data != 0)]

                # Check if we have any valid non-zero pixels
                if nonzero.size > 0:
                    # Check for infinite values
                    if np.any(np.isinf(nonzero)):
                        raise ValueError(f"Channel {i} contains infinite values")

                    try:
                        stats.mean = np.mean(nonzero)
                        stats.median = np.median(nonzero)
                        stats.sigma = np.std(nonzero)
                        stats.min = np.min(nonzero)
                        stats.max = np.max(nonzero)

                        # Verify results are finite
                        if not all(np.isfinite(x) for x in [stats.mean, stats.median, stats.sigma, stats.min, stats.max]):
                            raise ValueError(f"Non-finite statistics computed for channel {i}")

                        stats.bgnoise = self._estimate_noise(channel_data)

                        # More complex statistics
                        deviations = np.abs(nonzero - stats.median)
                        stats.mad = np.median(deviations)
                        stats.avgDev = np.mean(deviations)

                        # Verify complex stats are finite
                        if not all(np.isfinite(x) for x in [stats.bgnoise, stats.mad, stats.avgDev]):
                            raise ValueError(f"Non-finite advanced statistics computed for channel {i}")

                    except (RuntimeWarning, RuntimeError) as e:
                        # Handle any numerical computation errors
                        raise ValueError(f"Error computing statistics for channel {i}: {str(e)}")

                else:
                    # Set all statistics to zero when there are no valid non-zero pixels
                    stats.mean = 0
                    stats.median = 0
                    stats.sigma = 0
                    stats.min = 0
                    stats.max = 0
                    stats.bgnoise = 0
                    stats.mad = 0
                    stats.avgDev = 0

                self.stats[i] = stats

            except Exception as e:
                # Log the error and set all statistics to zero for this channel
                logging.error(f"Error processing channel {i}: {str(e)}")
                stats = ImageStats()
                stats.total = channel_data.size if 'channel_data' in locals() else 0
                stats.ngoodpix = 0
                stats.mean = 0
                stats.median = 0
                stats.sigma = 0
                stats.min = 0
                stats.max = 0
                stats.bgnoise = 0
                stats.mad = 0
                stats.avgDev = 0
                self.stats[i] = stats

@dataclass
class Homography:
    """
    Python equivalent of the Siril Homography structure. Contains coefficients
    for the Homography matrix that maps a sequence frame onto the reference
    frame.
    """
    h00: float = 0.0 #: Homography matrix H00
    h01: float = 0.0 #: Homography matrix H01
    h02: float = 0.0 #: Homography matrix H02
    h10: float = 0.0 #: Homography matrix H10
    h11: float = 0.0 #: Homography matrix H11
    h12: float = 0.0 #: Homography matrix H12
    h20: float = 0.0 #: Homography matrix H20
    h21: float = 0.0 #: Homography matrix H21
    h22: float = 0.0 #: Homography matrix H22
    pair_matched: int = 0 #: number of pairs matched
    Inliers: int = 0 #: number of inliers kept after RANSAC step

@dataclass
class BGSample:
    """
    Python equivalent of the Siril background_sample struct. Used to hold
    background sample data obtained from Siril, or to generate or modify
    background sample data to set in Siril.
    A BGSample can be constructed as:
    - s1 = BGSample(x=1.0, y=2.0)
    - s2 = BGSample(position=(1.0, 2.0))
    - s3 = BGSample(x=1.0, y=2.0, mean=0.5, size=31)
    """
    median: Tuple[float, float, float] = (0.0, 0.0, 0.0) #: Median values for R, G and B channels. For mono images only median[0] is used.
    mean: float = 0.0
    min: float = 0.0
    max: float = 0.0
    size: int = field(default=25, init=False)  #: The default size matches the size of Siril bg samples.
    valid: bool = True  #: Samples default to being valid
    position: Optional[Tuple[float, float]] = field(default=None, init=False)  #: Position in (x, y) image coordinates

    def __init__(self, *args, x=None, y=None, position=None, size=25, **kwargs):
        """
        Custom constructor to handle both (x, y) and position arguments while allowing other attributes.
        Ensures `size`, if specified, is an odd number.
        """
        if (x is not None or y is not None) and position is not None:
            raise ValueError("Cannot provide both position tuple and x,y coordinates")
        if (x is not None) ^ (y is not None):  # XOR check
            raise ValueError("Must provide both x and y coordinates")
        if position is None and x is None and y is None:
            raise ValueError("Must provide either both x and y coordinates or a position tuple")

        # Assign position
        self.position = position if position is not None else (float(x), float(y))

        # Validate and assign size
        if size % 2 == 0:
            raise ValueError("Size must be an odd number")
        if size < 0:
            raise ValueError("Size must be positive")
        self.size = size

        # Manually initialize other dataclass fields from kwargs
        for field_name in self.__dataclass_fields__:
            if field_name not in {"position", "size"}:  # Already set manually
                setattr(self, field_name, kwargs.get(field_name, getattr(self.__class__, field_name)))

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
    x0: float = 0.0            #: x coordinate of the peak
    y0: float = 0.0            #: y coordinate of the peak
    sx: float = 0.0            #: Size of the fitted function on the x axis in PSF coordinates
    sy: float = 0.0            #: Size of the fitted function on the y axis in PSF coordinates
    fwhmx: float = 0.0         #: FWHM in x axis in pixels
    fwhmy: float = 0.0         #: FWHM in y axis in pixels
    fwhmx_arcsec: float = 0.0  #: FWHM in x axis in arc seconds
    fwhmy_arcsec: float = 0.0  #: FWHM in y axis in arc seconds
    angle: float = 0.0         #: angle of the x and yaxes with respect to the image x and y axes
    rmse: float = 0.0          #: RMSE of the minimization
    sat: float = 0.0           #: Level above which pixels have satured
    R: int = 0                 #: Optimized box size to enclose sufficient pixels in the background
    has_saturated: bool = False #: Shows whether the star is saturated or not

    # Moffat parameters
    beta: float = 0.0          #: Moffat equation beta parameter
    profile: StarProfile = StarProfile.GAUSSIAN  # Whether profile is Gaussian or Moffat

    xpos: float = 0.0          #: x position of the star in the image
    ypos: float = 0.0          #: y position of the star in the image

    # photometry data
    mag: float = 0.0           #: (V) magnitude, approximate or accurate
    Bmag: float = 0.0          #: B magnitude
    s_mag: float = 999.99      #: error on the (V) magnitude
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
    units: Optional[str] = None #: Units
    ra: float = 0.0            #: Right Ascension
    dec: float = 0.0           #: Declination

@dataclass
class RegData:
    """Python equivalent of Siril regdata structure"""
    fwhm: float = 0.0                    #: copy of fwhm->fwhmx, used as quality indicator
    weighted_fwhm: np.float32 = 0.0      #: used to exclude spurious images
    roundness: np.float32 = 0.0          #: fwhm->fwhmy / fwhm->fwhmx, 0 when uninit, ]0, 1] when set
    quality: float = 0.0                 #: measure of image quality
    background_lvl: np.float32 = 0.0     #: background level
    number_of_stars: int = 0             #: number of stars detected in the image
    H: Homography = field(default_factory=Homography)   #: Stores a homography matrix describing the affine transform from this frame to the reference frame

    def __repr__(self):
        attrs = [f"    {k}={getattr(self, k)}" for k in self.__dataclass_fields__]
        return f"{self.__class__.__name__}(\n" + ",\n".join(attrs) + "\n)"

@dataclass
class ImgData:
    """Python equivalent of Siril imgdata structure"""
    filenum: int = 0              #: real file index in the sequence
    incl: bool = False           #: selected in the sequence
    date_obs: Optional[datetime] = None  #: date of the observation
    airmass: float = 0.0         #: airmass of the image
    rx: int = 0                 #: width
    ry: int = 0                 #: height

    def __repr__(self):
        attrs = [f"    {k}={getattr(self, k)}" for k in self.__dataclass_fields__]
        return f"{self.__class__.__name__}(\n" + ",\n".join(attrs) + "\n)"

class DistoType(IntEnum):
    """Python equivalent of the Siril disto_source enum"""
    DISTO_UNDEF = 0      #: No distortion
    DISTO_IMAGE = 1      #: Distortion from current image
    DISTO_FILE = 2       #: Distortion from given file
    DISTO_MASTER = 3     #: Distortion from master files
    DISTO_FILES = 4      #: Distortion stored in each file (true only from seq platesolve, even with no distortion, it will be checked upon reloading)
    DISTO_FILE_COMET = 5 #: special for cometary alignement, to be detected by apply reg

    def __str__(self):
        if self == DistoType.DISTO_UNDEF:
            return "No distortion"
        if self == DistoType.DISTO_IMAGE:
            return "Distortion from current image"
        if self == DistoType.DISTO_FILE:
            return "Distortion from given file"
        if self == DistoType.DISTO_MASTER:
            return "Distortion from master files"
        if self == DistoType.DISTO_FILES:
            return "Distortion stored in each file"
        if self == DistoType.DISTO_FILE_COMET:
            return "Cometary alignement"
        return "Unknown distortion type"


@dataclass
class DistoData:
    """Python equivalent of Siril disto_params structure"""
    index: DistoType = DistoType.DISTO_UNDEF #: Specifies the distrosion type
    filename: str = ""                     #: filename if DISTO_FILE or DISTO_MASTER (and optional for DISTO_FILE_COMET)
    velocity: Tuple[float, float] = (0, 0) #: shift velocity if DISTO_FILE_COMET

    def __str__(self):
        """For pretty-printing distortion information"""
        pretty = f'{DistoType(self.index)}'
        if len(self.filename) > 0:
            pretty += f'\nDistorsion file: {self.filename}'
        if self.index == DistoType.DISTO_FILE_COMET:
            pretty += f'\nVelocity X/Y: {self.velocity[0]:.2f} {self.velocity[1]:.2f}'
        return pretty

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
    imgparam: List[ImgData] = None       #: a structure for each image of the sequence [number]
    regparam: List[List[RegData]] = None #: registration parameters for each layer [nb_layers][number]
    stats: List[List[ImageStats]] = None #: statistics of the images for each layer [nb_layers][number]
    distoparam: List[DistoData] = None   #: distortion data for the sequence [nb_layers]
    beg: int = 0                         #: imgparam[0]->filenum
    end: int = 0                         #: imgparam[number-1]->filenum
    exposure: float = 0.0                #: exposure of frames
    fz: bool = False                     #: whether the file is compressed
    type: SequenceType = None            #: the type of sequence
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
        if self.distoparam is None:
            self.distoparam = []
    
    def __str__(self):
        """For pretty-printing sequence information"""
        pretty = f'Sequence: {self.seqname}'
        pretty += f'\nImages [selected/total]: {self.selnum} / {self.number}'
        pretty += f'\nNumber of layers: {self.nb_layers}'
        pretty += f'\nBitdepth: {self.bitpix}'
        pretty += f'\nReference image: {self.reference_image + 1}'
        if not self.is_variable:
            pretty += f'\nImage size: {self.rx}x{self.ry}'
        else:
            pretty += '\nImages have variable sizes'
        for i, r in enumerate(self.regparam):
            if any(rr is not None for rr in r):
                pretty += '\nSequence has registration data'
                pretty += f' from layer {i}' if self.nb_layers > 1 else ''
                if self.distoparam is not None and self.distoparam[i] is not None and self.distoparam[i].index != DistoType.DISTO_UNDEF:
                    pretty += f'\nDistortion found in this layer: {self.distoparam[i]}'
        return pretty

@dataclass
class SirilPoint:
    """
    Represents a 2D point in the Siril image with x and y coordinates.
    """
    x: float
    y: float

MAX_POINTS_PER_POLYGON = 100

@dataclass
class UserPolygon:
    """
    Represents a user-defined polygon.

    Attributes:
        polygon_id (int): A unique identifier for the polygon.
        points (List[Point]): List of points defining the polygon's shape.
        color (int): Packed RGBA color (32-bit integer).
        fill (bool): If True, the polygon should be filled when drawn.
        legend (str): Optional legend for the polygon.
    """
    points: List[SirilPoint] #: List of points defining the polygon's shape
    polygon_id: int = 0 #: unique identifier
    color: int = 0xFFFFFFFF #: 32-bit RGBA color (packed, uint_8 per component)
    fill: bool = False #: whether or not the polygon should be filled when drawn
    legend: str = None #: an optional legend

    def __str__(self):
        """For pretty-printing polygon information"""
        pretty = f'User-defined overlay polygon: {self.legend}'
        pretty += f'\nID: {self.polygon_id}'
        pretty += f'\nColor (RGBA): 0x{self.color:08X}'
        pretty += f'\nFill: {self.fill}'
        for i, point in enumerate(self.points):
            pretty += f'\nPoint {i}: {point.x}, {point.y}'
        return pretty

    def serialize(self) -> bytes:
        """
        Serializes a single UserPolygon object into a byte array.
        Returns:
            bytes: A byte array representing the serialized polygon data.
        Raises:
            ValueError: If the number of points exceeds the allowed limit.
        """
        if len(self.points) > MAX_POINTS_PER_POLYGON:
            raise ValueError(f"Too many points in polygon {self.polygon_id}: max allowed is {MAX_POINTS_PER_POLYGON}")

        # Pack ID, number of points, color, and fill flag
        # Use 'I' for unsigned int and '?' for boolean
        buffer = bytearray()
        buffer.extend(struct.pack('!iiI?', self.polygon_id, len(self.points), self.color, self.fill))

        # Pack each point as float
        for point in self.points:
            buffer.extend(struct.pack('!dd', point.x, point.y))

        # Pack the legend (if it exists)
        if self.legend is not None:
            legend_bytes = self.legend.encode('utf-8')
            # Pack the length of the string first, then the string itself
            buffer.extend(struct.pack('!i', len(legend_bytes)))
            buffer.extend(legend_bytes)
        else:
            # If legend is None, pack a length of 0
            buffer.extend(struct.pack('!i', 0))

        return bytes(buffer)

    @classmethod
    def deserialize_polygon(cls, data: bytes) -> Tuple['UserPolygon', bytes]:
        if len(data) < 13:
            raise ValueError("Invalid data size for polygon")

        polygon_id, n_points, color, fill = struct.unpack('!iiI?', data[:13])
        data = data[13:]

        if n_points < 0 or n_points > MAX_POINTS_PER_POLYGON:
            raise ValueError(f"Invalid number of points: {n_points}")

        points = []
        for _ in range(n_points):
            if len(data) < 16:
                raise ValueError("Not enough data for points")

            x, y = struct.unpack('!dd', data[:16])
            data = data[16:]
            points.append(SirilPoint(x, y))

        # Read legend length
        if len(data) < 4:
            raise ValueError("Not enough data for legend length")

        legend_length = struct.unpack('!i', data[:4])[0]
        data = data[4:]

        if legend_length > 0:
            if len(data) < legend_length:
                raise ValueError("Not enough data for legend string")

            legend = data[:legend_length].decode('utf-8')
            data = data[legend_length:]
        else:
            legend = None

        polygon = cls(points, polygon_id, color, fill, legend)
        return polygon, data

    @classmethod
    def deserialize_polygon_list(cls, data: bytes) -> List['UserPolygon']:
        if len(data) < 4:
            raise ValueError("Invalid data size for polygon list")

        num_polygons = struct.unpack('!I', data[:4])[0]
        data = data[4:]

        polygons = []
        for _ in range(num_polygons):
            # Create each polygon from the remaining data
            polygon, data = cls.deserialize_polygon(data)
            polygons.append(polygon)

        return polygons
