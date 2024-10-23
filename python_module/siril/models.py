from dataclasses import dataclass, field
from datetime import datetime
from typing import List, Optional, Tuple, Union
import numpy as np
from enum import IntEnum

class DataType(IntEnum):
    """Mimics the Siril data_type enum"""
    BYTE_IMG = 8
    SHORT_IMG = 16
    USHORT_IMG = 16
    LONG_IMG = 32
    FLOAT_IMG = 32
    DOUBLE_IMG = 64

@dataclass
class ImageStats:
    """Python equivalent of Siril imstats structure"""
    total: int = 0           # number of pixels
    ngoodpix: int = 0        # number of non-zero pixels
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

    # Getters
    def get_total(self) -> int:
        return self.total

    def get_ngoodpix(self) -> int:
        return self.ngoodpix

    def get_mean(self) -> float:
        return self.mean

    def get_median(self) -> float:
        return self.median

    def get_sigma(self) -> float:
        return self.sigma

    def get_avgDev(self) -> float:
        return self.avgDev

    def get_mad(self) -> float:
        return self.mad

    def get_sqrtbwmv(self) -> float:
        return self.sqrtbwmv

    def get_location(self) -> float:
        return self.location

    def get_scale(self) -> float:
        return self.scale

    def get_min(self) -> float:
        return self.min

    def get_max(self) -> float:
        return self.max

    def get_normValue(self) -> float:
        return self.normValue

    def get_bgnoise(self) -> float:
        return self.bgnoise

@dataclass
class FKeywords:
    """Python equivalent of fkeywords structure"""
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

    # DateTime fields
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

    # Getters
    def get_bscale(self) -> float:
        return self.bscale

    def get_bzero(self) -> float:
        return self.bzero

    def get_lo(self) -> int:
        return self.lo

    def get_hi(self) -> int:
        return self.hi

    def get_flo(self) -> float:
        return self.flo

    def get_fhi(self) -> float:
        return self.fhi

    def get_program(self) -> str:
        return self.program

    def get_filename(self) -> str:
        return self.filename

    def get_data_max(self) -> float:
        return self.data_max

    def get_data_min(self) -> float:
        return self.data_min

    def get_pixel_size_x(self) -> float:
        return self.pixel_size_x

    def get_pixel_size_y(self) -> float:
        return self.pixel_size_y

    def get_binning_x(self) -> int:
        return self.binning_x

    def get_binning_y(self) -> int:
        return self.binning_y

    def get_row_order(self) -> str:
        return self.row_order

    def get_date(self) -> Optional[datetime]:
        return self.date

    def get_date_obs(self) -> Optional[datetime]:
        return self.date_obs

    def get_expstart(self) -> float:
        return self.expstart

    def get_expend(self) -> float:
        return self.expend

    def get_filter(self) -> str:
        return self.filter

    def get_image_type(self) -> str:
        return self.image_type

    def get_object(self) -> str:
        return self.object

    def get_instrume(self) -> str:
        return self.instrume

    def get_telescop(self) -> str:
        return self.telescop

    def get_observer(self) -> str:
        return self.observer

    def get_centalt(self) -> float:
        return self.centalt

    def get_centaz(self) -> float:
        return self.centaz

    def get_sitelat(self) -> float:
        return self.sitelat

    def get_sitelong(self) -> float:
        return self.sitelong

    def get_sitelat_str(self) -> str:
        return self.sitelat_str

    def get_sitelong_str(self) -> str:
        return self.sitelong_str

    def get_siteelev(self) -> float:
        return self.siteelev

    def get_bayer_pattern(self) -> str:
        return self.bayer_pattern

    def get_bayer_xoffset(self) -> int:
        return self.bayer_xoffset

    def get_bayer_yoffset(self) -> int:
        return self.bayer_yoffset

    def get_airmass(self) -> float:
        return self.airmass

    def get_focal_length(self) -> float:
        return self.focal_length

    def get_flength(self) -> float:
        return self.flength

    def get_iso_speed(self) -> float:
        return self.iso_speed

    def get_exposure(self) -> float:
        return self.exposure

    def get_aperture(self) -> float:
        return self.aperture

    def get_ccd_temp(self) -> float:
        return self.ccd_temp

    def get_set_temp(self) -> float:
        return self.set_temp

    def get_livetime(self) -> float:
        return self.livetime

    def get_stackcnt(self) -> int:
        return self.stackcnt

    def get_cvf(self) -> float:
        return self.cvf

    def get_key_gain(self) -> int:
        return self.key_gain

    def get_key_offset(self) -> int:
        return self.key_offset

    def get_focname(self) -> str:
        return self.focname

    def get_focuspos(self) -> int:
        return self.focuspos

    def get_focussz(self) -> int:
        return self.focussz

    def get_foctemp(self) -> float:
        return self.foctemp

    # Setters
    def set_bscale(self, value: float) -> None:
        self.bscale = value

    def set_bzero(self, value: float) -> None:
        self.bzero = value

    def set_lo(self, value: int) -> None:
        self.lo = value

    def set_hi(self, value: int) -> None:
        self.hi = value

    def set_flo(self, value: float) -> None:
        self.flo = value

    def set_fhi(self, value: float) -> None:
        self.fhi = value

    def set_program(self, value: str) -> None:
        self.program = value

    def set_filename(self, value: str) -> None:
        self.filename = value

    def set_data_max(self, value: float) -> None:
        self.data_max = value

    def set_data_min(self, value: float) -> None:
        self.data_min = value

    def set_pixel_size_x(self, value: float) -> None:
        if value > 0:
            self.pixel_size_x = value
        else:
            raise ValueError("pixel_size_x must be greater than 0")

    def set_pixel_size_y(self, value: float) -> None:
        if value > 0:
            self.pixel_size_y = value
        else:
            raise ValueError("pixel_size_y must be greater than 0")

    def set_binning_x(self, value: int) -> None:
        if value >= 1:
            self.binning_x = value
        else:
            raise ValueError("binning_x must be greater than or equal to 1")

    def set_binning_y(self, value: int) -> None:
        if value >= 1:
            self.binning_y = value
        else:
            raise ValueError("binning_y must be greater than or equal to 1")

    def set_row_order(self, value: str) -> None:
        self.row_order = value

    def set_date(self, value: Optional[datetime]) -> None:
        self.date = value

    def set_date_obs(self, value: Optional[datetime]) -> None:
        self.date_obs = value

    def set_expstart(self, value: float) -> None:
        self.expstart = value

    def set_expend(self, value: float) -> None:
        self.expend = value

    def set_filter(self, value: str) -> None:
        self.filter = value

    def set_image_type(self, value: str) -> None:
        self.image_type = value

    def set_object(self, value: str) -> None:
        self.object = value

    def set_instrume(self, value: str) -> None:
        self.instrume = value

    def set_telescop(self, value: str) -> None:
        self.telescop = value

    def set_observer(self, value: str) -> None:
        self.observer = value

    def set_centalt(self, value: float) -> None:
        if value <= 90:
            self.centalt = value
        else:
            raise ValueError("centalt must be less than or equal to 90 degrees")

    def set_centaz(self, value: float) -> None:
        if 0 <= value < 360:
            self.centaz = value
        else:
            raise ValueError("centaz must be between 0 and 360 degrees (exclusive)")

    def set_sitelat(self, value: float) -> None:
        if -90 <= value <= 90:
            self.sitelat = value
        else:
            raise ValueError("sitelat must be between -90 and 90 degrees")

    def set_sitelong(self, value: float) -> None:
        if 0 <= value < 360:
            self.sitelong = value
        else:
            raise ValueError("sitelong must be between 0 and 360 degrees (exclusive)")

    def set_sitelat_str(self, value: str) -> None:
        self.sitelat_str = value

    def set_sitelong_str(self, value: str) -> None:
        self.sitelong_str = value

    def set_siteelev(self, value: float) -> None:
        self.siteelev = value

    def set_bayer_pattern(self, value: str) -> None:
        self.bayer_pattern = value

    def set_bayer_xoffset(self, value: int) -> None:
        self.bayer_xoffset = value

    def set_bayer_yoffset(self, value: int) -> None:
        self.bayer_yoffset = value

    def set_airmass(self, value: float) -> None:
        if value >= 1:
            self.airmass = value
        else:
            raise ValueError("airmass must be greater than or equal to 1")

    def set_focal_length(self, value: float) -> None:
        self.focal_length = value

    def set_flength(self, value: float) -> None:
        self.flength = value

    def set_iso_speed(self, value: float) -> None:
        self.iso_speed = value

    def set_exposure(self, value: float) -> None:
        self.exposure = value

    def set_aperture(self, value: float) -> None:
        self.aperture = value

    def set_ccd_temp(self, value: float) -> None:
        self.ccd_temp = value

    def set_set_temp(self, value: float) -> None:
        self.set_temp = value

    def set_livetime(self, value: float) -> None:
        self.livetime = value

    def set_stackcnt(self, value: int) -> None:
        self.stackcnt = value

    def set_cvf(self, value: float) -> None:
        self.cvf = value

    def set_key_gain(self, value: int) -> None:
        self.key_gain = value

    def set_key_offset(self, value: int) -> None:
        self.key_offset = value

    def set_focname(self, value: str) -> None:
        self.focname = value

    def set_focuspos(self, value: int) -> None:
        self.focuspos = value

    def set_focussz(self, value: int) -> None:
        self.focussz = value

    def set_foctemp(self, value: float) -> None:
        self.foctemp = value

@dataclass
class FFit:
    """Python equivalent of ffit structure with automated dimension handling"""
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

    history: List[str] = field(default_factory=list)

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
        """Get the data array"""
        return self._data

    @data.setter
    def data(self, value: Optional[np.ndarray]):
        """Set the data array and update dimensions"""
        if value is not None:
            shape = value.shape
            if len(shape) == 2:
                self._naxes = (shape[1], shape[0], 1)  # width, height, channels
            elif len(shape) == 3:
                if shape[2] not in (1, 3):
                    raise ValueError("Third dimension must be 1 or 3")
                self._naxes = (shape[1], shape[0], shape[2])  # width, height, channels
            else:
                raise ValueError("Data must be 2D or 3D")
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
    def rx(self) -> int:
        """Get image width"""
        return self._naxes[0]

    @property
    def ry(self) -> int:
        """Get image height"""
        return self._naxes[1]

    @property
    def icc_profile(self) -> Optional[bytes]:
        """Get the ICC profile"""
        return self._icc_profile

    @icc_profile.setter
    def icc_profile(self, value: Optional[bytes]):
        """Set ICC profile and update color_managed flag"""
        self._icc_profile = value
        self.color_managed = value is not None

    @property
    def dtype(self) -> np.dtype:
        """Get the numpy dtype based on the current type"""
        return np.uint16 if self.type == DataType.USHORT_IMG else np.float32

    def allocate_data(self):
        """Allocate memory for image data with appropriate type"""
        shape = (self.ry, self.rx) if self.naxis == 2 else (self.ry, self.rx, self.naxes[2])
        self.data = np.zeros(shape, dtype=self.dtype)

    def ensure_data_type(self):
        """Ensure data is in the correct type with proper scaling"""
        if self.data is not None and self.data.dtype != self.dtype:
            if self.type == DataType.FLOAT_IMG and self.data.dtype == np.uint16:
                # Convert from USHORT (0-65535) to FLOAT (0.0-1.0)
                self._data = self.data.astype(np.float32) / 65535.0
            elif self.type == DataType.USHORT_IMG and self.data.dtype == np.float32:
                # Convert from FLOAT (0.0-1.0) to USHORT (0-65535)
                self._data = (self.data * 65535.0).clip(0, 65535).astype(np.uint16)
            else:
                raise ValueError(f"Unsupported type conversion from {self.data.dtype} to {self.dtype}")

    def get_channel(self, channel: int) -> np.ndarray:
        """Get a specific channel of the data"""
        if self.data is None:
            raise ValueError("No data allocated")
        if self.naxis == 2:
            if channel != 0:
                raise ValueError("Cannot get channel > 0 for 2D data")
            return self.data
        return self.data[:, :, channel]

    def set_channel(self, channel: int, data: np.ndarray):
        """Set a specific channel of the data"""
        if self.data is None:
            raise ValueError("No data allocated")
        if self.naxis == 2:
            if channel != 0:
                raise ValueError("Cannot set channel > 0 for 2D data")
            self.data[:] = data
        else:
            self.data[:, :, channel] = data

    def update_stats(self):
        """Update image statistics for all channels"""
        if self.data is None:
            raise ValueError("No data allocated")

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
