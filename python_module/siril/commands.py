from .models import FKeywords, ImageStats, DataType, FFit

def get_ffit() -> FFit:
    """
    Get complete information about the currently loaded image as an FFit structure.

    Returns:
        FFit: A complete FFit object containing all image information

    Raises:
        CommandError: If the command fails or returns invalid data
    """
    conn = _ensure_connection()
    cmd, data = conn._send_command(_Command.GET_FFIT)

    if cmd == _Command.ERROR or not data:
        raise CommandError("Failed to get FFit data")

    # Create a BytesIO object to read the packed data
    buffer = io.BytesIO(data)

    # Initialize FFit object
    ffit = FFit()

    # Unpack basic image information
    (ffit.bitpix, ffit.orig_bitpix, ffit.naxis) = struct.unpack('!iii', buffer.read(12))

    # Unpack dimensions (naxes)
    naxes = struct.unpack('!iii', buffer.read(12))
    ffit._naxes = naxes

    # Unpack simple numeric fields from keywords
    (
        ffit.keywords.bscale, ffit.keywords.bzero,
        ffit.keywords.lo, ffit.keywords.hi,
        ffit.keywords.flo, ffit.keywords.fhi,
        ffit.keywords.data_max, ffit.keywords.data_min,
        ffit.keywords.pixel_size_x, ffit.keywords.pixel_size_y,
    ) = struct.unpack('!ffiiffffff', buffer.read(40))

    # Unpack integer fields
    (
        ffit.keywords.binning_x, ffit.keywords.binning_y,
        ffit.keywords.bayer_xoffset, ffit.keywords.bayer_yoffset,
        ffit.keywords.focuspos, ffit.keywords.focussz,
        ffit.keywords.stackcnt, ffit.keywords.key_gain,
        ffit.keywords.key_offset
    ) = struct.unpack('!iiiiiiiii', buffer.read(36))

    # Unpack float fields
    (
        ffit.keywords.expstart, ffit.keywords.expend,
        ffit.keywords.centalt, ffit.keywords.centaz,
        ffit.keywords.sitelat, ffit.keywords.sitelong,
        ffit.keywords.siteelev, ffit.keywords.airmass,
        ffit.keywords.focal_length, ffit.keywords.flength,
        ffit.keywords.iso_speed, ffit.keywords.exposure,
        ffit.keywords.aperture, ffit.keywords.ccd_temp,
        ffit.keywords.set_temp, ffit.keywords.livetime,
        ffit.keywords.cvf, ffit.keywords.foctemp
    ) = struct.unpack('!ffffffffffffffffff', buffer.read(72))

    # Function to read null-terminated string
    def read_string():
        chars = []
        while True:
            c = buffer.read(1)
            if c == b'\0' or not c:
                break
            chars.append(c.decode('utf-8'))
        return ''.join(chars)

    # Read string fields
    ffit.keywords.program = read_string()
    ffit.keywords.filename = read_string()
    ffit.keywords.row_order = read_string()
    ffit.keywords.filter = read_string()
    ffit.keywords.image_type = read_string()
    ffit.keywords.object = read_string()
    ffit.keywords.instrume = read_string()
    ffit.keywords.telescop = read_string()
    ffit.keywords.observer = read_string()
    ffit.keywords.sitelat_str = read_string()
    ffit.keywords.sitelong_str = read_string()
    ffit.keywords.bayer_pattern = read_string()
    ffit.keywords.focname = read_string()

    # Read date fields - assuming they're packed as timestamps
    timestamp = struct.unpack('!Q', buffer.read(8))[0]
    ffit.keywords.date = datetime.fromtimestamp(timestamp) if timestamp else None

    timestamp = struct.unpack('!Q', buffer.read(8))[0]
    ffit.keywords.date_obs = datetime.fromtimestamp(timestamp) if timestamp else None

    # Read boolean flags
    flags = struct.unpack('!I', buffer.read(4))[0]
    ffit.checksum = bool(flags & 1)
    ffit.top_down = bool(flags & 2)
    ffit.focalkey = bool(flags & 4)
    ffit.pixelkey = bool(flags & 8)
    ffit.color_managed = bool(flags & 16)

    # Read header and unknown keys if present
    header_size = struct.unpack('!I', buffer.read(4))[0]
    ffit.header = buffer.read(header_size).decode('utf-8') if header_size > 0 else None

    unknown_size = struct.unpack('!I', buffer.read(4))[0]
    ffit.unknown_keys = buffer.read(unknown_size).decode('utf-8') if unknown_size > 0 else None

    # Read statistics for each channel
    for i in range(3):
        stats = ImageStats()
        (
            stats.total, stats.ngoodpix,
            stats.mean, stats.median, stats.sigma,
            stats.avgDev, stats.mad, stats.sqrtbwmv,
            stats.location, stats.scale,
            stats.min, stats.max,
            stats.normValue, stats.bgnoise
        ) = struct.unpack('!iifffffffffffff', buffer.read(60))
        ffit.stats[i] = stats

    # Read remaining float values
    ffit.mini, ffit.maxi, ffit.neg_ratio = struct.unpack('!fff', buffer.read(12))

    # Read data type
    dtype = struct.unpack('!I', buffer.read(4))[0]
    ffit.type = DataType(dtype)

    # Read image data if present
    data_size = struct.unpack('!Q', buffer.read(8))[0]
    if data_size > 0:
        raw_data = buffer.read(data_size)
        # Convert raw data to numpy array based on data type
        dtype_map = {
            DataType.BYTE_IMG: np.uint8,
            DataType.SHORT_IMG: np.int16,
            DataType.USHORT_IMG: np.uint16,
            DataType.LONG_IMG: np.int32,
            DataType.FLOAT_IMG: np.float32,
            DataType.DOUBLE_IMG: np.float64
        }
        np_dtype = dtype_map.get(ffit.type, np.float32)
        ffit._data = np.frombuffer(raw_data, dtype=np_dtype).reshape(naxes)

    # Read ICC profile if present
    icc_size = struct.unpack('!I', buffer.read(4))[0]
    ffit.icc_profile = buffer.read(icc_size) if icc_size > 0 else None

    # Read history
    history_count = struct.unpack('!I', buffer.read(4))[0]
    ffit.history = []
    for _ in range(history_count):
        ffit.history.append(read_string())

    return ffit
