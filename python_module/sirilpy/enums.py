# Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
# Copyright (C) 2012-2025 team free-astro (see more in AUTHORS file)
# Reference site is https://siril.org
# SPDX-License-Identifier: GPL-3.0-or-later

from enum import IntEnum, unique

@unique
class _Status(IntEnum):
    """
    Returns the status of a command. NONE is for commands that
    may legitimately fail to return data but which should not be
    regarded as an error, instead this triggers the command processor
    to return the special python value None
    Internal class: this is not intended for use in scripts.
    """

    OK = 0
    NONE = 1
    ERROR = 0xFF

@unique
class _Command(IntEnum):
    """
    Enumerates the commands. This enum MUST match the one in
    siril_pythonmodule.h. Internal class: this is not intended for
    use in scripts.
    """
    SEND_COMMAND = 1
    LOG_MESSAGE = 2
    UPDATE_PROGRESS = 3
    GET_WORKING_DIRECTORY = 4
    GET_FILENAME = 5
    GET_DIMENSIONS = 6
    GET_PIXELDATA = 7
    GET_PIXELDATA_REGION = 8
    RELEASE_SHM = 9
    SET_PIXELDATA = 10
    GET_IMAGE_STATS = 11
    GET_KEYWORDS = 12
    GET_ICC_PROFILE = 13
    GET_FITS_HEADER = 14
    GET_FITS_HISTORY = 15
    GET_FITS_UNKNOWN_KEYS = 16
    GET_IMAGE = 17
    GET_PSFSTARS = 18
    GET_SEQ_STATS = 19
    GET_SEQ_REGDATA = 20
    GET_SEQ_IMGDATA = 21
    GET_SEQ_PIXELDATA = 22
    GET_SEQ_IMAGE = 23
    GET_SEQ = 24
    GET_CONFIG = 25
    GET_USERCONFIGDIR = 26
    GET_IS_IMAGE_LOADED = 27
    GET_IS_SEQUENCE_LOADED = 28
    GET_SELECTION = 29
    SET_SELECTION = 30
    GET_ACTIVE_VPORT = 31
    GET_STAR_IN_SELECTION = 32
    GET_STATS_FOR_SELECTION = 33
    PIX2WCS = 34
    WCS2PIX = 35
    UNDO_SAVE_STATE = 36
    GET_BUNDLE_PATH = 37
    ERROR_MESSAGEBOX = 38
    ERROR_MESSAGEBOX_MODAL = 39
    SIRIL_PLOT = 40
    CLAIM_THREAD = 41
    RELEASE_THREAD = 42
    SET_SEQ_FRAME_PIXELDATA = 43
    REQUEST_SHM = 44
    SET_SEQ_FRAME_INCL = 45
    GET_USERDATADIR = 46
    GET_SYSTEMDATADIR = 47
    GET_BGSAMPLES = 48
    SET_BGSAMPLES = 49
    GET_SEQ_FRAME_FILENAME = 50
    INFO_MESSAGEBOX = 51
    INFO_MESSAGEBOX_MODAL = 52
    WARNING_MESSAGEBOX = 53
    WARNING_MESSAGEBOX_MODAL = 54
    GET_SEQ_DISTODATA = 55
    SET_IMAGE_HEADER = 56
    ADD_USER_POLYGON = 57
    DELETE_USER_POLYGON = 58
    CLEAR_USER_POLYGONS = 59
    GET_USER_POLYGON = 60
    GET_USER_POLYGON_LIST = 61
    CONFIRM_MESSAGEBOX = 62
    ERROR = 0xFF

@unique
class _CommandStatus(IntEnum):
    CMD_NOT_FOUND = 1
    CMD_NO_WAIT = 1 << 1
    CMD_NO_CWD = 1 << 2
    CMD_NOT_SCRIPTABLE = 1 << 3
    CMD_WRONG_N_ARG = 1 << 4
    CMD_ARG_ERROR = 1 << 5
    CMD_SELECTION_ERROR = 1 << 6
    CMD_OK = 0
    CMD_GENERIC_ERROR = 1 << 7
    CMD_IMAGE_NOT_FOUND = 1 << 8
    CMD_SEQUENCE_NOT_FOUND = 1 << 9
    CMD_INVALID_IMAGE = 1 << 10
    CMD_LOAD_IMAGE_FIRST = 1 << 11
    CMD_ONLY_SINGLE_IMAGE = 1 << 12
    CMD_NOT_FOR_SINGLE = 1 << 13
    CMD_NOT_FOR_MONO = 1 << 14
    CMD_NOT_FOR_RGB = 1 << 15
    CMD_FOR_CFA_IMAGE = 1 << 16
    CMD_FILE_NOT_FOUND = 1 << 17
    CMD_FOR_PLATE_SOLVED = 1 << 18
    CMD_NEED_INIT_FIRST = 1 << 19
    CMD_ALLOC_ERROR = 1 << 20
    CMD_THREAD_RUNNING = 1 << 21
    CMD_DIR_NOT_FOUND = 1 << 22

class _Defaults:
    """
    Contains default values for different datatypes, matching Siril
    """
    DEFAULT_DOUBLE_VALUE = -999.0
    DEFAULT_FLOAT_VALUE = -999.0
    DEFAULT_INT_VALUE = -2147483647
    DEFAULT_UINT_VALUE = 2147483647
    VALUES = {DEFAULT_DOUBLE_VALUE, DEFAULT_FLOAT_VALUE, DEFAULT_INT_VALUE, DEFAULT_UINT_VALUE}

@unique
class _ConfigType(IntEnum):
    """
    Enumerates config variable types for use with the
    ``get_siril_config()`` method. Internal class: this is not intended
    for use in scripts.
    """
    BOOL = 0
    INT = 1
    DOUBLE = 2
    STR = 3
    STRDIR = 4
    STRLIST = 5