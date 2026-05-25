# Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
# Copyright (C) 2012-2026 team free-astro (see more in AUTHORS file)
# Reference site is https://siril.org
# SPDX-License-Identifier: GPL-3.0-or-later

import os
import struct
import ctypes
from multiprocessing import shared_memory
from multiprocessing.resource_tracker import unregister

class _SharedMemoryInfo(ctypes.Structure):
    """
    Structure matching the C-side shared memory info. Internal class:
    this is not intended for use in scripts.
    """
    _fields_ = [
        ("size", ctypes.c_size_t),
        ("data_type", ctypes.c_int),  # 0 for WORD, 1 for float
        ("width", ctypes.c_int),
        ("height", ctypes.c_int),
        ("channels", ctypes.c_int),
        ("shm_name", ctypes.c_char * 256)
    ]

    # Wire format of the shm-info struct emitted by Siril. Matches the
    # ctypes layout above field-for-field. Use this from connection.py
    # rather than copy-pasting `struct.unpack_from('=Qiiii256s', ...)`.
    _STRUCT_FORMAT = '=Qiiii256s'
    _STRUCT_SIZE = struct.calcsize('=Qiiii256s')

    def __len__(self):
        """Return the size in bytes of this structure."""
        return ctypes.sizeof(self)

    @classmethod
    def from_wire_buffer(cls, data: bytes, offset: int = 0) -> '_SharedMemoryInfo':
        """
        Decode the on-wire shm-info struct (`=Qiiii256s`, 280 bytes
        starting at `offset`) into a _SharedMemoryInfo instance.
        Retires the copy-pasted unpack pattern in connection.py.
        """
        size, data_type, width, height, channels, shm_name = struct.unpack_from(
            cls._STRUCT_FORMAT, data, offset
        )
        return cls(
            size=size,
            data_type=data_type,
            width=width,
            height=height,
            channels=channels,
            shm_name=shm_name,
        )

class SharedMemoryWrapper:
    """
    Wrapper class to handle shared memory creation and cleanup across platforms.
    """
    def __init__(self, name: str, size: int):
        if os.name != "nt":
            name = name.lstrip('/') # Remove leading '/' on POSIX systems
                                    # because SharedMemory.__init__ will add it back
        self.name = name
        self.size = size  # Store intended size separately
        self._shm = None
        try:
            # First try to attach to existing shared memory
            self._shm = shared_memory.SharedMemory(name=self.name)
            if os.name != "nt":
                # Unregister from the resource tracker as the shm is cleaned by Siril
                # Don't do this on Windows, it doesn't need it and doesn't work!
                unregister(self._shm._name, "shared_memory")
            if self._shm.size < self.size:
                raise ValueError(f"Existing shared memory size {self._shm.size} does not match expected {self.size}")
        except FileNotFoundError:
            # If it doesn't exist, create new shared memory
            print("Existing SHM not found, creating a new one...")
            self._shm = shared_memory.SharedMemory(name=self.name, create=True, size=self.size)

    @property
    def buf(self):
        # Return only the intended size of the buffer to handle MacOS page size rounding
        return memoryview(self._shm.buf)[:self.size]

    def close(self):
        if self._shm is not None:
            try:
                self._shm.close()
            except Exception:
                pass

    def unlink(self):
        if self._shm is not None:
            try:
                self._shm.unlink()
            except Exception:
                pass
