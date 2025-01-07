# Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
# Copyright (C) 2012-2024 team free-astro (see more in AUTHORS file)
# Reference site is https://siril.org
# SPDX-License-Identifier: GPL-3.0-or-later

from multiprocessing import shared_memory
from multiprocessing.resource_tracker import unregister

class SharedMemoryWrapper:
    """
    Wrapper class to handle shared memory creation and cleanup across platforms.
    """
    def __init__(self, name: str, size: int):
        self.name = name
        self.size = size  # Store intended size separately
        self._shm = None
        try:
            # First try to attach to existing shared memory
            self._shm = shared_memory.SharedMemory(name=self.name)
            unregister(self._shm._name, "shared_memory")
            if self._shm.size < self.size:
                raise ValueError(f"Existing shared memory size {self._shm.size} does not match expected {self.size}")
        except FileNotFoundError:
            # If it doesn't exist, create new shared memory
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
