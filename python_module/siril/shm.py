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
        self.size = size
        self._shm = None

        # Create or attach to the shared memory block
        try:
            # If the shm allocation doesn't exist, create one
            self._shm = shared_memory.SharedMemory(name=self.name, create=True, size=self.size)
            # If we create the shm (e.g. for set_image_pixeldata) we must clean it up, so no call
            # to unregister() here
        except FileExistsError:
            # Otherwise, open the existing one created by Siril
            self._shm = shared_memory.SharedMemory(name=self.name)
            # All shm allocations are unlinked by Siril after use, we do not need to resource track them
            unregister(self._shm._name, "shared_memory")
            if self._shm.size < self.size:
                raise ValueError(f"Existing shared memory size {self._shm.size} does not match expected {self.size}")

    @property
    def buf(self):
        return self._shm.buf

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
