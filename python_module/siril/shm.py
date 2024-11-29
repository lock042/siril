# Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
# Copyright (C) 2012-2024 team free-astro (see more in AUTHORS file)
# Reference site is https://siril.org
# SPDX-License-Identifier: GPL-3.0-or-later

from multiprocessing.shared_memory import SharedMemory
import os
import sys
import mmap
from .translations import _

class SharedMemoryWrapper:
    """
    Wrapper class to handle shared memory creation and cleanup across platforms.
    This is an internal class and is not intended for use by script authors
    """
    def __init__(self, name: str, size: int):
        self.name = name
        self.size = size
        self._shm = None
        self._fd = None
        self._mmap = None

        if sys.platform == 'win32':
            self._create_windows()
        else:
            self._create_unix()

    def _create_windows(self):
        """
        Create shared memory (Windows version). Not intended as a function for use
        directly in scripts: this is an internal method.
        """

        try:
            # First try to create
            self._shm = SharedMemory(name=self.name, create=True, size=self.size)
        except FileExistsError:
            # If it exists, try to open the existing segment instead
            try:
                self._shm = SharedMemory(name=self.name, create=False)
                # Verify size matches expected size
                if self._shm.buf.nbytes < self.size:
                    # Check that the actual buffer size is at least as big as expected
                    # (it may be bigger owing to rounding up to the next page boundary)
                    raise ValueError(f"Existing shared memory size {self._shm.buf.nbytes} does not match expected {self.size}")
            except Exception as e:
                raise RuntimeError(f"Failed to open existing shared memory: {e}")

    def _create_unix(self):
        """
        Create shared memory (POSIX version). Not intended as a function for use
        directly in scripts: this is an internal method.
        """
        try:
            # Create a temporary file as fallback
            fd = os.open(f"/dev/shm/{self.name}", os.O_CREAT | os.O_RDWR, 0o600)
            os.ftruncate(fd, self.size)
            self._fd = fd
            self._mmap = mmap.mmap(fd, self.size)
        except Exception as e:
            if self._fd is not None:
                try:
                    os.close(self._fd)
                except:
                    pass
            raise RuntimeError(f"Failed to create shared memory on Unix: {e}")

    @property
    def buf(self):
        if sys.platform == 'win32':
            return self._shm.buf
        return self._mmap

    def close(self):
        if sys.platform == 'win32':
            if self._shm is not None:
                try:
                    self._shm.close()
                except:
                    pass
        else:
            if self._mmap is not None:
                try:
                    self._mmap.close()
                except:
                    pass
            if self._fd is not None:
                try:
                    os.close(self._fd)
                except:
                    pass

    def unlink(self):
        if sys.platform == 'win32':
            if self._shm is not None:
                try:
                    self._shm.unlink()
                except:
                    pass
        else:
            try:
                os.unlink(f"/dev/shm/{self.name}")
            except:
                pass
