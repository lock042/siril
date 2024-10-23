from multiprocessing.shared_memory import SharedMemory
import os
import sys
import struct
from typing import Optional
import tempfile
import mmap

class SharedMemoryWrapper:
    """Wrapper class to handle shared memory creation and cleanup across platforms."""
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
        try:
            self._shm = SharedMemory(name=self.name, create=True, size=self.size)
        except FileExistsError:
            # If it exists, try to remove it and create again
            try:
                SharedMemory(name=self.name).unlink()
                self._shm = SharedMemory(name=self.name, create=True, size=self.size)
            except Exception as e:
                raise RuntimeError(f"Failed to create shared memory on Windows: {e}")

    def _create_unix(self):
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
