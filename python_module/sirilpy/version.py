"""Version information for sirilpy package."""

try:  # import from the packaging specification
    from importlib.metadata import metadata, PackageNotFoundError
    meta = metadata("sirilpy")
    __version__ = meta.get("version", "unknown")
    __author__ = meta.get("author", "unknown")
    __license__ = meta.get("license", "unknown")
except (ImportError, PackageNotFoundError):
    # Specific exceptions rather than general Exception
    __version__ = "unknown"
    __author__ = "unknown"
    __license__ = "unknown"

__copyright__ = " (c) Team free-astro 2024-2025"  # not a standard metadata
