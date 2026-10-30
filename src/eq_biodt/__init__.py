"""Empty Quarter Biodiversity Digital Twin (EQ-BioDT)."""

from importlib.metadata import PackageNotFoundError, version as _version

try:
    __version__ = _version("eq-biodt")
except PackageNotFoundError:
    __version__ = "0.6.0.dev0"
