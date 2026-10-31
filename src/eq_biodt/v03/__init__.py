"""EQ-BioDT causal-ontological twin (Python and HTTP paths keep the historical ``v03``/``/v3`` names)."""
from importlib.metadata import PackageNotFoundError, version as _version

try:
    VERSION = _version("eq-biodt")
except PackageNotFoundError:  # running from an uninstalled source tree
    VERSION = "0.7.0"
