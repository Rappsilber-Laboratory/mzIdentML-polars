from ._mzidentml_polars import *

try:
    from ._mzidentml_polars import __version__
except ImportError:
    try:
        from ._version import version as __version__
    except ImportError:
        __version__ = "unknown"
