# Author: Thomas-Otavio Peulen <thomas@peulen.xyz>

from __future__ import annotations
import sys

from . import restraints  # noqa: F401
from . import fret  # noqa: F401

# Take care of typing
try:
    if sys.version_info >= (3, 8):
        import typing
    elif sys.version_info >= (3, 7):
        # monkey patch the 3.7 typing system as TypedDict etc. is missing
        import typing
        try:
            import typing_extensions
            for key in typing_extensions.__dict__.keys():
                f = typing_extensions.__dict__[key]
                if callable(f):
                    typing.__dict__[key] = f
        except ModuleNotFoundError:
            print("WARNING typing_extensions not found", file=sys.stderr)
    else:
        import typing_extensions as typing
except ModuleNotFoundError:
    print("WARNING typing not found", file=sys.stderr)
    typing = None

__version__ = "0.14.0"
