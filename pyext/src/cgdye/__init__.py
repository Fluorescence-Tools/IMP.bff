"""Explicit (coarse-grained) dye modelling for IMP.bff -- the *source* tree.

Users reach every public name flat as ``IMP.bff.<Name>`` (see
``IMP.bff.api``); the sub-packages here (``labeling``, ``topology``,
``sampling``, ``rotamer``, ``io``, ``analysis``, ``sim``, ``scripts``) are how
the code is organised. Importing this package is cheap and needs no click.
"""

from __future__ import annotations

from IMP.bff.api import EXPORTS as _EXPORTS

#: the flat public names that live in cgdye (subset of IMP.bff.api.EXPORTS)
__all__ = sorted(n for n, m in _EXPORTS.items() if m.startswith("IMP.bff.cgdye."))


def __getattr__(name):
    """Resolve a flat public name lazily (PEP 562)."""
    if name in __all__:
        from IMP.bff.api import resolve
        value = resolve(name)
        globals()[name] = value
        return value
    raise AttributeError(f"module 'IMP.bff.cgdye' has no attribute {name!r}")


def __dir__():
    return sorted(set(globals()) | set(__all__))
