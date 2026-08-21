"""Simulation runners for cgdye.

The MD/MC runner (run_dye_simulation) is in :file:`pyext/IMP_bff.sim.i` as
%pythoncode. Uses IMP.pmi and IMP.rmf (deferred inside functions). This module
re-exports the surface.
"""

from IMP.bff import (
    run_dye_simulation,
)

__all__ = [
    "run_dye_simulation",
]
