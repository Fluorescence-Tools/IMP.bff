"""IMP + IMP.bff FRET-restrained rigid-body docking engine.

The docking engine (build_assembly, score, dock, dock_minimize, refine, screen)
is in :file:`pyext/IMP_bff.docking.i` as %pythoncode. Uses IMP.pmi (deferred
inside functions). This module re-exports the surface.

Port status (PRD-117): 0 of 13 names are C++-covered; all 13 (``build_assembly``,
``dock*``, ``refine``, ``screen``, the dock result/parameters) are %pythoncode
in ``pyext/IMP_bff.docking.i`` and still to port. Pure re-export.
"""

from IMP.bff import (
    DockingParameters,
    DockingResult,
    MeanDistanceRestraint,
    PairDistance,
    build_assembly,
    capture_poses,
    apply_poses,
    dock,
    dock_minimize,
    refine,
    screen,
    score,
    estimate_errors,
)

__all__ = [
    "DockingParameters",
    "DockingResult",
    "MeanDistanceRestraint",
    "PairDistance",
    "build_assembly",
    "capture_poses",
    "apply_poses",
    "dock",
    "dock_minimize",
    "refine",
    "screen",
    "score",
    "estimate_errors",
]
