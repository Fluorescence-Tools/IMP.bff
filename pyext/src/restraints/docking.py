"""IMP + IMP.bff FRET-restrained rigid-body docking engine.

The docking engine (build_assembly, score, dock, dock_minimize, refine, screen)
is in :file:`pyext/IMP_bff.docking.i` as %pythoncode. Uses IMP.pmi (deferred
inside functions). This module re-exports the surface.
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
