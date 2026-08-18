"""FRET-restrained docking, accessible volumes, and the fps.json format.

This package is the home of the FRET structural-modelling algorithms and of
the fps.json data schema (PRD-97): the docking engine (:mod:`.imp_engine`),
the AV backend (:mod:`.av`),
P(R_DA) distributions (:mod:`.distributions`), Olga-style informative pair
selection (:mod:`.greedy_olga`), docking-precision estimation
(:mod:`.uncertainty`) and PMI stat-file reading (:mod:`.stat`).

Model-distance calculation moved to :mod:`IMP.bff.representation.distance` in
stage 4b and the 28-line re-export shim left behind was retired in the cleanup:
a distance between two labels is a property of what *represents* them -- an
accessible volume, a rotamer library, a Gaussian -- not of the FRET engine.

File I/O and the fps.json schema moved to :mod:`IMP.bff.io` in PRD-113 stage 7
and are **not** re-exported here: ``IMP.bff.fret.read_fps_json`` is gone, and
``IMP.bff.io.read_fps_json`` (or the flat ``IMP.bff.read_fps_json``) replaces
it. A FRET model should not be the thing that owns a file format.

Everything here is code/script level — no GUI. Applications (e.g. ChiSurf)
build their workflow and views on top of these functions.
"""

from . import av
from . import engine
from . import distributions
from . import greedy_olga
from . import stat
from . import uncertainty
from . import imp_engine

from .av import (
    AccessibleVolume,
    compute_av,
    compute_avs_for_structure,
    load_structure_with_vdw,
)
from .engine import RigidBody, DistanceRestraint, SpringParameters
from .greedy_olga import select_informative_pairs
from .uncertainty import estimate_position_uncertainty
from .distributions import compute_distance_distributions
from .imp_engine import (
    DockingParameters,
    DockingResult,
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
    ensure_fps_json,
)

__all__ = [
    "av", "engine", "distributions",
    "greedy_olga", "stat", "uncertainty", "imp_engine",
    # av
    "AccessibleVolume", "compute_av", "compute_avs_for_structure",
    "load_structure_with_vdw",
    # engine
    "RigidBody", "DistanceRestraint", "SpringParameters",
    # algorithms
    "select_informative_pairs", "estimate_position_uncertainty",
    "compute_distance_distributions",
    # docking
    "DockingParameters", "DockingResult", "PairDistance", "build_assembly",
    "capture_poses", "apply_poses", "dock", "dock_minimize", "refine",
    "screen", "score", "estimate_errors", "ensure_fps_json",
]
