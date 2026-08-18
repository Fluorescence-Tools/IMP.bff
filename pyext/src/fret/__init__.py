"""FRET-restrained docking, accessible volumes, and the fps.json format.

This package is the home of the FRET structural-modelling algorithms and of
the fps.json data schema (PRD-97): the docking engine (:mod:`.imp_engine`),
the AV backend (:mod:`.av`), model-distance calculation (:mod:`.distance`),
P(R_DA) distributions (:mod:`.distributions`), Olga-style informative pair
selection (:mod:`.olga_greedy`), docking-precision estimation
(:mod:`.uncertainty`) and PMI stat-file reading (:mod:`.stat`).

File I/O and the fps.json schema moved to :mod:`IMP.bff.io` in PRD-113 stage 7
and are **not** re-exported here: ``IMP.bff.fret.read_fps_json`` is gone, and
``IMP.bff.io.read_fps_json`` (or the flat ``IMP.bff.read_fps_json``) replaces
it. A FRET model should not be the thing that owns a file format.

Everything here is code/script level — no GUI. Applications (e.g. ChiSurf)
build their workflow and views on top of these functions.
"""

from . import av
from . import distance
from . import engine
from . import distributions
from . import olga_greedy
from . import stat
from . import uncertainty
from . import imp_engine

from .av import (
    AccessibleVolume,
    compute_av,
    compute_avs_for_structure,
    load_structure_with_vdw,
)
from .distance import (
    average_distance,
    mean_fret_distance,
    distance_between_mean_positions,
    model_distance,
    chi2_score,
)
from .engine import RigidBody, DistanceRestraint, SpringParameters
from .olga_greedy import select_informative_pairs
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
    "av", "distance", "engine", "distributions",
    "olga_greedy", "stat", "uncertainty", "imp_engine",
    # av
    "AccessibleVolume", "compute_av", "compute_avs_for_structure",
    "load_structure_with_vdw",
    # distance
    "average_distance", "mean_fret_distance",
    "distance_between_mean_positions", "model_distance", "chi2_score",
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
