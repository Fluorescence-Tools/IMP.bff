"""FRET-restrained docking, accessible volumes, and the fps.json format.

This package is the home of the FRET structural-modelling algorithms and of
the fps.json data schema (PRD-97): the docking engine (:mod:`.imp_engine`),
the AV backend (:mod:`.av`), model-distance calculation (:mod:`.distance`),
P(R_DA) distributions (:mod:`.distributions`), Olga-style informative pair
selection (:mod:`.olga_greedy`), docking-precision estimation
(:mod:`.uncertainty`), PMI stat-file reading (:mod:`.stat`), file I/O
(:mod:`.io`) and the authored fps.json schema (:mod:`.fps_schema`).

Everything here is code/script level — no GUI. Applications (e.g. ChiSurf)
build their workflow and views on top of these functions.
"""

from . import fps_schema
from . import io
from . import av
from . import distance
from . import engine
from . import distributions
from . import olga_greedy
from . import stat
from . import uncertainty
from . import imp_engine

from .io import (
    read_fps_json,
    write_fps_json,
    read_evaluators_json,
    write_evaluators_json,
    read_old_lps_txt,
    read_old_distances_txt,
    load_structure,
    write_pdb,
    compute_rmsd,
)
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
    "fps_schema", "io", "av", "distance", "engine", "distributions",
    "olga_greedy", "stat", "uncertainty", "imp_engine",
    # io
    "read_fps_json", "write_fps_json", "read_evaluators_json",
    "write_evaluators_json", "read_old_lps_txt", "read_old_distances_txt",
    "load_structure", "write_pdb", "compute_rmsd",
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
