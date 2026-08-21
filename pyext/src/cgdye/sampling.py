"""Samplers for explicit dyes: linker Metropolis, rotamer libraries, RRT, kinetics.

The C++ kernels (Clustering, LinkerGeometry) are in their own headers. The
samplers (LangevinDyeSampler, LinkerSampler, RRT) and kinetic helpers are in
:file:`pyext/IMP_bff.sampling.i` as %pythoncode. This module re-exports the surface.

Port status (PRD-117): 1 of 24 names is C++-covered; the remaining 23 (all the
samplers and RRT/kinetic helpers) are %pythoncode here and still to port. Pure
re-export -- no logic of its own to move.
"""

from IMP.bff import (
    KB_KCAL,
    IMPRRTTree,
    LangevinDyeSampler,
    LangevinTrajectory,
    LinkerSampler,
    TorsionRRTTree,
    assign_frames_to_clusters,
    cluster_frames_leader,
    generate_linker_rotamers,
    is_collision_sphere,
    make_langevin_simulator,
    make_transform,
    reconstruct_rotamer_trajectory,
    rrt_grow_step,
    rotamer_correlation_times,
    rotamer_rotational_correlation_time,
    rotamer_transition_matrix,
    rmsd_no_align,
    run_rigid_body_rrt,
    run_torsion_rrt,
    sample_transform,
    steer_transform,
    torsion_distance,
    transform_distance,
)

__all__ = [
    "KB_KCAL",
    "IMPRRTTree",
    "LangevinDyeSampler",
    "LangevinTrajectory",
    "LinkerSampler",
    "TorsionRRTTree",
    "assign_frames_to_clusters",
    "cluster_frames_leader",
    "generate_linker_rotamers",
    "is_collision_sphere",
    "make_langevin_simulator",
    "make_transform",
    "reconstruct_rotamer_trajectory",
    "rrt_grow_step",
    "rotamer_correlation_times",
    "rotamer_rotational_correlation_time",
    "rotamer_transition_matrix",
    "rmsd_no_align",
    "run_rigid_body_rrt",
    "run_torsion_rrt",
    "sample_transform",
    "steer_transform",
    "torsion_distance",
    "transform_distance",
]
