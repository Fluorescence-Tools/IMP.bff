"""Rotamer-based cgdye utilities.

The C++ kernels (FRETPair, RotamerEnergy, States, RotamerStatistics) are in
their own headers. The rotamer library IO, ensemble, and FRET driver are in
:file:`pyext/IMP_bff.rotamer.i` as %pythoncode. This module re-exports the surface.

Port status (PRD-117): 1 of 20 names is C++-covered; the remaining 19 (the
rotamer library IO, ``RotamerEnsemble``/``RotamerFRET``, ``*_fps`` and
``rotamer_*`` payload/registry helpers) are %pythoncode in
``rotamer.i``/``rotamer_ensemble.i`` and still to port. Pure re-export.
"""

from IMP.bff import (
    SIMULATION_TYPE_R1,
    FRETFrameResult,
    RotamerDistance,
    RotamerEnsemble,
    RotamerFRET,
    RotamerPosition,
    distances_from_ensembles,
    load_protein_frames,
    load_rotamer_library,
    read_rotamer_fps,
    resolve_rotamer_library_path,
    rotamer_ensemble_payload,
    rotamer_ensembles_from_fps,
    rotamer_fret_from_fps,
    rotamer_library_metadata,
    rotamer_library_registry,
    rotamer_position_payload,
    selector_atom_indices,
    transform_library_to_site,
    write_rotamer_fps,
)

__all__ = [
    "FRETFrameResult",
    "RotamerDistance",
    "RotamerEnsemble",
    "RotamerFRET",
    "RotamerPosition",
    "distances_from_ensembles",
    "load_protein_frames",
    "load_rotamer_library",
    "read_rotamer_fps",
    "resolve_rotamer_library_path",
    "rotamer_ensemble_payload",
    "rotamer_ensembles_from_fps",
    "rotamer_fret_from_fps",
    "rotamer_library_metadata",
    "rotamer_library_registry",
    "rotamer_position_payload",
    "selector_atom_indices",
    "transform_library_to_site",
    "write_rotamer_fps",
]
