"""Reading and writing the formats this package speaks.

Serialisation, kept apart from the models that produce and consume it. A format
is a contract with everything outside the process -- other tools, older files,
other groups -- and it changes for reasons that have nothing to do with the
physics; leaving readers next to the model they happen to feed is how two
dialects of one format come to exist without anyone deciding.

* ``fps.json`` -- the authored definition (every field mapped to its flrCIF item
  where the dictionary defines one), the one reader and writer, and the legacy
  C# FPS ``.txt`` formats, **read only**. All C++ now:
  `include/IMP/bff/FPSSchema.h` and `FPSIO.h`. The C++ functions speak JSON
  text -- a caller passes a dict as `json.dumps(...)` and parses the returned
  text with `json.loads`.
* :mod:`IMP.bff.io.cif` -- BinaryCIF and mmCIF for dye templates and force
  fields.
* PDB, MOL2, mmCIF, DCD and coordinate comparison -- C++, in
  `include/IMP/bff/StructureIO.h` and `TrajectoryIO.h`. The **RMF** door stays
  Python, lazily: writing RMF needs `IMP.rmf`, which is not one of this module's
  `required_modules`.

Moved out of ``fret/`` by PRD-113 stage 7 and split three ways there. One module
held all of it, which is how ``fret/io.py`` came to be the thing that writes
PDBs -- and why ``IMP.bff.fret.read_fps_json`` no longer exists.

Port status (PRD-117): 15 of 17 names are C++-covered; the remaining 2
(``write_rmf`` and the lazy RMF writers) cannot be C++ without widening
`required_modules`. Pure re-export -- holds no logic of its own to move.
"""

from IMP.bff import (  # noqa: F401
    fps_av_simulation_types,
    fps_distance_types,
    fps_positions_for_docking,
    fps_simulation_types,
    read_evaluators_json,
    read_fps_json,
    read_old_distances_txt,
    read_old_lps_txt,
    write_evaluators_json,
    write_fps_json,
)
from IMP.bff import (  # noqa: F401
    compute_rmsd,
    load_structure,
    load_structure_with_particles,
    read_dcd,
    read_trajectory,
    write_pdb,
    write_rmf,
)

__all__ = [
    "fps_av_simulation_types", "fps_distance_types", "fps_positions_for_docking",
    "fps_simulation_types", "read_evaluators_json", "read_fps_json",
    "read_old_distances_txt", "read_old_lps_txt", "write_evaluators_json",
    "write_fps_json",
    "compute_rmsd", "load_structure", "load_structure_with_particles",
    "read_dcd", "read_trajectory", "write_pdb", "write_rmf",
]
