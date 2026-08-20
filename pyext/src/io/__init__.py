"""Reading and writing the formats this package speaks.

Serialisation, kept apart from the models that produce and consume it. A format
is a contract with everything outside the process -- other tools, older files,
other groups -- and it changes for reasons that have nothing to do with the
physics; leaving readers next to the model they happen to feed is how two
dialects of one format come to exist without anyone deciding.

* ``fps.json`` -- the authored definition (every field mapped to its flrCIF item
  where the dictionary defines one), the one reader and writer, and the legacy
  C# FPS ``.txt`` formats, **read only**. All C++ now:
  `include/IMP/bff/FPSSchema.h` and `FPSIO.h`, with the marshalling in
  `pyext/IMP_bff.fps.i`. A format nobody can still read is data that has been
  lost, and a decade of measurements live in those ``.txt`` files; nothing
  should write them again.
* :mod:`IMP.bff.io.cif` -- BinaryCIF and mmCIF for dye templates and force
  fields.
* PDB, MOL2, mmCIF, DCD and coordinate comparison -- C++, in
  `include/IMP/bff/StructureIO.h` and `TrajectoryIO.h`. The **RMF** door stays
  Python, lazily: writing RMF needs `IMP.rmf`, which is not one of this module's
  `required_modules`.

Moved out of ``fret/`` by PRD-113 stage 7 and split three ways there. One module
held all of it, which is how ``fret/io.py`` came to be the thing that writes
PDBs -- and why ``IMP.bff.fret.read_fps_json`` no longer exists.
"""

from IMP.bff import (  # noqa: F401
    AV_SIMULATION_TYPES,
    DISTANCE_TYPES,
    SIMULATION_TYPES,
    fps_positions_for_docking,
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
    "AV_SIMULATION_TYPES", "DISTANCE_TYPES", "SIMULATION_TYPES",
    "fps_positions_for_docking", "read_evaluators_json", "read_fps_json",
    "read_old_distances_txt", "read_old_lps_txt", "write_evaluators_json",
    "write_fps_json",
    "compute_rmsd", "load_structure", "load_structure_with_particles",
    "read_dcd", "read_trajectory", "write_pdb", "write_rmf",
]
