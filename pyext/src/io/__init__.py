"""Reading and writing the formats this package speaks.

Serialisation, kept apart from the models that produce and consume it. A format
is a contract with everything outside the process -- other tools, older files,
other groups -- and it changes for reasons that have nothing to do with the
physics; leaving readers next to the model they happen to feed is how two
dialects of one format come to exist without anyone deciding.

* :mod:`IMP.bff.io.fps` -- the authored definition of ``fps.json``, every
  field mapped to its flrCIF item where the dictionary defines one.
* :mod:`IMP.bff.io.fps` -- the one reader and writer for ``fps.json``.
* :mod:`IMP.bff.io.fps` -- the legacy C# FPS ``.txt`` formats,
  **read only**. A format nobody can still read is data that has been lost, and
  a decade of measurements live in these files; nothing should write them again.
* :mod:`IMP.bff.io.structure` -- PDB and RMF, and coordinate comparison.

Moved out of ``fret/`` by PRD-113 stage 7 and split three ways there. One module
held all of it, which is how ``fret/io.py`` came to be the thing that writes
PDBs -- and why ``IMP.bff.fret.read_fps_json`` no longer exists. The dye
templates and spectra readers belong here too and have not moved yet; they are
still in :mod:`IMP.bff.dye`, next to the ``Dye`` they build.
"""

from IMP.bff.io.fps import (  # noqa: F401
    AV_SIMULATION_TYPES,
    fps_positions_for_docking,
    read_evaluators_json,
    read_fps_json,
    write_evaluators_json,
    write_fps_json,
)
from IMP.bff.io.fps import (  # noqa: F401
    read_old_distances_txt,
    read_old_lps_txt,
)
from IMP.bff.io.structure import (  # noqa: F401
    compute_rmsd,
    load_structure,
    load_structure_with_particles,
    write_pdb,
    write_rmf,
)

__all__ = [
    "AV_SIMULATION_TYPES", "fps_positions_for_docking", "read_evaluators_json",
    "read_fps_json", "write_evaluators_json", "write_fps_json",
    "read_old_distances_txt", "read_old_lps_txt",
    "compute_rmsd", "load_structure", "load_structure_with_particles",
    "write_pdb", "write_rmf",
]
