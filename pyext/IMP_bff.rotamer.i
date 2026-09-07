/*
 * The rotamer layer, part one: placing a library on a residue.
 *
 *  - The backbone frame, the atom selectors and the library registry are
 *    `RotamerSite.h` (included here): nine multiplications a caller should not
 *    re-type, a selector grammar with its own ambiguity rules, and a JSON
 *    lookup with a cutoff-suffix grammar.
 *  - `load_rotamer_library` and `load_protein_frames` are `RotamerSite.h` and
 *    `HierarchyFrame.h`. The registry answers in JSON text, which is what the
 *    rest of this module's file layer speaks -- a caller that wants dicts and
 *    `Path`s makes them.
 *  - `RotamerEnsemble`, the fps payload layer and the `RotamerFRET` driver are
 *    in `IMP_bff.rotamer_ensemble.i`, which comes after `States`.
 *  - The RMF readers are `RmfIO.h`: `rmf` is one of this module's
 *    `required_modules`, so they are ordinary C++ like every other reader.
 */

%include "IMP/bff/RotamerSite.h"
