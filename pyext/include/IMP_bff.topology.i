/*
 * The probe topology: MOL2 in, a typed force-field system out.
 *
 *  - `create_forcefield_system` takes `FFComponentSpec` values, so a caller
 *    describes components directly rather than marshalling a dict.
 *    `create_probe_protein_system` and `probe_forcefield_system` are two- and
 *    one-component calls of it, the second adding the probe's `N`/`CA`/`C`/`O`
 *    anchor group. All three are in `TopologyBuild.h`. The impropers are the
 *    templates': a builder that fills `impropers = []` costs 81 of them on the
 *    shipped components (okf/validation/impropers_are_dropped.md).
 *  - Building, writing and reporting a system from `name=X,mol2=Y` strings is
 *    CLI, and lives in `bin/imp_bff`.
 *  - Node numbering belongs with the graph, not with the caller: a
 *    first-appearance numbering and a sorted one give different orders out,
 *    so `MolecularGraph` keys by the caller's int, `LabelledGraph` by a string
 *    site id, and the rotor and ring questions are methods on them.
 *  - A distance is `IMP.algebra.get_distance` and an angle is
 *    `bond_angle_rad`; the typed `Mol2Component` the reader returns is what
 *    they act on, with no dict in between.
 *
 * The CHARMM36 LJ table is `charmm36_lj(element)` in `Scoring.h` -- one table,
 * which is what keeps the two scorers from drifting apart.
 */

%include "IMP/bff/TopologyBuild.h"

%template(FFComponentSpecList) std::vector<IMP::bff::FFComponentSpec>;
