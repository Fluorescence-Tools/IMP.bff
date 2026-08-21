"""Stage 2 -- is this dye configuration allowed, and how heavily does it count?

<!-- port-status: PRD-117. 1 of 34 re-exported names are C++-covered already
     (``compute_rotamer_score`` et al. -- the other 33 are %pythoncode in
     ``pyext/IMP_bff.scoring.i`` and remain to be ported). This module is pure
     re-export; the port just deletes more of the surface it copies. -->

The score of a *configuration*, from the structure alone. A rotamer that clashes
with the backbone is not a rotamer the dye adopts; a conformer at high internal
energy is one it adopts rarely. This is what turns a raw set of candidate
positions into a weighted ensemble, and it is the same question whichever
representation produced the candidates:

======================  ==============================================
representation          what scoring it means
======================  ==============================================
accessible volume       occupancy: a voxel is reachable or it is not
rotamer library         clash energy of each conformer against the site
coarse-grained dye      the simplified force field's internal energy
======================  ==============================================

**Not the same thing as** :mod:`IMP.bff.restraints`. That answers *does this
model agree with a measurement* -- it takes experimental distances and returns
a chi-squared. This package never sees an experiment: it takes coordinates and
returns an energy or a weight. Two different questions that both get called
"scoring", which is why they are two packages and why this paragraph exists.

The C++ kernels and the orchestration around them live in
:file:`include/IMP/bff/Scoring.h` and :file:`include/IMP/bff/RotamerEnergy.h`,
wrapped through :file:`pyext/IMP_bff.scoring.i`. This module re-exports the
Python surface so ``from IMP.bff.scoring import ...`` keeps working.
"""

from IMP.bff import (
    CHARMM36_LJ,
    BoundingBoxFilter,
    DyeInternalEnergyEvaluator,
    RotamerScoreResultPy,
    _aabb_overlap,
    _atom_type,
    _build_cross_lj_params,
    _hydrogen_mask,
    _lj_energy_pairs,
    _pair_energy_matrix,
    _protein_charge_mask,
    _rotamer_charge_mask,
    _scaled_parameters,
    _selector_atom_names,
    _selector_matches,
    _site_mask,
    boltzmann_weights,
    build_dye_restraints,
    build_lj_type_table,
    compute_lj_pair_sites,
    compute_rotamer_score,
    dye_internal_system,
    lj_cross,
    lj_cross_params,
    lj_energy,
    lj_parameter_arrays,
    lj_params,
    lj_score,
    rotamer_cluster_weights,
    rotamer_mean_field_weights,
    rotamer_mean_field_weights_multi_dye,
    selector_resnames,
    site_element_map,
    torsion_cosine,
)

__all__ = [
    'CHARMM36_LJ',
    'BoundingBoxFilter',
    'DyeInternalEnergyEvaluator',
    'RotamerScoreResultPy',
    'boltzmann_weights',
    'build_dye_restraints',
    'build_lj_type_table',
    'compute_lj_pair_sites',
    'compute_rotamer_score',
    'dye_internal_system',
    'lj_cross',
    'lj_cross_params',
    'lj_energy',
    'lj_parameter_arrays',
    'lj_params',
    'lj_score',
    'rotamer_cluster_weights',
    'rotamer_mean_field_weights',
    'rotamer_mean_field_weights_multi_dye',
    'selector_resnames',
    'site_element_map',
    'torsion_cosine',
]
