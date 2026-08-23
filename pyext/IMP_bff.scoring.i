/*
 * Stage-2 scoring, from C++ to the Python surface it had as a module.
 *
 * The inner kernels are in RotamerEnergy.h, wrapped directly in swig.i-in.
 * The orchestration -- the CHARMM36 table, the Lorentz-Berthelot rules, the
 * Boltzmann weight, the AABB pre-filter, the masks and selectors, the
 * end-to-end rotamer score, the single-dye mean-field update and the
 * typed-system walkers -- is C++ now, in Scoring.h. Three things remain here:
 *
 * - `rotamer_mean_field_weights_multi_dye`: numpy iteration over the C++
 *   pair-energy matrix, one matrix per dye pair. The single-dye update is
 *   C++ (`rotamer_mean_field_weights`).
 * - `dye_internal_system`: assembles a typed system from `parse_dye_mol2`
 *   output, which is Python (topology.i) -- it moves with that batch.
 * - `torsion_cosine` / `build_dye_restraints`: they create IMP.core objects
 *   (Cosine, DistanceRestraint), which is IMP API glue and stays Python for
 *   the same reason every IMP API caller does.
 */

%include "IMP/bff/Scoring.h"

// RotamerScoreResult is a value type with vector members.
IMP_SWIG_VALUE(IMP::bff, RotamerScoreResult, RotamerScoreResults);

%attribute_np(IMP::bff::RotamerScoreResult, std::vector<double>, weights, get_weights);
%attribute_np(IMP::bff::RotamerScoreResult, std::vector<double>, energies, get_energies);
%attribute_py(IMP::bff::RotamerScoreResult, double, partition, get_partition);

%template(LJSitePairList) std::vector<IMP::bff::LJSitePair>;

%pythoncode %{
import math as _math
import numpy as _np


def rotamer_mean_field_weights_multi_dye(
    rotamer_coords_list, initial_weights_list, protein_coords,
    dye_elements_list, protein_elements, K=1.0, n_iter=10,
    aabb_pad=3.5, r_cutoff=12.0,
):
    """Multi-dye mean-field weight update including dye-dye cross-interactions.

    Iterates in numpy over the C++ conformer-pair energy matrix
    (:func:`pair_energy_matrix`), one matrix per dye pair; the weight update
    itself is the same log-sum-exp form the single-dye C++
    :func:`rotamer_mean_field_weights` performs.
    """
    n_dyes = len(rotamer_coords_list)
    q_list = [_np.asarray(w, dtype=_np.float64).copy()
              for w in initial_weights_list]
    E_bb_all = [
        _np.asarray(pair_energy_matrix(
            rotamer_coords_list[d], protein_coords,
            dye_elements_list[d], protein_elements,
            len(q_list[d]), 1, r_cutoff, aabb_pad),
            dtype=_np.float64).reshape(len(q_list[d]))
        for d in range(n_dyes)
    ]
    E_sc = {}
    for d1 in range(n_dyes):
        for d2 in range(d1 + 1, n_dyes):
            m = _np.asarray(pair_energy_matrix(
                rotamer_coords_list[d1], rotamer_coords_list[d2],
                dye_elements_list[d1], dye_elements_list[d2],
                len(q_list[d1]), len(q_list[d2]), r_cutoff, aabb_pad),
                dtype=_np.float64).reshape(len(q_list[d1]), len(q_list[d2]))
            E_sc[(d1, d2)] = m
            E_sc[(d2, d1)] = m.T
    for _iteration in range(n_iter):
        E_P_all = [E_bb_all[d].copy() for d in range(n_dyes)]
        for d1 in range(n_dyes):
            for d2 in range(n_dyes):
                if d1 == d2:
                    continue
                E_P_all[d1] += E_sc[(d1, d2)] @ q_list[d2]
        for d in range(n_dyes):
            log_q = _np.log(_np.maximum(q_list[d], 1e-300)) - K * E_P_all[d]
            log_q -= log_q.max()
            q_new = _np.exp(log_q)
            denom = q_new.sum()
            q_list[d] = (q_new / denom if denom > 0
                         else _np.ones_like(q_new) / len(q_new))
    return q_list


def dye_internal_system(atoms, bonds):
    """A typed system for a lone dye from `parse_dye_mol2` output.

    Python because its input is: `parse_dye_mol2` (topology.i) hands back
    dicts, and this bridges them into the typed value. The graph helpers it
    calls (build_graph / build_angles / build_dihedrals) are C++ already.
    """
    import json as _json
    ordered = sorted(atoms.values(), key=lambda x: x["serial"])
    sid = {a["serial"]: f"dye:{a['serial']}:{a['atom_name']}" for a in ordered}
    sites = [{"id": sid[a["serial"]], "atom_name": a["atom_name"]}
             for a in ordered]
    graph = build_graph(bonds)
    return forcefield_system_from_json(_json.dumps({
        "sites": sites,
        "bonds": [[sid[a], sid[b], 0.0, None] for a, b in sorted(bonds)],
        "angles": [[sid[a], sid[b], sid[c], 0.0, None]
                   for a, b, c in build_angles(graph)],
        "dihedrals": [[sid[a], sid[b], sid[c], sid[d], None]
                      for a, b, c, d in build_dihedrals(graph)],
    }))


def torsion_cosine(torsion_type):
    """The ``IMP.core.Cosine`` for a torsion type stored in the CHARMM convention."""
    if isinstance(torsion_type, dict):
        k = torsion_type["k"]
        n = torsion_type["periodicity"]
        phase = torsion_type["phase_rad"]
    else:
        k, n, phase = torsion_type.k, torsion_type.periodicity, torsion_type.phase
    return IMP.core.Cosine(float(k), int(n), float(phase) + _math.pi)


def build_dye_restraints(model, system, site_particles):
    """Build bonded + element-aware nonbonded restraints."""
    system = as_forcefield_system(system)
    restraints = []
    bt = system.bond_types
    at = system.angle_types
    tt = system.torsion_types
    it = system.improper_types
    for bd in system.bonds:
        if bd.site_a not in site_particles or bd.site_b not in site_particles:
            continue
        p1, p2 = site_particles[bd.site_a], site_particles[bd.site_b]
        d0 = float(bd.length)
        k = float(bt[bd.type_id])
        restraints.append(
            IMP.core.DistanceRestraint(model, IMP.core.Harmonic(d0, k), p1, p2))
    for an in system.angles:
        if (an.site_a not in site_particles or an.site_b not in site_particles
                or an.site_c not in site_particles):
            continue
        p1, p2, p3 = (site_particles[an.site_a], site_particles[an.site_b],
                      site_particles[an.site_c])
        k = float(at[an.type_id])
        restraints.append(
            IMP.core.AngleRestraint(
                model, IMP.core.Harmonic(float(an.theta), k), p1, p2, p3))
    for to in system.dihedrals:
        if any(x not in site_particles
               for x in (to.site_a, to.site_b, to.site_c, to.site_d)):
            continue
        p1, p2, p3, p4 = (site_particles[to.site_a], site_particles[to.site_b],
                          site_particles[to.site_c], site_particles[to.site_d])
        restraints.append(IMP.core.DihedralRestraint(
            model, torsion_cosine(tt[to.type_id]), p1, p2, p3, p4))
    for to in system.impropers:
        if any(x not in site_particles
               for x in (to.site_a, to.site_b, to.site_c, to.site_d)):
            continue
        p1, p2, p3, p4 = (site_particles[to.site_a], site_particles[to.site_b],
                          site_particles[to.site_c], site_particles[to.site_d])
        t = it[to.type_id]
        k = t["k"] if isinstance(t, dict) else t.k
        theta0 = IMP.core.get_dihedral(
            IMP.core.XYZ(p1), IMP.core.XYZ(p2), IMP.core.XYZ(p3),
            IMP.core.XYZ(p4))
        fun = IMP.core.Harmonic(theta0, float(k))
        restraints.append(IMP.core.DihedralRestraint(model, fun, p1, p2, p3, p4))
    excluded = {frozenset(p) for p in system.exclusions()}
    for pair in compute_lj_pair_sites(system):
        if frozenset({pair.site_a, pair.site_b}) in excluded:
            continue    # compute_lj_pair_sites already excluded these; kept
                        # for parity with the bonded loops above
        if pair.site_a not in site_particles or pair.site_b not in site_particles:
            continue
        p1, p2 = site_particles[pair.site_a], site_particles[pair.site_b]
        if pair.eps <= 0:
            continue
        restraints.append(IMP.core.DistanceRestraint(
            model,
            IMP.core.HarmonicLowerBound(float(pair.rmin), float(pair.eps)),
            p1, p2))
    return restraints
%}