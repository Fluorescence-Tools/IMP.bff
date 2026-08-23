/*
 * Stage-2 scoring, from C++ to the Python surface it had as a module.
 *
 * The inner kernels (all-pairs steric/electrostatic, the pair energy matrix,
 * the per-frame LJ over an explicit pair list) are in RotamerEnergy.h and
 * wrapped directly in swig.i-in. What lives here is the orchestration around
 * them: the CHARMM36 table, the Lorentz-Berthelot combining rules, the
 * Boltzmann weight, the AABB pre-filter, and the end-to-end rotamer score.
 *
 * The mask-building helpers (_atom_type, _selector_matches, _site_mask,
 * _hydrogen_mask, _protein_charge_mask, _rotamer_charge_mask) operate on
 * Python lists of atom names and residue names. They are called once per
 * structure, not per frame, and they stay as %pythoncode because their input
 * is Python lists -- porting them to C++ would mean designing a C++ API for
 * FRETpredict-style selector string parsing, which is glue, not a kernel.
 *
 * torsion_cosine and build_dye_restraints are IMP API glue (they create
 * IMP.core.Cosine / IMP.core.DistanceRestraint objects) and stay Python for
 * the same reason every IMP API caller does.
 */

// charmm36_lj and lj_cross return std::array<double, 2>, which SWIG does not
// wrap by default. They are reimplemented in %pythoncode below (returning
// Python tuples/dicts), so the C++ versions are internal only.
%ignore IMP::bff::charmm36_lj;
%ignore IMP::bff::lj_cross;
%ignore IMP::bff::lj_parameter_arrays;

%include "IMP/bff/Scoring.h"

// RotamerScoreResult is a value type with vector members.
IMP_SWIG_VALUE(IMP::bff, RotamerScoreResult, RotamerScoreResults);

%attribute_np(IMP::bff::RotamerScoreResult, std::vector<double>, weights, get_weights);
%attribute_np(IMP::bff::RotamerScoreResult, std::vector<double>, energies, get_energies);
%attribute_py(IMP::bff::RotamerScoreResult, double, partition, get_partition);

%pythoncode %{
import math as _math
import re as _re
from functools import lru_cache as _lru_cache
from dataclasses import dataclass as _dataclass

CHARMM36_LJ = {
    "C": {"rmin_half": 2.02446316, "epsilon": -0.06394724},
    "N": {"rmin_half": 1.89285714, "epsilon": -0.15428571},
    "O": {"rmin_half": 1.693, "epsilon": -0.12642017},
    "S": {"rmin_half": 2.1, "epsilon": -0.47},
    "H": {"rmin_half": 0.98357778, "epsilon": -0.03466645},
}


def lj_params(element):
    return CHARMM36_LJ.get(element, CHARMM36_LJ["C"])


def lj_cross(elem_i, elem_j):
    """Lorentz-Berthelot cross parameters ``(rmin, eps)`` for two elements."""
    pi = lj_params(elem_i)
    pj = lj_params(elem_j)
    rmin = pi["rmin_half"] + pj["rmin_half"]
    eps = (pi["epsilon"] * pj["epsilon"]) ** 0.5
    return rmin, eps


def lj_cross_params(elem_i, elem_j):
    """Lorentz-Berthelot combining rules for two elements. Returns (rmin_ij, epsilon_ij)."""
    return lj_cross(elem_i, elem_j)


def lj_parameter_arrays(elements):
    """``(rmin_half, epsilon)`` arrays for a sequence of element symbols."""
    import numpy as _np
    params = [lj_params(e) for e in elements]
    rmin_half = _np.array([p["rmin_half"] for p in params], dtype=_np.float64)
    epsilon = _np.array([p["epsilon"] for p in params], dtype=_np.float64)
    return rmin_half, epsilon


def lj_energy(r, rmin, eps, *, repulsive_only=False, cutoff=None, r_floor=0.01):
    """12-6 Lennard-Jones energy ``eps * ((rmin/r)^12 - 2 (rmin/r)^6)``."""
    import numpy as _np
    r_arr = _np.asarray(r, dtype=_np.float64)
    r_safe = _np.maximum(r_arr, r_floor)
    ratio6 = _np.power(_np.asarray(rmin, dtype=_np.float64) / r_safe, 6)
    energy = _np.asarray(eps, dtype=_np.float64) * (ratio6 * ratio6 - 2.0 * ratio6)
    if repulsive_only:
        energy = _np.where(r_arr < rmin, energy, 0.0)
    if cutoff is not None:
        energy = _np.where(r_arr < cutoff, energy, 0.0)
    if _np.ndim(energy) == 0:
        return float(energy)
    return energy


def lj_score(r, rmin, epsilon):
    """Repulsive-only LJ energy for one pair (scalar)."""
    if rmin == 0 or epsilon == 0:
        return 0.0
    return lj_energy(r, rmin, epsilon, repulsive_only=True)


def boltzmann_weights(energies, temperature=298.15):
    """Normalised Boltzmann weights from energies (log-sum-exp stabilised)."""
    import numpy as _np
    KB = 0.0019872041
    kt = KB * temperature
    e = _np.asarray(energies, dtype=_np.float64)
    e_min = _np.min(e)
    w = _np.exp(-(e - e_min) / kt)
    return w / _np.sum(w)


def rotamer_cluster_weights(assignments, frame_weights, n_clusters):
    """Aggregate frame weights into cluster weights."""
    import numpy as _np
    a = _np.asarray(assignments)
    fw = _np.asarray(frame_weights, dtype=_np.float64)
    w = _np.zeros(n_clusters)
    for i in range(len(a)):
        w[a[i]] += fw[i]
    return w / _np.sum(w)


class BoundingBoxFilter:
    """AABB clash pre-filter for rotamer coordinate batches."""

    def __init__(self, pad=3.5):
        self.pad = float(pad)

    def build(self, coords):
        import numpy as _np
        c = _np.asarray(coords, dtype=_np.float64)
        mins = c.min(axis=1) - self.pad
        maxs = c.max(axis=1) + self.pad
        return _np.concatenate([mins, maxs], axis=1)

    def build_single(self, coords):
        import numpy as _np
        c = _np.asarray(coords, dtype=_np.float64)
        return _np.concatenate([
            c.min(axis=0) - self.pad,
            c.max(axis=0) + self.pad,
        ])

    @staticmethod
    def intersects(box_a, box_b):
        return (
            box_a[3] >= box_b[0] and box_a[0] <= box_b[3] and
            box_a[4] >= box_b[1] and box_a[1] <= box_b[4] and
            box_a[5] >= box_b[2] and box_a[2] <= box_b[5]
        )

    def intersects_reference(self, rot_boxes, ref_box):
        import numpy as _np
        rb = _np.asarray(rot_boxes, dtype=_np.float64)
        ref = _np.asarray(ref_box, dtype=_np.float64)
        overlap_x = (rb[:, 3] >= ref[0]) & (rb[:, 0] <= ref[3])
        overlap_y = (rb[:, 4] >= ref[1]) & (rb[:, 1] <= ref[4])
        overlap_z = (rb[:, 5] >= ref[2]) & (rb[:, 2] <= ref[5])
        return overlap_x & overlap_y & overlap_z

    def filter_frames(self, coords, reference_coords):
        import numpy as _np
        c = _np.asarray(coords, dtype=_np.float64)
        rot_boxes = self.build(c)
        ref_box = self.build_single(reference_coords)
        mask = self.intersects_reference(rot_boxes, ref_box)
        return c[mask], mask


# -- Atom-type and selector helpers (string parsing, once per structure) ----

@_lru_cache(maxsize=4096)
def _atom_type(atom_name):
    name = str(atom_name).strip().upper()
    for atom in ("S", "N", "O", "H", "C"):
        if name.startswith(atom):
            return atom
    return "C"


def _selector_matches(selector, atom_name, resname=None):
    parts = [part.strip() for part in str(selector).split(" and ")]
    atom = parts[0].strip().upper() if parts else str(selector).strip().upper()
    if atom.startswith("name "):
        atom = atom.split(None, 1)[1].upper()
    if atom != str(atom_name).upper():
        return False
    selector_resname = None
    for part in parts[1:]:
        tokens = part.split()
        if len(tokens) >= 2 and tokens[0].lower() == "resname":
            selector_resname = tokens[1].upper()
    return selector_resname is None or (resname is not None and str(resname).upper() == selector_resname)


def selector_resnames(selector):
    if selector is None:
        return set()
    if isinstance(selector, str):
        selectors = [selector]
    else:
        selectors = list(selector)
    resnames = set()
    for item in selectors:
        for part in str(item).split(" and "):
            tokens = part.strip().split()
            if len(tokens) >= 2 and tokens[0].lower() == "resname":
                resnames.add(tokens[1].upper())
    return resnames


def _selector_atom_names(selector, atom_names=None, resnames=None):
    if selector is None:
        return set()
    if isinstance(selector, str):
        selectors = [selector]
    else:
        selectors = list(selector)
    wanted_resnames = selector_resnames(selectors)
    if wanted_resnames and len(wanted_resnames) > 1:
        return set()
    names = set()
    if resnames is None:
        for item in selectors:
            parts = [part.strip() for part in str(item).split(" and ")]
            atom = parts[0].strip().upper() if parts else str(item).strip().upper()
            if atom.startswith("name "):
                atom = atom.split(None, 1)[1].upper()
            names.add(atom)
    else:
        if atom_names is None:
            raise ValueError("atom_names is required when resnames is provided")
        for atom, resname in zip(atom_names, resnames):
            if any(_selector_matches(item, atom, resname) for item in selectors):
                names.add(str(atom).upper())
    return names


def _site_mask(atom_names, residue_indices=None, site_residue=None,
               site_chain=None, chain_ids=None, mask_backbone=True):
    import numpy as _np
    mask = _np.zeros(len(atom_names), dtype=bool)
    if residue_indices is None:
        residue_indices = [-1] * len(atom_names)
    wanted = set(site_residue) if isinstance(site_residue, (list, tuple, set)) else {site_residue}
    for i, (name, residue_index) in enumerate(zip(atom_names, residue_indices)):
        frame_chain = chain_ids[i] if chain_ids is not None and i < len(chain_ids) else None
        same_chain = site_chain is None or frame_chain is None or str(frame_chain).upper() == str(site_chain).upper()
        if site_residue is not None and same_chain and residue_index in wanted:
            mask[i] = True
            continue
        if mask_backbone and name.upper() in {"CA", "C", "N", "O"}:
            mask[i] = True
    return mask


def _hydrogen_mask(atom_names):
    import numpy as _np
    return _np.array([_atom_type(name) == "H" for name in atom_names], dtype=bool)


def _protein_charge_mask(atom_names, resnames):
    import numpy as _np
    charges = _np.zeros(len(atom_names), dtype=float)
    for i, (atom, resname) in enumerate(zip(atom_names, resnames)):
        name = atom.upper()
        res = resname.upper()
        if (name == "CZ" and res == "ARG") or (name == "NZ" and res == "LYS"):
            charges[i] = 1.0
        elif (name == "CG" and res == "ASP") or (name == "CD" and res == "GLU"):
            charges[i] = -1.0
        elif name in {"ND1", "NE2"} and res == "HIS":
            charges[i] = 0.25
    return charges


def _rotamer_charge_mask(atom_names, positive=None, negative=None, resnames=None):
    import numpy as _np
    positive = [] if positive is None else list(positive)
    negative = [] if negative is None else list(negative)
    if len(selector_resnames(positive)) > 1:
        positive = []
    if len(selector_resnames(negative)) > 1:
        negative = []
    _selector_atom_names(positive, atom_names=atom_names, resnames=resnames)
    _selector_atom_names(negative, atom_names=atom_names, resnames=resnames)
    charges = _np.zeros(len(atom_names), dtype=float)
    for i, name in enumerate(atom_names):
        resname = None if resnames is None or i >= len(resnames) else resnames[i]
        is_negative = any(_selector_matches(selector, name, resname) for selector in negative)
        is_positive = any(_selector_matches(selector, name, resname) for selector in positive)
        if is_negative:
            charges[i] = -1.0
        elif is_positive:
            charges[i] = 0.5
    return charges


def _scaled_parameters(atom_names, sigma_scaling, epsilon_scaling):
    import numpy as _np
    rmin2, eps = lj_parameter_arrays([_atom_type(name) for name in atom_names])
    return rmin2 * sigma_scaling, eps * epsilon_scaling


# -- The end-to-end rotamer score -------------------------------------------

_GAS_CONSTANT = 1.9858775e-3


@_dataclass
class RotamerScoreResultPy:
    """Result of rotamer-protein scoring."""
    weights: object
    partition: float
    energies: object


def compute_rotamer_score(
    rotamer_coords, protein_coords, protein_atom_names, protein_resnames,
    rotamer_atom_names, rotamer_metadata=None, rotamer_resnames=None,
    protein_residue_indices=None, protein_chain_ids=None,
    site_residue=None, site_chain=None, rotamer_weights=None,
    temperature=300.0, ignore_h=True, electrostatic=False,
    potential="lj", sigma_scaling=0.5, epsilon_scaling=1.0,
):
    """Compute Boltzmann weights for a rotamer ensemble."""
    import numpy as _np

    if potential not in ("lj", "gauss"):
        raise ValueError(f"Unknown potential {potential!r}")

    metadata = rotamer_metadata or {}
    protein_coords = _np.asarray(protein_coords, dtype=_np.float64)
    rotamer_coords = _np.asarray(rotamer_coords, dtype=_np.float64)

    protein_h = _hydrogen_mask(protein_atom_names)
    protein_mask = ~(protein_h if ignore_h else _np.zeros(len(protein_atom_names), dtype=bool))
    protein_mask &= ~_site_mask(
        protein_atom_names, protein_residue_indices,
        site_residue=site_residue, site_chain=site_chain,
        chain_ids=protein_chain_ids, mask_backbone=False)

    rotamer_backbone = _site_mask(rotamer_atom_names)
    rotamer_h = _hydrogen_mask(rotamer_atom_names)
    rotamer_mask = ~(rotamer_backbone | (rotamer_h if ignore_h else _np.zeros(len(rotamer_atom_names), dtype=bool)))

    protein_idx = _np.flatnonzero(protein_mask)
    rotamer_idx = _np.flatnonzero(rotamer_mask)
    if protein_idx.size == 0 or rotamer_idx.size == 0:
        weights = _np.ones(rotamer_coords.shape[0], dtype=_np.float64)
        weights /= weights.sum()
        return RotamerScoreResultPy(
            weights=weights, partition=1.0,
            energies=_np.zeros(rotamer_coords.shape[0]))

    protein_types = [_atom_type(name) for name in _np.asarray(protein_atom_names)[protein_idx]]
    rotamer_types = [_atom_type(name) for name in _np.asarray(rotamer_atom_names)[rotamer_idx]]
    protein_rmin2, protein_eps = _scaled_parameters(protein_types, sigma_scaling, epsilon_scaling)
    rotamer_rmin2, rotamer_eps = _scaled_parameters(rotamer_types, sigma_scaling, epsilon_scaling)

    eps_ij = _np.sqrt(_np.multiply.outer(rotamer_eps, protein_eps))
    rmin_ij = _np.add.outer(rotamer_rmin2, protein_rmin2)

    protein_positions = protein_coords[protein_idx]
    q_rotamer = _np.zeros(rotamer_idx.size, dtype=float)
    q_protein = _np.zeros(protein_idx.size, dtype=float)

    if electrostatic:
        q_rotamer = _rotamer_charge_mask(
            [rotamer_atom_names[i] for i in rotamer_idx],
            metadata.get("positive"), metadata.get("negative"),
            None if rotamer_resnames is None else [rotamer_resnames[i] for i in rotamer_idx])
        q_protein = _protein_charge_mask(
            [protein_atom_names[i] for i in protein_idx],
            [protein_resnames[i] for i in protein_idx])

    selected = _np.ascontiguousarray(
        _np.asarray(rotamer_coords, dtype=_np.float64)[:, rotamer_idx, :])
    energies = _np.asarray(IMP.bff.rotamer_interaction_energies(
        selected.ravel(),
        _np.ascontiguousarray(protein_positions, dtype=_np.float64).ravel(),
        _np.ascontiguousarray(rmin_ij, dtype=_np.float64).ravel(),
        _np.ascontiguousarray(eps_ij, dtype=_np.float64).ravel(),
        _np.ascontiguousarray(q_rotamer, dtype=_np.float64) if electrostatic else _np.empty(0),
        _np.ascontiguousarray(q_protein, dtype=_np.float64) if electrostatic else _np.empty(0),
        int(rotamer_coords.shape[0]), int(rotamer_idx.size), int(protein_idx.size),
        IMP.bff.ROTAMER_POTENTIAL_GAUSS if potential == "gauss"
        else IMP.bff.ROTAMER_POTENTIAL_LJ,
    ), dtype=_np.float64).reshape(-1, 2)
    pot_energy = _np.ascontiguousarray(energies[:, 0])
    dh_energy = _np.ascontiguousarray(energies[:, 1])

    boltzmann = _np.exp(-pot_energy / (_GAS_CONSTANT * temperature) - dh_energy)
    library_weights = _np.asarray(
        rotamer_weights if rotamer_weights is not None
        else metadata.get("weights", _np.ones(rotamer_coords.shape[0])),
        dtype=_np.float64)
    if library_weights.shape != boltzmann.shape:
        library_weights = _np.ones_like(boltzmann)
    boltzmann *= library_weights
    boltzmann = _np.nan_to_num(boltzmann, nan=0.0, posinf=0.0, neginf=0.0)
    partition = float(_np.sum(boltzmann))
    if partition <= 0.0:
        weights = _np.ones(rotamer_coords.shape[0], dtype=_np.float64) / rotamer_coords.shape[0]
        return RotamerScoreResultPy(weights=weights, partition=0.0,
                                    energies=pot_energy + dh_energy)
    weights = boltzmann / partition
    return RotamerScoreResultPy(weights=weights, partition=partition,
                                energies=pot_energy + dh_energy)


# -- Mean-field weights (Eq. 39-42 from IMP RotamerCalculator) ---------------

def _lj_energy_pairs(coords_a, coords_b, rmin, eps, r_cutoff=12.0):
    import numpy as _np
    diff = coords_a[:, None, :] - coords_b[None, :, :]
    r = _np.linalg.norm(diff, axis=-1)
    r_flat = r.ravel()
    if rmin is None or eps is None:
        return 0.0
    lj = lj_energy(r_flat, rmin, eps, repulsive_only=True, cutoff=r_cutoff)
    return float(lj.sum())


def _build_cross_lj_params(dye_elements, ref_elements):
    import numpy as _np
    n_d = len(dye_elements)
    n_r = len(ref_elements)
    rmin_arr = _np.zeros(n_d * n_r, dtype=_np.float64)
    eps_arr = _np.zeros(n_d * n_r, dtype=_np.float64)
    k = 0
    for ea in dye_elements:
        for eb in ref_elements:
            r, e = lj_cross(ea, eb)
            rmin_arr[k] = r
            eps_arr[k] = e
            k += 1
    return rmin_arr, eps_arr


def _aabb_overlap(dye_box, ref_box):
    return (
        dye_box[3] >= ref_box[0] and dye_box[0] <= ref_box[3] and
        dye_box[4] >= ref_box[1] and dye_box[1] <= ref_box[4] and
        dye_box[5] >= ref_box[2] and dye_box[2] <= ref_box[5]
    )


def _pair_energy_matrix(coords_a, coords_b, elements_a, elements_b,
                        r_cutoff, aabb_pad):
    import numpy as _np
    a = _np.asarray(coords_a, dtype=_np.float64)
    b = _np.asarray(coords_b, dtype=_np.float64)
    if b.ndim == 2:
        b = b[None, :, :]
    rmin, eps = _build_cross_lj_params(list(elements_a), list(elements_b))
    flat = IMP.bff.rotamer_pair_energy_matrix(
        _np.ascontiguousarray(a).ravel(), _np.ascontiguousarray(b).ravel(),
        rmin, eps,
        int(a.shape[0]), int(a.shape[1]), int(b.shape[0]), int(b.shape[1]),
        float(r_cutoff), float(aabb_pad))
    return _np.asarray(flat, dtype=_np.float64).reshape(a.shape[0], b.shape[0])


def rotamer_mean_field_weights(
    rotamer_coords, initial_weights, protein_coords,
    dye_elements, protein_elements, K=1.0, n_iter=10,
    aabb_pad=3.5, r_cutoff=12.0,
):
    """Iterative mean-field update of dye rotamer cluster weights."""
    import numpy as _np
    n_clusters = rotamer_coords.shape[0]
    q = initial_weights.copy().astype(_np.float64)
    E_bb = _pair_energy_matrix(
        rotamer_coords, protein_coords, dye_elements, protein_elements,
        r_cutoff, aabb_pad)[:, 0]
    for _iteration in range(n_iter):
        E_P = E_bb.copy()
        log_q = _np.log(_np.maximum(q, 1e-300)) - K * E_P
        log_q -= log_q.max()
        q_new = _np.exp(log_q)
        denom = q_new.sum()
        if denom > 0:
            q = q_new / denom
        else:
            q = _np.ones(n_clusters, dtype=_np.float64) / n_clusters
    return q


def rotamer_mean_field_weights_multi_dye(
    rotamer_coords_list, initial_weights_list, protein_coords,
    dye_elements_list, protein_elements, K=1.0, n_iter=10,
    aabb_pad=3.5, r_cutoff=12.0,
):
    """Multi-dye mean-field weight update including dye-dye cross-interactions."""
    import numpy as _np
    n_dyes = len(rotamer_coords_list)
    q_list = [w.copy().astype(_np.float64) for w in initial_weights_list]
    E_bb_all = [
        _pair_energy_matrix(rotamer_coords_list[d], protein_coords,
                            dye_elements_list[d], protein_elements,
                            r_cutoff, aabb_pad)[:, 0]
        for d in range(n_dyes)
    ]
    E_sc = {}
    for d1 in range(n_dyes):
        for d2 in range(d1 + 1, n_dyes):
            m = _pair_energy_matrix(
                rotamer_coords_list[d1], rotamer_coords_list[d2],
                dye_elements_list[d1], dye_elements_list[d2],
                r_cutoff, aabb_pad)
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
            q_list[d] = q_new / denom if denom > 0 else _np.ones_like(q_new) / len(q_new)
    return q_list


# -- IMP API glue: stays Python because it creates IMP objects --------------

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
        p1, p2, p3 = site_particles[an.site_a], site_particles[an.site_b], site_particles[an.site_c]
        k = float(at[an.type_id])
        restraints.append(
            IMP.core.AngleRestraint(model, IMP.core.Harmonic(float(an.theta), k), p1, p2, p3))
    for to in system.dihedrals:
        if any(x not in site_particles for x in (to.site_a, to.site_b, to.site_c, to.site_d)):
            continue
        p1, p2, p3, p4 = (site_particles[to.site_a], site_particles[to.site_b],
                          site_particles[to.site_c], site_particles[to.site_d])
        restraints.append(IMP.core.DihedralRestraint(model, torsion_cosine(tt[to.type_id]), p1, p2, p3, p4))
    for to in system.impropers:
        if any(x not in site_particles for x in (to.site_a, to.site_b, to.site_c, to.site_d)):
            continue
        p1, p2, p3, p4 = (site_particles[to.site_a], site_particles[to.site_b],
                          site_particles[to.site_c], site_particles[to.site_d])
        t = it[to.type_id]
        k = t["k"] if isinstance(t, dict) else t.k
        theta0 = IMP.core.get_dihedral(
            IMP.core.XYZ(p1), IMP.core.XYZ(p2), IMP.core.XYZ(p3), IMP.core.XYZ(p4))
        fun = IMP.core.Harmonic(theta0, float(k))
        restraints.append(IMP.core.DihedralRestraint(model, fun, p1, p2, p3, p4))
    excluded = {frozenset(p) for p in system.exclusions()}
    for sa, sb, rmin, eps in compute_lj_pair_sites(system, excluded=excluded):
        if sa not in site_particles or sb not in site_particles:
            continue
        p1, p2 = site_particles[sa], site_particles[sb]
        k = float(eps)
        if k <= 0:
            continue
        restraints.append(
            IMP.core.DistanceRestraint(model, IMP.core.HarmonicLowerBound(float(rmin), k), p1, p2))
    return restraints


def compute_lj_pair_sites(system, excluded=None):
    """Compute all non-excluded site pairs with LJ cross parameters."""
    system = as_forcefield_system(system)
    if excluded is None:
        excluded = {frozenset(p) for p in system.exclusions()}
    elem_map = site_element_map(system)
    sites = system.sites
    pairs = []
    for i in range(len(sites)):
        for j in range(i + 1, len(sites)):
            sa = sites[i].id
            sb = sites[j].id
            pair = frozenset({sa, sb})
            if pair in excluded:
                continue
            ea = elem_map.get(sa, "C")
            eb = elem_map.get(sb, "C")
            rmin, eps = lj_cross(ea, eb)
            pairs.append((sa, sb, rmin, eps))
    return pairs


def site_element_map(system):
    """Extract {site_id: element} from system sites."""
    system = as_forcefield_system(system)
    result = {}
    for s in system.sites:
        sid = s.id
        aname = s.atom_name
        m = _re.match(r"([A-Za-z])", aname)
        elem = m.group(1).upper() if m else "C"
        result[sid] = elem
    return result


def build_lj_type_table(elements):
    """Build _ff_lj_type entries for a set of elements."""
    lj_types = {}
    for elem in sorted(elements):
        params = CHARMM36_LJ.get(elem, CHARMM36_LJ["C"])
        lj_types[f"LJ_{elem}"] = {
            "element": elem,
            "rmin_half": params["rmin_half"],
            "epsilon": params["epsilon"],
        }
    return lj_types


def dye_internal_system(atoms_dict, bonds):
    """Minimal force-field system for a lone dye from parse_dye_mol2 output."""
    from IMP.bff import build_angles, build_dihedrals, build_graph
    ordered = sorted(atoms_dict.values(), key=lambda x: x["serial"])
    sid = {a["serial"]: f"dye:{a['serial']}:{a['atom_name']}" for a in ordered}
    sites = [{"id": sid[a["serial"]], "atom_name": a["atom_name"]} for a in ordered]
    graph = build_graph(bonds)
    return {
        "sites": sites,
        "bonds": [(sid[a], sid[b], 0.0, None) for a, b in sorted(bonds)],
        "angles": [(sid[a], sid[b], sid[c], 0.0, None) for a, b, c in build_angles(graph)],
        "dihedrals": [(sid[a], sid[b], sid[c], sid[d], None) for a, b, c, d in build_dihedrals(graph)],
    }


class DyeInternalEnergyEvaluator:
    """Evaluates the internal LJ energy of dye conformations."""

    def __init__(self, system, excluded=None):
        system = as_forcefield_system(system)
        self.pairs = compute_lj_pair_sites(system, excluded=excluded)
        self.id_to_idx = {s.id: i for i, s in enumerate(system.sites)}
        if self.pairs:
            import numpy as _np
            self._idx_a = _np.array([self.id_to_idx[p[0]] for p in self.pairs], dtype=_np.intp)
            self._idx_b = _np.array([self.id_to_idx[p[1]] for p in self.pairs], dtype=_np.intp)
            self._rmin = _np.array([p[2] for p in self.pairs], dtype=_np.float64)
            self._eps = _np.array([p[3] for p in self.pairs], dtype=_np.float64)
        else:
            import numpy as _np
            self._idx_a = _np.array([], dtype=_np.intp)
            self._idx_b = _np.array([], dtype=_np.intp)
            self._rmin = _np.array([], dtype=_np.float64)
            self._eps = _np.array([], dtype=_np.float64)

    def evaluate(self, coords):
        import numpy as _np
        return float(self.evaluate_batch(_np.asarray(coords, dtype=_np.float64)[None, :, :])[0])

    def evaluate_batch(self, coords_batch):
        import numpy as _np
        batch = _np.asarray(coords_batch, dtype=_np.float64)
        if len(self.pairs) == 0:
            return _np.zeros(batch.shape[0])
        return _np.asarray(IMP.bff.lj_pair_energies(
            _np.ascontiguousarray(batch).ravel(),
            self._idx_a.astype(_np.int32), self._idx_b.astype(_np.int32),
            self._rmin, self._eps,
            int(batch.shape[0]), int(batch.shape[1]), int(self._rmin.size), True),
            dtype=_np.float64)

    def evaluate_batch_filtered(self, coords_batch, reference_coords, pad=3.5):
        import numpy as _np
        bbf = BoundingBoxFilter(pad=pad)
        surviving, clash_mask = bbf.filter_frames(coords_batch, reference_coords)
        energies = _np.zeros(coords_batch.shape[0])
        if surviving.shape[0] > 0:
            energies[clash_mask] = self.evaluate_batch(surviving)
        return energies, clash_mask
%}
