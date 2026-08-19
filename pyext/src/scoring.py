"""Stage 2 — is this dye configuration allowed, and how heavily does it count?

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

* :mod:`~IMP.bff.scoring` -- sterics, and the CHARMM36 table.
  Shared: the rotamer scorer and the explicit dye's force field read the same
  parameters, so they cannot drift apart.
* :mod:`~IMP.bff.scoring` -- a conformer against its site.
* :mod:`~IMP.bff.scoring` -- LJ pair scoring and exclusion lists for an
  explicit dye.
* :mod:`~IMP.bff.scoring` -- the same, as IMP restraints.
* :mod:`~IMP.bff.scoring` -- weights from energies.
* :mod:`~IMP.bff.scoring` -- weights solved self-consistently, for
  when two dyes see each other.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from typing import Any
import math

import numpy as np

from IMP.bff.photophysics import kappa2_from_dipoles
import IMP
import IMP.bff
import IMP.container
import IMP.core

__all__ = [
    'CHARMM36_LJ',
    'lj_cross',
    'lj_energy',
    'lj_parameter_arrays',
    'lj_params',
    'torsion_cosine',
]

# --------------------------------------------------------------------------
# boltzmann
# --------------------------------------------------------------------------
"""Boltzmann scoring and statistical weight calculation."""

def boltzmann_weights(energies: np.ndarray, temperature: float = 298.15) -> np.ndarray:
    """Compute normalized Boltzmann weights from energies.
    
    Args:
        energies: Array of energies in kcal/mol (or consistent units with kT)
        temperature: Temperature in Kelvin
        
    Returns:
        Normalized weights summing to 1.0.
    """
    # kB in kcal/(mol*K)
    KB = 0.0019872041
    kt = KB * temperature
    
    # Use log-sum-exp trick for numerical stability
    e_min = np.min(energies)
    shifted_energies = energies - e_min
    unnormalized_weights = np.exp(-shifted_energies / kt)
    
    return unnormalized_weights / np.sum(unnormalized_weights)


def rotamer_cluster_weights(assignments: np.ndarray, frame_weights: np.ndarray, n_clusters: int) -> np.ndarray:
    """Aggregate frame weights into cluster weights.
    
    Args:
        assignments: Cluster index for each frame (n_frames,)
        frame_weights: Normalized weight for each frame (n_frames,)
        n_clusters: Total number of clusters
        
    Returns:
        Array of shape (n_clusters,) with normalized cluster weights.
    """
    weights = np.zeros(n_clusters)
    for i in range(len(assignments)):
        weights[assignments[i]] += frame_weights[i]
    
    return weights / np.sum(weights)


# --------------------------------------------------------------------------
# lennard_jones
# --------------------------------------------------------------------------
"""The Lennard-Jones term, and the CHARMM36 parameters it reads.

Sterics: how hard two atoms push each other apart. Nothing here is about
fluorescence, which is why it is a tool rather than a domain -- and why it sits
here rather than inside either of the two things that need it. A rotamer's score
against a protein is a clash energy
(:mod:`IMP.bff.scoring`); so is an explicit dye's
internal energy (:mod:`IMP.bff.cgdye.topology`). Those are a core
representation and a legacy package respectively, and the first must not depend
on the second to compute a steric term.

Lifted out of ``cgdye/topology/dye.py`` in the 2026-08-18 cleanup, when
``rotamer/`` moved into ``representation/`` and left that edge pointing the
wrong way.
"""

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
    """Lorentz-Berthelot cross parameters ``(rmin, eps)`` for two elements.

    ``rmin = rmin_half_i + rmin_half_j`` (Angstrom), ``eps = sqrt(eps_i eps_j)``
    (kcal/mol, positive well depth).
    """
    pi = lj_params(elem_i)
    pj = lj_params(elem_j)
    rmin = pi["rmin_half"] + pj["rmin_half"]
    eps = (pi["epsilon"] * pj["epsilon"]) ** 0.5
    return rmin, eps


def lj_parameter_arrays(elements):
    """``(rmin_half, epsilon)`` arrays for a sequence of element symbols.

    Unknown elements fall back to carbon, as ``lj_params`` does.
    """
    params = [lj_params(e) for e in elements]
    rmin_half = np.array([p["rmin_half"] for p in params], dtype=np.float64)
    epsilon = np.array([p["epsilon"] for p in params], dtype=np.float64)
    return rmin_half, epsilon


def lj_energy(r, rmin, eps, *, repulsive_only=False, cutoff=None, r_floor=0.01):
    """12-6 Lennard-Jones energy ``eps * ((rmin/r)^12 - 2 (rmin/r)^6)``.

    Vectorised over numpy arrays (``r``, ``rmin``, ``eps`` broadcast against
    each other); scalars work too. ``r`` is clamped to ``r_floor`` to keep the
    energy finite. ``repulsive_only=True`` zeroes the attractive tail
    (``r >= rmin``), ``cutoff`` zeroes pairs beyond that distance. With
    positive ``eps`` (the ``lj_cross`` convention) the well depth at ``rmin``
    is ``-eps``.
    """
    r_arr = np.asarray(r, dtype=np.float64)
    r_safe = np.maximum(r_arr, r_floor)
    ratio6 = np.power(np.asarray(rmin, dtype=np.float64) / r_safe, 6)
    energy = np.asarray(eps, dtype=np.float64) * (ratio6 * ratio6 - 2.0 * ratio6)
    if repulsive_only:
        energy = np.where(r_arr < rmin, energy, 0.0)
    if cutoff is not None:
        energy = np.where(r_arr < cutoff, energy, 0.0)
    if np.ndim(energy) == 0:
        return float(energy)
    return energy


# --------------------------------------------------------------------------
# dye_lj
# --------------------------------------------------------------------------
"""LJ pair scoring with CHARMM36 parameters and exclusion list generation.

Tasks 2.5/2.6: per-element LJ scoring, 1-2/1-3/1-4 exclusions.

Algorithmic tricks ported from IMP RotamerCalculator.cpp:
  - Trick 1: Axis-aligned bounding-box (AABB) clash pre-filter, analogous to
    ResidueRotamer::create_bounding_boxes / intersect().  Skips expensive pairwise
    LJ evaluation for rotamer frames whose atom cloud cannot possibly clash with
    a reference (protein) coordinate set.
  - Trick 3: Vectorized batch LJ evaluation using NumPy broadcasting — replaces
    the Python-level for-loop over frames with a single (n_frames, n_pairs)
    tensor operation.
"""

#!/usr/bin/env python






# ---------------------------------------------------------------------------
# Scalar helpers (kept for single-frame Metropolis acceptance in LinkerSampler)
# ---------------------------------------------------------------------------

def lj_score(r, rmin, epsilon):
    """Repulsive-only Lennard-Jones energy for one pair (scalar).

    ``epsilon * ((rmin/r)**12 - 2*(rmin/r)**6)`` for ``r < rmin``, 0 beyond
    (no attractive tail in the coarse-grained linker sampler); ``r`` is clamped
    to 0.01 A. Thin scalar front to :func:`IMP.bff.cgdye.topology.dye.lj_energy`.
    """
    if rmin == 0 or epsilon == 0:
        return 0.0
    return lj_energy(r, rmin, epsilon, repulsive_only=True)


# ---------------------------------------------------------------------------
# Exclusion / topology helpers
# ---------------------------------------------------------------------------

def compute_exclusions(system):
    """Compute 1-2, 1-3, and 1-4 excluded pairs from system topology.

    Returns a set of frozenset({site_id_a, site_id_b}) pairs.
    """
    excluded = set()

    bond_pairs = set()
    for a, b, _length, _tid in system.get("bonds", []):
        bond_pairs.add(frozenset({a, b}))
    excluded.update(bond_pairs)

    angle_ends = set()
    for a, b, c, _theta, _tid in system.get("angles", []):
        pair = frozenset({a, c})
        if pair not in excluded:
            angle_ends.add(pair)
    excluded.update(angle_ends)

    dihedral_ends = set()
    for a, b, c, d, _tid in system.get("dihedrals", []):
        pair = frozenset({a, d})
        if pair not in excluded:
            dihedral_ends.add(pair)
    excluded.update(dihedral_ends)

    for a, b, c, d, _tid in system.get("impropers", []):
        for x, y in [(a, c), (a, d), (b, d)]:
            excluded.add(frozenset({x, y}))

    return excluded


def dye_internal_system(atoms_dict, bonds):
    """Minimal force-field ``system`` for a lone dye from ``parse_dye_mol2`` output.

    Sites are the MOL2 atoms in serial order (ids ``dye:<atom_name>``, the
    order ``LinkerSampler.get_coords`` uses); bonds, angles and dihedrals are
    derived from the MOL2 connectivity so that :func:`compute_exclusions`
    removes the 1-2, 1-3 and 1-4 pairs from the LJ pair list. Without them the
    linker sampler scored bonded neighbours with the repulsive LJ, which is
    unphysical and biases torsion/angle sampling. Site ids are
    ``dye:<serial>:<atom_name>`` (unique; see below).
    """
    # Deferred deliberately. An exclusion list is derived from connectivity, so
    # scoring genuinely needs the molecular graph -- but importing `scoring`
    # must not drag in the explicit-dye package, which is off the domain
    # layout. Keeping the edge inside the function keeps load order clean.
    from IMP.bff.cgdye.topology import build_angles, build_dihedrals, build_graph

    ordered = sorted(atoms_dict.values(), key=lambda x: x["serial"])
    # Site ids carry the serial: MOL2 atom names are not unique across the
    # dye + linker residues (alexa488_r48 has 83 atoms, 68 distinct names), and
    # name-only ids made same-named atoms one site, i.e. self-pairs at r = 0
    # with 1e31 energies in every frame.
    sid = {a["serial"]: f"dye:{a['serial']}:{a['atom_name']}" for a in ordered}
    sites = [{"id": sid[a["serial"]], "atom_name": a["atom_name"]} for a in ordered]
    graph = build_graph(bonds)
    return {
        "sites": sites,
        "bonds": [(sid[a], sid[b], 0.0, None) for a, b in sorted(bonds)],
        "angles": [(sid[a], sid[b], sid[c], 0.0, None) for a, b, c in build_angles(graph)],
        "dihedrals": [(sid[a], sid[b], sid[c], sid[d], None) for a, b, c, d in build_dihedrals(graph)],
    }


def build_lj_type_table(elements):
    """Build _ff_lj_type entries for a set of elements.

    Returns dict {type_id: {element, rmin_half, epsilon}}.
    """

    lj_types = {}
    for elem in sorted(elements):
        params = CHARMM36_LJ.get(elem, CHARMM36_LJ["C"])
        lj_types[f"LJ_{elem}"] = {
            "element": elem,
            "rmin_half": params["rmin_half"],
            "epsilon": params["epsilon"],
        }
    return lj_types


def site_element_map(system):
    """Extract {site_id: element} from system sites.

    Derives element from atom_name prefix (first letter capitalized).
    """
    import re

    result = {}
    for s in system.get("sites", []):
        sid = s["id"]
        aname = s.get("atom_name", "")
        m = re.match(r"([A-Za-z])", aname)
        elem = m.group(1).upper() if m else "C"
        result[sid] = elem
    return result


def lj_cross_params(elem_i, elem_j):
    """Lorentz-Berthelot combining rules for two elements.

    Returns (rmin_ij, epsilon_ij).
    """
    return lj_cross(elem_i, elem_j)


def compute_lj_pair_sites(system, excluded=None):
    """Compute all non-excluded site pairs with LJ cross parameters.

    Returns list of (site_id_a, site_id_b, rmin_ij, epsilon_ij).
    """
    if excluded is None:
        excluded = compute_exclusions(system)
    elem_map = site_element_map(system)
    sites = system.get("sites", [])
    pairs = []
    for i in range(len(sites)):
        for j in range(i + 1, len(sites)):
            sa = sites[i]["id"]
            sb = sites[j]["id"]
            pair = frozenset({sa, sb})
            if pair in excluded:
                continue
            ea = elem_map.get(sa, "C")
            eb = elem_map.get(sb, "C")
            rmin, eps = lj_cross_params(ea, eb)
            pairs.append((sa, sb, rmin, eps))
    return pairs


# ---------------------------------------------------------------------------
# Trick 1 — Axis-Aligned Bounding Box (AABB) pre-filter
# Mirrors ResidueRotamer::create_bounding_boxes / intersect() from IMP.
# Boxes are padded by half the maximum LJ rmin (≈ 3.5 Å default) so that
# any overlapping box pair *could* have a clashing atom pair.
# ---------------------------------------------------------------------------

class BoundingBoxFilter:
    """AABB clash pre-filter for rotamer coordinate batches.

    Algorithmic source: IMP RotamerCalculator.cpp,
    ResidueRotamer::create_bounding_boxes() and intersect().

    Usage::

        bbf = BoundingBoxFilter(pad=3.5)
        # boxes for all rotamer frames: shape (n_frames, 6)  [xmin,ymin,zmin,xmax,ymax,zmax]
        rot_boxes = bbf.build(rotamer_coords)       # (n_frames, n_atoms, 3)
        prot_box  = bbf.build_single(protein_coords)  # (n_prot_atoms, 3)
        mask = bbf.intersects_reference(rot_boxes, prot_box)  # bool (n_frames,)
        surviving = rotamer_coords[mask]
    """

    def __init__(self, pad: float = 3.5):
        """
        Args:
            pad: Padding in Angstroms added to each face of the AABB.
                 IMP uses 3.5 Å (≈ C–C LJ rmin/2).
        """
        self.pad = float(pad)

    def build(self, coords: np.ndarray) -> np.ndarray:
        """Build AABBs for a batch of coordinate frames.

        Args:
            coords: shape (n_frames, n_atoms, 3)

        Returns:
            boxes: shape (n_frames, 6) — columns are
                   [xmin, ymin, zmin, xmax, ymax, zmax] with padding applied.
        """
        mins = coords.min(axis=1) - self.pad   # (n_frames, 3)
        maxs = coords.max(axis=1) + self.pad   # (n_frames, 3)
        return np.concatenate([mins, maxs], axis=1)  # (n_frames, 6)

    def build_single(self, coords: np.ndarray) -> np.ndarray:
        """Build a single AABB for a static coordinate set.

        Args:
            coords: shape (n_atoms, 3)

        Returns:
            box: shape (6,) — [xmin, ymin, zmin, xmax, ymax, zmax]
        """
        return np.concatenate([
            coords.min(axis=0) - self.pad,
            coords.max(axis=0) + self.pad,
        ])

    @staticmethod
    def intersects(box_a: np.ndarray, box_b: np.ndarray) -> bool:
        """Test AABB intersection between two single boxes (shape (6,)).

        Direct port of ResidueRotamer::intersect() from IMP.
        """
        return (
            box_a[3] >= box_b[0] and box_a[0] <= box_b[3] and
            box_a[4] >= box_b[1] and box_a[1] <= box_b[4] and
            box_a[5] >= box_b[2] and box_a[2] <= box_b[5]
        )

    def intersects_reference(
        self,
        rot_boxes: np.ndarray,
        ref_box: np.ndarray,
    ) -> np.ndarray:
        """Vectorized AABB intersection test of many rotamer frames vs one reference.

        Args:
            rot_boxes: shape (n_frames, 6)
            ref_box:   shape (6,)

        Returns:
            mask: bool array (n_frames,) — True if frame *could* clash.
        """
        # rot_boxes columns: [xmin0, ymin1, zmin2, xmax3, ymax4, zmax5]
        # ref_box:           [xmin0, ymin1, zmin2, xmax3, ymax4, zmax5]
        overlap_x = (rot_boxes[:, 3] >= ref_box[0]) & (rot_boxes[:, 0] <= ref_box[3])
        overlap_y = (rot_boxes[:, 4] >= ref_box[1]) & (rot_boxes[:, 1] <= ref_box[4])
        overlap_z = (rot_boxes[:, 5] >= ref_box[2]) & (rot_boxes[:, 2] <= ref_box[5])
        return overlap_x & overlap_y & overlap_z

    def filter_frames(
        self,
        coords: np.ndarray,
        reference_coords: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Filter rotamer frames that cannot clash with a reference coordinate set.

        Args:
            coords: shape (n_frames, n_atoms, 3) — rotamer frames
            reference_coords: shape (n_ref_atoms, 3) — static reference (e.g. protein)

        Returns:
            surviving_coords: subset of coords where AABB overlaps reference
            mask:             bool (n_frames,) survivor mask
        """
        rot_boxes = self.build(coords)
        ref_box = self.build_single(reference_coords)
        mask = self.intersects_reference(rot_boxes, ref_box)
        return coords[mask], mask


# ---------------------------------------------------------------------------
# Trick 3 — Vectorized batch LJ scoring
# Replaces the Python for-loop in the old evaluate_batch with NumPy
# (n_frames, n_pairs) broadcasting — O(1) Python loops.
# ---------------------------------------------------------------------------

class DyeInternalEnergyEvaluator:
    """Evaluates the internal LJ energy of dye conformations.

    evaluate()       — scalar, used inside Metropolis loop (unchanged API).
    evaluate_batch() — vectorized NumPy over all frames simultaneously
                       (Trick 3 from IMP RotamerCalculator inner-loop pattern).
    evaluate_batch_filtered() — AABB pre-filter (Trick 1) then vectorized LJ.
    """

    def __init__(self, system, excluded=None):
        self.pairs = compute_lj_pair_sites(system, excluded=excluded)
        # Map site_id to index in the coordinate array
        self.id_to_idx = {s["id"]: i for i, s in enumerate(system["sites"])}

        # Pre-extract index arrays and parameter arrays for vectorized eval
        if self.pairs:
            self._idx_a = np.array([self.id_to_idx[p[0]] for p in self.pairs], dtype=np.intp)
            self._idx_b = np.array([self.id_to_idx[p[1]] for p in self.pairs], dtype=np.intp)
            self._rmin  = np.array([p[2] for p in self.pairs], dtype=np.float64)
            self._eps   = np.array([p[3] for p in self.pairs], dtype=np.float64)
        else:
            self._idx_a = np.array([], dtype=np.intp)
            self._idx_b = np.array([], dtype=np.intp)
            self._rmin  = np.array([], dtype=np.float64)
            self._eps   = np.array([], dtype=np.float64)

    # ------------------------------------------------------------------
    # Scalar evaluate — kept for Metropolis single-step acceptance
    # ------------------------------------------------------------------
    def evaluate(self, coords):
        """Evaluate LJ energy for a single frame.

        Args:
            coords: Array of shape (n_atoms, 3)
        """
        # One frame through the vectorised kernel: the same arithmetic as
        # evaluate_batch (regression-tested equal to 1e-10), no Python pair loop.
        return float(self.evaluate_batch(np.asarray(coords, dtype=np.float64)[None, :, :])[0])

    # ------------------------------------------------------------------
    # Trick 3 — vectorized batch LJ evaluation
    # ------------------------------------------------------------------
    def evaluate_batch(self, coords_batch: np.ndarray) -> np.ndarray:
        """Vectorized LJ evaluation for all frames simultaneously.

        Trick 3: instead of a Python for-loop over frames, we broadcast
        the full (n_frames, n_pairs) distance tensor in one NumPy call.

        Args:
            coords_batch: shape (n_frames, n_atoms, 3)

        Returns:
            energies: shape (n_frames,)
        """
        batch = np.asarray(coords_batch, dtype=np.float64)
        if len(self.pairs) == 0:
            return np.zeros(batch.shape[0])
        n_frames, n_atoms = batch.shape[0], batch.shape[1]
        return np.asarray(IMP.bff.lj_pair_energies(
            np.ascontiguousarray(batch).ravel(),
            self._idx_a.astype(np.int32), self._idx_b.astype(np.int32),
            self._rmin, self._eps,
            int(n_frames), int(n_atoms), int(self._rmin.size), True),
            dtype=np.float64)

    # ------------------------------------------------------------------
    # Combined Trick 1 + 3: filter then score
    # ------------------------------------------------------------------
    def evaluate_batch_filtered(
        self,
        coords_batch: np.ndarray,
        reference_coords: np.ndarray,
        pad: float = 3.5,
    ) -> tuple[np.ndarray, np.ndarray]:
        """AABB pre-filter (Trick 1) then vectorized LJ (Trick 3).

        Frames whose bounding box does not overlap the reference (protein)
        bounding box are assigned energy=0 without any pairwise computation.

        Args:
            coords_batch:     shape (n_frames, n_atoms, 3) — dye rotamers
            reference_coords: shape (n_ref_atoms, 3)       — protein atoms
            pad:              AABB padding in Å (default 3.5, same as IMP)

        Returns:
            energies: shape (n_frames,) — 0 for non-overlapping frames
            clash_mask: bool (n_frames,) — True where AABB overlapped
        """
        bbf = BoundingBoxFilter(pad=pad)
        surviving, clash_mask = bbf.filter_frames(coords_batch, reference_coords)

        energies = np.zeros(coords_batch.shape[0])
        if surviving.shape[0] > 0:
            energies[clash_mask] = self.evaluate_batch(surviving)

        return energies, clash_mask


# --------------------------------------------------------------------------
# torsion
# --------------------------------------------------------------------------
"""The torsion potential, and the sign convention that made it wrong once.

CHARMM writes a dihedral term as ``k(1 + cos(n*phi - delta))``; IMP's
``IMP.core.Cosine`` writes it as ``k(1 + cos(n*phi + delta))``. The two differ
by the sign of the phase, and getting it backwards makes a planar torsion
*minimal* at 90 degrees instead of 0 -- which is what PRD-108's Langevin sampler
found by way of a thermodynamic pin, not by inspection.

Moved out of ``cgdye/topology/dye.py`` in the four-stage restructure: a torsion
term is an energy, so it is scoring, and leaving it in the topology builder made
``scoring`` depend on a representation.
"""

def torsion_cosine(torsion_type):
    """The ``IMP.core.Cosine`` for a torsion type stored in the CHARMM convention.

    cgdye's ``torsion_types`` (``k``, ``periodicity`` n, ``phase_rad`` δ) mean
    the CHARMM/AMBER form ``V = k (1 + cos(n φ − δ))``: ``T_PI`` (n = 2,
    δ = π) is minimal at the planar 0°/180°, ``T_LINK`` (n = 3, δ = 0) at
    the staggered ±60°/180°. ``IMP.core.Cosine(k, n, δ')`` scores
    ``k (1 − cos(n φ − δ'))`` — the *opposite* sign — so δ' = δ + π. Passing
    δ straight through (as the code did) put every conjugated torsion's
    minimum at 90° and every linker torsion's at the eclipsed 0°/±120°.
    """
    return IMP.core.Cosine(float(torsion_type["k"]), int(torsion_type["periodicity"]),
                           float(torsion_type["phase_rad"]) + math.pi)


# --------------------------------------------------------------------------
# dye_restraints
# --------------------------------------------------------------------------
"""IMP restraint builder with element-aware dye nonbonded terms."""

def build_dye_restraints(model, system, site_particles):
    """Build bonded + element-aware nonbonded restraints.

    Nonbonded uses per-pair HarmonicLowerBound at LJ rmin.
    """
    restraints = []

    bt = system.get("bond_types", {})
    at = system.get("angle_types", {})
    tt = system.get("torsion_types", {})
    it = system.get("improper_types", {})

    for a, b, length, tid in system.get("bonds", []):
        if a not in site_particles or b not in site_particles:
            continue
        p1, p2 = site_particles[a], site_particles[b]
        d0 = float(length)
        k = float(bt[tid]["k"])
        restraints.append(
            IMP.core.DistanceRestraint(model, IMP.core.Harmonic(d0, k), p1, p2)
        )

    for a, b, c, theta, tid in system.get("angles", []):
        if (
            a not in site_particles
            or b not in site_particles
            or c not in site_particles
        ):
            continue
        p1, p2, p3 = site_particles[a], site_particles[b], site_particles[c]
        k = float(at[tid]["k"])
        restraints.append(
            IMP.core.AngleRestraint(
                model, IMP.core.Harmonic(float(theta), k), p1, p2, p3
            )
        )

    for a, b, c, d, tid in system.get("dihedrals", []):
        if any(x not in site_particles for x in (a, b, c, d)):
            continue
        p1, p2, p3, p4 = (
            site_particles[a],
            site_particles[b],
            site_particles[c],
            site_particles[d],
        )
        # CHARMM-convention type -> IMP.core.Cosine (sign flip, see torsion_cosine)
        restraints.append(IMP.core.DihedralRestraint(model, torsion_cosine(tt[tid]), p1, p2, p3, p4))

    for a, b, c, d, tid in system.get("impropers", []):
        if any(x not in site_particles for x in (a, b, c, d)):
            continue
        p1, p2, p3, p4 = (
            site_particles[a],
            site_particles[b],
            site_particles[c],
            site_particles[d],
        )
        t = it[tid]
        theta0 = IMP.core.get_dihedral(
            IMP.core.XYZ(p1), IMP.core.XYZ(p2), IMP.core.XYZ(p3), IMP.core.XYZ(p4)
        )
        fun = IMP.core.Harmonic(theta0, float(t["k"]))
        restraints.append(IMP.core.DihedralRestraint(model, fun, p1, p2, p3, p4))

    excluded = compute_exclusions(system)
    for sa, sb, rmin, eps in compute_lj_pair_sites(system, excluded=excluded):
        if sa not in site_particles or sb not in site_particles:
            continue
        p1, p2 = site_particles[sa], site_particles[sb]
        k = float(eps)
        if k <= 0:
            continue
        restraints.append(
            IMP.core.DistanceRestraint(
                model,
                IMP.core.HarmonicLowerBound(float(rmin), k),
                p1,
                p2,
            )
        )

    return restraints


# --------------------------------------------------------------------------
# mean_field
# --------------------------------------------------------------------------
"""Mean-field rotamer probability update for dye rotamer libraries.

Algorithmic source: IMP RotamerCalculator.cpp, RotamerCalculator::transform(),
equations (39)–(42) from Bauer et al. / Dunbrack mean-field algorithm.

The IMP C++ implementation operates on amino-acid side chains using a hardcoded
topology table (chi-axis atoms: N-CA-CB-CG etc.).  For the dye linker — a
non-canonical flexible chain — we cannot reuse that table, so we reimplement
the *algorithm* (equations 39–42) in Python/NumPy against our own coordinate
arrays and LJ scoring backend.

What this module does
---------------------
Given:
  - n_clusters rotamer cluster coordinate sets for the dye  (n_clusters, n_atoms, 3)
  - initial Boltzmann weights  q[i]  for each cluster  (n_clusters,)
  - static protein coordinates                          (n_prot_atoms, 3)
  - LJ pair parameters between dye atoms and protein atoms

It iterates (Eq. 41/42):

  E_P[i] = E_bb[i] + Σ_{j≠i} q[j] * E_SC[ij]          (Eq. 41)
  q[i]   ∝ P_boltz[i] * exp(-K * E_P[i])               (Eq. 42)

until convergence or max_iter, returning updated cluster weights.

For a *single* dye (no dye–dye cross-terms) this reduces to a pure
protein–dye interaction reweighting:

  E_P[i] = E_bb[i]  (Eq. 39 only)
  q[i]   ∝ P_boltz[i] * exp(-K * E_bb[i])

For *multi-dye* systems the E_SC cross term between the two dyes' rotamers
is computed as the pairwise LJ interaction between rotamer i of dye A and
rotamer j of dye B, then averaged with the current q of the partner dye.
"""

# ---------------------------------------------------------------------------
# LJ helpers (mirror IMP's PairScore, but in NumPy)
# ---------------------------------------------------------------------------

def _lj_energy_pairs(
    coords_a: np.ndarray,   # (n_a, 3)  dye atoms in one rotamer
    coords_b: np.ndarray,   # (n_b, 3)  reference atoms (protein or other dye)
    rmin: np.ndarray,       # (n_a * n_b,) or None for cut-off only
    eps: np.ndarray,        # (n_a * n_b,)
    r_cutoff: float = 12.0,
) -> float:
    """Compute total LJ interaction energy between two atom sets.

    Vectorized over all n_a × n_b pairs.
    """
    # (n_a, n_b, 3)
    diff = coords_a[:, None, :] - coords_b[None, :, :]
    r = np.linalg.norm(diff, axis=-1)           # (n_a, n_b)
    r_flat = r.ravel()

    # rmin / eps broadcast: if None assume uniform zeros → no interaction
    if rmin is None or eps is None:
        return 0.0

    # Repulsive-only LJ within the cutoff, through the shared kernel
    lj = lj_energy(r_flat, rmin, eps, repulsive_only=True, cutoff=r_cutoff)
    return float(lj.sum())


def _build_cross_lj_params(
    dye_elements: list[str],
    ref_elements: list[str],
) -> tuple[np.ndarray, np.ndarray]:
    """Build flattened LJ rmin / eps arrays for all dye × ref atom pairs.

    Returns (rmin, eps) arrays of shape (n_dye * n_ref,).
    """
    n_d = len(dye_elements)
    n_r = len(ref_elements)
    rmin_arr = np.zeros(n_d * n_r, dtype=np.float64)
    eps_arr  = np.zeros(n_d * n_r, dtype=np.float64)

    k = 0
    for ea in dye_elements:
        for eb in ref_elements:
            r, e = lj_cross(ea, eb)
            rmin_arr[k] = r
            eps_arr[k]  = e
            k += 1
    return rmin_arr, eps_arr


# ---------------------------------------------------------------------------
# Trick 1 (AABB) used here too — import from scoring
# ---------------------------------------------------------------------------

def _aabb_overlap(
    dye_box: np.ndarray,   # (6,)  dye rotamer box
    ref_box: np.ndarray,   # (6,)  protein box
) -> bool:
    """Return True if the two AABBs overlap.  Port of IMP::intersect()."""
    return (
        dye_box[3] >= ref_box[0] and dye_box[0] <= ref_box[3] and
        dye_box[4] >= ref_box[1] and dye_box[1] <= ref_box[4] and
        dye_box[5] >= ref_box[2] and dye_box[2] <= ref_box[5]
    )



def _pair_energy_matrix(coords_a, coords_b, elements_a, elements_b,
                        r_cutoff, aabb_pad):
    """Interaction energy of every conformer of *a* against every one of *b*.

    ``coords_a`` is ``(n_a_conf, n_a_atoms, 3)``; ``coords_b`` likewise, or
    ``(n_b_atoms, 3)`` for a single static set such as a protein. Returns
    ``(n_a_conf, n_b_conf)``.

    The parameters are built for **this** atom ordering. That matters: the
    version this replaces stored one matrix per unordered dye pair and reused it
    for both orders, so with ``d1 > d2`` it multiplied a ``(n_1, n_2)`` distance
    matrix by ``(n_2, n_1)`` parameters. The flattened lengths match, so numpy
    broadcast it silently and paired every atom with the wrong partner's
    parameters. Measured on a 3-atom against a 2-atom dye: 8366.65 where the
    right answer is 8448.92.
    """
    a = np.asarray(coords_a, dtype=np.float64)
    b = np.asarray(coords_b, dtype=np.float64)
    if b.ndim == 2:
        b = b[None, :, :]
    rmin, eps = _build_cross_lj_params(list(elements_a), list(elements_b))
    flat = IMP.bff.rotamer_pair_energy_matrix(
        np.ascontiguousarray(a).ravel(), np.ascontiguousarray(b).ravel(),
        rmin, eps,
        int(a.shape[0]), int(a.shape[1]), int(b.shape[0]), int(b.shape[1]),
        float(r_cutoff), float(aabb_pad))
    return np.asarray(flat, dtype=np.float64).reshape(a.shape[0], b.shape[0])


# ---------------------------------------------------------------------------
# Main public function — single-dye mean-field reweighting (Eq. 39 + 42)
# ---------------------------------------------------------------------------

def rotamer_mean_field_weights(
    rotamer_coords: np.ndarray,         # (n_clusters, n_dye_atoms, 3)
    initial_weights: np.ndarray,        # (n_clusters,)  Boltzmann prior
    protein_coords: np.ndarray,         # (n_prot_atoms, 3)
    dye_elements: list[str],            # element per dye atom, len = n_dye_atoms
    protein_elements: list[str],        # element per protein atom, len = n_prot_atoms
    K: float = 1.0,
    n_iter: int = 10,
    aabb_pad: float = 3.5,
    r_cutoff: float = 12.0,
) -> np.ndarray:
    """Iterative mean-field update of dye rotamer cluster weights.

    Implements equations (39) and (42) from IMP RotamerCalculator.cpp for a
    *single* dye (no dye–dye cross-terms).  The protein acts as the fixed
    backbone environment (E_bb), and q[i] is updated each iteration.

    For a single dye E_P[i] = E_bb[i] (the Σ q_lm E_SC term vanishes because
    there is only one dye).  The loop still converges in one iteration but is
    kept general for the multi-dye extension.

    Args:
        rotamer_coords:   Cluster centre coordinates, shape (n_clusters, n_dye_atoms, 3).
        initial_weights:  Initial Boltzmann weights, shape (n_clusters,), sum ≈ 1.
        protein_coords:   Static protein atom positions, shape (n_prot_atoms, 3).
        dye_elements:     Element symbol for each dye atom.
        protein_elements: Element symbol for each protein atom.
        K:                Boltzmann scaling factor (IMP's ``K`` in Eq. 42).
                          K=1 ≈ kT-scale energies; increase to sharpen distribution.
        n_iter:           Number of mean-field iterations (default 10, usually
                          converges in 3–5).
        aabb_pad:         AABB padding in Å (same 3.5 Å default as IMP).
        r_cutoff:         Atom-pair distance cut-off in Å.

    Returns:
        weights: shape (n_clusters,), updated and renormalized.
    """
    n_clusters = rotamer_coords.shape[0]
    q = initial_weights.copy().astype(np.float64)

    # E_bb[i], the protein interaction of each cluster (Eq. 39). Constant for a
    # static protein, so it is computed once -- the all-pairs work is C++
    # (:file:`include/IMP/bff/RotamerEnergy.h`), bounding-box pre-filter and all.
    E_bb = _pair_energy_matrix(
        rotamer_coords, protein_coords, dye_elements, protein_elements,
        r_cutoff, aabb_pad)[:, 0]

    # Iterative mean-field update — Eq. 41 / 42
    # For single dye: E_P[i] = E_bb[i]  (no SC-SC cross term)
    for _iteration in range(n_iter):
        E_P = E_bb.copy()  # Eq. 41 simplified for single dye

        # Eq. 42: q[i] ∝ P_boltz[i] * exp(-K * E_P[i])
        log_q = np.log(np.maximum(q, 1e-300)) - K * E_P
        log_q -= log_q.max()               # numerical stability (log-sum-exp)
        q_new = np.exp(log_q)
        denom = q_new.sum()
        if denom > 0:
            q = q_new / denom
        else:
            # Degenerate: fall back to uniform
            q = np.ones(n_clusters, dtype=np.float64) / n_clusters

    return q


# ---------------------------------------------------------------------------
# Multi-dye extension — Eq. 40/41/42 with dye–dye SC cross term
# ---------------------------------------------------------------------------

def rotamer_mean_field_weights_multi_dye(
    rotamer_coords_list: list[np.ndarray],    # list of (n_i, n_atoms_i, 3)
    initial_weights_list: list[np.ndarray],   # list of (n_i,)
    protein_coords: np.ndarray,               # (n_prot_atoms, 3)
    dye_elements_list: list[list[str]],
    protein_elements: list[str],
    K: float = 1.0,
    n_iter: int = 10,
    aabb_pad: float = 3.5,
    r_cutoff: float = 12.0,
) -> list[np.ndarray]:
    """Multi-dye mean-field weight update including dye–dye cross-interactions.

    Implements equations (39)–(42) from IMP RotamerCalculator.cpp for
    multiple dyes whose rotamer clouds may interact with each other.

    Args:
        rotamer_coords_list:  One array per dye, shape (n_i, n_atoms_i, 3).
        initial_weights_list: Initial weights per dye, each shape (n_i,).
        protein_coords:       Static protein atoms, shape (n_prot_atoms, 3).
        dye_elements_list:    Element lists per dye.
        protein_elements:     Element list for protein.
        K:                    Boltzmann scaling (Eq. 42).
        n_iter:               Mean-field iterations.
        aabb_pad:             AABB padding in Å.
        r_cutoff:             Atom-pair distance cut-off in Å.

    Returns:
        Updated weight arrays, one per dye.
    """
    n_dyes = len(rotamer_coords_list)
    q_list = [w.copy().astype(np.float64) for w in initial_weights_list]

    # Everything below the iteration is independent of the weights, so it is
    # computed once. The version this replaces rebuilt every dye-dye energy
    # inside the loop -- n_iter times the work for the same answer.
    E_bb_all = [
        _pair_energy_matrix(rotamer_coords_list[d], protein_coords,
                            dye_elements_list[d], protein_elements,
                            r_cutoff, aabb_pad)[:, 0]
        for d in range(n_dyes)
    ]

    # One matrix per *ordered* pair, each built with its own atom ordering.
    E_sc = {}
    for d1 in range(n_dyes):
        for d2 in range(d1 + 1, n_dyes):
            m = _pair_energy_matrix(
                rotamer_coords_list[d1], rotamer_coords_list[d2],
                dye_elements_list[d1], dye_elements_list[d2],
                r_cutoff, aabb_pad)
            E_sc[(d1, d2)] = m
            # The interaction is symmetric, so the reverse order is the
            # transpose -- not the same matrix read the other way round, which
            # is the mistake the old code made.
            E_sc[(d2, d1)] = m.T

    # Mean-field iteration — Eq. 41 / 42
    for _iteration in range(n_iter):
        E_P_all = [E_bb_all[d].copy() for d in range(n_dyes)]

        # Eq. 41: each partner's energies weighted by its current probability.
        # A matrix-vector product now that the energies are precomputed.
        for d1 in range(n_dyes):
            for d2 in range(n_dyes):
                if d1 == d2:
                    continue
                E_P_all[d1] += E_sc[(d1, d2)] @ q_list[d2]

        # Eq. 42: renormalize each dye independently
        for d in range(n_dyes):
            log_q = np.log(np.maximum(q_list[d], 1e-300)) - K * E_P_all[d]
            log_q -= log_q.max()
            q_new = np.exp(log_q)
            denom = q_new.sum()
            q_list[d] = q_new / denom if denom > 0 else np.ones_like(q_new) / len(q_new)

    return q_list


# --------------------------------------------------------------------------
# rotamer
# --------------------------------------------------------------------------
"""Rotamer-protein interaction scoring."""

_GAS_CONSTANT = 1.9858775e-3


@dataclass
class RotamerScoreResult:
    """Result of rotamer-protein scoring.

    Attributes
    ----------
    weights : numpy.ndarray
        Normalized Boltzmann weights.
    partition : float
        Partition function before normalization.
    energies : numpy.ndarray
        Raw potential plus electrostatic energies.
    """

    weights: np.ndarray
    partition: float
    energies: np.ndarray


@lru_cache(maxsize=4096)
def _atom_type(atom_name: str) -> str:
    """Return a coarse atom type from an atom name.

    Memoised: a structure has thousands of atoms and a handful of distinct
    names, so the same strings arrive over and over. Profiling one FRETpredict
    comparison counted 292 872 calls, all of them resolving one of a few dozen
    names to one of five letters.

    Parameters
    ----------
    atom_name : str
        Atom name.

    Returns
    -------
    str
        One-letter atom type.
    """
    name = str(atom_name).strip().upper()
    for atom in ("S", "N", "O", "H", "C"):
        if name.startswith(atom):
            return atom
    return "C"


def _selector_matches(selector: str, atom_name: str, resname: str | None = None) -> bool:
    """Return whether one selector matches an atom.

    Parameters
    ----------
    selector : str
        FRETpredict-style selector.
    atom_name : str
        Atom name.
    resname : str, optional
        Atom residue name.

    Returns
    -------
    bool
        True if the selector matches the atom.
    """
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


def selector_resnames(selector: str | list[str] | None) -> set[str]:
    """Extract residue names from FRETpredict-style selectors.

    Parameters
    ----------
    selector : str or list of str or None
        Selector such as ``C7 and resname T48``.

    Returns
    -------
    set[str]
        Residue names appearing in the selector.
    """
    if selector is None:
        return set()
    if isinstance(selector, str):
        selectors = [selector]
    else:
        selectors = list(selector)
    resnames: set[str] = set()
    for item in selectors:
        for part in str(item).split(" and "):
            tokens = part.strip().split()
            if len(tokens) >= 2 and tokens[0].lower() == "resname":
                resnames.add(tokens[1].upper())
    return resnames


def _selector_atom_names(selector: str | list[str] | None, atom_names: list[str] | None = None, resnames: list[str] | None = None) -> set[str]:
    """Extract atom names from FRETpredict-style selectors.

    Parameters
    ----------
    selector : str or list of str or None
        Selector such as ``C7 and resname T48``.
    atom_names : list of str, optional
        Atom names aligned with ``resnames``.
    resnames : list of str, optional
        Residue names aligned with atom names. When provided, selectors are
        matched against their ``resname`` terms.

    Returns
    -------
    set[str]
        Atom names.
    """
    if selector is None:
        return set()
    if isinstance(selector, str):
        selectors = [selector]
    else:
        selectors = list(selector)
    wanted_resnames = selector_resnames(selectors)
    if wanted_resnames and len(wanted_resnames) > 1:
        return set()
    names: set[str] = set()
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


def _site_mask(
    atom_names: list[str],
    residue_indices: list[int] | None = None,
    site_residue: list[int] | int | None = None,
    site_chain: str | None = None,
    chain_ids: list[str] | None = None,
    mask_backbone: bool = True,
) -> np.ndarray:
    """Return a mask for placement-residue atoms.

    Parameters
    ----------
    atom_names : list of str
        Atom names.
    residue_indices : list of int, optional
        Residue indices aligned with atom names.
    site_residue : list of int, int, or None
        Placement residue indices used to exclude site atoms.
    site_chain : str, optional
        Placement chain ID used with ``site_residue``.
    chain_ids : list of str, optional
        Chain IDs aligned with atom names.
    mask_backbone : bool
        Also mask backbone atoms by name.

    Returns
    -------
    numpy.ndarray
        Boolean mask.
    """
    mask = np.zeros(len(atom_names), dtype=bool)
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


def _hydrogen_mask(atom_names: list[str]) -> np.ndarray:
    """Return a mask for hydrogen atoms.

    Parameters
    ----------
    atom_names : list[str]
        Atom names.

    Returns
    -------
    numpy.ndarray
        Boolean mask.
    """
    return np.array([_atom_type(name) == "H" for name in atom_names], dtype=bool)


def _protein_charge_mask(atom_names: list[str], resnames: list[str]) -> np.ndarray:
    """Return protein charge codes for Debye-Huckel scoring.

    Parameters
    ----------
    atom_names : list[str]
        Atom names.
    resnames : list[str]
        Residue names.

    Returns
    -------
    numpy.ndarray
        Charge code array.
    """
    charges = np.zeros(len(atom_names), dtype=float)
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


def _rotamer_charge_mask(
    atom_names: list[str],
    positive: list[str] | None = None,
    negative: list[str] | None = None,
    resnames: list[str] | None = None,
) -> np.ndarray:
    """Return rotamer charge codes for Debye-Huckel scoring.

    Parameters
    ----------
    atom_names : list[str]
        Atom names.
    positive : list of str, optional
        Positively charged atom selectors.
    negative : list of str, optional
        Negatively charged atom selectors.
    resnames : list of str, optional
        Residue names aligned with atom names.

    Returns
    -------
    numpy.ndarray
        Charge code array.
    """
    positive = [] if positive is None else list(positive)
    negative = [] if negative is None else list(negative)
    if len(selector_resnames(positive)) > 1:
        positive = []
    if len(selector_resnames(negative)) > 1:
        negative = []
    _selector_atom_names(positive, atom_names=atom_names, resnames=resnames)
    _selector_atom_names(negative, atom_names=atom_names, resnames=resnames)
    charges = np.zeros(len(atom_names), dtype=float)
    for i, name in enumerate(atom_names):
        resname = None if resnames is None or i >= len(resnames) else resnames[i]
        is_negative = any(_selector_matches(selector, name, resname) for selector in negative)
        is_positive = any(_selector_matches(selector, name, resname) for selector in positive)
        if is_negative:
            charges[i] = -1.0
        elif is_positive:
            charges[i] = 0.5
    return charges


def _scaled_parameters(
    atom_names: list[str],
    sigma_scaling: float,
    epsilon_scaling: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Return scaled LJ parameters for atom names.

    Parameters
    ----------
    atom_names : list[str]
        Atom names.
    sigma_scaling : float
        Sigma scaling factor.
    epsilon_scaling : float
        Epsilon scaling factor.

    Returns
    -------
    tuple[numpy.ndarray, numpy.ndarray]
        Scaled ``Rmin2`` and ``epsilon`` arrays.
    """
    rmin2, eps = lj_parameter_arrays(_atom_type(name) for name in atom_names)
    return rmin2 * sigma_scaling, eps * epsilon_scaling


def compute_rotamer_score(
    rotamer_coords: np.ndarray,
    protein_coords: np.ndarray,
    protein_atom_names: list[str],
    protein_resnames: list[str],
    rotamer_atom_names: list[str],
    rotamer_metadata: dict[str, Any] | None = None,
    rotamer_resnames: list[str] | None = None,
    protein_residue_indices: list[int] | None = None,
    protein_chain_ids: list[str] | None = None,
    site_residue: list[int] | int | None = None,
    site_chain: str | None = None,
    rotamer_weights: np.ndarray | None = None,
    temperature: float = 300.0,
    ignore_h: bool = True,
    electrostatic: bool = False,
    potential: str = "lj",
    sigma_scaling: float = 0.5,
    epsilon_scaling: float = 1.0,
) -> RotamerScoreResult:
    """Compute Boltzmann weights for a rotamer ensemble.

    Parameters
    ----------
    rotamer_coords : numpy.ndarray
        Rotamer coordinates with shape ``(n_rotamers, n_atoms, 3)`` in Å.
    protein_coords : numpy.ndarray
        Protein coordinates with shape ``(n_atoms, 3)`` in Å.
    protein_atom_names : list of str
        Protein atom names.
    protein_resnames : list of str
        Protein residue names.
    rotamer_atom_names : list of str
        Rotamer atom names.
    rotamer_metadata : dict, optional
        Rotamer metadata containing charged atom selectors.
    rotamer_resnames : list of str, optional
        Rotamer residue names aligned with atom names.
    protein_residue_indices : list of int, optional
        Protein residue numbers aligned with atom names.
    protein_chain_ids : list of str, optional
        Protein chain IDs aligned with atom names.
    rotamer_weights : numpy.ndarray, optional
        External rotamer population weights.
    temperature : float
        Temperature in K.
    ignore_h : bool
        Ignore hydrogen atoms.
    electrostatic : bool
        Include Debye-Huckel electrostatics.
    potential : {'lj', 'gauss'}
        Potential type.
    sigma_scaling : float
        LJ sigma scaling.
    epsilon_scaling : float
        LJ epsilon scaling.
    site_residue : list of int, int, or None
        Placement residue indices used to exclude site atoms.
    site_chain : str, optional
        Placement chain ID used with ``site_residue``.

    Returns
    -------
    RotamerScoreResult
        Normalized weights, partition function, and raw energies.
    """
    # Validated first, before any of the early returns below: an unknown
    # potential is a caller error whatever the geometry happens to be, and it
    # used to slip through whenever the atom selection came out empty.
    if potential not in ("lj", "gauss"):
        raise ValueError(f"Unknown potential {potential!r}")

    metadata = rotamer_metadata or {}
    protein_coords = np.asarray(protein_coords, dtype=np.float64)
    rotamer_coords = np.asarray(rotamer_coords, dtype=np.float64)

    protein_h = _hydrogen_mask(protein_atom_names)
    protein_mask = ~(protein_h if ignore_h else np.zeros(len(protein_atom_names), dtype=bool))
    protein_mask &= ~_site_mask(
        protein_atom_names,
        protein_residue_indices,
        site_residue=site_residue,
        site_chain=site_chain,
        chain_ids=protein_chain_ids,
        mask_backbone=False,
    )

    rotamer_backbone = _site_mask(rotamer_atom_names)
    rotamer_h = _hydrogen_mask(rotamer_atom_names)
    rotamer_mask = ~(rotamer_backbone | (rotamer_h if ignore_h else np.zeros(len(rotamer_atom_names), dtype=bool)))

    protein_idx = np.flatnonzero(protein_mask)
    rotamer_idx = np.flatnonzero(rotamer_mask)
    if protein_idx.size == 0 or rotamer_idx.size == 0:
        weights = np.ones(rotamer_coords.shape[0], dtype=np.float64)
        weights /= weights.sum()
        return RotamerScoreResult(weights=weights, partition=1.0, energies=np.zeros(rotamer_coords.shape[0]))

    protein_types = [_atom_type(name) for name in np.asarray(protein_atom_names)[protein_idx]]
    rotamer_types = [_atom_type(name) for name in np.asarray(rotamer_atom_names)[rotamer_idx]]
    protein_rmin2, protein_eps = _scaled_parameters(protein_types, sigma_scaling, epsilon_scaling)
    rotamer_rmin2, rotamer_eps = _scaled_parameters(rotamer_types, sigma_scaling, epsilon_scaling)

    eps_ij = np.sqrt(np.multiply.outer(rotamer_eps, protein_eps))
    rmin_ij = np.add.outer(rotamer_rmin2, protein_rmin2)

    protein_positions = protein_coords[protein_idx]
    q_rotamer = np.zeros(rotamer_idx.size, dtype=float)
    q_protein = np.zeros(protein_idx.size, dtype=float)

    if electrostatic:
        q_rotamer = _rotamer_charge_mask(
            [rotamer_atom_names[i] for i in rotamer_idx],
            metadata.get("positive"),
            metadata.get("negative"),
            None if rotamer_resnames is None else [rotamer_resnames[i] for i in rotamer_idx],
        )
        q_protein = _protein_charge_mask(
            [protein_atom_names[i] for i in protein_idx],
            [protein_resnames[i] for i in protein_idx],
        )

    # The all-pairs inner loop is C++ (:file:`include/IMP/bff/RotamerEnergy.h`).
    # Every dye atom against every protein atom, for every conformer -- and the
    # Python it replaces built an (n_dye x n_protein) distance matrix per
    # conformer, so it allocated once per rotamer and took the square root of
    # every pair before deciding almost all of them were beyond the cutoff.
    selected = np.ascontiguousarray(
        np.asarray(rotamer_coords, dtype=np.float64)[:, rotamer_idx, :])
    energies = np.asarray(IMP.bff.rotamer_interaction_energies(
        selected.ravel(),
        np.ascontiguousarray(protein_positions, dtype=np.float64).ravel(),
        np.ascontiguousarray(rmin_ij, dtype=np.float64).ravel(),
        np.ascontiguousarray(eps_ij, dtype=np.float64).ravel(),
        np.ascontiguousarray(q_rotamer, dtype=np.float64) if electrostatic else np.empty(0),
        np.ascontiguousarray(q_protein, dtype=np.float64) if electrostatic else np.empty(0),
        int(rotamer_coords.shape[0]), int(rotamer_idx.size), int(protein_idx.size),
        IMP.bff.ROTAMER_POTENTIAL_GAUSS if potential == "gauss"
        else IMP.bff.ROTAMER_POTENTIAL_LJ,
    ), dtype=np.float64).reshape(-1, 2)
    pot_energy = np.ascontiguousarray(energies[:, 0])
    dh_energy = np.ascontiguousarray(energies[:, 1])

    boltzmann = np.exp(-pot_energy / (_GAS_CONSTANT * temperature) - dh_energy)
    library_weights = np.asarray(rotamer_weights if rotamer_weights is not None else metadata.get("weights", np.ones(rotamer_coords.shape[0])), dtype=np.float64)
    if library_weights.shape != boltzmann.shape:
        library_weights = np.ones_like(boltzmann)
    boltzmann *= library_weights
    boltzmann = np.nan_to_num(boltzmann, nan=0.0, posinf=0.0, neginf=0.0)
    partition = float(np.sum(boltzmann))
    if partition <= 0.0:
        weights = np.ones(rotamer_coords.shape[0], dtype=np.float64) / rotamer_coords.shape[0]
        return RotamerScoreResult(weights=weights, partition=0.0, energies=pot_energy + dh_energy)
    weights = boltzmann / partition
    return RotamerScoreResult(weights=weights, partition=partition, energies=pot_energy + dh_energy)
