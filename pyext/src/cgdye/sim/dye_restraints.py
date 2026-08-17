"""IMP restraint builder with element-aware dye nonbonded terms."""

import IMP
import IMP.container
import IMP.core

from IMP.bff.cgdye.sampling.scoring import compute_exclusions, compute_lj_pair_sites
from IMP.bff.cgdye.topology.dye import torsion_cosine


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
