"""The cgdye MD runner's NMR distance-restraint path, retired 2026-08-19.

The only consumer of ``read_nmr_restraints``. It turned parsed restraints
into ``IMP.core.DistanceRestraint`` pairs with harmonic upper and lower
bounds, scaled by ``nmr_strength / error``.
"""

def _build_nmr_restraints(
    model,
    component_hiers,
    nmr_file,
    guest_component,
    nmr_strength,
    restraint_set,
    system=None,
    system_cif=None,
    use_minimum_restraints=True,
):
    if not nmr_file or not os.path.exists(nmr_file):
        return []

    atom_maps = _build_component_atom_maps(
        component_hiers, system=system, system_cif=system_cif
    )
    host_candidates = [c for c in atom_maps.keys() if c != guest_component]
    host_component = host_candidates[0] if len(host_candidates) == 1 else None

    data = read_nmr_restraints(nmr_file, fixed_component_name=host_component)
    distances = data.get("Distances", {})
    positions = data.get("Positions", {})
    sets = data.get("restraint_sets", {})

    # Only filter if the user explicitly requested a set. 
    # Otherwise, we include everything (CSP + all NOE groups).
    chosen = restraint_set
    prefix = sets.get(f"{chosen}_prefix", "") if chosen else ""
    common_prefix = sets.get("common_prefix", "")
    restraints = []

    noe_groups = defaultdict(list)

    for dname, d in distances.items():
        if prefix:
            if not dname.startswith(prefix):
                if not common_prefix or not dname.startswith(common_prefix):
                    continue

        p1d = positions.get(d.get("position1_name"))
        p2d = positions.get(d.get("position2_name"))
        if not p1d or not p2d:
            continue

        p1 = _resolve_position_particle(
            d.get("position1_name"),
            p1d,
            guest_component,
            atom_maps,
            host_component=host_component,
        )
        p2 = _resolve_position_particle(
            d.get("position2_name"),
            p2d,
            guest_component,
            atom_maps,
            host_component=host_component,
        )
        if p1 is None or p2 is None:
            continue

        target = float(d.get("distance", 0.0))
        err_pos = float(d.get("error_pos", 1.0))
        err_neg = float(d.get("error_neg", 1.0))
        k_pos = float(nmr_strength) / max(err_pos, 1e-6)
        k_neg = float(nmr_strength) / max(err_neg, 1e-6)
        rtype = d.get("distance_type", "AtomDistance")

        # Group NOESY (non-CSP) restraints by their set prefix if requested
        if use_minimum_restraints and "_RID" in dname and not dname.startswith(common_prefix):
            peak_id = dname.split("_RID")[0]
            noe_groups[peak_id].append((p1, p2, target, k_pos, k_neg, rtype, dname))
            continue

        if rtype in ("AtomDistance", "Atom", "AtomUpperBound"):
            r = IMP.core.DistanceRestraint(
                model, IMP.core.HarmonicUpperBound(target, k_pos), p1, p2
            )
            r.set_name(f"{dname}_upper")
            restraints.append(r)
        if rtype in ("AtomDistance", "Atom", "AtomLowerBound"):
            r = IMP.core.DistanceRestraint(
                model, IMP.core.HarmonicLowerBound(target, k_neg), p1, p2
            )
            r.set_name(f"{dname}_lower")
            restraints.append(r)

    # Process grouped NOESY ambiguous restraints
    for peak_id, items in noe_groups.items():
        rs = []
        for (p1, p2, target, k_pos, k_neg, rtype, dname) in items:
            if rtype in ("AtomDistance", "Atom", "AtomUpperBound"):
                r = IMP.core.DistanceRestraint(
                    model, IMP.core.HarmonicUpperBound(target, k_pos), p1, p2
                )
                r.set_name(f"{dname}_upper")
                rs.append(r)
            if rtype in ("AtomDistance", "Atom", "AtomLowerBound"):
                r = IMP.core.DistanceRestraint(
                    model, IMP.core.HarmonicLowerBound(target, k_neg), p1, p2
                )
                r.set_name(f"{dname}_lower")
                rs.append(r)
        
        if len(rs) == 1:
            restraints.append(rs[0])
        elif len(rs) > 1:
            mr = IMP.core.MinimumRestraint(1, rs)
            mr.set_name(peak_id)
            restraints.append(mr)

    return restraints
