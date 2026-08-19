"""Stage 4 — what the sampled ensemble looks like, as a structural quantity.

The dye's occupancy in space, and quantities read off it. Distinct from
:mod:`IMP.bff.observables`, deliberately: that package owns the projection onto
what an *experiment* would see, with a contract that keeps convolution, pileup,
binning and counting noise outside. This one answers structural questions — how
the dye density is distributed, how a kinetic ensemble transfers — without
committing to a measurement.

Both are stage 4. Keeping them apart means the experiment-neutral contract has a
package to live in, rather than being a rule about part of a larger one.

* :mod:`~IMP.bff.analysis` — occupancy grids for a mobile component,
  region-resolved, and their projections.
* :mod:`~IMP.bff.analysis` — FRET over a kinetic rotamer ensemble,
  where the transfer is resolved per conformer rather than averaged first.
"""

from __future__ import annotations

from pathlib import Path
import RMF
import math

import numpy as np

from IMP.bff.io.cif import read_dye_forcefield_cif
from IMP.bff.io.cif import read_component_template_cif, region_features
import IMP
import IMP.algebra
import IMP.atom
import IMP.core
import IMP.em
import IMP.rmf

# --------------------------------------------------------------------------
# density
# --------------------------------------------------------------------------
"""Build region-specific mobile-component occupancy EM grids using IMP.em."""

#!/usr/bin/env python






SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_TRAJ_ROOT = SCRIPT_DIR / "traj_latest"
DEFAULT_OUTPUT_DIR = SCRIPT_DIR / "analysis" / "latest"
REGION_ORDER = ("linker", "top", "middle", "bottom")


def atom_name(atom, site_name_lookup=None):
    # Priority 1: Use IMP.atom.Atom input index to look up formal name from system CIF
    if IMP.atom.Atom.get_is_setup(atom):
        idx = IMP.atom.Atom(atom).get_input_index()
        if site_name_lookup and idx in site_name_lookup:
            return site_name_lookup[idx]
        elif site_name_lookup and (idx + 1) in site_name_lookup:
            return site_name_lookup[idx + 1]

    # Priority 2: Use established AtomType if available
    if IMP.atom.Atom.get_is_setup(atom):
        at = IMP.atom.Atom(atom).get_atom_type().get_string()
        nm = at.replace("HET: ", "").replace("HET:", "").strip()
    else:
        nm = atom.get_name()

    # Always take the last part of a slash-separated name
    if "/" in nm:
        nm = nm.split("/")[-1]
    return nm


def read_pdb_quiet(path, model):
    old = IMP.get_log_level()
    try:
        IMP.set_log_level(IMP.SILENT)
        return IMP.atom.read_pdb(str(path), model, IMP.atom.AllPDBSelector())
    finally:
        IMP.set_log_level(old)


def vector_length(v):
    return math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2])


def center_of_mass(atoms):
    c = IMP.algebra.Vector3D(0.0, 0.0, 0.0)
    for a in atoms:
        c += IMP.core.XYZ(a).get_coordinates()
    return c / max(1, len(atoms))


def _vec_to_np(v):
    return np.array([float(v[0]), float(v[1]), float(v[2])], dtype=float)


def _compute_fixed_axis(
    fixed_atoms, fixed_name_lookup=None, ref_axis=None, axis_element=None
):
    """Compute the principal (minimal-variance) axis of the fixed-component atoms.

    Parameters
    ----------
    axis_element : str or None
        One-letter element symbol used to orient the axis sign (e.g. ``"S"``
        for sulfonyl atoms).  When provided the axis is flipped so that it
        points *toward* the centroid of atoms whose name starts with this
        letter.  When ``None`` the sign is determined by ``ref_axis`` only.
    """
    coords = np.array(
        [_vec_to_np(IMP.core.XYZ(a).get_coordinates()) for a in fixed_atoms],
        dtype=float,
    )
    center = coords.mean(axis=0)
    centered = coords - center
    cov = centered.T @ centered / max(1, centered.shape[0])
    evals, evecs = np.linalg.eigh(cov)
    axis = evecs[:, int(np.argmin(evals))]
    axis /= max(np.linalg.norm(axis), 1e-12)

    if axis_element:
        marker_coords = []
        for a in fixed_atoms:
            nm = atom_name(a, fixed_name_lookup)
            if nm.startswith(axis_element):
                marker_coords.append(_vec_to_np(IMP.core.XYZ(a).get_coordinates()))
        if marker_coords:
            marker_center = np.array(marker_coords, dtype=float).mean(axis=0)
            v = marker_center - center
            if np.dot(axis, v) < 0:
                axis = -axis

    if ref_axis is not None and np.dot(axis, ref_axis) < 0:
        axis = -axis

    return center, axis


def _compute_long_axis(atoms, ref_axis=None):
    if not atoms:
        return None
    coords = np.array(
        [_vec_to_np(IMP.core.XYZ(a).get_coordinates()) for a in atoms], dtype=float
    )
    if coords.shape[0] < 2:
        return None
    center = coords.mean(axis=0)
    centered = coords - center
    cov = centered.T @ centered / max(1, centered.shape[0])
    evals, evecs = np.linalg.eigh(cov)
    axis = evecs[:, int(np.argmax(evals))]
    axis /= max(np.linalg.norm(axis), 1e-12)
    if ref_axis is not None and np.dot(axis, ref_axis) < 0:
        axis = -axis
    return axis


def _write_binned_profile(
    values_by_region, region_order, out_csv, bin_width, region_colors,
    origin=None,
):
    """A binned profile per region, as CSV.

    One function for both profiles the analysis writes. They differed in one
    line: the axis profile bins from ``min(all_values)``, the radial one from
    zero, because a distance from an axis has an origin and a position along it
    does not. That is *origin* -- ``None`` means "from the smallest value".
    """
    all_values = []
    for region in region_order:
        all_values.extend(values_by_region.get(region, []))

    col_header = (
        ",".join(f"{r}_count" for r in region_order)
        + ","
        + ",".join(f"{r}_density" for r in region_order)
    )

    if not all_values:
        with open(out_csv, "w") as out:
            out.write(f"bin_start_A,bin_end_A,center_A,{col_header}\n")
        return

    zmin = min(all_values) if origin is None else float(origin)
    zmax = max(all_values)
    n_bins = max(1, int(math.ceil((zmax - zmin) / bin_width)))
    edges = [zmin + i * bin_width for i in range(n_bins + 1)]
    centers = [0.5 * (edges[i] + edges[i + 1]) for i in range(n_bins)]

    counts = {r: [0] * n_bins for r in region_order}
    for region in region_order:
        vals = values_by_region.get(region, [])
        for z in vals:
            idx = int((z - zmin) // bin_width)
            idx = max(0, min(n_bins - 1, idx))
            counts[region][idx] += 1

    density = {}
    for region in region_order:
        total = max(1, sum(counts[region]))
        density[region] = [c / (total * bin_width) for c in counts[region]]

    with open(out_csv, "w") as out:
        out.write(f"bin_start_A,bin_end_A,center_A,{col_header}\n")
        for i in range(n_bins):
            count_vals = ",".join(str(counts[r][i]) for r in region_order)
            density_vals = ",".join(f"{density[r][i]:.8f}" for r in region_order)
            out.write(
                f"{edges[i]:.4f},{edges[i + 1]:.4f},{centers[i]:.4f},{count_vals},{density_vals}\n"
            )


    # Rendering deliberately omitted: the CSV above carries every
    # number the plot showed, and IMP.bff carries no dependency
    # beyond what IMP itself brings. Plot it in the application.


    # Rendering deliberately omitted: the CSV above carries every
    # number the plot showed, and IMP.bff carries no dependency
    # beyond what IMP itself brings. Plot it in the application.


def _write_orientation_profile(frame_rows, out_csv):
    with open(out_csv, "w") as out:
        out.write("frame_index,angle_deg,abs_cos_theta,cos_theta\n")
        for row in frame_rows:
            out.write(
                f"{row['frame_index']},{row['angle_deg']:.6f},{row['abs_cos_theta']:.8f},{row['cos_theta']:.8f}\n"
            )

    if not frame_rows:
        return

    frames = [r["frame_index"] for r in frame_rows]
    angles = [r["angle_deg"] for r in frame_rows]

    # Rendering deliberately omitted: the CSV above carries every
    # number the plot showed, and IMP.bff carries no dependency
    # beyond what IMP itself brings. Plot it in the application.


def resolve_rmf_path(traj_root, system_name):
    """Locate the RMF file for *system_name* under *traj_root*."""
    base = traj_root / system_name
    candidates = [
        base / "rmfs" / "0.rmf3",
        base / "traj.rmf3",
        base / "traj_init.rmf3",
    ]
    for p in candidates:
        if p.exists():
            return p
    raise FileNotFoundError(f"No RMF found for {system_name} under {traj_root}")


def discover_run_roots(runs_root):
    roots = []
    for p in sorted(runs_root.glob("run*")):
        if p.is_dir():
            roots.append(p)
    if not roots:
        raise FileNotFoundError(f"No run directories found under {runs_root}")
    return roots


def load_site_name_lookup(traj_root, system_name, component):
    """Load site-name lookup indexed by serial for *component* from the system CIF."""
    run_dir = traj_root / system_name
    system_cif = run_dir / "system.cif"

    if not system_cif.exists():
        return {}

    system = read_dye_forcefield_cif(str(system_cif))
    lookup = {}
    for site in system.sites:
        if site.get("component") == component:
            serial = site.get("site_serial")
            if serial is not None:
                lookup[int(serial)] = site.get("atom_name", site.get("id"))
    return lookup


def _load_regions_from_template(template_cif_path):
    """Return (regions_dict, region_colors_dict, region_order) from a template CIF.

    regions_dict : {region_name: [{"name": atom_name, "occurrence": int}, ...]}
    region_colors_dict : {region_name: color_str}
    region_order : list of region names in template definition order
    """
    template = read_component_template_cif(str(template_cif_path))
    rf = region_features(template)  # {feature_id: color}
    features = template.get("features", {})

    region_order = [fid for fid in features if fid in rf]
    regions = {}
    colors = {}
    for fid in region_order:
        spec = features[fid]
        regions[fid] = spec.get("atoms", [])
        colors[fid] = rf[fid]
    return regions, colors, region_order


def build_name_map(atoms, site_name_lookup=None):
    use_enum = False
    if site_name_lookup:
        # If all particles have lost their input_index (-1), fallback to enumeration
        use_enum = all(
            (not IMP.atom.Atom.get_is_setup(a))
            or (IMP.atom.Atom(a).get_input_index() == -1)
            for a in atoms
        )

    out = {}
    for i, a in enumerate(atoms):
        if use_enum:
            nm = site_name_lookup.get(i + 1, atom_name(a))
        else:
            nm = atom_name(a, site_name_lookup)
        out.setdefault(nm, []).append(a)
    return out


def collect_region_atoms(mobile_atoms, regions, site_name_lookup=None):
    by_name = build_name_map(mobile_atoms, site_name_lookup)
    selected = {}
    for region_name, atom_specs in regions.items():
        atoms = []
        missing = []
        for spec in atom_specs:
            name = spec.get("name") if isinstance(spec, dict) else str(spec)
            occ = int(spec.get("occurrence", 1)) if isinstance(spec, dict) else 1
            matches = by_name.get(name, [])
            if not matches:
                # Fallback: look for generic name + suffix (e.g. O3 -> O3A, O3B)
                # Sort to ensure consistent occurrence mapping
                candidates = sorted(
                    [
                        k
                        for k in by_name.keys()
                        if k.startswith(name) and k[len(name) :].isalpha()
                    ]
                )
                if candidates:
                    for cand in candidates:
                        matches.extend(by_name[cand])

            if not matches:
                missing.append(name)
                continue
            idx = occ - 1
            if 0 <= idx < len(matches):
                atoms.append(matches[idx])
            else:
                atoms.append(matches[0])
        selected[region_name] = atoms
        if missing:
            print(
                f"    Warning: {region_name} missing {len(missing)} atoms: "
                + ", ".join(missing)
            )
    return selected


def write_region_density(points, out_path, resolution, voxel_size):
    model = IMP.Model()
    particles = []
    for x, y, z in points:
        p = IMP.Particle(model)
        xyzr = IMP.core.XYZR.setup_particle(p)
        xyzr.set_coordinates(IMP.algebra.Vector3D(x, y, z))
        xyzr.set_radius(1.0)
        IMP.atom.Mass.setup_particle(p, 1.0)
        particles.append(p)

    if not particles:
        raise ValueError(f"No points available for {out_path}")

    density = IMP.em.particles2density(particles, resolution, voxel_size)
    density.std_normalize()
    IMP.em.write_map(density, str(out_path), IMP.em.MRCReaderWriter())


def write_radial_histogram(distances, out_csv, bin_width):
    if not distances:
        with open(out_csv, "w") as out:
            out.write("bin_start_A,bin_end_A,count\n")
        return

    max_d = max(distances)
    n_bins = max(1, int(math.ceil(max_d / bin_width)))
    counts = [0] * n_bins
    for d in distances:
        idx = min(n_bins - 1, int(d // bin_width))
        counts[idx] += 1

    with open(out_csv, "w") as out:
        out.write("bin_start_A,bin_end_A,count\n")
        for i, c in enumerate(counts):
            lo = i * bin_width
            hi = lo + bin_width
            out.write(f"{lo:.3f},{hi:.3f},{c}\n")


def analyze_dye_density(
    traj_roots,
    output_dir,
    mobile,
    mobile_template_cif,
    resolution,
    voxel_size,
    bin_width,
    max_frames,
    fixed_name=None,
    run_dir_name=None,
    axis_element=None,
):
    """Analyze occupancy for one mobile component.

    Parameters
    ----------
    mobile : str
        Component name of the mobile molecule in the RMF hierarchy.
    mobile_template_cif : str or Path
        Path to the template CIF for the mobile component.  Regions and their
        colors are read from this file.
    fixed_name : str, optional
        Component name of the fixed component in the RMF hierarchy.
        Inferred from non-mobile children if not provided.
    run_dir_name : str, optional
        Name of the run output subdirectory (e.g. ``"CX4_atto655_imp"``).
        This is the directory name under each traj_root that was written by
        the runner.  Defaults to ``"{fixed_name}_{mobile}_imp"`` when
        ``fixed_name`` is known, or auto-detected from the traj_root contents.
    axis_element : str or None
        Element symbol used to orient the fixed-component axis sign.
        ``None`` means no orientation correction.
    """
    print(f"\nAnalyzing {mobile}")
    print(f"  Template: {mobile_template_cif}")

    regions, region_colors, region_order = _load_regions_from_template(
        mobile_template_cif
    )

    total_used = 0
    total_available = 0

    all_points = {k: [] for k in region_order}
    radial_distances = {k: [] for k in region_order}
    axis_distances = {k: [] for k in region_order}
    xy_distances = {k: [] for k in region_order}

    axis_ref = None
    pi_axis_ref = None
    orientation_rows = []
    global_frame_index = 0

    for traj_root in traj_roots:
        _fixed_name = fixed_name
        # Determine the run subdirectory name.
        # Priority: explicit run_dir_name > "{fixed_name}_{mobile}_imp" > auto-detect.
        if run_dir_name:
            _run_dir = run_dir_name
        elif _fixed_name:
            _run_dir = f"{_fixed_name}_{mobile}_imp"
        else:
            # Auto-detect: find the first subdirectory ending with _{mobile}_imp
            candidates = sorted(
                p.name
                for p in traj_root.iterdir()
                if p.is_dir() and p.name.endswith(f"_{mobile}_imp")
            )
            if not candidates:
                raise FileNotFoundError(
                    f"No run directory ending with '_{mobile}_imp' found under {traj_root}"
                )
            _run_dir = candidates[0]

        rmf_path = resolve_rmf_path(traj_root, _run_dir)
        print(f"  RMF: {rmf_path}")

        fh = RMF.open_rmf_file_read_only(str(rmf_path))
        model = IMP.Model()
        IMP.set_check_level(IMP.NONE)
        roots = IMP.rmf.create_hierarchies(fh, model)
        if not roots:
            raise RuntimeError(f"No hierarchies found in {rmf_path}")
        root = roots[0]

        # Identify fixed and mobile hierarchies by name.
        fixed_hier = None
        mobile_hier = None
        for child in root.get_children():
            nm = child.get_name()
            if nm == mobile:
                mobile_hier = child
            elif _fixed_name is None or nm == _fixed_name:
                fixed_hier = child
        if fixed_hier is None or mobile_hier is None:
            raise RuntimeError(
                f"Could not find fixed/mobile hierarchies in {rmf_path}; "
                f"found {[c.get_name() for c in root.get_children()]}"
            )
        if _fixed_name is None:
            _fixed_name = fixed_hier.get_name()

        fixed_atoms = IMP.atom.get_by_type(fixed_hier, IMP.atom.ATOM_TYPE)
        mobile_atoms = IMP.atom.get_by_type(mobile_hier, IMP.atom.ATOM_TYPE)

        fixed_name_lookup = load_site_name_lookup(traj_root, _run_dir, _fixed_name)
        mobile_name_lookup = load_site_name_lookup(traj_root, _run_dir, mobile)
        if not fixed_atoms:
            fixed_atoms = IMP.core.get_leaves(fixed_hier)
        if not mobile_atoms:
            mobile_atoms = IMP.core.get_leaves(mobile_hier)

        region_atoms = collect_region_atoms(mobile_atoms, regions, mobile_name_lookup)

        # Build pi-system atoms from all region atoms except linker
        pi_atoms = []
        seen_ids = set()
        for key in region_order:
            if key == "linker":
                continue
            for a in region_atoms.get(key, []):
                aid = a.get_index() if hasattr(a, "get_index") else id(a)
                if aid in seen_ids:
                    continue
                seen_ids.add(aid)
                pi_atoms.append(a)

        n_frames_total = fh.get_number_of_frames()
        n_frames = (
            n_frames_total if max_frames <= 0 else min(n_frames_total, max_frames)
        )
        total_available += n_frames_total
        total_used += n_frames

        for fi in range(n_frames):
            IMP.rmf.load_frame(fh, RMF.FrameID(fi))
            fixed_center, axis_vec = _compute_fixed_axis(
                fixed_atoms, fixed_name_lookup, axis_ref, axis_element=axis_element
            )
            if axis_ref is None:
                axis_ref = axis_vec.copy()

            pi_axis = _compute_long_axis(pi_atoms, pi_axis_ref)
            if pi_axis is not None:
                if pi_axis_ref is None:
                    pi_axis_ref = pi_axis.copy()
                cos_theta = float(np.dot(pi_axis, axis_vec))
                cos_theta = max(-1.0, min(1.0, cos_theta))
                abs_cos = abs(cos_theta)
                angle_deg = math.degrees(math.acos(max(-1.0, min(1.0, abs_cos))))
                orientation_rows.append(
                    {
                        "frame_index": global_frame_index,
                        "angle_deg": angle_deg,
                        "abs_cos_theta": abs_cos,
                        "cos_theta": cos_theta,
                    }
                )
            global_frame_index += 1

            for region_name in region_order:
                for atom_p in region_atoms[region_name]:
                    c = _vec_to_np(IMP.core.XYZ(atom_p).get_coordinates())
                    rel_np = c - fixed_center
                    all_points[region_name].append((rel_np[0], rel_np[1], rel_np[2]))
                    radial_distances[region_name].append(float(np.linalg.norm(rel_np)))
                    zproj = float(np.dot(rel_np, axis_vec))
                    axis_distances[region_name].append(zproj)
                    perp = rel_np - zproj * axis_vec
                    xy_distances[region_name].append(float(np.linalg.norm(perp)))

        del fh

    print(f"  Frames used: {total_used}/{total_available}")

    mobile_out = output_dir / mobile
    mobile_out.mkdir(parents=True, exist_ok=True)

    for region_name in region_order:
        mrc_out = mobile_out / f"occupancy_{region_name}.mrc"
        write_region_density(
            all_points[region_name],
            mrc_out,
            resolution=resolution,
            voxel_size=voxel_size,
        )
        print(f"  Wrote grid: {mrc_out}")

        hist_out = mobile_out / f"radial_{region_name}.csv"
        write_radial_histogram(radial_distances[region_name], hist_out, bin_width)
        print(f"  Wrote radial histogram: {hist_out}")

    axis_csv = mobile_out / "axis_z_profile_regions.csv"
    _write_binned_profile(
        axis_distances, region_order, axis_csv, bin_width, region_colors
    )
    print(f"  Wrote axis profile CSV: {axis_csv}")

    xy_csv = mobile_out / "axis_xy_profile_regions.csv"
    # origin=0.0: a distance from the axis is measured from the axis, so the
    # first bin starts at zero rather than at the smallest distance observed.
    _write_binned_profile(
        xy_distances, region_order, xy_csv, bin_width, region_colors, origin=0.0
    )
    print(f"  Wrote xy profile CSV: {xy_csv}")

    orientation_csv = mobile_out / "axis_mobile_vs_fixed_orientation.csv"
    _write_orientation_profile(orientation_rows, orientation_csv)
    print(f"  Wrote orientation CSV: {orientation_csv}")

    if axis_ref is not None:
        import json

        axis_def = mobile_out / "fixed_axis_definition.json"
        summary = {}
        if orientation_rows:
            angles = np.array([r["angle_deg"] for r in orientation_rows], dtype=float)
            summary = {
                "n_frames": int(len(angles)),
                "mean_angle_deg": float(np.mean(angles)),
                "median_angle_deg": float(np.median(angles)),
                "std_angle_deg": float(np.std(angles)),
                "min_angle_deg": float(np.min(angles)),
                "max_angle_deg": float(np.max(angles)),
            }
        with open(axis_def, "w") as out:
            json.dump(
                {
                    "axis_unit_vector": axis_ref.tolist(),
                    "note": "Axis from smallest-variance PCA direction of fixed-component atoms; oriented toward axis_element centroid.",
                    "mobile_vs_fixed_orientation": summary,
                },
                out,
                indent=2,
            )
        print(f"  Wrote axis definition: {axis_def}")


# The click decorators for this moved to `bin/imp_bff` as
# `analyze-trajectories`; what is left is the function they called, which is
# library code and now importable without click.
def analyze_dye_trajectories(
    traj_root,
    runs_root,
    combine_runs,
    output_dir,
    mobiles,
    mobile_template_cifs,
    fixed_name,
    system_name,
    axis_element,
    resolution,
    voxel_size,
    bin_width,
    max_frames,
):
    if combine_runs:
        roots = discover_run_roots(runs_root if runs_root else (traj_root / "runs"))
    else:
        roots = [traj_root]

    # Build template CIF path list — one per mobile, or fall back to default location
    template_map = {}
    tpl_list = list(mobile_template_cifs) if mobile_template_cifs else []
    for i, mob in enumerate(mobiles):
        if i < len(tpl_list):
            tpl = Path(tpl_list[i])
        else:
            # Default: the bundled template, which is IMP module data
            from IMP.bff.tools import get_template_dir
            tpl = get_template_dir() / f"{mob}.template.cif"
        template_map[mob] = tpl

    for mob in mobiles:
        tpl = template_map[mob]
        if not tpl.exists():
            raise ValueError(
                f"Template CIF not found for mobile component '{mob}': {tpl}\n"
                "Pass --mobile-template-cif explicitly."
            )
        analyze_dye_density(
            traj_roots=roots,
            output_dir=output_dir,
            mobile=mob,
            mobile_template_cif=tpl,
            resolution=resolution,
            voxel_size=voxel_size,
            bin_width=bin_width,
            max_frames=max_frames,
            fixed_name=fixed_name,
            run_dir_name=system_name,
            axis_element=axis_element,
        )




# --------------------------------------------------------------------------
# kinetic_fret
# --------------------------------------------------------------------------
"""Physically exact FRET analysis using kinetic rotamer ensembles."""

def fret_efficiency_exact_kinetic(
    p_matrix: np.ndarray, 
    fret_rates: np.ndarray, 
    tau0: float, 
    dt: float = 1.0,
    weights: np.ndarray | None = None
) -> float:
    """Compute exact FRET efficiency by solving the Master Equation.
    
    This avoids approximations by accounting for the competition between 
    fluorescence decay and conformational transitions.
    
    Args:
        p_matrix: Transition probability matrix P (n x n) for time step dt.
        fret_rates: FRET rates for each state (n,). Units: 1/ns.
        tau0: Donor lifetime in absence of acceptor (ns).
        dt: Time step of the transitions in p_matrix (ns).
        weights: Stationary distribution (n,). If None, calculated as first eigenvector.
        
    Returns:
        Exact FRET efficiency.
    """
    n = p_matrix.shape[0]
    k_rad = 1.0 / tau0
    
    # 1. Convert Probability Matrix P to Rate Matrix M
    # For small dt: P = exp(M*dt) approx I + M*dt => M = (P - I) / dt
    # This assumes the transitions recorded are for a physical time interval.
    M = (p_matrix - np.eye(n)) / dt
    
    # 2. Setup the sink matrix (Radiative decay + FRET rates)
    # K is diagonal matrix of total decay rates
    K_fret = np.diag(fret_rates)
    
    # 3. Stationary distribution (initial population)
    if weights is None:
        # Find eigenvector with eigenvalue 1 for P (or 0 for M)
        evals, evecs = np.linalg.eig(p_matrix.T)
        idx = np.argmin(np.abs(evals - 1.0))
        weights = np.real(evecs[:, idx])
        weights /= weights.sum()
        
    # 4. Integrated populations G_j = int_0^inf p_j(t) dt.
    # With P[i, j] the probability of i -> j (rows sum to 1) the populations
    # evolve as a row vector, dp/dt = p M - p (k_rad I + K_fret), so
    # p0 = G (k_rad I + K_fret - M), i.e. G solves the *transposed* system
    # A^T G = w. Solving A G = w instead is only right for symmetric M; for a
    # non-symmetric P it gave the wrong fast-exchange limit (test_physics_invariants).
    A = k_rad * np.eye(n) + K_fret - M
    G = np.linalg.solve(A.T, weights)
    
    # 5. E = sum_i k_fret_i G_i (the fraction of excitations that leave via FRET)
    efficiency = np.sum(fret_rates * G)
    
    return float(efficiency)


def fret_efficiency_exact_kinetic_pair(
    dist_matrix: np.ndarray,
    kappa2_matrix: np.ndarray,
    p_d: np.ndarray,
    p_a: np.ndarray,
    weights_d: np.ndarray,
    weights_a: np.ndarray,
    R0: float = 52.0,
    tau0: float = 4.0,
    dt: float = 0.1
) -> float:
    """Calculate exact FRET efficiency for two kinetic ensembles.
    
    Args:
        dist_matrix: (nd, na) distances
        kappa2_matrix: (nd, na) orientation factors
        p_d: Donor transition matrix (nd, nd)
        p_a: Acceptor transition matrix (na, na)
        weights_d: Donor stationary weights
        weights_a: Acceptor stationary weights
        R0: Förster radius (A)
        tau0: Lifetime (ns)
        dt: Time step (ns)
    """
    # 1. Product space transition matrix P_total = P_d kron P_a
    # This represents the joint kinetics of both dyes.
    # Note: kron is expensive for large libraries.
    # P_total[i*na + j, k*na + l] = P_d[i,k] * P_a[j,l]
    P_total = np.kron(p_d, p_a)
    
    # 2. Product space weights and FRET rates
    w_total = np.outer(weights_d, weights_a).flatten()
    
    k_rad = 1.0 / tau0
    # k_fret_ij = k_rad * (R0/r_ij)^6 * 1.5 * kappa2_ij
    rate_ratios = (R0 / dist_matrix)**6 * (1.5 * kappa2_matrix)
    fret_rates = (k_rad * rate_ratios).flatten()
    
    return fret_efficiency_exact_kinetic(P_total, fret_rates, tau0, dt, weights=w_total)


def fret_efficiency_regimes(
    dist_matrix: np.ndarray,
    kappa2_matrix: np.ndarray,
    weights_d: np.ndarray,
    weights_a: np.ndarray,
    R0: float = 52.0
) -> dict:
    """Calculate FRET efficiency under different kinetic regimes.
    
    Args:
        dist_matrix: (n_donor, n_acceptor) distances
        kappa2_matrix: (n_donor, n_acceptor) orientation factors
        weights_d: (n_donor,) donor rotamer weights
        weights_a: (n_acceptor,) acceptor rotamer weights
        R0: Förster radius for kappa2=2/3 (A)
        
    Returns:
        Dict with 'static', 'dynamic', 'dynamic_plus' efficiencies.
    """
    # Combined weights
    w_ij = np.outer(weights_d, weights_a)
    
    # Rate constant k_fret / k_rad = (R0/r)^6 * (k2 / (2/3))
    # We use (1.5 * k2) to normalize kappa2 relative to 2/3
    rate_ratio = (R0 / dist_matrix)**6 * (1.5 * kappa2_matrix)
    
    # 1. Static Regime: Average of efficiencies
    eff_static = rate_ratio / (1 + rate_ratio)
    e_static = np.sum(w_ij * eff_static)
    
    # 2. Dynamic Regime: Use <kappa^2>
    k2_avg = np.sum(w_ij * kappa2_matrix)
    rate_ratio_dyn = (R0 / dist_matrix)**6 * (1.5 * k2_avg)
    eff_dyn = rate_ratio_dyn / (1 + rate_ratio_dyn)
    e_dynamic = np.sum(w_ij * eff_dyn)
    
    # 3. Dynamic+ Regime: Average of rates
    rate_avg = np.sum(w_ij * rate_ratio)
    e_dynamic_plus = rate_avg / (1 + rate_avg)
    
    return {
        "static": float(e_static),
        "dynamic": float(e_dynamic),
        "dynamic_plus": float(e_dynamic_plus),
        "kappa2_avg": float(k2_avg)
    }
