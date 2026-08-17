#!/usr/bin/env python
"""Write a PyMOL scene script for mobile-component region coloring and MRC overlays.

Reads region definitions and colors from the mobile component's template CIF,
removing all system-specific hardcoding.
"""

from __future__ import annotations

import sys
from pathlib import Path

import click



from IMP.bff.cgdye.io.template_cif import read_component_template_cif, region_features


def parse_mol2_serials_by_name(mol2_path: Path) -> dict[str, list[int]]:
    """Return {atom_name: [serial, ...]} from the @<TRIPOS>ATOM section of a MOL2 file."""
    by_name: dict[str, list[int]] = {}
    in_atom = False
    with mol2_path.open() as fh:
        for line in fh:
            if line.startswith("@<TRIPOS>ATOM"):
                in_atom = True
                continue
            if line.startswith("@<TRIPOS>"):
                in_atom = False
                continue
            if not in_atom:
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            serial = int(parts[0])
            atom_name = parts[1]
            by_name.setdefault(atom_name, []).append(serial)
    return by_name


def region_serials(atom_specs: list, by_name: dict[str, list[int]]) -> list[int]:
    """Return sorted MOL2 serials for the given atom specs."""
    serials: list[int] = []
    for item in atom_specs:
        if isinstance(item, str):
            serials.extend(by_name.get(item, []))
            continue
        if isinstance(item, dict):
            name = str(item.get("name", "")).strip()
            occ = int(item.get("occurrence", 1))
            vals = by_name.get(name, [])
            idx = occ - 1
            if 0 <= idx < len(vals):
                serials.append(vals[idx])
    return sorted(set(serials))


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option(
    "--system-cif",
    type=click.Path(path_type=Path),
    required=True,
    help="Built system mmCIF (loaded into PyMOL as 'system').",
)
@click.option(
    "--mobile-mol2",
    type=click.Path(path_type=Path),
    required=True,
    help="MOL2 file of the mobile component.",
)
@click.option(
    "--mobile-name",
    type=str,
    required=True,
    help="Name of the mobile component (used as PyMOL object name).",
)
@click.option(
    "--mobile-template-cif",
    type=click.Path(path_type=Path),
    required=True,
    help="Template CIF for the mobile component; regions and colors are read from here.",
)
@click.option(
    "--analysis-dir",
    type=click.Path(path_type=Path),
    required=True,
    help="Directory containing occupancy_<region>.mrc files.",
)
@click.option(
    "--output-pml",
    type=click.Path(path_type=Path),
    default=None,
    help="Output PML path. Defaults to <analysis-dir>/view_regions.pml.",
)
def main(
    system_cif: Path,
    mobile_mol2: Path,
    mobile_name: str,
    mobile_template_cif: Path,
    analysis_dir: Path,
    output_pml: Path | None,
):
    analysis_dir = analysis_dir.resolve()
    system_cif = system_cif.resolve()
    mobile_mol2 = mobile_mol2.resolve()
    mobile_template_cif = mobile_template_cif.resolve()

    if output_pml is None:
        output_pml = analysis_dir / "view_regions.pml"
    output_pml = output_pml.resolve()
    output_pml.parent.mkdir(parents=True, exist_ok=True)

    # Load region definitions and colors from template CIF
    template = read_component_template_cif(str(mobile_template_cif))
    rf = region_features(template)  # {feature_id: color}
    features = template.get("features", {})
    region_order = [fid for fid in features if fid in rf]
    regions = {fid: features[fid].get("atoms", []) for fid in region_order}
    region_colors = rf

    by_name = parse_mol2_serials_by_name(mobile_mol2)

    lines: list[str] = []
    lines.append("reinitialize")
    lines.append(f'load "{system_cif}", system')
    lines.append(f'load "{mobile_mol2}", {mobile_name}')
    lines.append("hide everything")
    lines.append("show sticks, system")
    lines.append(f"show sticks, {mobile_name}")
    lines.append("color gray70, system")
    lines.append("set stick_radius, 0.18")

    for region in region_order:
        color = region_colors.get(region, "white")
        if color.startswith("#"):
            rgb = [int(color[i:i+2], 16) / 255.0 for i in (1, 3, 5)]
            pymol_color = f"color_{region}"
            lines.append(f"set_color {pymol_color}, {list(rgb)}")
            color = pymol_color

        serials = region_serials(regions.get(region, []), by_name)
        if not serials:
            continue
        sel = f"{mobile_name}_{region}"
        serial_expr = "+".join(str(s) for s in serials)
        lines.append(f"select {sel}, {mobile_name} and id {serial_expr}")
        lines.append(f"color {color}, {sel}")

    lines.append("set transparency, 0.20, system")
    lines.append(f"set transparency, 0.05, {mobile_name}")

    volume_regions = []
    for region in region_order:
        mrc_path = analysis_dir / f"occupancy_{region}.mrc"
        if not mrc_path.exists():
            continue
        color = region_colors.get(region, "white")
        if color.startswith("#"):
            color = f"color_{region}"
            
        volume_regions.append((region, color))
        map_obj = f"map_{region}"
        vol_obj = f"vol_{region}"
        lines.append(f'load "{mrc_path}", {map_obj}')
        lines.append(f"volume {vol_obj}, {map_obj}")
        lines.append(
            f"volume_color {vol_obj}, 0.00 {color} 0.00 0.02 {color} 0.30 0.08 {color} 0.95"
        )
        lines.append(f"set volume_mode, 2, {vol_obj}")
        lines.append(f"set volume_sigma, 0.9, {vol_obj}")
        lines.append(f"set volume_panel, 0, {vol_obj}")

    lines.append("python")
    lines.append("from pymol import cmd")
    for region, color in volume_regions:
        vol_obj = f"vol_{region}"
        lines.append(
            f"cmd.volume_color('{vol_obj}', '0.00 {color} 0.00 0.02 {color} 0.30 0.08 {color} 0.95')"
        )
    lines.append("python end")

    # Build axis vectors in PyMOL via embedded Python
    # Non-linker regions are used for mobile pi-axis
    non_linker_regions = [r for r in region_order if r != "linker"]
    pi_sel_parts = " or ".join(f"{mobile_name}_{r}" for r in non_linker_regions)
    if not pi_sel_parts:
        pi_sel_parts = mobile_name

    lines.extend(
        [
            "python",
            "import numpy as np",
            "from pymol import cmd",
            "",
            "def _coords(selection):",
            "    m = cmd.get_model(selection)",
            "    if not m.atom:",
            "        return None",
            "    return np.array([[a.coord[0], a.coord[1], a.coord[2]] for a in m.atom], dtype=float)",
            "",
            "def _principal_axis(coords, mode='long'):",
            "    center = coords.mean(axis=0)",
            "    x = coords - center",
            "    cov = x.T @ x / max(1, x.shape[0])",
            "    evals, evecs = np.linalg.eigh(cov)",
            "    idx = int(np.argmax(evals)) if mode == 'long' else int(np.argmin(evals))",
            "    axis = evecs[:, idx]",
            "    n = np.linalg.norm(axis)",
            "    if n < 1e-12:",
            "        return center, None",
            "    return center, axis / n",
            "",
            "fixed_coords = _coords('system and not elem H')",
            f"pi_coords = _coords('({pi_sel_parts}) and not elem H')",
            f"if pi_coords is None:\n    pi_coords = _coords('{mobile_name} and not elem H')",
            "",
            "if fixed_coords is not None and fixed_coords.shape[0] >= 3:",
            "    f_center, f_axis = _principal_axis(fixed_coords, mode='short')",
            "    if f_axis is not None:",
            "        f_len = 18.0",
            "        f1 = f_center - f_axis * f_len",
            "        f2 = f_center + f_axis * f_len",
            "        cmd.pseudoatom('fixed_axis_start', pos=list(f1))",
            "        cmd.pseudoatom('fixed_axis_end', pos=list(f2))",
            "        cmd.distance('fixed_axis_vec', 'fixed_axis_start', 'fixed_axis_end')",
            "        cmd.color('orange', 'fixed_axis_vec')",
            "        cmd.set('dash_width', 4.0, 'fixed_axis_vec')",
            "        cmd.set('dash_gap', 0.0, 'fixed_axis_vec')",
            "",
            "if pi_coords is not None and pi_coords.shape[0] >= 3:",
            "    p_center, p_axis = _principal_axis(pi_coords, mode='long')",
            "    if p_axis is not None:",
            "        p_len = 14.0",
            "        p1 = p_center - p_axis * p_len",
            "        p2 = p_center + p_axis * p_len",
            "        cmd.pseudoatom('mobile_axis_start', pos=list(p1))",
            "        cmd.pseudoatom('mobile_axis_end', pos=list(p2))",
            "        cmd.distance('mobile_axis_vec', 'mobile_axis_start', 'mobile_axis_end')",
            "        cmd.color('red', 'mobile_axis_vec')",
            "        cmd.set('dash_width', 4.0, 'mobile_axis_vec')",
            "        cmd.set('dash_gap', 0.0, 'mobile_axis_vec')",
            "",
            "python end",
        ]
    )

    lines.append("orient")
    lines.append(
        f'png "{analysis_dir / "view_regions.png"}", width=1800, height=1200, dpi=220, ray=1'
    )

    output_pml.write_text("\n".join(lines) + "\n")
    print(f"Wrote {output_pml}")


if __name__ == "__main__":
    main()
