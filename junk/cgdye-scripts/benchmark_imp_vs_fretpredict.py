#!/usr/bin/env python
"""Performance benchmark: IMP placement vs reference numpy transform."""

import time
import json
import os
import sys
from pathlib import Path
import numpy as np

import click
import IMP
import IMP.atom
import IMP.core
import IMP.algebra



from IMP.bff.cgdye.labeling.attachment import attach_dyes, place_dye_from_coords

def _structure(name):
    """A bundled input structure, wherever IMP.bff is installed.

    These paths used to be relative to the working directory, so a script only
    ran from one place -- and stopped running at all once the data became IMP
    module data under data/cgdye.
    """
    from IMP.bff.cgdye.utils import get_structure_dir
    return str(get_structure_dir(name))



def _numpy_transform(ca, n, c, coords):
    """Reference numpy-based transform (backbone frame)."""
    origin = ca
    v1 = n - origin
    v1 /= np.linalg.norm(v1)
    v2 = c - origin
    v2 -= np.dot(v2, v1) * v1
    v2 /= np.linalg.norm(v2)
    v3 = np.cross(v1, v2)
    # Rotation matrix
    R = np.vstack([v1, v2, v3]).T
    return np.dot(coords, R.T) + origin


@click.command()
@click.option("--n-trials", default=1000, show_default=True)
@click.option("--output-json", default="output/benchmarks/placement_speed.json")
def main(n_trials, output_json):
    model = IMP.Model()
    protein = IMP.atom.read_pdb(_structure("1DG3.pdb"), model, IMP.atom.NonWaterPDBSelector())
    dye = IMP.atom.read_pdb(_structure("alexa488_r48.pdb"), model, IMP.atom.AllPDBSelector())

    # Get CA, N, C coordinates for hGBP1 residue 481
    sel = IMP.atom.Selection(protein, chain_id="A", residue_index=481)
    ca_p = IMP.atom.Selection(sel, atom_type=IMP.atom.AtomType("CA")).get_selected_particles()[0]
    n_p = IMP.atom.Selection(sel, atom_type=IMP.atom.AtomType("N")).get_selected_particles()[0]
    c_p = IMP.atom.Selection(sel, atom_type=IMP.atom.AtomType("C")).get_selected_particles()[0]

    ca_v = IMP.core.XYZ(ca_p).get_coordinates()
    n_v = IMP.core.XYZ(n_p).get_coordinates()
    c_v = IMP.core.XYZ(c_p).get_coordinates()

    ca_np = np.array([ca_v[0], ca_v[1], ca_v[2]])
    n_np = np.array([n_v[0], n_v[1], n_v[2]])
    c_np = np.array([c_v[0], c_v[1], c_v[2]])

    dye_atoms = IMP.atom.get_by_type(dye, IMP.atom.ATOM_TYPE)
    base_coords = np.array([IMP.core.XYZ(a).get_coordinates() for a in dye_atoms])

    # Benchmark IMP placement
    t0 = time.time()
    for _ in range(n_trials):
        place_dye_from_coords(dye, ca_v, n_v, c_v)
    t_imp = (time.time() - t0) / n_trials

    # Benchmark Numpy placement
    t0 = time.time()
    for _ in range(n_trials):
        _ = _numpy_transform(ca_np, n_np, c_np, base_coords)
    t_np = (time.time() - t0) / n_trials

    results = {
        "n_trials": n_trials,
        "imp_placement_sec": t_imp,
        "numpy_placement_sec": t_np,
        "speedup_factor": t_np / t_imp if t_imp > 0 else 0
    }

    os.makedirs(os.path.dirname(os.path.abspath(output_json)), exist_ok=True)
    with open(output_json, "w") as fh:
        json.dump(results, fh, indent=2)

    click.echo(f"Benchmark finished. IMP: {t_imp:.6f}s, Numpy: {t_np:.6f}s")
    click.echo(f"Wrote {output_json}")


if __name__ == "__main__":
    main()
