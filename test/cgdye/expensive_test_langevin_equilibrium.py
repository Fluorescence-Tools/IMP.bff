"""md and bd sample the same equilibrium of the attached dye at hGBP1 481 (PRD-108 stage 3).

Both integrators run on the same energy function; the distribution of the
dye's heavy-atom centroid distance from the site CA must agree between them
(loose bounds: the linker equilibrates on the ns scale).
"""

import numpy as np
import pytest

import IMP
import IMP.atom
import IMP.core

from IMP.bff.label.attachment import attach_dyes
from IMP.bff.cgdye.sampling.langevin import LangevinDyeSampler
from IMP.bff.tools.paths import get_structure_dir


def _centroid_distances(integrator, n_steps, write_every, seed):
    m = IMP.Model()
    prot = IMP.atom.read_pdb(str(get_structure_dir("1DG3.pdb")), m, IMP.atom.NonWaterPDBSelector())
    dye = IMP.atom.read_mol2(str(get_structure_dir("alexa488_r48.mol2")), m)
    attach_dyes(prot, [(dye, "A", 481)], strip_site_sidechain=True)
    s = LangevinDyeSampler(prot, dye, str(get_structure_dir("alexa488_r48.mol2")), "A", 481, integrator=integrator, seed=seed)
    s.minimize(300)
    traj = s.run(n_steps, write_every=write_every)
    heavy = [i for i, n in enumerate(s.atom_names) if not n.upper().startswith("H")]
    centroids = traj.coordinates[:, heavy, :].mean(axis=1)
    return np.linalg.norm(centroids - s.site_ca, axis=1), traj


def test_md_and_bd_equilibrium_of_the_dye_centroid():
    d_md, traj_md = _centroid_distances("md", 600000, 200, 1)     # 1.2 ns
    d_bd, traj_bd = _centroid_distances("bd", 800000, 400, 2)     # 0.4 ns
    # discard the first quarter as equilibration
    d_md, d_bd = d_md[len(d_md) // 4:], d_bd[len(d_bd) // 4:]
    assert np.isfinite(traj_md.coordinates).all() and np.isfinite(traj_bd.coordinates).all()
    assert abs(d_md.mean() - d_bd.mean()) < 3.0, (d_md.mean(), d_bd.mean())
    edges = np.linspace(0, 30, 16)
    h_md, _ = np.histogram(d_md, bins=edges, density=True)
    h_bd, _ = np.histogram(d_bd, bins=edges, density=True)
    ks = np.abs(np.cumsum(h_md) - np.cumsum(h_bd)).max() * (edges[1] - edges[0])
    assert ks < 0.3, ks


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-q", "-p", "no:cacheprovider"]))
