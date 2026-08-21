"""The AV's grid resolution comes from ``disc_step``, and a disagreeing
``source_info`` is refused rather than discarded.

``compute_av_from_structure`` takes a position definition and honours its AV
fields -- ``allowed_sphere_radius``, ``strip_mask``, ``contact_volume_*``. It does
**not** honour ``simulation_grid_resolution``: that field is *written* into the
particle from the ``disc_step`` argument, so a caller who states
the resolution in ``source_info`` and leaves ``disc_step`` at its default
silently gets 1.5 A.

That is not a hypothetical. Both quenching identifiability benchmarks
(PRD-110, PRD-111) passed the resolution this way, so their ``--resolution`` flag
was inert and every run was made at 1.5 A while the recorded validation pages
said 2.5 A. The numbers were internally consistent -- only the label was wrong --
but a resolution being wrong is invisible in the result, which is exactly what
makes it worth refusing loudly.
"""

import numpy as np
import pytest

import IMP.bff

_SOURCE = dict(
    chain_identifier="A", residue_seq_number=132, atom_name="CB",
    simulation_type="AV1", linker_length=20.0, linker_width=0.5, radius1=3.5,
    allowed_sphere_radius=2.1,
)

@pytest.fixture(scope="module")
def pdb_path():
    return IMP.bff.get_example_path("structure/T4L/3GUN.pdb")


def _density(av):
    ng = av.get_ng()
    return np.asarray(av.get_density()).reshape(ng, ng, ng)


@pytest.mark.parametrize("disc_step", [1.5, 2.0, 2.5])
def test_disc_step_sets_the_grid_spacing(pdb_path, disc_step):
    av = IMP.bff.compute_av_from_structure(pdb_path, dict(_SOURCE),
                                           disc_step=disc_step)
    assert av.get_grid_step() == pytest.approx(disc_step)


def test_a_coarser_grid_is_a_smaller_array(pdb_path):
    fine = IMP.bff.compute_av_from_structure(pdb_path, dict(_SOURCE), disc_step=1.5)
    coarse = IMP.bff.compute_av_from_structure(pdb_path, dict(_SOURCE), disc_step=2.5)
    assert np.prod(_density(coarse).shape) < np.prod(_density(fine).shape)
    # Same cloud, so the physical extent survives the coarsening.
    fine_extent = np.array(_density(fine).shape) * fine.get_grid_step()
    coarse_extent = np.array(_density(coarse).shape) * coarse.get_grid_step()
    assert coarse_extent == pytest.approx(fine_extent, rel=0.25)


def test_an_agreeing_declaration_is_accepted(pdb_path):
    av = IMP.bff.compute_av_from_structure(
        pdb_path, dict(_SOURCE, simulation_grid_resolution=2.0), disc_step=2.0)
    assert av.get_grid_step() == pytest.approx(2.0)


def test_a_disagreeing_declaration_raises(pdb_path):
    with pytest.raises(ValueError, match="simulation_grid_resolution"):
        IMP.bff.compute_av_from_structure(
            pdb_path, dict(_SOURCE, simulation_grid_resolution=2.5), disc_step=1.5)


def test_the_default_is_not_silently_imposed_on_a_declared_resolution(pdb_path):
    """The exact shape of the bug: declare 2.5, pass no disc_step, get 1.5."""
    with pytest.raises(ValueError, match="disc_step=1.5"):
        IMP.bff.compute_av_from_structure(
            pdb_path, dict(_SOURCE, simulation_grid_resolution=2.5))


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q", "-p", "no:cacheprovider"]))
