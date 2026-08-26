"""`attach_probes` takes `ProbeAttachment` values, and `bin/imp_bff` now does.

Three commands called it with `(hierarchy, chain, residue)` tuples and one with
four parallel lists -- both shapes older than the C++ port, and both a
`TypeError` at the binding. `dye sample-rotamer` then read `attached[0]["site"]`
from a return that is a list of values, not of dicts. None of it had a test, so
none of it was known to be dead.
"""

import pytest

import IMP
import IMP.atom
import IMP.core
import IMP.bff


@pytest.fixture(scope="module")
def protein_and_probe():
    m = IMP.Model()
    protein = IMP.atom.read_pdb(
        IMP.bff.get_example_path("structure/T4L/3GUN.pdb"), m,
        IMP.atom.NonWaterNonHydrogenPDBSelector())
    probe = IMP.atom.read_mol2(
        str(IMP.bff.get_structure_dir("alexa488_r48.mol2")), m)
    return m, protein, probe


def test_attach_probes_takes_values_not_tuples(protein_and_probe):
    m, protein, probe = protein_and_probe
    with pytest.raises(TypeError):
        IMP.bff.attach_probes(protein, [(probe, "A", 132)], True)

    out = IMP.bff.attach_probes(
        protein, [IMP.bff.ProbeAttachment(probe, "A", 132)], True)
    assert len(out) == 1
    assert out[0].get_chain() == "A"
    assert out[0].get_residue() == 132
    assert out[0].get_n_stripped() > 0, "the side chain was to be stripped"


def test_the_site_comes_from_resolve_probe_site(protein_and_probe):
    """What the command needed out of the attachment, and where it is."""
    m, protein, probe = protein_and_probe
    site = IMP.bff.resolve_probe_site(protein, "A", 132)
    assert [IMP.atom.Atom(p).get_atom_type().get_string() for p in site] == \
        ["CA", "N", "C"]


def test_dye_sample_rotamer_runs(imp_bff_program, tmp_path):
    """The command was dead: it raised at `attach_probes` before it could get
    to reading a dict out of a list of values."""
    from click.testing import CliRunner
    out = tmp_path / "sampled.rmf3"
    result = CliRunner().invoke(imp_bff_program.dye, [
        "sample-rotamer",
        "--protein-pdb", IMP.bff.get_example_path("structure/T4L/3GUN.pdb"),
        "--chain", "A", "--residue", "132", "--dye", "alexa488",
        "--n-samples", "3", "--output-rmf", str(out)])
    assert result.exit_code == 0, result.output
    assert out.is_file()

    import RMF
    fh = RMF.open_rmf_file_read_only(str(out))
    assert fh.get_number_of_frames() == 3
