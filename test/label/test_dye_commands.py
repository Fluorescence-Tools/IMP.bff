"""The `dye` commands, run.

Every one of these was dead, and none of them had a test -- which is the only
reason a command can be dead. What each was:

* `sample-dof-walk` imported `LinkerSampler`, a name from before the linker
  sampler became C++; then it read `resolve_probe_site`'s answer as a dict;
  then, once running, it rejected every proposal, because it counted the bond
  between the probe and its neighbouring residues as a clash.
* `sample-langevin` passed `timestep_fs=None` into a `double`, and asked
  `run()` to write an RMF -- which it stopped doing when the sampler became
  C++, because writing files is a program's business.
* `label-fusion` used `sys.exit` in a module that never imported `sys`.

They are smoke tests: a few steps each, checking the command runs and writes
what it says it wrote. What the samplers *compute* is pinned elsewhere
(`test/cgdye/test_langevin_sampler.py`).
"""

import pytest
from click.testing import CliRunner

import IMP.bff

PDB = IMP.bff.get_example_path("structure/T4L/3GUN.pdb")


@pytest.fixture(scope="module")
def dye(imp_bff_program):
    return imp_bff_program.dye


@pytest.fixture(scope="module")
def rotamer(imp_bff_program):
    return imp_bff_program.rotamer


def test_sample_dof_walk_accepts_something(dye, tmp_path):
    """A walk that accepts nothing is not a sampler.

    The probe is placed *bonded into* the site, so its first atoms sit a bond
    length from residues i-1 and i+1; counting those as clashes rejected every
    proposal. They are excluded now, and a proposal that does not make the
    clash count worse is accepted -- otherwise a walk that starts inside the
    protein could never leave.
    """
    out = tmp_path / "walk.rmf3"
    r = CliRunner().invoke(dye, [
        "sample-dof-walk", "--protein-pdb", PDB, "--chain", "A",
        "--residue", "132", "--dye", "alexa488", "--n-steps", "40",
        "--output-rmf", str(out)])
    assert r.exit_code == 0, r.output
    assert out.is_file()
    line = [l for l in r.output.splitlines() if "Finished" in l][-1]
    accepted = int(line.split()[1].split("/")[0])
    assert accepted > 0, line


def test_sample_langevin_runs_and_writes_frames(dye, tmp_path):
    out = tmp_path / "md.rmf3"
    r = CliRunner().invoke(dye, [
        "sample-langevin", "--protein-pdb", PDB, "--chain", "A",
        "--residue", "132", "--dye", "alexa488", "--n-steps", "100",
        "--write-every", "50", "--minimize-steps", "10",
        "--output-rmf", str(out)])
    assert r.exit_code == 0, r.output
    assert out.is_file()

    import RMF
    fh = RMF.open_rmf_file_read_only(str(out))
    assert fh.get_number_of_frames() == 2


def test_the_timestep_is_chosen_when_it_is_not_given(dye, tmp_path):
    """A negative timestep means "one per integrator": 2 fs for md, 0.5 for
    bd. It used to be `None`, which is not a `double`."""
    out = tmp_path / "bd.rmf3"
    r = CliRunner().invoke(dye, [
        "sample-langevin", "--protein-pdb", PDB, "--chain", "A",
        "--residue", "132", "--dye", "alexa488", "--integrator", "bd",
        "--n-steps", "50", "--write-every", "50", "--minimize-steps", "0",
        "--output-rmf", str(out)])
    assert r.exit_code == 0, r.output
    assert "0.5 fs" in r.output, r.output


def test_label_writes_a_labelled_structure(dye, tmp_path):
    out = tmp_path / "labelled.pdb"
    r = CliRunner().invoke(dye, [
        "label", PDB, "--chain", "A", "--residue", "132",
        "--dye", "alexa488", "--output", str(out)])
    assert r.exit_code == 0, r.output
    assert out.is_file()
    assert out.stat().st_size > 0


def test_label_fusion_runs(dye, tmp_path):
    """`sys.exit` in a module that never imported `sys` is a `NameError` on
    the way out, not an exit."""
    out = tmp_path / "fusion.pdb"
    r = CliRunner().invoke(dye, [
        "label-fusion", PDB, "--chain", "A", "--output", str(out)])
    assert r.exit_code == 0, r.output


def test_rotamer_r0_reports_a_forster_radius(rotamer):
    r = CliRunner().invoke(rotamer, [
        "r0", "--donor", "AlexaFluor 488", "--acceptor", "AlexaFluor 594",
        "--k2", "0.6667"])
    assert r.exit_code == 0, r.output
    value = float(r.output.split()[0])
    assert 40.0 < value < 70.0, r.output
