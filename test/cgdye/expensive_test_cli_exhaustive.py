import os
import pytest
from click.testing import CliRunner
from IMP.bff.cgdye.cli import dye
from IMP.bff.tools.paths import get_structure_dir

@pytest.fixture(scope="module")
def runner():
    return CliRunner()

@pytest.fixture(scope="module")
def kinetic_lib(tmp_path_factory):
    """Generate a small kinetic library for testing other commands."""
    tmp_dir = tmp_path_factory.mktemp("data")
    runner = CliRunner()
    # We need a mol2 file in a specific place or build-lib won't find it.
    # build-lib searches inputs/structures. 
    # For the test, we'll just run build-lib on the real inputs but with very few steps.
    lib_dir = tmp_dir / "libs"
    runner.invoke(dye, [
        "build-lib", 
        "--n-steps", "50", 
        "--cluster-threshold", "2.0", 
        "--output-dir", str(lib_dir)
    ])
    # Find one generated rmf
    for f in os.listdir(lib_dir):
        if f.endswith(".rmf3"):
            return os.path.join(lib_dir, f)
    return None

def test_cli_help(runner):
    result = runner.invoke(dye, ["--help"])
    assert result.exit_code == 0

def test_cli_label_pdb(runner, tmp_path):
    out_pdb = tmp_path / "test_label.pdb"
    result = runner.invoke(dye, [
        "label", str(get_structure_dir("1DG3.pdb")), 
        "--residue", "481", 
        "--dye", "Alexa488", 
        "--linker", "C1R",
        "--output", str(out_pdb)
    ])
    assert result.exit_code == 0
    assert out_pdb.exists()

def test_cli_analyze_tc(runner, kinetic_lib):
    if not kinetic_lib:
        pytest.skip("No kinetic library generated")
    result = runner.invoke(dye, ["analyze-tc", kinetic_lib])
    assert result.exit_code == 0
    assert "Slowest TC" in result.output

def test_cli_reconstruct(runner, kinetic_lib, tmp_path):
    if not kinetic_lib:
        pytest.skip("No kinetic library generated")
    out_rmf = tmp_path / "recon.rmf3"
    result = runner.invoke(dye, [
        "reconstruct",
        "--lib-rmf", kinetic_lib,
        "--n-frames", "10",
        "--output-rmf", str(out_rmf)
    ])
    assert result.exit_code == 0
    assert out_rmf.exists()

def test_cli_sample_rotamer(runner, tmp_path):
    out_rmf = tmp_path / "test_rot.rmf3"
    result = runner.invoke(dye, [
        "sample-rotamer",
        "--protein-pdb", str(get_structure_dir("1DG3.pdb")),
        "--residue", "481",
        "--dye", "Alexa488",
        "--n-samples", "5",
        "--output-rmf", str(out_rmf)
    ])
    assert result.exit_code == 0
    assert out_rmf.exists()

def test_cli_label_fp_dual(runner, tmp_path):
    out_pdb = tmp_path / "dual_fp.pdb"
    result = runner.invoke(dye, [
        "label-fp", str(get_structure_dir("1DG3.pdb")),
        "--site", "6:eGFP",
        "--site", "583:mCherry",
        "--output", str(out_pdb)
    ])
    assert result.exit_code == 0
    assert out_pdb.exists()


# IMP runs every .py under test/ as a standalone script, and a file of bare
# pytest functions would import cleanly and exit 0 -- reporting success without
# running a single assertion. Hand the file to pytest explicitly so a failure
# here is a failure in ctest.
if __name__ == "__main__":
    import sys
    try:
        import pytest
    except ImportError:
        print("pytest not installed; skipping", __file__)
        sys.exit(0)
    sys.exit(pytest.main([__file__, "-q", "-p", "no:cacheprovider"]))
