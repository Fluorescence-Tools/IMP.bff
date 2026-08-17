"""The flat ``IMP.bff.*`` surface: every export resolves, is documented, and says what it is."""

import re
import subprocess
import sys

import pytest

import IMP.bff
import IMP.bff.api as api

# a name is self-explanatory when it carries one of these tokens
_FAMILY = re.compile(
    r"(Rotamer|Dye|Linker|Langevin|AccessibleVolume|CHARMM36|SITE_KEEP"
    r"|rotamer|dye|linker|langevin|fret|fps|forster|kappa2|strip|backbone|_av\b|_avs_|^av_"
    r"|boltzmann|cluster|rrt|SIMULATION_TYPE|hierarch|mol2|cif|rmf|dcd|protein_frames|nmr|lj_|select_atoms|attach)"
)


@pytest.mark.parametrize("name", api.public_names())
def test_export_resolves_and_is_documented(name):
    value = api.resolve(name)
    assert value is not None
    if not isinstance(value, (dict, tuple, list)):
        assert getattr(value, "__doc__", None), f"{name} has no docstring"
    assert _FAMILY.search(name), f"{name} does not say what it is (naming families)"
    # reachable flat, and cached after first access
    assert getattr(IMP.bff, name) is value
    assert name in dir(IMP.bff)


def test_no_retired_names_survive():
    retired = [
        "System", "read_ff_system", "write_ff_system", "read_cgdye_template",
        "write_cgdye_template", "resolve_site", "apply_rotamer_coords",
        "generate_rotamers", "compute_exact_efficiency", "calculate_fret_exact",
        "calculate_fret_regimes", "calculate_r0", "kappa2_from_vectors", "run_imp_rrt",
        "compute_boltzmann_weights", "mean_field_weights", "build_transition_probability_matrix",
        "load_reference_rotamers", "InternalEnergyEvaluator", "build_combined_system",
        "parse_mol2", "analyze_mobile", "LJ_PARAMETERS",
    ]
    for name in retired:
        assert name not in api.EXPORTS
        with pytest.raises(AttributeError):
            getattr(IMP.bff, name)
    import IMP.bff.cgdye as cgdye
    for name in ("System", "resolve_site", "generate_rotamers"):
        assert not hasattr(cgdye, name)


def test_import_is_lazy_and_click_free():
    code = (
        "import sys; sys.modules['click'] = None\n"
        "import IMP.bff\n"
        "assert 'IMP.bff.cgdye' not in sys.modules, 'import IMP.bff must not import cgdye'\n"
        "IMP.bff.RotamerFRET; IMP.bff.forster_radius_from_spectra; IMP.bff.strip_hierarchy\n"
        "import IMP.bff.cgdye as c; c.attach_dyes\n"
        "print('ok')\n"
    )
    result = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, timeout=180)
    assert result.returncode == 0, result.stderr
    assert "ok" in result.stdout


def test_cgdye_mirror_matches_api():
    import IMP.bff.cgdye as cgdye
    expected = sorted(n for n, m in api.EXPORTS.items() if m.startswith("IMP.bff.cgdye."))
    assert cgdye.__all__ == expected
    for name in expected[:5]:
        assert getattr(cgdye, name) is api.resolve(name)


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q", "-p", "no:cacheprovider"]))
