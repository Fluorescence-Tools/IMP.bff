"""The public surface: domain-scoped, resolvable, documented, and correctly filed.

This replaced a naming-families regex. That test asserted every export's *name*
matched one of ``Rotamer|Dye|Linker|Langevin|...`` -- a list that had to grow
with every feature, and which could only ever say that a name looked plausible.
It could not tell a misfiled name from a well-filed one, and it grew by four
tokens in the six months before this was written.

What is checked instead is structural, and does not grow:

* every export resolves and carries a docstring;
* every export's module lies **inside the domain it is filed under**, so the
  grouping in ``api.py`` is enforced rather than decorative;
* no name is exported by two domains;
* every domain in the map is a real subpackage.
"""

import subprocess
import sys

import pytest

import IMP.bff
import IMP.bff.api as api


@pytest.mark.parametrize("name", api.public_names())
def test_export_resolves_and_is_documented(name):
    value = api.resolve(name)
    assert value is not None
    if not isinstance(value, (dict, tuple, list)):
        assert getattr(value, "__doc__", None), f"{name} has no docstring"
    # reachable flat, and cached after first access
    assert getattr(IMP.bff, name) is value
    assert name in dir(IMP.bff)


@pytest.mark.parametrize("domain", sorted(api.BY_DOMAIN))
def test_every_export_lives_in_the_domain_it_is_filed_under(domain):
    """The check the naming regex could not make.

    A name filed under ``fret`` whose module is ``IMP.bff.quenching.model`` is
    a mistake, and it is exactly the kind of mistake that accumulates when the
    grouping is a comment.
    """
    prefix = f"IMP.bff.{domain}."
    wrong = {name: module for name, module in api.BY_DOMAIN[domain].items()
             if not module.startswith(prefix)}
    assert not wrong, wrong


def test_the_flat_view_is_derived_and_complete():
    flat = {}
    for members in api.BY_DOMAIN.values():
        flat.update(members)
    assert flat == api.EXPORTS


def test_no_name_is_exported_by_two_domains():
    seen = {}
    clashes = {}
    for domain, members in api.BY_DOMAIN.items():
        for name in members:
            if name in seen:
                clashes[name] = (seen[name], domain)
            seen[name] = domain
    assert not clashes, clashes


def test_every_domain_is_a_real_subpackage():
    import importlib
    for domain in api.BY_DOMAIN:
        importlib.import_module(f"IMP.bff.{domain}")


def test_domain_of_agrees_with_the_module_path():
    for name, module in api.EXPORTS.items():
        assert module.split(".")[2] == api.domain_of(name), name


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
