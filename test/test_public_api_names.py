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
    # `IMP.bff.observables` and `IMP.bff.io.fps` are both "in the observables
    # / io domain"; a domain that is one flat module *is* the module, so an
    # exact match counts as much as a prefix. Before the consolidation only the
    # prefix form existed and this read `startswith`.
    module_form = f"IMP.bff.{domain}"
    package_form = f"IMP.bff.{domain}."
    wrong = {name: module for name, module in api.BY_DOMAIN[domain].items()
             if module != module_form and not module.startswith(package_form)}
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
        "import IMP.bff.label as l; l.attach_dyes\n"
        "print('ok')\n"
    )
    result = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, timeout=180)
    assert result.returncode == 0, result.stderr
    assert "ok" in result.stdout


def test_cgdye_is_off_the_public_surface():
    """It is explicit-dye molecular mechanics, kept but not a domain.

    ``cgdye`` used to mirror its slice of the flat map. As of 2026-08-18 it has
    no flat names at all: propagating an all-atom dye under a force field is
    IMP's territory, not a spectroscopy library's, so the package stays
    importable by module path and off the surface. What was *not* molecular
    mechanics was taken out of it first -- the rotamer library into
    ``representation``, attachment into ``label``, its readers into ``io``.
    """
    import IMP.bff.cgdye as cgdye

    assert cgdye.__all__ == []
    assert not [m for m in api.EXPORTS.values() if m.startswith("IMP.bff.cgdye.")]
    assert "cgdye" not in api.BY_DOMAIN

    # still reachable by module path -- kept, not deleted
    import importlib
    assert importlib.import_module("IMP.bff.cgdye.sampling").LangevinDyeSampler

    # and the harvested pieces landed where they were supposed to
    assert api.domain_of("RotamerFRET") == "representation"
    assert api.domain_of("attach_dyes") == "label"
    assert api.domain_of("read_rotamer_library_rmf") == "io"


def test_nothing_in_the_package_imports_cgdye_at_module_scope():
    """No domain may *load* cgdye. If one grows a module-level edge into it,
    that is the signal something was filed in the wrong place -- as the LJ term
    and the rotamer library both were.

    Module scope, not every import: ``scoring/dye_lj.py`` reaches into
    ``cgdye.topology.dye`` from inside a function for the molecular graph, since
    an exclusion list is derived from connectivity. That is a genuine
    dependency of scoring on topology and it is deferred, so importing a domain
    never pulls cgdye in. A module-level edge would.
    """
    import ast
    from pathlib import Path

    src = Path(__file__).resolve().parent.parent / "pyext" / "src"
    offenders = []
    for path in sorted(src.rglob("*.py")):
        if "cgdye" in path.parts:
            continue
        tree = ast.parse(path.read_text())
        names = set()
        for node in tree.body:                       # body, not walk
            if isinstance(node, ast.ImportFrom) and node.module:
                names.add(node.module)
            elif isinstance(node, ast.Import):
                names |= {a.name for a in node.names}
        if any(n.startswith("IMP.bff.cgdye") for n in names):
            offenders.append(str(path.relative_to(src)))
    assert not offenders, offenders


def test_importing_a_domain_does_not_load_cgdye():
    """The property the test above is a proxy for, checked directly."""
    code = (
        "import sys\n"
        "import IMP.bff.scoring, IMP.bff.representation, IMP.bff.label\n"
        "import IMP.bff.observables, IMP.bff.io, IMP.bff.sampling\n"
        "loaded = [m for m in sys.modules if m.startswith('IMP.bff.cgdye')]\n"
        "assert not loaded, loaded\n"
        "print('ok')\n"
    )
    result = subprocess.run([sys.executable, "-c", code], capture_output=True,
                            text=True, timeout=300)
    assert result.returncode == 0, result.stderr[-2000:]
    assert "ok" in result.stdout


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q", "-p", "no:cacheprovider"]))
