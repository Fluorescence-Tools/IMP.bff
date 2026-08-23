"""The public surface: flat, resolvable, documented, and shim-free.

This replaced a naming-families regex, then a domain-map audit; both were about
a structure that no longer exists. Every public name is an attribute of
``IMP.bff`` itself now -- a SWIG wrapper or a `%pythoncode` def in the one
generated ``__init__.py`` -- so what is left to check is structural:

* every public name resolves, carries a docstring, and is in ``dir()``;
* the lazy names (IMP.pmi/IMP.rmf wrappers) build on first use and cache;
* the retired names stay retired;
* ``import IMP.bff`` stays click-free and pulls no optional dependency;
* no ``IMP.bff.<submodule>`` package survives -- the shims are gone and the
  flat namespace is the only surface.
"""

import subprocess
import sys

import pytest

import IMP.bff

#: Names built on first use because defining them at import time would make
#: `import IMP.bff` require a module it does not depend on.
LAZY_NAMES = [
    "AVNetworkRestraintWrapper",
    "write_rmf",
    "write_rotamer_library_rmf",
    "read_rotamer_library_rmf",
]

#: A sample of the always-defined flat surface, one per source `.i` group.
FLAT_NAMES = [
    "compute_av",                    # avbuilder
    "AccessibleVolume",              # avmodel
    "read_fps_json",                 # fps
    "read_component_template_cif",   # cif
    "forcefield_system_from_json",   # forcefield / DyeForceField
    "compute_rotamer_score",         # scoring (pythoncode, flat)
    "attach_dyes",                   # label (pythoncode, flat)
    "build_dye_protein_system",      # topology (pythoncode, flat)
    "LangevinDyeSampler",            # sampling (pythoncode, flat)
    "RotamerEnsemble",               # rotamer_ensemble (pythoncode, flat)
    "build_system_from_specs",       # sim/topology (pythoncode, flat)
    "simulate_dye_diffusion",        # dyesampling
    "fret_rate_trace",               # quenching / FRETRateTrace
]

#: The module paths that used to be the domain packages. None may exist.
RETIRED_SUBMODULES = [
    "IMP.bff.api",
    "IMP.bff.io",
    "IMP.bff.io.cif",
    "IMP.bff.label",
    "IMP.bff.scoring",
    "IMP.bff.representation",
    "IMP.bff.representation.rotamer",
    "IMP.bff.restraints",
    "IMP.bff.cgdye",
    "IMP.bff.cgdye.topology",
    "IMP.bff.cgdye.sampling",
    "IMP.bff.cgdye.sim",
]


@pytest.mark.parametrize("name", FLAT_NAMES)
def test_flat_name_resolves_and_is_documented(name):
    value = getattr(IMP.bff, name)
    assert value is not None
    if not isinstance(value, (dict, tuple, list)):
        assert getattr(value, "__doc__", None), f"{name} has no docstring"
    assert name in dir(IMP.bff)


@pytest.mark.parametrize("name", LAZY_NAMES)
def test_lazy_name_builds_and_caches(name):
    value = getattr(IMP.bff, name)
    assert value is not None
    # cached: the second access is the same object, not a rebuild
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
        with pytest.raises(AttributeError):
            getattr(IMP.bff, name)


def test_import_is_lazy_and_click_free():
    code = (
        "import sys; sys.modules['click'] = None\n"
        "import IMP.bff\n"
        "IMP.bff.RotamerFRET; IMP.bff.forster_radius_from_spectra; IMP.bff.strip_hierarchy\n"
        "IMP.bff.attach_dyes\n"
        "print('ok')\n"
    )
    result = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, timeout=180)
    assert result.returncode == 0, result.stderr
    assert "ok" in result.stdout


def test_the_domain_packages_are_gone():
    """The shims under pyext/src are deleted; the flat namespace is the surface.

    A submodule reappearing is a regression: it forks the name space again,
    which is exactly what the domain layout was retired to end.
    """
    import importlib.util
    for module in RETIRED_SUBMODULES:
        try:
            spec = importlib.util.find_spec(module)
        except ModuleNotFoundError:      # a missing parent is the same miss
            spec = None
        assert spec is None, f"{module} exists again"


def test_import_works_with_the_optional_dependencies_blocked():
    """RMF / IMP.pmi are optional: blocking them must not break the import.

    sim.i guards them in try/except at module scope, so they load *when
    present* -- the property that must hold either way is that the import
    survives without them.
    """
    code = (
        "import sys\n"
        "for m in ('RMF', 'IMP.pmi', 'IMP.rmf'):\n"
        "    sys.modules[m] = None\n"
        "import IMP.bff\n"
        "print('ok')\n"
    )
    result = subprocess.run([sys.executable, "-c", code], capture_output=True,
                            text=True, timeout=300)
    assert result.returncode == 0, result.stderr[-2000:]
    assert "ok" in result.stdout


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q", "-p no:cacheprovider"]))
