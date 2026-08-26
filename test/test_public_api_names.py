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

#: There are no lazily-built names any more. These four were built on first
#: use, through a `_LAZY` table and a PEP 562 `__getattr__`, because naming
#: them at import time would have made `import IMP.bff` require `IMP.rmf`.
#: `rmf` is one of this module's modules now and they are C++ (`RmfIO.h`), so
#: what this list checks is that they are ordinary attributes like any other.
FORMERLY_LAZY_NAMES = [
    "write_rmf",
    "write_rotamer_library_rmf",
    "read_rotamer_library_rmf",
    "protein_frames_from_rmf",
    "get_anchor_cb_position",
]

#: A sample of the always-defined flat surface, one per source `.i` group.
FLAT_NAMES = [
    "compute_av",                    # avbuilder
    "AccessibleVolume",              # avmodel
    "read_fps_json",                 # fps
    "read_component_template_cif",   # cif
    "forcefield_system_from_json",   # forcefield / DyeForceField
    "compute_rotamer_score",         # scoring (pythoncode, flat)
    "attach_probes",                   # label (pythoncode, flat)
    "build_dye_protein_system",      # topology (C++, TopologyBuild.h)
    "probe_forcefield_system",         # topology (C++, TopologyBuild.h)
    "AttachedProbeDynamics",            # sampling (pythoncode, flat)
    "RotamerEnsemble",               # rotamer_ensemble (C++, RotamerEnsemble.h)
    "RotamerFRET",                   # rotamer_ensemble (C++, RotamerFret.h)
    "load_rotamer_library",          # rotamer (C++, RotamerSite.h)
    # `build_system_from_specs` was here. It was the `build-system` command's
    # body -- build a system, write the CIF, print what it holds -- so it is
    # in `bin/imp_bff` with the command, not in the library. What the library
    # offers is `build_forcefield_system`, which takes `FFComponentSpec`
    # values rather than the CLI's `name=X,mol2=Y` strings.
    "build_forcefield_system",       # topology (C++, TopologyBuild.h)
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


@pytest.mark.parametrize("name", FORMERLY_LAZY_NAMES)
def test_the_formerly_lazy_names_are_ordinary_attributes(name):
    value = getattr(IMP.bff, name)
    assert value is not None
    assert getattr(IMP.bff, name) is value
    assert name in dir(IMP.bff)


def test_there_is_no_lazy_machinery_left():
    """No `_LAZY` table, and no module `__getattr__` to consult it."""
    assert not hasattr(IMP.bff, "_LAZY")
    assert "__getattr__" not in vars(IMP.bff)


def test_no_retired_names_survive():
    retired = [
        "System", "read_ff_system", "write_ff_system", "read_cgdye_template",
        "write_cgdye_template", "resolve_site", "apply_rotamer_coords",
        "generate_rotamers", "compute_exact_efficiency", "calculate_fret_exact",
        "calculate_fret_regimes", "calculate_r0", "kappa2_from_vectors", "run_imp_rrt",
        "compute_boltzmann_weights", "mean_field_weights", "build_transition_probability_matrix",
        "load_reference_rotamers", "InternalEnergyEvaluator", "build_combined_system",
        "parse_mol2", "analyze_mobile", "LJ_PARAMETERS",
        # `AVNetworkRestraint` is `ProbeNetworkRestraint`, and the PMI wrapper
        # around it is gone entirely -- a `RestraintBase` subclass has no C++
        # spelling, and `probe_network_restraint_set` is what it did besides
        # PMI's bookkeeping.
        "AVNetworkRestraint", "AVNetworkRestraintWrapper",
        "SimpleAVNetworkRestraint", "av_network_restraint_set",
    ]
    for name in retired:
        with pytest.raises(AttributeError):
            getattr(IMP.bff, name)


def test_import_is_lazy_and_click_free():
    code = (
        "import sys; sys.modules['click'] = None\n"
        "import IMP.bff\n"
        "IMP.bff.RotamerFRET; IMP.bff.forster_radius_from_spectra; IMP.bff.strip_hierarchy\n"
        "IMP.bff.attach_probes\n"
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
    """`IMP.pmi` is the optional one: blocking it must not break the import.

    `RMF` and `IMP.rmf` are no longer on that list -- `rmf` is one of this
    module's `required_modules`, because `RmfIO.h` reads and writes RMF in
    C++ where four `%pythoncode` builders used to. IMP.pmi is still optional,
    and stays optional: a `RestraintBase` subclass has no C++ spelling, so
    what used to be `AVNetworkRestraintWrapper` is a program's business now.
    """
    code = (
        "import sys\n"
        "sys.modules['IMP.pmi'] = None\n"
        "import IMP.bff\n"
        "print('ok')\n"
    )
    result = subprocess.run([sys.executable, "-c", code], capture_output=True,
                            text=True, timeout=300)
    assert result.returncode == 0, result.stderr[-2000:]
    assert "ok" in result.stdout


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q", "-p no:cacheprovider"]))
