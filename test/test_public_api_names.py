"""The public surface: flat, resolvable, documented, and shim-free.

Every public name is an attribute of ``IMP.bff`` itself, wrapped from a C++
declaration, so what there is to check is structural:

* every public name resolves, carries a docstring, and is in ``dir()``;
* nothing is built lazily -- a missing attribute is a missing attribute;
* the retired names stay retired;
* ``import IMP.bff`` stays click-free and pulls no optional dependency;
* no ``IMP.bff.<submodule>`` package survives: the flat namespace is the only
  surface.
"""

import subprocess
import sys

import pytest

import IMP.bff

#: The RMF readers and writers. `rmf` is one of this module's modules
#: (`dependencies.py`), so these are ordinary C++ attributes -- nothing here
#: is built on first use, and this list checks exactly that.
RMF_NAMES = [
    "write_rmf",
    "write_rotamer_library_rmf",
    "read_rotamer_library_rmf",
    "protein_frames_from_rmf",
    "get_anchor_cb_position",
]

#: A sample of the always-defined flat surface, one per source `.i` group.
FLAT_NAMES = [
    "GraphNode",                  # graph runtime
    "get_av",                    # avbuilder
    "AccessibleVolume",              # avmodel
    "read_fps_json",                 # fps
    "read_component_template_cif",   # cif
    "forcefield_system_from_json",   # forcefield / ProbeForceField
    "get_rotamer_score",         # scoring (pythoncode, flat)
    "attach_probes",                   # label (pythoncode, flat)
    "create_probe_protein_system",      # topology (C++, TopologyBuild.h)
    "probe_forcefield_system",         # topology (C++, TopologyBuild.h)
    "AttachedProbeDynamics",            # sampling (pythoncode, flat)
    "RotamerEnsemble",               # rotamer_ensemble (C++, Rotamer.h)
    "RotamerFRET",                   # rotamer_ensemble (C++, Rotamer.h)
    "load_rotamer_library",          # rotamer (C++, RotamerSite.h)
    # `build_system_from_specs` was here. It was the `build-system` command's
    # body -- build a system, write the CIF, print what it holds -- so it is
    # in `bin/imp_bff` with the command, not in the library. What the library
    # offers is `create_forcefield_system`, which takes `FFComponentSpec`
    # values rather than the CLI's `name=X,mol2=Y` strings.
    "create_forcefield_system",       # topology (C++, TopologyBuild.h)
    "simulate_probe_diffusion",        # probesampling
    "fret_rate_trace",               # quenching / FRETRateTrace
]

#: Submodule paths that must not exist: the namespace is flat.
RETIRED_SUBMODULES = [
    "IMP.bff.api",
    "IMP.bff.io",
    "IMP.bff.io.cif",
    "IMP.bff.label",
    "IMP.bff.scoring",
    "IMP.bff.representation",
    "IMP.bff.representation.rotamer",
    "IMP.bff.restraints",
    "IMP.bff.cgprobe",
    "IMP.bff.cgprobe.topology",
    "IMP.bff.cgprobe.sampling",
    "IMP.bff.cgprobe.sim",
]


@pytest.mark.parametrize("name", FLAT_NAMES)
def test_flat_name_resolves_and_is_documented(name):
    value = getattr(IMP.bff, name)
    assert value is not None
    if not isinstance(value, (dict, tuple, list)):
        assert getattr(value, "__doc__", None), f"{name} has no docstring"
    assert name in dir(IMP.bff)


@pytest.mark.parametrize("name", RMF_NAMES)
def test_the_rmf_names_are_ordinary_attributes(name):
    value = getattr(IMP.bff, name)
    assert value is not None
    assert getattr(IMP.bff, name) is value
    assert name in dir(IMP.bff)


def test_there_is_no_lazy_machinery():
    """No lazy-name table, and no module `__getattr__` to consult one."""
    assert not hasattr(IMP.bff, "_LAZY")
    assert "__getattr__" not in vars(IMP.bff)


def test_no_retired_names_survive():
    retired = [
        # The graph runtime moved out of chinet's generic vocabulary.
        "Node", "Port", "Session", "EvaluationGraph", "Expression",
        "System", "read_ff_system", "write_ff_system", "read_cgprobe_template",
        "write_cgprobe_template", "resolve_site", "apply_rotamer_coords",
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
        # Writers with no caller. Templates and numpy/text rotamer libraries
        # are read here and authored elsewhere -- by hand, by `.drot`, or by
        # RMF -- so writing them back was API that only its own tests used.
        "write_component_template_cif", "write_probe_template_cif",
        "write_rotamer_library",
        # The probe vocabulary is flrCIF's: `_flr_probe_list`, not "dye".
        "DYE_PAIR_DISTANCE_MEAN", "DYE_PAIR_EFFICIENCY", "LL_DYE_CBETA",
        "LL_DYE_ACCESSIBLE_VOLUME", "DEFAULT_DYE_RADIUS",
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

    `RMF` and `IMP.rmf` are not optional -- `rmf` is one of this module's
    `required_modules`, because `RmfIO.h` reads and writes RMF in C++. IMP.pmi
    is optional and stays optional: a `RestraintBase` subclass has no C++
    spelling, so a PMI wrapper is a program's business.
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
