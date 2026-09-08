"""The plugin's copy of the compute contract must match the library's.

`gpu/imp_bff_wgpu.c` restates `ImpBffComputeBackend` and the propagate
signature because that header is C++ and the plugin is C -- a backend may be
built by a different compiler, and a std::vector across that line is a promise
neither side can keep. Two declarations of one ABI drift, and the failure mode
is not a compile error: it is a plugin that loads and reads the wrong
arguments off the stack.
"""

import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
HEADER = ROOT / "include" / "Compute.h"
PLUGIN = ROOT / "gpu" / "imp_bff_wgpu.c"

pytestmark = pytest.mark.skipif(not PLUGIN.is_file(), reason="no gpu/ in this tree")


def _abi(path):
    m = re.search(r"#define\s+IMPBFF_COMPUTE_BACKEND_ABI\s+(\d+)",
                  path.read_text(encoding="utf-8"))
    assert m, "no IMPBFF_COMPUTE_BACKEND_ABI in %s" % path.name
    return int(m.group(1))


def _fields(text, struct):
    body = re.search(r"struct\s+%s\s*\{(.*?)\}" % struct, text, re.S)
    assert body, struct
    return [ln.split("//")[0].strip().rstrip(";").split()[-1].lstrip("*")
            for ln in body.group(1).splitlines()
            if ln.strip() and not ln.strip().startswith(("//", "/*", "*", "!"))]


def test_the_abi_numbers_agree():
    assert _abi(HEADER) == _abi(PLUGIN)


def test_the_struct_has_the_same_fields_in_the_same_order():
    a = _fields(HEADER.read_text(encoding="utf-8"), "ImpBffComputeBackend")
    b = _fields(PLUGIN.read_text(encoding="utf-8"), "ImpBffComputeBackend")
    assert a == b, "%s vs %s" % (a, b)


def test_the_propagate_signature_has_the_same_parameters():
    def params(text):
        m = re.search(r"ImpBffPropagateFn\)\((.*?)\);", text, re.S)
        assert m
        return [p.strip().split()[-1].lstrip("*") for p in m.group(1).split(",")]
    assert params(HEADER.read_text(encoding="utf-8")) == \
           params(PLUGIN.read_text(encoding="utf-8"))


def test_both_module_recipes_build_the_plugin():
    """The plugin is one C file that links nothing and is found by path at run
    time, so IMP's module tooling never has to build it -- but somebody must.
    For the wheel and the `bff` conda package that is
    `standalone/CMakeLists.txt`; for the `imp.bff` module package it is the
    recipe, in one compiler command. Dropping that line would leave the module
    package silently CPU-only, which is exactly the kind of absence a plugin
    is designed to make invisible."""
    for script in ("conda-recipe/build.sh", "conda-recipe/bld.bat"):
        text = (ROOT / script).read_text(encoding="utf-8")
        assert "imp_bff_wgpu.c" in text, "%s no longer builds the backend" % script
    cmake = (ROOT / "standalone" / "CMakeLists.txt").read_text(encoding="utf-8")
    assert "IMPBFF_WITH_GPU" in cmake and 'add_subdirectory("${IMPBFF_ROOT}/gpu"' in cmake


def test_the_build_says_whether_it_can_load_a_backend_at_all():
    """`IMPBFF_WITH_GPU=0` compiles the loader out -- no dlopen anywhere in the
    library. A build that lost the door must say so rather than look like a
    machine with no GPU, because the two want different things done about
    them."""
    import IMP.bff
    if IMP.bff.built_with_gpu_support():
        return                                   # the other tests cover this
    assert IMP.bff.get_compute_backend_name() == "cpu"
    assert not IMP.bff.load_compute_backend("/nonexistent", "")
    assert "IMPBFF_WITH_GPU=0" in IMP.bff.get_compute_backend_error()
    assert IMP.bff.enable_gpu() == "cpu"
