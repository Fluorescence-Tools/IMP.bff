"""The committed shader header must be what the .wgsl says.

`gpu/diffusion_wgsl.h` carries the compute shader as a C string literal so
that a built `libimp_bff_wgpu` needs no data file beside it, but the shader is
*authored* as `gpu/diffusion.wgsl` so it can be read and edited as WGSL. Two
copies of the same thing drift, so this fails when they do -- the same
arrangement, and the same reason, as `test_vendored_headers.py`.

Regenerate with:

    python utility/wgsl_to_header.py gpu/diffusion.wgsl gpu/diffusion_wgsl.h
"""

import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
GENERATOR = ROOT / "utility" / "wgsl_to_header.py"
# (shader, generated header, symbol)
SHADERS = [("diffusion", "kDiffusionWgsl"), ("mlp", "kMlpWgsl")]
WGSL = ROOT / "gpu" / "diffusion.wgsl"
HEADER = ROOT / "gpu" / "diffusion_wgsl.h"


def _generator():
    spec = importlib.util.spec_from_file_location("wgsl_to_header", GENERATOR)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.mark.skipif(not WGSL.is_file(), reason="no gpu/ in this tree")
@pytest.mark.parametrize("stem,symbol", SHADERS)
def test_the_header_is_what_the_shader_says(stem, symbol):
    wgsl = ROOT / "gpu" / ("%s.wgsl" % stem)
    header = ROOT / "gpu" / ("%s_wgsl.h" % stem)
    expected = _generator().generate(str(wgsl), symbol)
    assert header.read_text(encoding="utf-8") == expected, (
        "gpu/%s_wgsl.h is out of date; regenerate it with "
        "python utility/wgsl_to_header.py gpu/%s.wgsl gpu/%s_wgsl.h"
        % (stem, stem, stem))


@pytest.mark.skipif(not WGSL.is_file(), reason="no gpu/ in this tree")
def test_the_shader_declares_what_the_plugin_binds():
    """The C sets eight bindings by number. A shader that renumbered them
    would still compile and would read the wrong buffers."""
    text = WGSL.read_text(encoding="utf-8")
    for i, name in enumerate(["cur", "nxt", "w", "self_", "part", "p", "trace", "idx"]):
        assert "@binding(%d) var<%s" % (i, "uniform" if name == "p" else "storage") in text \
            or "@binding(%d) var<storage" % i in text, "binding %d moved" % i
        assert name in text
    for entry in ("fn sweep", "fn reduce_a", "fn reduce_b"):
        assert entry in text, "%s is named in the plugin" % entry


@pytest.mark.skipif(not WGSL.is_file(), reason="no gpu/ in this tree")
def test_the_network_shader_covers_every_activation():
    """`Activation` has seven members and the shader dispatches on their
    declaration order. One missing would return the identity silently."""
    text = (ROOT / "gpu" / "mlp.wgsl").read_text(encoding="utf-8")
    for code in range(7):
        assert "act == %du" % code in text, "activation %d is not handled" % code
    header = (ROOT / "include" / "internal" / "MlpCore.h").read_text(encoding="utf-8")
    declared = header.split("enum class Activation {")[1].split("}")[0]
    assert len([x for x in declared.split(",") if x.strip()]) == 7, (
        "MlpCore.h grew an activation; gpu/mlp.wgsl dispatches on seven")
