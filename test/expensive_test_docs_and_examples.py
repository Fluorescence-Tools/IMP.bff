"""Every example and every manual notebook, executed.

Neither is run by the ordinary suite, and both broke quietly as the C++ port
changed shapes under them: `attach_probes` took values instead of tuples,
`resolve_probe_site` answered a list instead of a dict, a keep-set became a
function, a sampler stopped writing files, decay convolution moved to tttrlib.
Nine dead call sites, none of which a unit test would have found, because what
was dead was the *documentation*.

It is an `expensive_test_` -- not collected by default, run deliberately:

    pytest test/expensive_test_docs_and_examples.py -q

A notebook runs with **its own directory as the working directory**, which is
what Jupyter does and what the relative paths in them assume.
"""

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
#: Long-running by nature -- a replica-exchange run over a dimer.
SLOW = {"mGBP2_dimer_fp.py"}
#: The TCSPC instrument layer was deleted from this package on purpose
#: (`PRD-113 stage 0`), and these still document it: `IMP.bff.DecayCurve`,
#: `DecayConvolution`, `DecayLifetimeHandler`, `DecayScale`, `DecayPattern`,
#: `DecayLinearization`. They are not broken call sites but pages describing a
#: layer that moved out; skipped until they are ported or dropped.
DELETED_DECAY_LAYER = {
    "decay_curves.ipynb",
    "decay_forward_model.ipynb",
    "decay_objective_function.ipynb",
    "programming_imp_decorator.ipynb",
}

CELL_RUNNER = r'''
import json, sys
nb = json.load(open(sys.argv[1]))
for i, c in enumerate(nb["cells"]):
    if c["cell_type"] != "code":
        continue
    src = "".join(c["source"])
    if src.strip().startswith("!") or src.startswith("%%"):
        continue
    exec(compile(src, "cell%d" % i, "exec"), globals())
'''


def _examples():
    return sorted(p for p in (ROOT / "examples").glob("*/*.py")
                  if p.name not in SLOW)


def _notebooks():
    return sorted(p for p in (ROOT / "doc").glob("**/*.ipynb")
                  if p.name not in DELETED_DECAY_LAYER)


def _run(args, cwd):
    env = dict(os.environ, MPLBACKEND="Agg")
    return subprocess.run(args, capture_output=True, text=True, timeout=600,
                          cwd=str(cwd), env=env)


@pytest.mark.parametrize("script", _examples(), ids=lambda p: p.name)
def test_example_runs(script):
    r = _run([sys.executable, str(script)], cwd=ROOT)
    assert r.returncode == 0, (r.stderr or "")[-2000:]


@pytest.mark.parametrize("notebook", _notebooks(), ids=lambda p: p.name)
def test_notebook_runs(notebook, tmp_path):
    runner = tmp_path / "run_cells.py"
    runner.write_text(CELL_RUNNER)
    r = _run([sys.executable, str(runner), notebook.name], cwd=notebook.parent)
    assert r.returncode == 0, (r.stderr or "")[-2000:]
