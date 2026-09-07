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
    # Every notebook under `doc/` runs. Nothing is skipped: a page this module
    # ships is a page this module can execute, or it does not ship.
    return sorted((ROOT / "doc").glob("**/*.ipynb"))


def _run(args, cwd):
    # MPLBACKEND: a figure window would block the run.
    # FI_PROVIDER: anything that pulls in IMP.mpi calls MPI_Finalize at exit, and
    # libfabric's default provider selection picks whichever NIC comes first --
    # on a machine with a VPN up that is a `utun`, and finalize aborts flushing
    # its send queue ("OFI poll failed"), killing an example that had already
    # finished its work. tcp over loopback is what a single-rank run wants
    # anyway.
    env = dict(os.environ, MPLBACKEND="Agg", FI_PROVIDER="tcp")
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
