"""Guards that outlived the source tree they audited.

This file's substance was an audit of ``pyext/src`` -- the subpackage import
graph, its cycles, its private reaching. That tree is gone: the shims are
deleted and every public name is an attribute of ``IMP.bff`` itself. Two tests
here never were about it, so they stayed; the rest went with the tree.
"""

import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pytest

import IMP.bff


def test_test_module_basenames_are_unique():
    """pytest imports test modules by basename, so a collision silently drops one.

    ``test/io/test_io.py`` collided with ``test/cgdye/test_io.py`` and was
    reported as a collection ERROR at the very bottom of a 780-test run -- easy
    to read past, and the file's assertions simply did not run.
    """
    root = Path(__file__).resolve().parent
    by_name = defaultdict(list)
    for path in sorted(root.rglob("test_*.py")):
        if "__pycache__" in path.parts:
            continue
        by_name[path.name].append(str(path.relative_to(root)))
    clashes = {name: paths for name, paths in by_name.items() if len(paths) > 1}
    assert not clashes, clashes


def test_the_array_door_computes_an_accessible_volume():
    """``compute_av`` end to end -- four atoms and a source between them is
    enough: what is pinned is that the flat C++ door is callable at all."""
    av = IMP.bff.compute_av(
        np.column_stack([
            np.array([[0.0, 0.0, 0.0], [6.0, 0.0, 0.0], [0.0, 6.0, 0.0],
                      [0.0, 0.0, 6.0]]),
            np.array([1.7, 1.7, 1.7, 1.7]),
        ]),
        np.array([2.0, 2.0, 2.0]),
    )
    assert av.get_density().reshape((av.get_ng(),) * 3).ndim == 3
    assert av.get_points().reshape(-1, 4).shape[0] > 0
    assert av.get_mean_position().shape == (3,)


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q", "-p no:cacheprovider"]))
