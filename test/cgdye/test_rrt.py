#!/usr/bin/env python3

import sys
from pathlib import Path
from IMP.bff.cgdye.utils import get_template_dir, get_structure_dir


from IMP.bff.cgdye.sampling.rrt import (
    IMPRRTTree,
    make_transform,
    transform_distance,
    is_collision_sphere,
    rrt_grow_step,
    run_imp_rrt,
    sample_transform,
    steer_transform,
)


def test_tree_init_and_add():
    t = IMPRRTTree(make_transform(0, 0, 0, 0, 0, 0))
    n1 = t.add_node(make_transform(1, 0, 0, 0, 0, 0), 0)
    assert n1 == 1
    assert t.nodes[1].parent_id == 0


def test_tree_nearest():
    t = IMPRRTTree(make_transform(0, 0, 0, 0, 0, 0))
    t.add_node(make_transform(5, 0, 0, 0, 0, 0), 0)
    nid = t.nearest(make_transform(4.9, 0, 0, 0, 0, 0))
    assert nid == 1


def test_config_distance_positive():
    d = transform_distance(
        make_transform(0, 0, 0, 0, 0, 0),
        make_transform(1, 2, 2, 0, 0, 0),
    )
    assert d > 0


def test_steer_limits_step_size():
    a = make_transform(0, 0, 0, 0, 0, 0)
    b = make_transform(10, 0, 0, 0, 0, 0)
    s = steer_transform(a, b, 1.0)
    tx = s.get_translation()[0]
    assert tx > 0
    assert tx < 10


def test_collision_checker():
    cfg = make_transform(0, 0, 0, 0, 0, 0)
    obstacles = [(0.5, 0, 0, 0.6)]
    assert is_collision_sphere(cfg, obstacles, 0.2)


def test_rrt_grow_step_noncolliding():
    t = IMPRRTTree(make_transform(0, 0, 0, 0, 0, 0))
    nid = rrt_grow_step(t, make_transform(1, 0, 0, 0, 0, 0), 0.5, lambda _c: False)
    assert nid is not None
    assert nid == 1


def test_run_rrt_reaches_goal_in_empty_space():
    start = make_transform(0, 0, 0, 0, 0, 0)
    goal = make_transform(1, 0, 0, 0, 0, 0)
    bounds = [(-2, 2), (-2, 2), (-2, 2), (-1, 1), (-1, 1), (-1, 1)]
    tree, goal_id = run_imp_rrt(
        start,
        bounds,
        n_iter=200,
        step_size=0.25,
        collision_fn=lambda _c: False,
        goal_transform=goal,
        goal_bias=0.3,
        goal_tolerance=0.35,
        seed=7,
    )
    assert goal_id is not None
    path = tree.path_to_root(goal_id)
    assert len(path) >= 2


def test_run_rrt_respects_collision():
    start = make_transform(0, 0, 0, 0, 0, 0)
    bounds = [(-2, 2), (-2, 2), (-2, 2), (-1, 1), (-1, 1), (-1, 1)]

    def coll(cfg):
        return cfg.get_translation()[0] > 0.2

    tree, _goal_id = run_imp_rrt(
        start,
        bounds,
        n_iter=100,
        step_size=0.5,
        collision_fn=coll,
        goal_transform=make_transform(1, 0, 0, 0, 0, 0),
        goal_bias=0.5,
        goal_tolerance=0.1,
        seed=4,
    )
    for n in tree.nodes.values():
        assert n.transform.get_translation()[0] <= 0.2 + 1e-9


def test_run_torsion_rrt_grows_collision_free_tree():
    from IMP.bff.cgdye.sampling.rrt import run_torsion_rrt, torsion_distance
    import math

    # forbid the half-space where the first torsion is > pi/2 (a "wall")
    def collides(cfg):
        return cfg[0] > math.pi / 2

    tree, goal_id = run_torsion_rrt(3, 200, 0.3, collides, seed=1)
    assert len(tree) > 1
    assert all(c[0] <= math.pi / 2 for c in tree.configs)
    # every non-root node is within one step of its parent
    for i in range(1, len(tree)):
        assert torsion_distance(tree.configs[i], tree.configs[tree.parents[i]]) <= 0.3 + 1e-9

    tree, goal_id = run_torsion_rrt(2, 500, 0.5, lambda c: False,
                                    goal_cfg=[1.0, -1.0], goal_bias=0.5, seed=2)
    assert goal_id is not None
    assert torsion_distance(tree.configs[goal_id], [1.0, -1.0]) <= 0.1



# IMP runs every .py under test/ as a standalone script, and a file of bare
# pytest functions would import cleanly and exit 0 -- reporting success without
# running a single assertion. Hand the file to pytest explicitly so a failure
# here is a failure in ctest.
if __name__ == "__main__":
    import sys
    try:
        import pytest
    except ImportError:
        print("pytest not installed; skipping", __file__)
        sys.exit(0)
    sys.exit(pytest.main([__file__, "-q", "-p", "no:cacheprovider"]))
