"""Rapidly-exploring random trees over rigid placements and torsions.

The planners are C++ (`IMP/bff/RRT.h`) and the collision test crosses back
into Python through a director: a caller subclasses `RRTCollision` and the
growth loop asks it about every candidate step.

**What these tests no longer cover, and why.** The Python had unit tests for
`steer_transform`, `is_collision_sphere`, `rrt_grow_step` and `tree.nearest` --
helpers that are now internal to the growth loop. What those tests were really
asserting is asserted here at the level that matters: no step is longer than
`step_size`, no node sits in a collision, and a reachable goal is reached.
A caller that wants a sphere-obstacle check writes it in its own
`RRTCollision`; it was eight lines of arithmetic, not an algorithm.
"""

import math

import pytest

import IMP.algebra
from IMP.bff import (
    RRTCollision,
    configuration_to_transformation,
    grow_rigid_body_rrt,
    grow_torsion_rrt,
    transformation_to_configuration,
)


class _Free(RRTCollision):
    """Nothing ever collides."""

    def is_collision(self, configuration):
        return False


class _Wall(RRTCollision):
    """Everything beyond `limit` along the first coordinate collides."""

    def __init__(self, limit):
        RRTCollision.__init__(self)
        self.limit = limit

    def is_collision(self, configuration):
        return configuration[0] > self.limit


def _placement(tx, ty, tz, rx=0.0, ry=0.0, rz=0.0):
    """A rigid placement as the six numbers the planner works in."""
    return [tx, ty, tz, rx, ry, rz]


def _torsion_distance(a, b):
    def wrap(x):
        return (x + math.pi) % (2.0 * math.pi) - math.pi
    return math.sqrt(sum(wrap(x - y) ** 2 for x, y in zip(a, b)))


def test_a_placement_survives_the_round_trip():
    """Six numbers in, the same placement out."""
    want = [1.0, -2.0, 3.0, 0.3, -0.2, 0.1]
    got = transformation_to_configuration(configuration_to_transformation(want))
    assert got == pytest.approx(want, abs=1e-9)


def test_torsion_tree_grows_and_never_enters_the_wall():
    tree = grow_torsion_rrt(3, 200, 0.3, _Wall(math.pi / 2), seed=1)
    assert tree.n_nodes > 1
    for i in range(tree.n_nodes):
        assert tree.get_configuration(i)[0] <= math.pi / 2

    # every node is within one step of the node it grew from
    for i in range(1, tree.n_nodes):
        parent = tree.parents[i]
        assert _torsion_distance(tree.get_configuration(i),
                                 tree.get_configuration(parent)) <= 0.3 + 1e-9


def test_torsion_tree_reaches_a_goal_in_open_space():
    goal = [1.0, -1.0]
    tree = grow_torsion_rrt(2, 500, 0.5, _Free(), goal=goal, goal_bias=0.5,
                            goal_tolerance=0.1, seed=2)
    assert tree.goal_node >= 0
    assert _torsion_distance(tree.get_configuration(tree.goal_node), goal) <= 0.1


def test_rigid_tree_reaches_a_goal_and_reports_the_path():
    start = _placement(0, 0, 0)
    goal = _placement(1, 0, 0)
    low, high = [-2, -2, -2, -1, -1, -1], [2, 2, 2, 1, 1, 1]
    tree = grow_rigid_body_rrt(start, low, high, 200, 0.25, _Free(), goal=goal,
                               goal_bias=0.3, goal_tolerance=0.35, seed=7)
    assert tree.goal_node >= 0
    path = list(tree.get_path_to_root(tree.goal_node))
    assert len(path) >= 2
    assert path[0] == 0                      # a path starts where the tree did


def test_rigid_tree_respects_the_collision_check():
    start = _placement(0, 0, 0)
    low, high = [-2, -2, -2, -1, -1, -1], [2, 2, 2, 1, 1, 1]
    tree = grow_rigid_body_rrt(start, low, high, 100, 0.5, _Wall(0.2),
                               goal=_placement(1, 0, 0),
                               goal_bias=0.5, goal_tolerance=0.1, seed=4)
    for i in range(tree.n_nodes):
        assert tree.get_configuration(i)[0] <= 0.2 + 1e-9


def test_a_goal_out_of_reach_is_reported_as_not_reached():
    """No goal node is an answer, not an error: the tree still explored."""
    start = _placement(0, 0, 0)
    low, high = [-1, -1, -1, 0, 0, 0], [1, 1, 1, 0, 0, 0]
    tree = grow_rigid_body_rrt(start, low, high, 20, 0.05, _Free(),
                               goal=_placement(50, 0, 0),
                               goal_tolerance=0.1, seed=3)
    assert tree.goal_node == -1
    assert tree.n_nodes > 1


# IMP runs every .py under test/ as a standalone script, and a file of bare
# pytest functions would import cleanly and exit 0 -- reporting success without
# running a single assertion.
if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-q", "-p", "no:cacheprovider"]))
