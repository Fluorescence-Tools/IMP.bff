"""IMP-native RRT sampling utilities."""

from dataclasses import dataclass
import math
import random

import IMP.algebra


@dataclass
class IMPNode:
    node_id: int
    transform: IMP.algebra.Transformation3D
    parent_id: int | None


class IMPRRTTree:
    def __init__(self, root_transform):
        self.nodes = {0: IMPNode(0, root_transform, None)}
        self.next_id = 1

    def add_node(self, transform, parent_id):
        nid = self.next_id
        self.nodes[nid] = IMPNode(nid, transform, parent_id)
        self.next_id += 1
        return nid

    def nearest(self, transform, rot_weight=0.25):
        best_id = None
        best_d = float("inf")
        for nid, n in self.nodes.items():
            d = transform_distance(n.transform, transform, rot_weight=rot_weight)
            if d < best_d:
                best_d = d
                best_id = nid
        return best_id

    def path_to_root(self, node_id):
        path = []
        cur = node_id
        while cur is not None:
            path.append(cur)
            cur = self.nodes[cur].parent_id
        return list(reversed(path))


def make_transform(tx, ty, tz, rx, ry, rz):
    rot = IMP.algebra.get_rotation_from_fixed_xyz(rx, ry, rz)
    tr = IMP.algebra.Vector3D(tx, ty, tz)
    return IMP.algebra.Transformation3D(rot, tr)


def transform_to_tuple(t):
    tr = t.get_translation()
    e = IMP.algebra.get_fixed_xyz_from_rotation(t.get_rotation())
    return (tr[0], tr[1], tr[2], e.get_x(), e.get_y(), e.get_z())


def transform_distance(a, b, rot_weight=0.25):
    at = a.get_translation()
    bt = b.get_translation()
    dt = math.sqrt((at[0] - bt[0]) ** 2 + (at[1] - bt[1]) ** 2 + (at[2] - bt[2]) ** 2)

    ae = IMP.algebra.get_fixed_xyz_from_rotation(a.get_rotation())
    be = IMP.algebra.get_fixed_xyz_from_rotation(b.get_rotation())
    dr = math.sqrt(
        (ae.get_x() - be.get_x()) ** 2
        + (ae.get_y() - be.get_y()) ** 2
        + (ae.get_z() - be.get_z()) ** 2
    )
    return dt + rot_weight * dr


def sample_transform(bounds, rng):
    vals = [rng.uniform(lo, hi) for lo, hi in bounds]
    return make_transform(*vals)


def steer_transform(from_t, to_t, step_size, rot_weight=0.25):
    d = transform_distance(from_t, to_t, rot_weight=rot_weight)
    if d == 0:
        return from_t
    if d <= step_size:
        return to_t
    alpha = step_size / d
    fa = transform_to_tuple(from_t)
    ta = transform_to_tuple(to_t)
    vals = [fa[i] + alpha * (ta[i] - fa[i]) for i in range(6)]
    return make_transform(*vals)


def is_collision_sphere(transform, obstacles, probe_radius):
    tr = transform.get_translation()
    x, y, z = tr[0], tr[1], tr[2]
    for ox, oy, oz, orad in obstacles:
        d2 = (x - ox) ** 2 + (y - oy) ** 2 + (z - oz) ** 2
        md = probe_radius + orad
        if d2 < md * md:
            return True
    return False


def rrt_grow_step(tree, target_t, step_size, collision_fn, rot_weight=0.25):
    near_id = tree.nearest(target_t, rot_weight=rot_weight)
    near_t = tree.nodes[near_id].transform
    new_t = steer_transform(near_t, target_t, step_size, rot_weight=rot_weight)
    if collision_fn(new_t):
        return None
    return tree.add_node(new_t, near_id)


def run_imp_rrt(
    start_transform,
    bounds,
    n_iter,
    step_size,
    collision_fn,
    goal_transform=None,
    goal_bias=0.1,
    goal_tolerance=1.0,
    rot_weight=0.25,
    seed=0,
):
    rng = random.Random(seed)
    tree = IMPRRTTree(start_transform)
    goal_node_id = None

    for _ in range(n_iter):
        if goal_transform is not None and rng.random() < goal_bias:
            target = goal_transform
        else:
            target = sample_transform(bounds, rng)

        new_id = rrt_grow_step(
            tree, target, step_size, collision_fn, rot_weight=rot_weight
        )
        if new_id is None:
            continue
        if goal_transform is not None:
            d = transform_distance(
                tree.nodes[new_id].transform, goal_transform, rot_weight=rot_weight
            )
            if d <= goal_tolerance:
                goal_node_id = new_id
                break

    return tree, goal_node_id
