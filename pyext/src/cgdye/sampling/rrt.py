"""RRT sampling core utilities."""

from dataclasses import dataclass
import math
import random


@dataclass
class RRTNode:
    node_id: int
    config: tuple
    parent_id: int | None


class RRTTree:
    def __init__(self, root_config):
        self.nodes = {0: RRTNode(0, tuple(root_config), None)}
        self.next_id = 1

    def add_node(self, config, parent_id):
        nid = self.next_id
        self.nodes[nid] = RRTNode(nid, tuple(config), parent_id)
        self.next_id += 1
        return nid

    def nearest(self, config):
        best_id = None
        best_d = float("inf")
        for nid, n in self.nodes.items():
            d = config_distance(n.config, config)
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


def config_distance(a, b, rot_weight=0.25):
    # config = (tx, ty, tz, rx, ry, rz)
    dt = math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2)
    dr = math.sqrt((a[3] - b[3]) ** 2 + (a[4] - b[4]) ** 2 + (a[5] - b[5]) ** 2)
    return dt + rot_weight * dr


def sample_config(bounds, rng):
    out = []
    for lo, hi in bounds:
        out.append(rng.uniform(lo, hi))
    return tuple(out)


def steer(from_cfg, to_cfg, step_size):
    d = config_distance(from_cfg, to_cfg)
    if d == 0:
        return tuple(from_cfg)
    if d <= step_size:
        return tuple(to_cfg)
    alpha = step_size / d
    return tuple(from_cfg[i] + alpha * (to_cfg[i] - from_cfg[i]) for i in range(6))


def is_collision_sphere(config, obstacles, probe_radius):
    # obstacles: list[(x,y,z,r)] ; only translation part for collision
    x, y, z = config[0], config[1], config[2]
    for ox, oy, oz, orad in obstacles:
        d2 = (x - ox) ** 2 + (y - oy) ** 2 + (z - oz) ** 2
        min_d = probe_radius + orad
        if d2 < min_d * min_d:
            return True
    return False


def rrt_grow_step(tree, target_cfg, step_size, collision_fn):
    near_id = tree.nearest(target_cfg)
    near_cfg = tree.nodes[near_id].config
    new_cfg = steer(near_cfg, target_cfg, step_size)
    if collision_fn(new_cfg):
        return None
    return tree.add_node(new_cfg, near_id)


def run_rrt(
    start_config,
    bounds,
    n_iter,
    step_size,
    collision_fn,
    goal_config=None,
    goal_bias=0.1,
    goal_tolerance=1.0,
    seed=0,
):
    rng = random.Random(seed)
    tree = RRTTree(start_config)
    goal_node_id = None

    for _ in range(n_iter):
        if goal_config is not None and rng.random() < goal_bias:
            target = goal_config
        else:
            target = sample_config(bounds, rng)
        new_id = rrt_grow_step(tree, target, step_size, collision_fn)
        if new_id is None:
            continue
        if goal_config is not None:
            d = config_distance(tree.nodes[new_id].config, goal_config)
            if d <= goal_tolerance:
                goal_node_id = new_id
                break

    return tree, goal_node_id
