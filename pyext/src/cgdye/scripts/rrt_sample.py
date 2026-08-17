#!/usr/bin/env python

import json
import os

import click

from IMP.bff.cgdye.sampling.rrt import make_transform, run_imp_rrt


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option("--n-iter", default=500, show_default=True, type=int)
@click.option("--step-size", default=0.5, show_default=True, type=float)
@click.option("--goal-bias", default=0.2, show_default=True, type=float)
@click.option("--seed", default=0, show_default=True, type=int)
@click.option("--output", default="output/rrt/rrt_result.json", show_default=True)
def main(n_iter, step_size, goal_bias, seed, output):
    start = make_transform(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    goal = make_transform(5.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    bounds = [
        (-10.0, 10.0),
        (-10.0, 10.0),
        (-10.0, 10.0),
        (-3.14, 3.14),
        (-3.14, 3.14),
        (-3.14, 3.14),
    ]

    tree, goal_id = run_imp_rrt(
        start,
        bounds,
        n_iter=n_iter,
        step_size=step_size,
        collision_fn=lambda _cfg: False,
        goal_transform=goal,
        goal_bias=goal_bias,
        goal_tolerance=0.5,
        seed=seed,
    )

    out = {
        "n_nodes": len(tree.nodes),
        "goal_reached": goal_id is not None,
        "goal_node_id": goal_id,
    }
    out_dir = os.path.dirname(os.path.abspath(output))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with open(output, "w") as fh:
        json.dump(out, fh, indent=2)
    click.echo(f"Wrote {output}")


if __name__ == "__main__":
    main()
