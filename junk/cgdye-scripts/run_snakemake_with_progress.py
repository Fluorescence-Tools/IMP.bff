#!/usr/bin/env python
"""Run Snakemake and render a compact progress bar."""

import re
import subprocess
import sys

import click


PROGRESS_RE = re.compile(r"(\d+) of (\d+) steps \(([0-9.]+)%\) done")


def _render_bar(done, total, width=36):
    if total <= 0:
        return "[" + ("-" * width) + "]"
    frac = max(0.0, min(1.0, float(done) / float(total)))
    fill = int(round(frac * width))
    return "[" + ("#" * fill) + ("-" * (width - fill)) + "]"


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option("-j", "--jobs", type=int, default=2, show_default=True)
@click.option(
    "--rerun-incomplete/--no-rerun-incomplete", default=True, show_default=True
)
@click.option(
    "--snakefile",
    type=click.Path(path_type=str),
    default="Snakefile",
    show_default=True,
)
@click.option(
    "--configfile",
    type=click.Path(path_type=str),
    default=None,
    help="Explicit Snakemake config file path.",
)
def main(jobs, rerun_incomplete, snakefile, configfile):
    cmd = ["snakemake", "-j", str(jobs), "--snakefile", snakefile]
    if configfile:
        cmd.extend(["--configfile", configfile])
    if rerun_incomplete:
        cmd.append("--rerun-incomplete")

    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
    )

    last_printed = None
    try:
        assert proc.stdout is not None
        for line in proc.stdout:
            line = line.rstrip("\n")
            m = PROGRESS_RE.search(line)
            if m:
                done = int(m.group(1))
                total = int(m.group(2))
                pct = float(m.group(3))
                state = f"{_render_bar(done, total)} {done}/{total} ({pct:.1f}%)"
                if state != last_printed:
                    print(state, flush=True)
                    last_printed = state
                continue
            print(line, flush=True)
    except KeyboardInterrupt:
        proc.terminate()
        raise

    rc = proc.wait()
    if rc != 0:
        sys.exit(rc)


if __name__ == "__main__":
    main()
