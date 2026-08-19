#!/usr/bin/env python
"""AV ↔ rotamer-ensemble comparison on the bundled hGBP1 and T4L systems (PRD-108).

Writes the authoritative table ``okf/validation/av_vs_rotamer.md`` and the
pin file ``test/references/cgdye_av_vs_rotamer_pins.json`` when run from the
repository (``--okf-dir``/``--pins`` override the paths); otherwise prints
the tables.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import click

from IMP.bff.representation.compare import (
    compare_pairs, compare_positions, markdown_table, summary_numbers)

DONOR_LIB = "AlexaFluor 488 C1R"      # + " cutoff<N>"
ACCEPTOR_LIB = "AlexaFluor 594 C1R"
R0_A488_A594 = 52.0   # Å, the value the T4L fps.json carries


def _libs(cutoff):
    return f"{DONOR_LIB} cutoff{cutoff}", f"{ACCEPTOR_LIB} cutoff{cutoff}"


def hgbp1_case(cutoff=30):
    """hGBP1 (1DG3, chain A): the fps.json positions of examples/structure/GBP/hGBP1.fps.json on chain A."""
    import IMP.bff
    from IMP.bff.tools.paths import get_structure_dir
    pdb = str(get_structure_dir("1DG3.pdb"))
    fps = json.load(open(IMP.bff.get_example_path("structure/GBP/hGBP1.fps.json")))
    # residue 254 is not resolved in 1DG3; the other chain-A sites are
    positions = {n: p for n, p in fps["Positions"].items()
                 if p.get("chain_identifier") == "A" and n.endswith("F") and p["residue_seq_number"] != 254}
    d_lib, a_lib = _libs(cutoff)
    libraries = {n: (d_lib if p["residue_seq_number"] == 481 else a_lib) for n, p in positions.items()}
    pairs = [("A481F", n) for n in sorted(positions) if n != "A481F"]
    return f"hGBP1 1DG3 chain A — donor Alexa488 C1R at 481, acceptor Alexa594 C1R elsewhere (cutoff{cutoff} libraries)", pdb, positions, libraries, pairs, {}


def t4l_case(cutoff=30):
    """T4 lysozyme (3GUN) with the shipped fret.fps.json: D positions Alexa488 C1R, A positions Alexa594 C1R."""
    import IMP.bff
    pdb = IMP.bff.get_example_path("structure/T4L/3GUN.pdb")
    fps = json.load(open(IMP.bff.get_example_path("structure/T4L/fret.fps.json")))
    positions = dict(fps["Positions"])
    d_lib, a_lib = _libs(cutoff)
    libraries = {n: (d_lib if n.endswith("D") else a_lib) for n in positions}
    experimental = {}
    pairs = []
    for name, d in fps["Distances"].items():
        key = (d["position1_name"], d["position2_name"])
        pairs.append(key)
        experimental[key] = {"distance": d["distance"], "error_neg": d["error_neg"], "error_pos": d["error_pos"],
                             "distance_type": d.get("distance_type", "RDAMean"), "name": name}
    return f"T4L 3GUN — the 99 fps.json distances (D = Alexa488 C1R, A = Alexa594 C1R, cutoff{cutoff} libraries)", pdb, positions, libraries, pairs, experimental


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option("--system", type=click.Choice(["hgbp1", "t4l", "both"]), default="both", show_default=True)
@click.option("--okf-dir", default=None, help="Write av_vs_rotamer.md here (default: repo okf/validation if found).")
@click.option("--pins", default=None, help="Write the pin JSON here (default: repo test/references).")
@click.option("--n-samples", default=50000, show_default=True, type=int, help="AV distance MC samples.")
@click.option("--temperature", default=298.15, show_default=True, type=float)
@click.option("--cutoffs", default="30,10", show_default=True, help="Rotamer-library cutoffs to compare (comma separated).")
def main(system, okf_dir, pins, n_samples, temperature, cutoffs):
    """Compare AV clouds and rotamer ensembles position by position and pair by pair."""
    repo = None
    for parent in Path(__file__).resolve().parents:
        if (parent / "okf" / "prds").exists():
            repo = parent
            break
    cases = {"hgbp1": hgbp1_case, "t4l": t4l_case}
    names = ["hgbp1", "t4l"] if system == "both" else [system]
    md = ["# AV ↔ rotamer-ensemble cross-validation (PRD-108)", "",
          "Recorded by `python -m IMP.bff.representation.compare_cli`. AVs from the fps.json positions "
          "(default FPS strip; the authored T4L `strip_mask` is outside the current dialect), rotamer ensembles "
          f"screened at T = {temperature} K without electrostatics; R0 = {R0_A488_A594} Å (κ² = 2/3). "
          "This table is authoritative; the test pins the numbers and asserts loose sanity bounds only.", ""]
    numbers = {"_note": "AV vs rotamer-ensemble numbers recorded by compare_av_rotamer (PRD-108 stage 2); drift pins, not physics gates.",
               "settings": {"n_samples": n_samples, "temperature": temperature, "forster_radius": R0_A488_A594,
                            "donor_library": DONOR_LIB, "acceptor_library": ACCEPTOR_LIB, "cutoffs": cutoffs}}
    for cutoff in [int(c) for c in cutoffs.split(",") if c.strip()]:
        for key in names:
            title, pdb, positions, libraries, pairs, experimental = cases[key](cutoff)
            per_pos = compare_positions(pdb, positions, libraries, n_samples=n_samples, temperature=temperature)
            rows = compare_pairs(per_pos, pairs, R0_A488_A594, experimental=experimental, n_samples=n_samples)
            md.append(markdown_table(per_pos, rows, title))
            numbers[f"{key}_cutoff{cutoff}"] = summary_numbers(per_pos, rows)
    text = "\n".join(md)
    click.echo(text)
    okf_path = Path(okf_dir) if okf_dir else (repo / "okf" / "validation" if repo else None)
    if okf_path is not None:
        okf_path.mkdir(parents=True, exist_ok=True)
        (okf_path / "av_vs_rotamer.md").write_text(text)
        click.echo(f"wrote {okf_path / 'av_vs_rotamer.md'}")
    pin_path = Path(pins) if pins else (repo / "test" / "references" / "cgdye_av_vs_rotamer_pins.json" if repo else None)
    if pin_path is not None:
        pin_path.write_text(json.dumps(numbers, indent=1))
        click.echo(f"wrote {pin_path}")


if __name__ == "__main__":
    main()
