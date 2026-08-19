"""Support for the ``imp_bff`` console script.

Not the whole command line. This package holds what ``bin/imp_bff`` imports --
the optimizer states that write frames during a flexible fit, and the angle-file
reader it starts from. The commands themselves (``flexfit``, ``rmsd``,
``auto_model``) live in ``bin/imp_bff``, because that is the file the console
script points at.

There is a second, unconnected command tree: ``python -m IMP.bff.cgdye.cli``
gives a ``dye`` group with its own subcommands, including ``rotamer``. It is not
reachable through ``imp_bff``. That split is a leftover of the PRD-107
migration, not a design; joining them is worth doing and has not been done.

Every command-line entry point that is *not* cgdye's now lives here, including
two that used to sit beside the code they drive. That is the rule, not tidiness:
a click command is a decorated function, so ``import click`` runs at module
scope, and a library module that carries one cannot be imported without click
installed.
"""

import click
import RMF
import os
import pathlib
import typing

import IMP.atom

# What the moved command sections call. They used to sit next to these names.
from IMP.bff.dye import forster_radius_from_spectra

__all__ = [
    'WritePDBFrame',
    'WriteRMFFrame',
    'read_angle_file',
]

# --------------------------------------------------------------------------
# flexfit
# --------------------------------------------------------------------------
def read_angle_file(
        hier: IMP.atom.Hierarchy,
        flex_dict: dict
) -> typing.Tuple[
    typing.List[IMP.atom.Residue],
    typing.List[IMP.atom.Bond]
]:
    flexible_residues = list()
    bonds = list()
    for fr in flex_dict["Flexible residues"]:
        chain_id = fr['chain_identifier']
        residue_id = fr['residue_seq_number']
        sel = IMP.atom.Selection(
            hierarchy=hier,
            chain_ids=[chain_id],
            residue_indexes=[residue_id]
        )
        flexible_residues.append(
            sel.get_selected_particles(False)[0]
        )
    for bnd in flex_dict["Bonds"]:
        a1 = bnd[0]
        a2 = bnd[1]
        sel1 = IMP.atom.Selection(
            hierarchy=hier,
            chain_ids=[a1['chain_identifier']],
            residue_indexes=[a1['residue_seq_number']],
            atom_type=IMP.atom.AtomType(a1['atom_name'])
        )
        sel2 = IMP.atom.Selection(
            hierarchy=hier,
            chain_ids=[a2['chain_identifier']],
            residue_indexes=[a2['residue_seq_number']],
            atom_type=IMP.atom.AtomType(a2['atom_name'])
        )
        p1 = sel1.get_selected_particles()[0]
        p2 = sel2.get_selected_particles()[0]
        b1 = IMP.atom.Bonded(p1)
        b2 = IMP.atom.Bonded(p2)
        bonds.append(
            IMP.atom.create_bond(b1, b2, IMP.atom.Bond.SINGLE)
        )
    return flexible_residues, bonds


class WriteRMFFrame(IMP.OptimizerState):

    def __init__(
            self,
            filename,
            root_hier: IMP.atom.Hierarchy,
            restraints: typing.List[IMP.Restraint],
            name: str = "WriteRMFFrame"
    ):
        model = root_hier.get_model()
        super().__init__(model, name)
        self.restraints = restraints
        fileName, fileExtension = os.path.splitext(filename)
        fn = pathlib.Path(fileName + ".0.rmf3")
        if pathlib.Path(fn).exists():
            for i in range(100):
                fn = pathlib.Path(fileName + ".{:d}.rmf3".format(i))
                if not fn.exists():
                    break
        self._rmf_filename = str(fn)
        self._rh = RMF.create_rmf_file(str(fn))
        IMP.rmf.add_hierarchies(self._rh, [root_hier])
        IMP.rmf.add_restraints(self._rh, restraints)
        IMP.rmf.save_frame(self._rh)

    def do_update(self, arg0):
        #print(*[r.evaluate(False) for r in self.restraints], sep="\t")
        IMP.rmf.save_frame(self._rh)


class WritePDBFrame(IMP.OptimizerState):

    def __init__(
            self,
            filename,
            root_hier: IMP.atom.Hierarchy,
            restraints: typing.List[IMP.Restraint],
            restraint_filename: str = None,
            name: str = "WriteDCDFrame",
            output_objects: typing.List = None,
            multi_state: bool = False
    ):
        model = root_hier.get_model()
        super().__init__(model, name)

        import os
        self._pdb_basename = filename
        if restraint_filename is None:
            restraint_filename = os.path.splitext(filename)[0] + ".rst.txt"
        self._restraint_filename = restraint_filename
        self.restraints = restraints
        self._hier = root_hier
        self.frame = 0
        self.output_objects = output_objects
        self.multi_state = multi_state

    def do_update(self, arg0):
        with open(self._restraint_filename, "a+") as fp:
            fp.write("%s\t" % self.frame)
            fp.write("\t".join(["{:.3f}".format(r.evaluate(False)) for r in self.restraints]))
            fp.write("\t")
            if isinstance(self.output_objects, list):
                for obj in self.output_objects:
                    fp.write(str(obj) + "\t")
            fp.write("\n")
        lead = os.path.splitext(self._pdb_basename)[0]
        if not self.multi_state:
            hiers = [self._hier]
        else:
            hiers = self._hier.get_children()
        for i, hier in enumerate(hiers):
            out_fn = lead + "_state_" + str(i) + "_" + "{:04d}".format(self.frame) + ".pdb"
            IMP.atom.write_pdb(hier, out=out_fn)
        self.frame += 1


# --------------------------------------------------------------------------
# where the command-line entry points live
# --------------------------------------------------------------------------
# They are gathered here rather than beside the code they drive, and that is a
# rule rather than a tidiness preference: a click command is a *decorated*
# function, so `import click` has to run at module scope, and merging
# `representation/rotamer/cli.py` into `representation/rotamer.py` made
# `import IMP.bff.representation.rotamer` require click. A library module must
# import cleanly without it -- `IMP.bff` is imported by code that will never
# see a terminal -- and `test_import_is_lazy_and_click_free` is what says so.
#
# A CLI is not one of the four stages. It is a way in.
# --------------------------------------------------------------------------
# rotamer commands (was representation/rotamer/cli.py)
# --------------------------------------------------------------------------
"""Click CLI for rotamer-based cgdye tools."""

@click.group()
def rotamer() -> None:
    """Rotamer-based cgdye tools."""


@rotamer.command()
@click.option("--protein", "protein_path", required=True, type=click.Path(exists=True), help="Protein PDB or RMF path.")
@click.option("--residue", "residues", multiple=True, required=True, type=int, help="Placement residue number. Repeat twice.")
@click.option("--chain", "chains", multiple=True, default=None, help="Placement chain ID. Repeat twice.")
@click.option("--donor", default="AlexaFluor 488", show_default=True, help="Donor dye name.")
@click.option("--acceptor", default="AlexaFluor 594", show_default=True, help="Acceptor dye name.")
@click.option("--libname-1", required=True, help="Donor rotamer library name or RMF path.")
@click.option("--libname-2", required=True, help="Acceptor rotamer library name or RMF path.")
@click.option("--temperature", default=300.0, show_default=True, type=float, help="Temperature in K.")
@click.option("--electrostatic", is_flag=True, help="Include Debye-Huckel electrostatics.")
@click.option("--fixed-r0", is_flag=True, help="Use a fixed Förster radius.")
@click.option("--r0", default=5.4, show_default=True, type=float, help="Fixed R0 in nm.")
@click.option("--output-prefix", default="rotamer_fret", show_default=True, help="Output prefix.")
@click.option("--z-cutoff", default=0.05, show_default=True, type=float, help="Partition-function cutoff.")
@click.option("--max-frames", default=None, type=int, help="Maximum protein frames to process.")
def predict(
    protein_path: str,
    residues: tuple[int, ...],
    chains: tuple[str, ...] | None,
    donor: str,
    acceptor: str,
    libname_1: str,
    libname_2: str,
    temperature: float,
    electrostatic: bool,
    fixed_r0: bool,
    r0: float,
    output_prefix: str,
    z_cutoff: float,
    max_frames: int | None,
) -> None:
    """Run rotamer-based FRET prediction."""
    if len(residues) != 2:
        raise click.UsageError("--residue must be provided exactly twice")
    if chains is None or len(chains) == 0:
        chain_values = [None, None]
    elif len(chains) == 2:
        chain_values = list(chains)
    else:
        raise click.UsageError("--chain must be omitted or provided exactly twice")

    fret = RotamerFRET(
        protein_path,
        list(residues),
        chains=chain_values,
        donor=donor,
        acceptor=acceptor,
        libname_1=libname_1,
        libname_2=libname_2,
        temperature=temperature,
        electrostatic=electrostatic,
        fixed_R0=fixed_r0,
        r0=r0,
        output_prefix=output_prefix,
        z_cutoff=z_cutoff,
        max_frames=max_frames,
    )
    fret.run()
    click.echo(f"Wrote FRET prediction files with prefix {output_prefix!r}")


@rotamer.command()
@click.option("--donor", required=True, help="Donor dye name.")
@click.option("--acceptor", required=True, help="Acceptor dye name.")
@click.option("--k2", required=True, type=float, help="Orientation factor.")
def r0(donor: str, acceptor: str, k2: float) -> None:
    """Calculate a Förster radius."""
    value = forster_radius_from_spectra(donor, acceptor, k2)
    click.echo(f"{value:.6f} nm")


if __name__ == "__main__":
    rotamer()


@rotamer.command("compare-av")
@click.option("--system", type=click.Choice(["hgbp1", "t4l", "both"]), default="both", show_default=True)
@click.option("--okf-dir", default=None, help="Write av_vs_rotamer.md here (default: the repo's okf/validation when run from a checkout).")
@click.option("--pins", default=None, help="Write the pin JSON here (default: the repo's test/references).")
@click.option("--n-samples", default=50000, show_default=True, type=int)
@click.option("--temperature", default=298.15, show_default=True, type=float)
@click.option("--cutoffs", default="30,10", show_default=True)
@click.pass_context
def compare_av(ctx, **kwargs) -> None:
    """AV ↔ rotamer-ensemble comparison on the bundled hGBP1 and T4L systems (PRD-108)."""
    from IMP.bff.representation.compare import main as _main
    ctx.invoke(_main, **kwargs)

# --------------------------------------------------------------------------
# av-vs-rotamer comparison (was representation/compare_cli.py)
# --------------------------------------------------------------------------
"""AV ↔ rotamer-ensemble comparison on the bundled hGBP1 and T4L systems (PRD-108).

Writes the authoritative table ``okf/validation/av_vs_rotamer.md`` and the
pin file ``test/references/cgdye_av_vs_rotamer_pins.json`` when run from the
repository (``--okf-dir``/``--pins`` override the paths); otherwise prints
the tables.
"""

#!/usr/bin/env python





DONOR_LIB = "AlexaFluor 488 C1R"      # + " cutoff<N>"
ACCEPTOR_LIB = "AlexaFluor 594 C1R"
R0_A488_A594 = 52.0   # Å, the value the T4L fps.json carries


def _libs(cutoff):
    return f"{DONOR_LIB} cutoff{cutoff}", f"{ACCEPTOR_LIB} cutoff{cutoff}"


def hgbp1_case(cutoff=30):
    """hGBP1 (1DG3, chain A): the fps.json positions of examples/structure/GBP/hGBP1.fps.json on chain A."""
    import IMP.bff
    from IMP.bff.tools import get_structure_dir
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


