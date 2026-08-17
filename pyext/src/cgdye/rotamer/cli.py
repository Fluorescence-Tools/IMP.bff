"""Click CLI for rotamer-based cgdye tools."""

from __future__ import annotations

import click
from IMP.bff.cgdye.rotamer.fret import RotamerFRET
from IMP.bff.fret.forster import forster_radius_from_spectra


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
