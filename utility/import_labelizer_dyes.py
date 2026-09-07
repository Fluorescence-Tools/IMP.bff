#!/usr/bin/env python3
"""Add the dyes Labelizer ships that this package's library does not have.

`data/rotamer_library/R0/probe_library.cif` came from the FRETpredict tables and
holds 38 dyes. The Labelizer package ships 15, and thirteen of them are already
here under the vendor's spelling with **identical** quantum yields and
extinction coefficients -- both sets descend from the same upstream data, and
`resolve_probe_name` bridges `Alexa488` to `AlexaFluor488`. Two are genuinely
absent:

    Atto532   QY 0.90   eps 115000
    Atto643   QY 0.62   eps 150000

This imports those two, and nothing else. Run it once; it is idempotent and
refuses to duplicate a dye that is already there.

Why the resampling matters
--------------------------
`spectral_overlap` requires that **every dye share one wavelength grid** and
throws otherwise, so the imported spectra are put on this library's grid --
300 to 970 nm at 1 nm. Labelizer's tables stop at 821 nm (Atto532) and 900 nm
(Atto643), and the gap is filled with zeros, which is what the reference does
too (`fluorophore.py:188` zero-pads outside the measured range).

Normalisation is not cosmetic here. The donor's emission is divided by its own
integral, so its scale cancels -- but the acceptor's excitation is multiplied
by the extinction coefficient to give a molar absorptivity, so it **must** be
peak-normalised to 1 or R0 is wrong by whatever factor it is out. Both curves
are therefore scaled to a peak of 1.

Usage
-----
::

    python utility/import_labelizer_dyes.py
    python utility/import_labelizer_dyes.py --dry-run
    python utility/import_labelizer_dyes.py --labelizer /path/to/resources
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

#: The library this package ships, relative to the repository root.
LIBRARY = Path("data/rotamer_library/R0/probe_library.cif")

#: Where Labelizer keeps its fluorophores, relative to the sibling checkout.
DEFAULT_RESOURCES = Path(
    "../labelizer-backend/src_labelizer/labelizer/resources/fluorophores"
)

#: The grid every dye in the library is on. `spectral_overlap` throws if two
#: dyes disagree, so this is not negotiable.
GRID_START, GRID_END, GRID_STEP = 300, 970, 1

#: What to import, and what to call it here. The library keys ATTO dyes in
#: caps with no space (`ATTO647N`), and carries the vendor and the number in
#: their own columns.
WANTED = {
    "Atto532": ("ATTO532", "ATTO", "532"),
    "Atto643": ("ATTO643", "ATTO", "643"),
}


def read_scalars(resources: Path, name: str) -> tuple[float, float]:
    """Quantum yield and extinction coefficient from a `.fluo` file."""
    values = {}
    for line in (resources / f"{name}.fluo").read_text().splitlines():
        if "," in line:
            key, value = line.split(",", 1)
            values[key.strip()] = float(value)
    return values["QY"], values["EC"]


def read_spectra(resources: Path, name: str) -> tuple[list[float], list[float]]:
    """Excitation and emission on the library grid, peak-normalised.

    The CSVs are ``Wavelength, <dye> EM, <dye> AB`` -- absorption is what the
    library calls excitation. A blank cell means the curve was not measured
    there and becomes zero, as does anything outside the file's range.
    """
    path = resources / f"{name}.csv"
    # utf-8-sig: the files carry a BOM, which otherwise ends up in the first
    # column name and makes the header unmatchable.
    with open(path, encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.reader(handle))

    header = [h.strip().strip('"') for h in rows[0]]
    try:
        i_em = next(i for i, h in enumerate(header) if h.endswith("EM"))
        i_ab = next(i for i, h in enumerate(header) if h.endswith(("AB", "EX")))
    except StopIteration:
        raise SystemExit(f"{path}: expected an EM and an AB/EX column, got {header}")

    by_wavelength: dict[int, tuple[float, float]] = {}
    for row in rows[1:]:
        if not row or not row[0].strip():
            continue
        wavelength = int(round(float(row[0])))

        def cell(index: int) -> float:
            if index >= len(row) or not row[index].strip():
                return 0.0
            return float(row[index])

        by_wavelength[wavelength] = (cell(i_ab), cell(i_em))

    excitation, emission = [], []
    for wavelength in range(GRID_START, GRID_END + 1, GRID_STEP):
        ex, em = by_wavelength.get(wavelength, (0.0, 0.0))
        excitation.append(ex)
        emission.append(em)

    for curve, what in ((excitation, "excitation"), (emission, "emission")):
        peak = max(curve)
        if peak <= 0.0:
            raise SystemExit(f"{path}: the {what} curve is empty")
        for i, value in enumerate(curve):
            curve[i] = value / peak
    return excitation, emission


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--labelizer", type=Path, default=DEFAULT_RESOURCES,
                    help="Labelizer's fluorophore resource directory")
    ap.add_argument("--library", type=Path, default=LIBRARY,
                    help="the dye library CIF to extend")
    ap.add_argument("--dry-run", action="store_true",
                    help="report what would be added and change nothing")
    args = ap.parse_args(argv)

    if not args.labelizer.is_dir():
        print(f"{args.labelizer} is not a directory; pass --labelizer",
              file=sys.stderr)
        return 1

    text = args.library.read_text()
    dye_rows, spectrum_rows = [], []

    for source, (key, vendor, number) in sorted(WANTED.items()):
        if f"\n{key} " in text:
            print(f"{key} is already in the library; skipping")
            continue
        quantum_yield, extinction = read_scalars(args.labelizer, source)
        excitation, emission = read_spectra(args.labelizer, source)

        # `dye` is the _bff_probe.probe_type column: the library is
        # probe-generic now, and everything Labelizer ships is an organic dye.
        dye_rows.append(
            f"{key} dye {vendor} {number} {extinction} {quantum_yield}")
        for i, wavelength in enumerate(
                range(GRID_START, GRID_END + 1, GRID_STEP)):
            spectrum_rows.append(
                f"{key} {float(wavelength)} {excitation[i]} {emission[i]}")
        print(f"{key}: QY {quantum_yield}, eps {extinction}, "
              f"{len(excitation)} points from {source}")

    if not dye_rows:
        print("nothing to add")
        return 0
    if args.dry_run:
        print(f"--dry-run: would add {len(dye_rows)} dyes, "
              f"{len(spectrum_rows)} spectrum rows")
        return 0

    # The scalar loop ends at the first `#` after it; the spectrum loop runs to
    # the end of the file. Insert at each, rather than rewriting the file, so
    # the 38 dyes already there are untouched byte for byte.
    marker = "\n#\n#\nloop_\n_bff_probe_spectrum.chromophore_name"
    if marker not in text:
        print("could not find the boundary between the two loops",
              file=sys.stderr)
        return 1
    text = text.replace(marker, "\n" + "\n".join(dye_rows) + marker, 1)

    tail = text.rstrip()
    if not tail.endswith("#"):
        print("the library does not end with the expected `#`", file=sys.stderr)
        return 1
    text = tail[:-1].rstrip() + "\n" + "\n".join(spectrum_rows) + "\n#\n"

    args.library.write_text(text)
    print(f"wrote {args.library}: +{len(dye_rows)} dyes, "
          f"+{len(spectrum_rows)} spectrum rows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
