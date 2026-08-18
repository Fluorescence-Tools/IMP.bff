"""Read the bundled dye library from mmCIF.

**CIF is the data format for this package** (PRD-113). The library was CSV --
one table of extinction coefficients and quantum yields plus one file of
excitation/emission curves per dye -- which carried no category names, no units
and no way to say where a number came from.

Two categories, both bff-native:

``_bff_dye``
    ``chromophore_name``, ``probe_type``, ``chromophore_number``,
    ``extinction_coefficient``, ``quantum_yield``
``_bff_dye_spectrum``
    ``chromophore_name``, ``wavelength``, ``excitation``, ``emission``

They are bff-native because **no dictionary in the stack defines an item for a
quantum yield, an extinction coefficient or a spectrum** -- checked across all
ten ``.dic`` files in ``../mmfdb/src/mmfdb/data`` (mmCIF std, PDBx v50 and
v5_next, DDL, MA, IHM, IHM-FLR, mmfdb FLR and workflow extensions); the only
matches are ``_em_detector.detective_quantum_efficiency`` and the NMR spectral
categories. They should be proposed for ``mmfdb_flr_ext.dic``.

Identifiers follow flrCIF: ``chromophore_name`` is
``_flr_probe_list.chromophore_name``.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Optional

import numpy as np
import ihm.format

from .species import Dye, Spectrum

__all__ = ["read_dye_library", "DYE_LIBRARY_CIF"]

#: Name of the bundled library inside the rotamer-library data directory.
DYE_LIBRARY_CIF = "dye_library.cif"


class _BaseHandler:
    """``ihm.format.CifReader`` reads the keywords from ``__call__``'s signature.

    So the parameter names below *are* the CIF item names -- a ``**kwargs``
    handler is handed nothing at all, silently.
    """

    not_in_file = object()
    omitted = object()
    unknown = object()

    def __init__(self):
        self.rows = []

    def end_save_frame(self):
        pass

    def _clean(self, **kwargs):
        self.rows.append({
            k: (None if v in (self.not_in_file, self.omitted, self.unknown) else v)
            for k, v in kwargs.items()})


class _DyeHandler(_BaseHandler):
    def __call__(self, chromophore_name, probe_type, chromophore_number,
                 extinction_coefficient, quantum_yield):
        self._clean(chromophore_name=chromophore_name, probe_type=probe_type,
                    chromophore_number=chromophore_number,
                    extinction_coefficient=extinction_coefficient,
                    quantum_yield=quantum_yield)


class _SpectrumHandler(_BaseHandler):
    def __call__(self, chromophore_name, wavelength, excitation, emission):
        self._clean(chromophore_name=chromophore_name, wavelength=wavelength,
                    excitation=excitation, emission=emission)


def _float(value) -> Optional[float]:
    if value is None or value in ("", ".", "?"):
        return None
    return float(value)


#: Parsed libraries, keyed by resolved path and modification time.
_LIBRARY_CACHE: Dict[tuple, Dict[str, "Dye"]] = {}


def read_dye_library(path: str | Path | None = None) -> Dict[str, Dye]:
    """Every dye in the bundled library, keyed by chromophore name.

    :param path: the CIF file; defaults to the bundled one.
    :returns: ``{chromophore_name: Dye}``. A dye with a table entry but no
        curves is returned without a spectrum -- it is still a usable species,
        it just cannot derive R0.

    **Cached on (path, mtime, size).** The bundled library is a shipped
    read-only file, and a spectrum lookup does not change it -- but the parse is
    not cheap: profiling one FRETpredict comparison found this called 21 times
    for the same file, 550 000 calls to the row cleaner and 1.65 million to the
    float converter, for 2.2 s of a 13 s run.

    The key includes mtime and size so editing the file during a session is
    picked up; only an edit that changes neither would be missed, which is not a
    thing that happens to a data file. Callers that mutate the returned ``Dye``
    objects would now be mutating the cached ones -- nothing does, and nothing
    should: a dye is a species, not a scratch pad.
    """
    if path is None:
        import IMP.bff
        path = Path(IMP.bff.get_data_path("rotamer_library")) / "R0" / DYE_LIBRARY_CIF
    resolved = Path(path)
    key = None
    if resolved.exists():
        stat = resolved.stat()
        key = (str(resolved.resolve()), stat.st_mtime_ns, stat.st_size)
        hit = _LIBRARY_CACHE.get(key)
        if hit is not None:
            # A fresh mapping each time, so a caller that adds or drops an entry
            # does not edit the library for everyone. The `Dye` values are
            # shared, and that is safe because `Dye` is frozen -- reaching round
            # that with `object.__setattr__` corrupts the species for every
            # later reader, which is exactly what a test did on the day this
            # cache landed.
            return dict(hit)
    path = resolved

    dyes = _DyeHandler()
    spectra = _SpectrumHandler()
    with path.open() as handle:
        reader = ihm.format.CifReader(handle, {"_bff_dye": dyes,
                                               "_bff_dye_spectrum": spectra})
        reader.read_file()

    curves: Dict[str, list] = {}
    for row in spectra.rows:
        name = row.get("chromophore_name")
        if name is None:
            continue
        curves.setdefault(name, []).append(
            (_float(row.get("wavelength")),
             _float(row.get("excitation")),
             _float(row.get("emission"))))

    out: Dict[str, Dye] = {}
    for row in dyes.rows:
        name = row.get("chromophore_name")
        if name is None:
            continue
        points = sorted(curves.get(name, []), key=lambda r: r[0])
        spectrum = None
        if points:
            spectrum = Spectrum(
                wavelength=np.array([p[0] for p in points], dtype=float),
                excitation=np.array([p[1] for p in points], dtype=float),
                emission=np.array([p[2] for p in points], dtype=float),
            )
        out[name] = Dye(
            name=name,
            spectrum=spectrum,
            extinction_coefficient=_float(row.get("extinction_coefficient")),
            quantum_yield=_float(row.get("quantum_yield")),
        )
    if key is not None:
        _LIBRARY_CACHE[key] = out
        return dict(out)
    return out
