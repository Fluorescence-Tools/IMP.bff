"""Read/write compact FF topology mmCIF tables.

The C++ FF writer, template CIF reader/writer, rotamer library IO and the
string utilities are in :file:`include/IMP/bff/CifIO.h`, wrapped through
:file:`pyext/IMP_bff.cif.i`. The FF reader is in
:file:`include/IMP/bff/ForceFieldCIF.h`; the dict-literal convenience bridge is
`forcefield_system_from_json` in C++ (nlohmann).

Port status (PRD-117): everything is C++ except `as_forcefield_system` and
`forcefield_system_from_dict`, which are duck typing on a Python object (does
the caller have a dict or a typed system?) and the one `json.dumps` beside it.
"""

import json

from IMP.bff import (
    RotamerLibraryData,
    ComponentTemplate,
    compress_int_ranges,
    expand_site_range,
    forcefield_system_from_json,
    normalize_weights,
    read_component_template_cif,
    read_forcefield_cif,
    read_dye_template_cif,
    read_rotamer_library,
    region_features,
    split_site_id,
    write_component_template_cif,
    write_dye_forcefield_cif as _write_system,
    write_dye_template_cif,
    write_rotamer_library,
)


def forcefield_system_from_dict(d):
    """A parsed or built dictionary as a typed DyeForceFieldSystem."""
    return forcefield_system_from_json(json.dumps(d))


def as_forcefield_system(system):
    """A system, whichever way it was given."""
    if isinstance(system, dict):
        return forcefield_system_from_dict(system)
    return system


def read_dye_forcefield_cif(path):
    """Read a force-field system from mmCIF."""
    return read_forcefield_cif(str(path))


def write_dye_forcefield_cif(path, system):
    """Write a force-field system to mmCIF. Accepts DyeForceFieldSystem or dict."""
    _write_system(str(path), as_forcefield_system(system))


__all__ = [
    "ComponentTemplate",
    "RotamerLibraryData",
    "as_forcefield_system",
    "compress_int_ranges",
    "expand_site_range",
    "forcefield_system_from_dict",
    "forcefield_system_from_json",
    "normalize_weights",
    "read_component_template_cif",
    "read_dye_forcefield_cif",
    "read_dye_template_cif",
    "read_rotamer_library",
    "region_features",
    "split_site_id",
    "write_component_template_cif",
    "write_dye_forcefield_cif",
    "write_dye_template_cif",
    "write_rotamer_library",
]
