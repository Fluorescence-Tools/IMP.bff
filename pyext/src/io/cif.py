"""Read/write compact FF topology mmCIF tables.

The C++ FF writer and rotamer library IO are in :file:`include/IMP/bff/CifIO.h`,
wrapped through :file:`pyext/IMP_bff.cif.i`. The FF reader is in
:file:`include/IMP/bff/ForceFieldCIF.h`. The template CIF reader/writer stays
as %pythoncode in the .i file (it returns dicts and uses ihm.format, which is
Python-only). This module re-exports the Python surface so
``from IMP.bff.io.cif import ...`` keeps working.
"""

from IMP.bff import (
    RotamerLibraryData,
    _base_path,
    _parse_mol2_atom_site_rows,
    _parse_pdb_atom_site_rows,
    _parse_struct_atom_site_rows,
    _write_cif_safe,
    as_forcefield_system,
    compress_int_ranges,
    expand_site_range,
    forcefield_system_from_dict,
    normalize_weights,
    read_component_template_cif,
    read_dye_forcefield_cif,
    read_dye_template_cif,
    read_rotamer_library,
    region_features,
    split_site_id,
    write_component_template_cif,
    write_dye_forcefield_cif,
    write_dye_template_cif,
    write_rotamer_library,
)

__all__ = [
    "RotamerLibraryData",
    "as_forcefield_system",
    "compress_int_ranges",
    "expand_site_range",
    "forcefield_system_from_dict",
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
