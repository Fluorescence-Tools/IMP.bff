/*
 * A coordinate file as flat columns, filled by IMP's readers.
 *
 * `StructureTable` belongs to the connection layer -- it exists only where
 * IMP is linked -- but nothing in its declaration names an IMP type, which
 * is what lets it be wrapped without IMP's own SWIG interfaces. A module
 * that wraps `IMP::atom::Hierarchy` has to `%import` IMP's kernel interface,
 * and the extension then calls `PyImport_ImportModule("_IMP_kernel")` at
 * load; behind columns, none of that follows (PRD-139).
 *
 * The string columns stay plain `std::vector<std::string>` -- SWIG hands
 * those out as sequences, and numpy builds its `<U5` column from one in a
 * single call. Only the numeric columns get out-views, because those are the
 * ones a caller wants as an array without a copy per atom.
 */

%feature("kwargs") IMP::bff::read_structure_table;

IMP_SWIG_VALUE(IMP::bff, StructureTable, StructureTables);

/* Hidden from SWIG's member wrapping and republished below as shaped arrays:
   a flat `xyz` and an `(n, 3)` `xyz` cannot both be called `xyz`. */
%ignore IMP::bff::StructureTable::xyz;
%ignore IMP::bff::StructureTable::radius;
%ignore IMP::bff::StructureTable::mass;
%ignore IMP::bff::StructureTable::bfactor;
%ignore IMP::bff::StructureTable::atom_id;
%ignore IMP::bff::StructureTable::res_id;

%include "IMP/bff/StructureTable.h"

%attribute_np2(IMP::bff::StructureTable, std::vector<double>, xyz,
               get_xyz, 3);
%attribute_np(IMP::bff::StructureTable, std::vector<double>, radius,
              get_radius);
%attribute_np(IMP::bff::StructureTable, std::vector<double>, mass,
              get_mass);
%attribute_np(IMP::bff::StructureTable, std::vector<double>, bfactor,
              get_bfactor);
%attribute_np(IMP::bff::StructureTable, std::vector<int>, atom_id,
              get_atom_id);
%attribute_np(IMP::bff::StructureTable, std::vector<int>, res_id,
              get_res_id);
%attribute(IMP::bff::StructureTable, int, n_atoms, get_n_atoms);
