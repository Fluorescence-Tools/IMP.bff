/*
 * A structure as flat columns -- the record a reader fills and a writer
 * takes. Core: it names neither an IMP nor an RMF type.
 *
 * The string columns stay plain `std::vector<std::string>` -- SWIG hands
 * those out as sequences, and numpy builds its `<U5` column from one in a
 * single call. Only the numeric columns get out-views, because those are the
 * ones a caller wants as an array without a copy per atom.
 */

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
               get_xyz, 3, set_xyz);
%attribute_np(IMP::bff::StructureTable, std::vector<double>, radius,
              get_radius, set_radius);
%attribute_np(IMP::bff::StructureTable, std::vector<double>, mass,
              get_mass, set_mass);
%attribute_np(IMP::bff::StructureTable, std::vector<double>, bfactor,
              get_bfactor, set_bfactor);
%attribute_np(IMP::bff::StructureTable, std::vector<int>, atom_id,
              get_atom_id, set_atom_id);
%attribute_np(IMP::bff::StructureTable, std::vector<int>, res_id,
              get_res_id, set_res_id);
%attribute(IMP::bff::StructureTable, int, n_atoms, get_n_atoms);
