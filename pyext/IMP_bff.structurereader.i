/*
 * IMP's PDB and mmCIF readers behind the flat table.
 *
 * `read_structure_table` belongs to the connection layer -- it exists only
 * where IMP is linked -- but nothing in its declaration names an IMP type,
 * which is what lets it be wrapped without IMP's own SWIG interfaces. A
 * module that wraps `IMP::atom::Hierarchy` has to `%import` IMP's kernel
 * interface, and the extension then calls
 * `PyImport_ImportModule("_IMP_kernel")` at load; behind a table, none of
 * that follows (PRD-139).
 */

%feature("kwargs") IMP::bff::read_structure_table;

%include "IMP/bff/StructureReader.h"
