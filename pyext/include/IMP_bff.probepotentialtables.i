/*
 * The parameter tables the coarse-grained potentials read (`ProbePotentialTables.h`,
 * core), before the scores that read them (potentials.i, connection layer).
 */

// The parameter tables, before the scores that read them.
IMP_SWIG_VALUE(IMP::bff, PotentialTable, PotentialTableList);
%include "IMP/bff/ProbePotentialTables.h"
%attribute_np(IMP::bff::PotentialTable, std::vector<double>, table,
              get_values);
%template(PotentialTableVector) std::vector<IMP::bff::PotentialTable>;
