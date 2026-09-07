/*
 * An atom selection written as text.
 *
 * `chain A and resid 10-20 and not name CA+CB` compiles to an
 * `IMP::atom::Selection`, so a selection a user types composes with everything
 * that takes one. Both spellings of the language are accepted -- `resi`/`resid`,
 * `name CA+CB`/`name CA CB`, `50-60`/`50 to 60`, `byres`/`same residue as`.
 */

IMP_SWIG_VALUE(IMP::bff, SelectionAtom, SelectionAtoms);
IMP_SWIG_VALUE(IMP::bff, SelectionExpression, SelectionExpressions);

%include "IMP/bff/SelectionExpression.h"

%template(SelectionAtomVector) std::vector<IMP::bff::SelectionAtom>;
