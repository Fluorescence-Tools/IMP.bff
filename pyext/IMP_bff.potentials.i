/*
 * The coarse-grained protein potentials, as IMP scores and restraints.
 *
 * The two typed contact scores are `IMP::core::StatisticalPairScore`
 * subclasses, so their template bases have to be instantiated before the
 * header is read -- the same three lines `IMP.atom` writes for
 * `DopePairScore` (modules/atom/pyext/swig.i-in:252-260). Everything else is
 * an ordinary object or a free function.
 */

// The residue types the contact tables are indexed by. A key family of its
// own, like `IMP::atom::DopeType` -- IMP keeps those two C++-only, but a
// caller wants to ask whether a residue is one the tables cover, and
// `ResidueContactType.get_key_exists(name)` is that question.
IMP_SWIG_VALUE_INSTANCE(IMP::bff, ResidueContactType, ResidueContactType,
                        ResidueContactTypes);

// The parameter tables, before the scores that read them.
IMP_SWIG_VALUE(IMP::bff, PotentialTable, PotentialTableList);
%include "IMP/bff/PotentialTables.h"
%attribute_np(IMP::bff::PotentialTable, std::vector<double>, table,
              get_values);
%template(PotentialTableVector) std::vector<IMP::bff::PotentialTable>;

IMP_SWIG_OBJECT(IMP::bff, MiyazawaJerniganPairScore,
                MiyazawaJerniganPairScores);
IMP_SWIG_OBJECT(IMP::bff, UNRESCentroidPairScore, UNRESCentroidPairScores);
IMP_SWIG_OBJECT(IMP::bff, LennardJonesBeadPairScore,
                LennardJonesBeadPairScores);
IMP_SWIG_OBJECT(IMP::bff, GoRestraint, GoRestraints);
IMP_SWIG_OBJECT(IMP::bff, HydrogenBondRestraint, HydrogenBondRestraints);
IMP_SWIG_OBJECT(IMP::bff, GeneralizedBornRestraint,
                GeneralizedBornRestraints);
IMP_SWIG_OBJECT(IMP::bff, RamachandranRestraint, RamachandranRestraints);
IMP_SWIG_VALUE(IMP::bff, GoContacts, GoContactsList);

// The template bases of the two statistical scores.
%template(_BffStatisticalResidueDistance)
        IMP::score_functor::DistancePairScore<
                IMP::score_functor::Statistical<IMP::bff::ResidueContactType,
                                                false, false, false> >;
%template(_BffStatisticalResidue)
        ::IMP::core::StatisticalPairScore<IMP::bff::ResidueContactType, false,
                                          false>;

%include "IMP/bff/Potentials.h"

%template(ResidueContactType) ::IMP::Key<8064531>;

// The contact list as numpy, so a caller can see which pairs a Go model took
// and how deep their wells are.
%attribute_np(IMP::bff::GoContacts, std::vector<double>, well_depths,
              get_energies);
%attribute_np(IMP::bff::GoContacts, std::vector<double>, well_minima,
              get_distances);

// The backbone dihedrals a Ramachandran restraint reads, radians.
%attribute_np(IMP::bff::RamachandranRestraint, std::vector<double>, phi,
              get_phi);
%attribute_np(IMP::bff::RamachandranRestraint, std::vector<double>, psi,
              get_psi);
