/*
 * The native Labelizer: per-residue label-site scoring, and the FRET pair
 * score built on it.
 *
 * All of it is C++ -- the structural features, the likelihood tables, the
 * weighted geometric mean and the pair layer. Nothing is left here but the
 * SWIG declarations the wrapper cannot infer: the value semantics of the
 * result structs, the containers they travel in, and the numpy array
 * properties.
 *
 * The `Labelizer` prefix and the `labelizer_` free functions mark what came from the
 * Labelizer package, so a reader can tell the ported model from the module's
 * own physics without consulting a PRD.
 */

/* Value semantics before the header that declares them: SWIG resolves a type
   at the point of use, and the containers below instantiate on these. */
IMP_SWIG_VALUE(IMP::bff, LabelizerResidue, LabelizerResidues);
IMP_SWIG_VALUE(IMP::bff, LabelizerTable, LabelizerTables);
IMP_SWIG_VALUE(IMP::bff, LabelizerParameter, LabelizerParameters);
IMP_SWIG_VALUE(IMP::bff, LabelizerScore, LabelizerScores);
IMP_SWIG_VALUE(IMP::bff, LabelizerFRETPairScore, LabelizerFRETPairScores);
/* The three Mfdb* value types are declared in swig.i-in now, ahead of
   ProbeRotamerLibrary.h, which needs MfdbTags for `write_probe_rotamer_drot_with_provenance`. */

%include "IMP/bff/LabelizerFeatures.h"
%include "IMP/bff/LabelizerScore.h"
%include "IMP/bff/LabelizerFRET.h"
%include "IMP/bff/LabelizerIO.h"

/* The containers the free functions above take and return. Instantiated after
   the headers, on types SWIG has already seen. */
%template(LabelizerResidueList) std::vector<IMP::bff::LabelizerResidue>;
%template(LabelizerParameterList) std::vector<IMP::bff::LabelizerParameter>;
%template(LabelizerScoreList) std::vector<IMP::bff::LabelizerScore>;
%template(LabelizerPairScoreList) std::vector<IMP::bff::LabelizerFRETPairScore>;
%template(MfdbTagList) std::vector<IMP::bff::MfdbTag>;
%template(MfdbColumnList) std::vector<IMP::bff::MfdbColumn>;
%template(MfdbAttributionList) std::vector<IMP::bff::MfdbAttribution>;
