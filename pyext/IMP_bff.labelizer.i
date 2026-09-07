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
 * The `Ll` prefix and the `ll_` free functions mark what came from the
 * Labelizer package, so a reader can tell the ported model from the module's
 * own physics without consulting a PRD.
 */

/* Value semantics before the header that declares them: SWIG resolves a type
   at the point of use, and the containers below instantiate on these. */
IMP_SWIG_VALUE(IMP::bff, LlResidue, LlResidues);
IMP_SWIG_VALUE(IMP::bff, LlTable, LlTables);
IMP_SWIG_VALUE(IMP::bff, LlParameter, LlParameters);
IMP_SWIG_VALUE(IMP::bff, LlScore, LlScores);
IMP_SWIG_VALUE(IMP::bff, LlPairScore, LlPairScores);
/* The three Mfdb* value types are declared in swig.i-in now, ahead of
   RotamerLibrary.h, which needs MfdbTags for `write_drot_with_provenance`. */

%include "IMP/bff/LabelizerFeatures.h"
%include "IMP/bff/LabelizerScore.h"
%include "IMP/bff/LabelizerFret.h"
%include "IMP/bff/LabelizerIO.h"
%include "IMP/bff/ProbeContainer.h"

/* The containers the free functions above take and return. Instantiated after
   the headers, on types SWIG has already seen. */
%template(LlResidueList) std::vector<IMP::bff::LlResidue>;
%template(LlParameterList) std::vector<IMP::bff::LlParameter>;
%template(LlScoreList) std::vector<IMP::bff::LlScore>;
%template(LlPairScoreList) std::vector<IMP::bff::LlPairScore>;
%template(MfdbTagList) std::vector<IMP::bff::MfdbTag>;
%template(MfdbColumnList) std::vector<IMP::bff::MfdbColumn>;
%template(MfdbAttributionList) std::vector<IMP::bff::MfdbAttribution>;
