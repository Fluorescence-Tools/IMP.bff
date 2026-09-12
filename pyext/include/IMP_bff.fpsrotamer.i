/* FPSRotamer public bindings. */
IMP_SWIG_VALUE(IMP::bff, FPSRotamerPosition, FPSRotamerPositions);
IMP_SWIG_VALUE(IMP::bff, FPSRotamerDistance, FPSRotamerDistances);
IMP_SWIG_VALUE(IMP::bff, FPSRotamerSelection, FPSRotamerSelections);

%feature("kwargs") IMP::bff::fps_rotamer_position_payload;
%feature("kwargs") IMP::bff::fps_rotamer_distances_from_ensembles;
%feature("kwargs") IMP::bff::write_fps_rotamer;


%feature("kwargs") IMP::bff::fps_rotamer_fret;

%include "IMP/bff/FPSRotamer.h"

%template(VectorPairStringString) std::vector<std::pair<std::string, std::string> >;
