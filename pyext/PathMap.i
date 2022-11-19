/* Make selected classes extensible in Python */
IMP_SWIG_OBJECT(IMP::bff, PathMap, PathMaps);

// Use numpy.i for outputs of AV densities
%ignore IMP::bff::PathMap::get_tile_values(
    int value_type,
    std::pair<float, float> bounds,
    const std::string &feature_name
);

%include "IMP/bff/PathMapHeader.h"
%include "IMP/bff/PathMap.h"
%include "IMP/bff/PathMapTile.h"

%template(VectorPathMapTile) std::vector<IMP::bff::PathMapTile>;
