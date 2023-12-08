/* Make selected classes extensible in Python */
IMP_SWIG_OBJECT(IMP::bff, PathMap, PathMaps);

// Use numpy.i for outputs of AV densities
%ignore IMP::bff::PathMap::get_tile_values(
    int value_type,
    std::pair<float, float> bounds,
    const std::string &feature_name
);
%ignore IMP::bff::PathMap::get_xyz_density();


%include "IMP/bff/PathMapHeader.h"
%include "IMP/bff/PathMap.h"
%include "IMP/bff/PathMapTile.h"
%include "IMP/bff/PathMapTileEdge.h"

%template(VectorPathMapTile) std::vector<IMP::bff::PathMapTile>;
%template(VectorPathMapTileEdge) std::vector<IMP::bff::PathMapTileEdge>;
// %template(VectorIMPVector4D) std::vector<IMP::algebra::Vector4D>;
