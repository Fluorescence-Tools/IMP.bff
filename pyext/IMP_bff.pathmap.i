/* Make selected classes extensible in Python */
IMP_SWIG_OBJECT(IMP::bff, PathMap, PathMaps);

// Use numpy.i for outputs of AV densities
%ignore IMP::bff::PathMap::get_tile_values(
    int value_type,
    std::pair<float, float> bounds,
    const std::string &feature_name
);
%ignore IMP::bff::PathMap::get_xyz_density();


// PathMapHeader is a value, and a value may be returned by value or const
// reference only -- never as a mutable reference, which would let Python edit
// the map's header behind its back. get_path_map_header_writable() stays in
// C++ (AV.cpp sets the path origin through it) but is not wrapped; from Python
// the value-correct pair is get_path_map_header() and set_path_map_header().
%ignore IMP::bff::PathMap::get_path_map_header_writable;

IMP_SWIG_VALUE(IMP::bff, PathMapHeader, PathMapHeaders)
IMP_SWIG_VALUE_SERIALIZE_IMPL(IMP::bff, PathMapHeader)

IMP_SWIG_VALUE(IMP::bff, PathMapTile, PathMapTiles)
IMP_SWIG_VALUE_SERIALIZE_IMPL(IMP::bff, PathMapTile)
IMP_SWIG_VALUE(IMP::bff, PathMapTileEdge, PathMapTileEdges)
IMP_SWIG_VALUE_SERIALIZE_IMPL(IMP::bff, PathMapTileEdge)
%include "IMP/bff/PathMap.h"

%template(VectorPathMapTile) std::vector<IMP::bff::PathMapTile>;
%template(VectorPathMapTileEdge) std::vector<IMP::bff::PathMapTileEdge>;
// %template(VectorIMPVector4D) std::vector<IMP::algebra::Vector4D>;
