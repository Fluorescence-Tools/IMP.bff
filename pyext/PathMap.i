/* Make selected classes extensible in Python */
IMP_SWIG_OBJECT(IMP::bff, PathMap, PathMaps);

// Use numpy.i for outputs of AV densities
%ignore IMP::bff::PathMap::get_tile_values(
    int value_type,
    std::pair<float, float> bounds,
    const std::string &feature_name
);
%ignore IMP::bff::PathMap::get_xyz_density();


// PathMapHeader is NOT declared IMP_SWIG_VALUE: PathMap returns it by
// pointer (get_path_map_header, get_path_map_header_writable) and takes it
// by non-const reference (set_path_map_header), which IMP's value machinery
// rejects at compile time. Making it a value is a public API change to
// PathMap, not a declaration -- see PRD-93.
IMP_SWIG_VALUE_SERIALIZE_IMPL(IMP::bff, PathMapHeader)

IMP_SWIG_VALUE(IMP::bff, PathMapTile, PathMapTiles)
IMP_SWIG_VALUE_SERIALIZE_IMPL(IMP::bff, PathMapTile)
IMP_SWIG_VALUE(IMP::bff, PathMapTileEdge, PathMapTileEdges)
IMP_SWIG_VALUE_SERIALIZE_IMPL(IMP::bff, PathMapTileEdge)
%include "IMP/bff/PathMapHeader.h"
%include "IMP/bff/PathMap.h"
%include "IMP/bff/PathMapTile.h"
%include "IMP/bff/PathMapTileEdge.h"

%template(VectorPathMapTile) std::vector<IMP::bff::PathMapTile>;
%template(VectorPathMapTileEdge) std::vector<IMP::bff::PathMapTileEdge>;
// %template(VectorIMPVector4D) std::vector<IMP::algebra::Vector4D>;
