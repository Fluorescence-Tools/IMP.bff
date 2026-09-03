IMP_SWIG_VALUE(IMP::bff, PathMapTileEdge, PathMapTileEdges);
IMP_SWIG_OBJECT(IMP::bff, PathMap, PathMaps);

// Use numpy.i for outputs of AV densities
%ignore IMP::bff::PathMap::get_tile_values(
    int value_type,
    std::pair<float, float> bounds,
    const std::string &feature_name
);
%ignore IMP::bff::PathMap::get_xyz_density();

/* Since PathMapHeaders are not values or Objects, we must ensure that whenever
   pointers to them are returned to Python, the object that manages the header's
   storage (usually a PathMap) must be kept alive so that the header is not
   prematurely freed. See modules/em/pyext/swig.i-in for a similar workaround
   for DensityHeader. */
namespace IMP {
 namespace bff {
  %feature("shadow") PathMap::get_path_map_header() const %{
    def get_path_map_header(self):
        h = _IMP_bff.DensityMap_get_path_map_header(self)
        h._owner = self
        return h
  %}
  %feature("shadow") PathMap::get_path_map_header_writable() %{
    def get_path_map_header_writable(self):
        h = _IMP_bff.PathMap_get_path_map_header_writable(self)
        h._owner = self
        return h
  %}
  }
}

%include "IMP/bff/PathMapHeader.h"
%include "IMP/bff/PathMap.h"
%include "IMP/bff/PathMapTile.h"
%include "IMP/bff/PathMapTileEdge.h"

%template(VectorPathMapTile) std::vector<IMP::bff::PathMapTile>;
%template(VectorPathMapTileEdge) std::vector<IMP::bff::PathMapTileEdge>;
// %template(VectorIMPVector4D) std::vector<IMP::algebra::Vector4D>;
