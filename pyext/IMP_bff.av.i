%ignore IMP::bff::AVOccupancyMap::read_window;
%ignore IMP::bff::AVOccupancyMap::read_window_strided;
%ignore IMP::bff::AVOccupancyMap::set_coordinate_snapshot;
IMP_SWIG_OBJECT(IMP::bff, AVOccupancyMap, AVOccupancyMaps);
IMP_SWIG_OBJECT(IMP::bff, AVOccupancyRegistry, AVOccupancyRegistries);
IMP_SWIG_DECORATOR(IMP::bff, AV, AVs);
IMP_SWIG_OBJECT_SERIALIZE(IMP::bff, ProbeNetworkRestraint, ProbeNetworkRestraints);

// Keyword arguments from SWIG rather than from a hand-written shadow. The
// default constructor exists for cereal to deserialise into and has no
// business in Python -- a restraint over no hierarchy is not a restraint
// -- and hiding it leaves one wrapped constructor, which is what SWIG
// needs before it will generate keywords at all.
%ignore IMP::bff::ProbeNetworkRestraint::ProbeNetworkRestraint();
%feature("kwargs") IMP::bff::ProbeNetworkRestraint::ProbeNetworkRestraint;

// A `std::map<std::string, ParticleIndex>` has no template here and would
// arrive in Python as an opaque handle. The two accessors that answer the same
// questions -- `get_point_position_names()` and
// `get_position_particle_index()` -- are wrapped instead; the map is for C++.
%ignore IMP::bff::ProbeNetworkRestraint::get_point_positions;

%template(MapStringAVPairDistanceMeasurement) std::map<std::string, IMP::bff::AVPairDistanceMeasurement>;
%attribute_py(IMP::bff::AV, IMP::bff::PathMap, map, get_map);

IMP_SWIG_VALUE(IMP::bff, AVPairDistanceMeasurement, AVPairDistanceMeasurements)
IMP_SWIG_VALUE_SERIALIZE_IMPL(IMP::bff, AVPairDistanceMeasurement)

%include "IMP/bff/AVOccupancyMap.h"
%include "IMP/bff/AV.h"
%include "IMP/bff/ProbeNetworkRestraint.h"
