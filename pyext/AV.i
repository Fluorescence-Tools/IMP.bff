%ignore IMP::bff::AVOccupancyMap::read_window;
%ignore IMP::bff::AVOccupancyMap::read_window_strided;
IMP_SWIG_OBJECT(IMP::bff, AVOccupancyMap, AVOccupancyMaps);
IMP_SWIG_OBJECT(IMP::bff, AVOccupancyRegistry, AVOccupancyRegistries);
IMP_SWIG_DECORATOR(IMP::bff, AV, AVs);
IMP_SWIG_OBJECT_SERIALIZE(IMP::bff, AVNetworkRestraint, AVNetworkRestraints);

// The restraint has two constructors (the default one exists for
// deserialization), and SWIG cannot generate keyword arguments for an
// overloaded function. This shadow restores keyword arguments for the real
// constructor -- the parameter order below is the C++ order.
%feature("shadow") IMP::bff::AVNetworkRestraint::AVNetworkRestraint %{
def __init__(self, *args, **kwargs):
    if not args and not kwargs:
        _IMP_bff.AVNetworkRestraint_swiginit(self, _IMP_bff.new_AVNetworkRestraint())
        return
    names = ("hier", "fps_json_fn", "name", "score_set", "n_samples",
             "space_fixed", "shared_map", "distance", "quad_k",
             "search_grid_factor", "search_stencil")
    defaults = {"name": "AVNetworkRestraint%1%", "score_set": "",
                "n_samples": 50000, "space_fixed": True, "shared_map": True,
                "distance": "quad", "quad_k": 50, "search_grid_factor": 1,
                "search_stencil": 26}
    if len(args) > len(names):
        raise TypeError("AVNetworkRestraint() takes at most %d positional "
                        "arguments (%d given)" % (len(names), len(args)))
    values = dict(zip(names, args))
    for k, v in kwargs.items():
        if k not in names:
            raise TypeError("AVNetworkRestraint() got an unexpected keyword "
                            "argument %r" % k)
        if k in values:
            raise TypeError("AVNetworkRestraint() got multiple values for "
                            "argument %r" % k)
        values[k] = v
    for k in names:
        if k not in values:
            if k not in defaults:
                raise TypeError("AVNetworkRestraint() missing required "
                                "argument %r" % k)
            values[k] = defaults[k]
    _IMP_bff.AVNetworkRestraint_swiginit(
        self, _IMP_bff.new_AVNetworkRestraint(*[values[k] for k in names]))
%}

%template(MapStringAVPairDistanceMeasurement) std::map<std::string, IMP::bff::AVPairDistanceMeasurement>;
%attribute_py(IMP::bff::AV, IMP::bff::PathMap, map, get_map);

IMP_SWIG_VALUE(IMP::bff, AVPairDistanceMeasurement, AVPairDistanceMeasurements)
IMP_SWIG_VALUE_SERIALIZE_IMPL(IMP::bff, AVPairDistanceMeasurement)

%include "IMP/bff/AVOccupancyMap.h"
%include "IMP/bff/AV.h"
%include "IMP/bff/AVNetworkRestraint.h"
