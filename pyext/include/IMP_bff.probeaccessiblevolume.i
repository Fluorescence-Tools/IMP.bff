/*
 * The accessible volume as a representation of a label's states:
 * `ProbeAccessibleVolume` -> `ACV`, one description of one grid-enumerated cloud,
 * and the label distributions that produce one on demand. The base class,
 * `States`, and everything representation-neutral are in states.i, just
 * above.
 *
 * The constructors carry everything as flat vectors (via the shared
 * `const std::vector<double>&` typemap, which accepts a contiguous array of
 * any shape or `None` for an absent grid); the only conversion left is `ACV`'s
 * scalar `slow_radius`, which C++ broadcasts over the centres.
 * `LabelDistribution` and its two subclasses are C++ values whose scalar
 * state becomes native attributes; the two string fields (position_name,
 * simulation_type) stay methods, since `%attribute` cannot compile a
 * std::string getter.
 */

IMP_SWIG_VALUE(IMP::bff, ProbeAccessibleVolume, ProbeAccessibleVolumes);
IMP_SWIG_VALUE(IMP::bff, ACV, ACVs);

%attribute(IMP::bff::LabelDistribution, double, simulation_grid_resolution,
          get_simulation_grid_resolution);
%attribute(IMP::bff::LabelDistribution, int, n_points, get_n_points);
%attribute(IMP::bff::ProbeDistributionNormal, double, width, get_width);

// The `*_vector` accessors are the C++ side of the same buffers the numpy
// views publish. Wrapping them would give Python a second, slower way to ask
// the same question, so they stay out of the binding.
%ignore IMP::bff::ProbeAccessibleVolume::get_density_vector;
%ignore IMP::bff::ProbeAccessibleVolume::get_grid_origin_vector;

%include "IMP/bff/ProbeAccessibleVolume.h"

%attribute(IMP::bff::ProbeAccessibleVolume, int, ng, get_ng);
%attribute(IMP::bff::ProbeAccessibleVolume, double, grid_step, get_grid_step);
%attribute(IMP::bff::ACV, double, trapped_fraction, get_trapped_fraction);

// The keyed set of volumes the structure doors return (get_avs_for_structure,
// the network restraint's used volumes); instantiated here, where the value
// type is known, so every later consumer sees a typed map.
%template(MapStringProbeAccessibleVolume) std::map<std::string, IMP::bff::ProbeAccessibleVolume>;
