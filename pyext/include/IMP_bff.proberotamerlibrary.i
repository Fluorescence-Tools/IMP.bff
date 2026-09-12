/* Probe conformers and their .drot representation. */
IMP_SWIG_VALUE(IMP::bff, ProbeRotamerLibrary, ProbeRotamerLibraries);
%ignore IMP::bff::ProbeRotamerLibrary::coords;
%ignore IMP::bff::ProbeRotamerLibrary::weights;
%include "IMP/bff/ProbeRotamerLibrary.h"
%attribute_np3(IMP::bff::ProbeRotamerLibrary, std::vector<double>, coords,
               get_coords, n_rotamers, 3, set_coords);
%attribute_np(IMP::bff::ProbeRotamerLibrary, std::vector<double>, weights,
              get_weights, set_weights);
