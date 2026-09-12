/* FRETRotamer public bindings. */
IMP_SWIG_VALUE(IMP::bff, FRETRotamerFrameResult, FRETRotamerFrameResults);

%feature("kwargs") IMP::bff::FRETRotamer::FRETRotamer;
%feature("kwargs") IMP::bff::FRETRotamer::from_frames;
%feature("kwargs") IMP::bff::FRETRotamer::reweight;
%feature("kwargs") IMP::bff::FRETRotamer::save;



%include "IMP/bff/FRETRotamer.h"

%attribute_np2(IMP::bff::FRETRotamer, std::vector<double>, z_values,
               get_z_values, 2);
%attribute_np(IMP::bff::FRETRotamer, std::vector<double>, k2_values,
              get_k2_values);
%attribute_np(IMP::bff::FRETRotamer, std::vector<double>, estatic_values,
              get_estatic_values);
%attribute_np(IMP::bff::FRETRotamer, std::vector<double>, edynamic1_values,
              get_edynamic1_values);
%attribute_np(IMP::bff::FRETRotamer, std::vector<double>, edynamic2_values,
              get_edynamic2_values);
%attribute(IMP::bff::FRETRotamer, double, r0, get_r0);
// What the driver was told to do: the sites, the dyes and their libraries.
%attribute_py(IMP::bff::FRETRotamer, VectorInt, residues, get_residues);
%attribute_py(IMP::bff::FRETRotamer, std::vector<std::string>, chains, get_chains);
%attributestring(IMP::bff::FRETRotamer, std::string, donor, get_donor);
%attributestring(IMP::bff::FRETRotamer, std::string, acceptor, get_acceptor);
%attributestring(IMP::bff::FRETRotamer, std::string, libname_1,
                 get_libname_1);
%attributestring(IMP::bff::FRETRotamer, std::string, libname_2,
                 get_libname_2);
%attribute(IMP::bff::FRETRotamer, int, n_frames, get_n_frames);
%attributestring(IMP::bff::FRETRotamer, std::string, output_prefix,
                 get_output_prefix);
