/* ProbeRotamer public bindings. */
IMP_SWIG_VALUE(IMP::bff, ProbeRotamerSiteOptions, ProbeRotamerSiteOptionsList);
IMP_SWIG_VALUE(IMP::bff, ProbeRotamerEnsemble, ProbeRotamerEnsembles);

// Keyword arguments for the two values a caller constructs by hand: an
// ensemble has thirteen fields and a site has six, and positional calls at
// that width are how a `sigma_scaling` ends up in `epsilon_scaling`.
%feature("kwargs") IMP::bff::ProbeRotamerEnsemble::ProbeRotamerEnsemble;
%feature("kwargs") IMP::bff::ProbeRotamerEnsemble::from_site;
%feature("kwargs") IMP::bff::ProbeRotamerSiteOptions::ProbeRotamerSiteOptions;


%include "IMP/bff/ProbeRotamer.h"

// The shapes, on the C++ getters (see the note in avmodel.i): the flat view is
// what the C++ callers want and the table shape is what a numpy caller cannot
// be left to guess. `points` and `mu` come from `States` and are declared
// there, for every representation at once.
%attribute_np3(IMP::bff::ProbeRotamerEnsemble, std::vector<double>, atoms,
               get_atoms, n_rotamers, 3);
%attribute_np2(IMP::bff::ProbeRotamerEnsemble, std::vector<double>, centres,
               get_centres, 3);
%attribute_np(IMP::bff::ProbeRotamerEnsemble, std::vector<double>, weights,
              get_weights);
%attribute_np(IMP::bff::ProbeRotamerEnsemble, std::vector<double>, energies,
              get_energies);

// The name lists come back as the wrapped `std::vector<std::string>`, which
// indexes and iterates like a tuple.
%attribute_py(IMP::bff::ProbeRotamerEnsemble, std::vector<std::string>, atom_names,
              get_atom_names);
%attribute_py(IMP::bff::ProbeRotamerEnsemble, std::vector<std::string>, resnames,
              get_resnames);

%attribute(IMP::bff::ProbeRotamerEnsemble, int, n_rotamers, get_n_rotamers);
%attribute(IMP::bff::ProbeRotamerEnsemble, int, n_atoms, get_n_atoms);
%attribute(IMP::bff::ProbeRotamerEnsemble, double, partition, get_partition);
%attribute(IMP::bff::ProbeRotamerEnsemble, double, effective_sample_size,
           get_effective_sample_size);
%attribute(IMP::bff::ProbeRotamerEnsemble, int, residue, get_residue);
%attributestring(IMP::bff::ProbeRotamerEnsemble, std::string, library, get_library);
%attributestring(IMP::bff::ProbeRotamerEnsemble, std::string, chain, get_chain);

%template(MapStringProbeRotamerEnsemble) std::map<std::string, IMP::bff::ProbeRotamerEnsemble>;
