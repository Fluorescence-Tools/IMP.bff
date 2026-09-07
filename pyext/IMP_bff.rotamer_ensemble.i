/*
 * RotamerEnsemble -- last in swig.i-in, because it derives from `States`
 * (avmodel.i) and returns `FRETPairGeometry` / `FRETPairEfficiencies`.
 *
 * `atoms`, `energies`, `partition`, `library`, `chain` and `residue` are part
 * of the C++ value, not attributes hung on it, so an ensemble built in C++ is
 * one and a consumer that wants the atoms need not be Python.
 * `RotamerEnsemble::from_site` -- place the library in the residue's backbone
 * frame, screen it, take the chromophore centre and the transition dipole --
 * keeps its four kernel calls and their reshapes on the C++ side of the
 * boundary, because the shapes are what the kernels themselves say.
 */

IMP_SWIG_VALUE(IMP::bff, RotamerSiteOptions, RotamerSiteOptionsList);
IMP_SWIG_VALUE(IMP::bff, RotamerEnsemble, RotamerEnsembles);

// Keyword arguments for the two values a caller constructs by hand: an
// ensemble has thirteen fields and a site has six, and positional calls at
// that width are how a `sigma_scaling` ends up in `epsilon_scaling`.
%feature("kwargs") IMP::bff::RotamerEnsemble::RotamerEnsemble;
%feature("kwargs") IMP::bff::RotamerEnsemble::from_site;
%feature("kwargs") IMP::bff::RotamerSiteOptions::RotamerSiteOptions;

%include "IMP/bff/RotamerEnsemble.h"

// The shapes, on the C++ getters (see the note in avmodel.i): the flat view is
// what the C++ callers want and the table shape is what a numpy caller cannot
// be left to guess. `points` and `mu` come from `States` and are declared
// there, for every representation at once.
%attribute_np3(IMP::bff::RotamerEnsemble, std::vector<double>, atoms,
               get_atoms, n_rotamers, 3);
%attribute_np2(IMP::bff::RotamerEnsemble, std::vector<double>, centres,
               get_centres, 3);
%attribute_np(IMP::bff::RotamerEnsemble, std::vector<double>, weights,
              get_weights);
%attribute_np(IMP::bff::RotamerEnsemble, std::vector<double>, energies,
              get_energies);

// The name lists come back as the wrapped `std::vector<std::string>`, which
// indexes and iterates like a tuple.
%attribute_py(IMP::bff::RotamerEnsemble, std::vector<std::string>, atom_names,
              get_atom_names);
%attribute_py(IMP::bff::RotamerEnsemble, std::vector<std::string>, resnames,
              get_resnames);

%attribute(IMP::bff::RotamerEnsemble, int, n_rotamers, get_n_rotamers);
%attribute(IMP::bff::RotamerEnsemble, int, n_atoms, get_n_atoms);
%attribute(IMP::bff::RotamerEnsemble, double, partition, get_partition);
%attribute(IMP::bff::RotamerEnsemble, double, effective_sample_size,
           get_effective_sample_size);
%attribute(IMP::bff::RotamerEnsemble, int, residue, get_residue);
%attributestring(IMP::bff::RotamerEnsemble, std::string, library, get_library);
%attributestring(IMP::bff::RotamerEnsemble, std::string, chain, get_chain);

%template(MapStringRotamerEnsemble) std::map<std::string, IMP::bff::RotamerEnsemble>;

/*
 * The fps.json layer over the ensembles: which position is labelled with
 * what, and what distance a pair of ensembles predicts. Entries cross as JSON
 * text, which is what the rest of the fps layer does -- an entry carries
 * whatever keys its writer put there, and the alias reading (a chain is
 * `chain_identifier`, `chain` or `segid`) is the format's business, done once
 * in C++ rather than in every caller.
 */
IMP_SWIG_VALUE(IMP::bff, RotamerPosition, RotamerPositions);
IMP_SWIG_VALUE(IMP::bff, RotamerDistance, RotamerDistances);
IMP_SWIG_VALUE(IMP::bff, RotamerFpsSelection, RotamerFpsSelections);

%feature("kwargs") IMP::bff::rotamer_position_payload;
%feature("kwargs") IMP::bff::distances_from_ensembles;
%feature("kwargs") IMP::bff::write_rotamer_fps;

%include "IMP/bff/RotamerFps.h"

%template(VectorPairStringString) std::vector<std::pair<std::string, std::string> >;

/*
 * The FRETpredict driver: two libraries, two sites, every frame. It was the
 * last of `rotamer.i`'s Python -- a class holding twenty options, five numpy
 * arrays and the file writing, over kernels that were all C++ already. Its
 * `distance_distributions` never held anything: `trajectory_analysis`
 * allocated the array when `calc_distr` was set and nothing ever wrote to it,
 * so the option, the array and the `rmin`/`rmax`/`dr` axis it was shaped
 * from are gone rather than ported.
 *
 * The parameter names are FRETpredict's (`fixed_R0`, `ign_H`, `libname_1`),
 * deliberately: the parity harness hands one keyword dictionary to both.
 */
IMP_SWIG_VALUE(IMP::bff, FRETFrameResult, FRETFrameResults);

%feature("kwargs") IMP::bff::RotamerFRET::RotamerFRET;
%feature("kwargs") IMP::bff::RotamerFRET::from_frames;
%feature("kwargs") IMP::bff::RotamerFRET::reweight;
%feature("kwargs") IMP::bff::RotamerFRET::save;
%feature("kwargs") IMP::bff::rotamer_fret_from_fps;

%include "IMP/bff/RotamerFret.h"

%attribute_np2(IMP::bff::RotamerFRET, std::vector<double>, z_values,
               get_z_values, 2);
%attribute_np(IMP::bff::RotamerFRET, std::vector<double>, k2_values,
              get_k2_values);
%attribute_np(IMP::bff::RotamerFRET, std::vector<double>, estatic_values,
              get_estatic_values);
%attribute_np(IMP::bff::RotamerFRET, std::vector<double>, edynamic1_values,
              get_edynamic1_values);
%attribute_np(IMP::bff::RotamerFRET, std::vector<double>, edynamic2_values,
              get_edynamic2_values);
%attribute(IMP::bff::RotamerFRET, double, r0, get_r0);
// What the driver was told to do: the sites, the dyes and their libraries.
%attribute_py(IMP::bff::RotamerFRET, VectorInt, residues, get_residues);
%attribute_py(IMP::bff::RotamerFRET, std::vector<std::string>, chains, get_chains);
%attributestring(IMP::bff::RotamerFRET, std::string, donor, get_donor);
%attributestring(IMP::bff::RotamerFRET, std::string, acceptor, get_acceptor);
%attributestring(IMP::bff::RotamerFRET, std::string, libname_1,
                 get_libname_1);
%attributestring(IMP::bff::RotamerFRET, std::string, libname_2,
                 get_libname_2);
%attribute(IMP::bff::RotamerFRET, int, n_frames, get_n_frames);
%attributestring(IMP::bff::RotamerFRET, std::string, output_prefix,
                 get_output_prefix);
