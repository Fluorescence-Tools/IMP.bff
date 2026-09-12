/** \file IMP/bff/ProbeRotamer.h
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#ifndef IMPBFF_PROBEROTAMER_H
#define IMPBFF_PROBEROTAMER_H


// -------- from RotamerEnergy.h --------
/**
 *  (formerly IMP/bff/RotamerEnergy.h, now a section of this file)
 *  \brief The interaction energy of each rotamer with its protein.
 *
 * A rotamer library is a list of candidate dye conformers; scoring turns it
 * into a weighted ensemble by asking how hard each one collides with the
 * structure it is attached to. That is an all-pairs problem — every dye atom
 * against every protein atom, for every conformer — and it is the inner loop
 * of `FRETRotamer`.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#include <IMP/bff/bff_config.h>

#include <vector>

IMPBFF_BEGIN_NAMESPACE

struct FRETPairGeometry;
struct FRETPairEfficiencies;

//! Which functional form the steric term takes.
enum ProbeRotamerPotential {
    PROBE_ROTAMER_POTENTIAL_LJ = 0,     //!< 12-6 Lennard-Jones
    PROBE_ROTAMER_POTENTIAL_GAUSS = 1   //!< a soft Gaussian well, for a coarse dye
};

//! Steric and electrostatic energy of every conformer against the structure.
/*!
    \param[in] rotamer_coords,n_rotamer_coords flat,
               `n_rotamers * n_probe_atoms * 3`. A raw buffer, so numpy's array
               passes straight through -- see #brownian_walk_in_volume on why.
    \param[in] protein_coords,n_protein_coords flat, `n_protein_atoms * 3`
    \param[in] rmin_ij combined \f$R_{min}\f$ per (dye atom, protein atom),
               flat `n_probe_atoms * n_protein_atoms`
    \param[in] eps_ij combined well depth, same shape
    \param[in] q_probe per dye atom, or **empty** to skip electrostatics
    \param[in] q_protein per protein atom, or empty
    \param[in] n_rotamers,n_probe_atoms,n_protein_atoms shapes
    \param[in] potential a #ProbeRotamerPotential
    \param[in] steric_cutoff pairs beyond this contribute nothing, Angstrom
    \param[in] coulomb_cutoff the same for the screened Coulomb term
    \param[in] debye_length screening length of the Debye-Huckel term, Angstrom
    \param[in] coulomb_prefactor the `q_i q_j / r` coefficient
    \return two values per rotamer: steric energy, then electrostatic energy.
            Separate rather than summed because the Boltzmann factor treats them
            differently -- the steric term is divided by RT and the screened
            Coulomb term is already in units of it.
*/
IMPBFFEXPORT std::vector<double> probe_rotamer_interaction_energies(
        double* rotamer_coords, int n_rotamer_coords,
        double* protein_coords, int n_protein_coords,
        double* rmin_ij, int n_rmin_ij,
        double* eps_ij, int n_eps_ij,
        const std::vector<double>& q_probe,
        const std::vector<double>& q_protein,
        int n_rotamers, int n_probe_atoms, int n_protein_atoms,
        int potential,
        double steric_cutoff = 10.0,
        double coulomb_cutoff = 20.0,
        double debye_length = 10.0,
        double coulomb_prefactor = 7.0);

//! Repulsive-only LJ energy between every pair of conformers of two sets.
/*!
    The mean-field weight update needs the interaction energy of conformer *i*
    of one dye against conformer *j* of another (or against a protein, which is
    one "conformer"). Those energies **do not depend on the weights**, so they
    are computed once here and the iteration becomes a matrix-vector product.
    Recomputing the matrix inside the iteration loop instead is ten times the
    work for the same answer.

    Repulsive-only: the attractive tail (\f$r \ge R_{min}\f$) is dropped, which
    is what the mean-field treatment wants -- it is asking what is *blocked*,
    not what is bound.

    Conformer pairs whose padded axis-aligned bounding boxes do not overlap are
    skipped entirely, which is most of them once two dyes are more than a linker
    apart.

    \param[in] coords_a flat, `n_a_conf * n_a_atoms * 3`
    \param[in] coords_b flat, `n_b_conf * n_b_atoms * 3`
    \param[in] rmin,eps combined parameters, flat `n_a_atoms * n_b_atoms`
    \param[in] r_cutoff pairs beyond this contribute nothing, Angstrom
    \param[in] aabb_pad padding on each bounding box, Angstrom
    \param[in] r_floor distance clamp keeping the energy finite
    \param[out] out_view,n_out_view flat `n_a_conf * n_b_conf`. Placed before
                the defaulted parameters, not after: a parameter without a
                default may not follow one that has it, and numpy's typemap
                binds on the pair's *names*, not its position.
*/
IMPBFFEXPORT void pair_energy_matrix_kernel(
        const std::vector<double>& coords_a,
        const std::vector<double>& coords_b,
        const std::vector<double>& rmin,
        const std::vector<double>& eps,
        int n_a_conf, int n_a_atoms, int n_b_conf, int n_b_atoms,
        double** out_view, int* n_out_view,
        double r_cutoff = 12.0,
        double aabb_pad = 3.5,
        double r_floor = 0.01);

//! Internal LJ energy of each frame, over an explicit list of atom pairs.
/*!
    The bonded exclusions have already decided which atom pairs interact, so
    this takes the list rather than rediscovering it. One energy per frame.

    The vectorised form it replaces gathered `(n_frames, n_pairs, 3)` twice --
    once for each end of every pair -- before taking a single difference. At ten
    thousand frames and twelve hundred pairs that is 288 MB per gather.

    \param[in] coords,n_coords flat, `n_frames * n_atoms * 3`, zero-copy
    \param[in] index_a,index_b atom indices of each pair
    \param[in] rmin,eps combined parameters per pair
    \param[in] n_frames,n_atoms,n_pairs shapes
    \param[in] repulsive_only drop the attractive tail beyond \f$R_{min}\f$
    \param[in] r_floor distance clamp keeping the energy finite
    \return one energy per frame
*/
IMPBFFEXPORT std::vector<double> lj_pair_energies(
        double* coords, int n_coords,
        const std::vector<int>& index_a,
        const std::vector<int>& index_b,
        const std::vector<double>& rmin,
        const std::vector<double>& eps,
        int n_frames, int n_atoms, int n_pairs,
        bool repulsive_only = true,
        double r_floor = 0.01);

IMPBFF_END_NAMESPACE

// -------- from RotamerSite.h --------
/**
 *  (formerly IMP/bff/RotamerSite.h, now a section of this file)
 *  \brief Placing a rotamer library on a residue: the backbone frame, the
 *         atom selectors, and the library registry.
 *
 * The transforms and the registry resolution were Python in
 * `representation/rotamer.py`. They are kernels and file logic, not glue:
 * the backbone frame is nine multiplications a caller should not re-type,
 * the selector matching has the same ambiguity rules as the scoring one, and
 * the registry is a JSON lookup with a cutoff-suffix grammar.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#include <IMP/bff/ProbeRotamerLibrary.h>

#include <string>

IMPBFF_BEGIN_NAMESPACE

// --------------------------------------------------------------------------
// The backbone frame
// --------------------------------------------------------------------------

//! CA, N, C coordinates of `(chain, residue)` in a protein frame.
/*!
    The first atom of each name that matches: same chain when both name one
    (case-blind), and the residue number when the frame carries one. A missing
    backbone atom raises -- a frame that cannot place a library is not a frame
    a caller can guess around.

    \param[in] coords flat, three per atom
    \param[in] atom_names,chain_ids one per atom; chain_ids may be empty
    \param[in] residue_indices one per atom, -1 where unknown; may be empty
    \param[in] chain the chain asked for; empty matches any
    \param[in] residue the residue number asked for
    \param[out] out_view,n_out_view nine values: CA, then N, then C
    \throw ValueException when CA, N or C is missing
*/
IMPBFFEXPORT void resolve_backbone_site(
        const std::vector<double>& coords,
        const std::vector<std::string>& atom_names,
        const std::vector<std::string>& chain_ids,
        const std::vector<int>& residue_indices,
        const std::string& chain, int residue,
        double** out_view, int* n_out_view);

//! The site frame's rotation: rows are x (CA->N), y (in the N-CA-C plane),
//! z = x cross y.
/*! \param[out] out_view,n_out_view nine values, row-major */
IMPBFFEXPORT void backbone_rotation(
        const std::vector<double>& ca, const std::vector<double>& n,
        const std::vector<double>& c,
        double** out_view, int* n_out_view);

//! Library coordinates into the backbone frame at CA.
/*!
    \param[in] coords flat, `n_rotamers * n_atoms * 3`, the library's own frame
    \param[in] ca,n,c the site's backbone atoms, three values each
    \param[out] out_view,n_out_view flat, the input's shape, in the protein's
               frame: `rotated + ca`
*/
IMPBFFEXPORT void transform_library_to_site(
        const std::vector<double>& coords,
        const std::vector<double>& ca, const std::vector<double>& n,
        const std::vector<double>& c,
        double** out_view, int* n_out_view);

// --------------------------------------------------------------------------
// Atom selectors
// --------------------------------------------------------------------------

//! Indices of the atoms named by FRETpredict selectors.
/*!
    A selector is `NAME` or `NAME and resname RES`. The resname clause is
    **honoured** when \p resnames is given: atom names repeat between the dye
    residue and its linker -- `C13` is in both `A48` and `C1R` -- and a
    selector that ignores the residue takes whichever comes first in the atom
    ordering, which is the dye today only by luck. Without \p resnames the
    clause cannot be checked and the first name match is returned.

    \param[in] atom_names one per atom, uppercased internally
    \param[in] selectors one entry per wanted atom
    \param[in] resnames one per atom, or empty to skip the clause check
    \return one index per selector entry, in order
    \throw ValueException naming the selector when one matches no atom
*/
IMPBFFEXPORT std::vector<int> selector_atom_indices(
        const std::vector<std::string>& atom_names,
        const std::vector<std::string>& selectors,
        const std::vector<std::string>& resnames = std::vector<std::string>());

// `selector_resnames(selectors)` lives in RotamerScoring.h beside the rest of the
// selector mini-language.

// --------------------------------------------------------------------------
// The bundled library registry
// --------------------------------------------------------------------------

//! The registry key of a library name: the name without its cutoff suffix.
/*! `'AlexaFluor 488 C1R cutoff30'` -> `'AlexaFluor 488 C1R'`. */
IMPBFFEXPORT std::string normalize_library_name(const std::string& name);

//! The cutoff of a library name, or -1 when the name carries none.
IMPBFFEXPORT int library_name_cutoff(const std::string& name);

//! The bundled library registry, as JSON text.
/*! `libraries.json` read once and cached. The whole table, for a caller that
    wants to list what is available rather than ask about one name. */
IMPBFFEXPORT std::string probe_rotamer_library_registry();

//! One registry entry of `libraries.json`, plus the resolved spelling.
/*! The returned JSON object is the registry's own entry with `name` (the
    key), `library_name` (the name as asked) and `cutoff` added. The registry
    is read once and cached. An unknown library raises. */
IMPBFFEXPORT std::string probe_rotamer_library_metadata(const std::string& name);

//! The registry entry a library *file* belongs to, or `{}` for none.
/*! The inverse of \c library_filename: a path names a library by its stem,
    so `A48_C1R_cutoff30.drot`, `A48_C1R_cutoff30.bcif` and the `A48_C1R.rmf3`
    template all resolve to Alexa488 C1R -- the first two carrying `cutoff`
    30, the third none. This is how an explicit path still knows its *dye*:
    the coordinates are in the file, but the transition-dipole and attachment
    selectors are in the registry, and without them an ensemble cannot orient
    itself.

    A library outside the registry -- a user's own dye -- returns the empty
    object `{}` rather than raising, because that is not an error: the caller
    supplies the selectors instead. */
IMPBFFEXPORT std::string probe_rotamer_library_metadata_for_path(
        const std::string& path);

//! The library file stem the metadata and cutoff select.
/*! The registry's `filename` carries FRETpredict's default cutoff; a name
    that asks for another cutoff replaces it. */
IMPBFFEXPORT std::string library_filename(const std::string& metadata_json,
                                          int cutoff);

//! Resolve a library name (or an explicit path) to a library file.
/*!
    The FRETpredict library files (module data, `data/rotamer_library`) are
    canonical: `<stem>.pdb` + `<stem>_cutoff<N>.bcif` (+ weights) per cutoff.
    They are tried first so the *requested cutoff* is the one loaded -- the
    RMF templates hold only the cutoff-30 clustering, so resolving every name
    to `<stem>.rmf3` silently returned the wrong library for cutoff10/20.

    \param[in] name a registry name, or a path that exists
    \param[in] lib_dir an extra directory to search after the canonical one;
               empty searches the module's template directory only
    \return the existing file path
    \throw IOException when nothing matches, or when only the cutoff-30 RMF
           template exists for another cutoff (that mismatch once answered
           quietly)
*/
IMPBFFEXPORT std::string resolve_probe_rotamer_library_path(
        const std::string& name, const std::string& lib_dir = "");

//! Read the library a name (or a path) resolves to, with its metadata.
/*!
    The one door to a rotamer library: it resolves the name, reads whichever
    container the resolution lands on, attaches the registry entry, and fills
    the residue names -- from the container when it carries them, from the
    `<stem>.pdb` beside it when it does not, and by inference from the dye and
    linker names as a last resort. Weights come back normalised whatever the
    file stored.

    \param[in] name a registry name (`"AlexaFluor 488 C1R cutoff10"`), a
               locator (`"dyes.drot.pto::A48_C1R_cutoff10"`) or a path
    \param[in] lib_dir an extra directory to search; empty searches the
               module's own
    \throw IOException when nothing resolves, or the container cannot be read
    \throw ValueException for a container this module does not read -- an
           `.rmf3` template is one, because reading it needs `IMP.rmf`
*/
IMPBFFEXPORT ProbeRotamerLibrary load_probe_rotamer_library(const std::string& name,
                                                 const std::string& lib_dir =
                                                         "");

//! Residue names per atom, inferred from a library's registry metadata.
/*!
    The fallback for a library whose container carries no residue names and
    has no `<stem>.pdb` beside it. The dye's own residue name comes from its
    selectors; the linker's from the registry name's `_C1R` suffix, whose
    atoms are a fixed set. Empty when the metadata names no residue at all --
    an inference with nothing to infer from is a guess, and a wrong residue
    name silently selects the wrong atom.
*/
IMPBFFEXPORT std::vector<std::string> infer_probe_rotamer_resnames(
        const std::vector<std::string>& atom_names,
        const std::string& metadata_json);

IMPBFF_END_NAMESPACE

// -------- from ProbeRotamerEnsemble.h --------
/**
 *  (formerly IMP/bff/ProbeRotamerEnsemble.h, now a section of this file)
 *  \brief A rotamer library placed and screened at one labelling site.
 *
 * The fps.json position type `R1` (PRD-108): a FRETpredict-style 1:1 rotamer
 * library transformed into a residue's backbone frame and Boltzmann-screened
 * against the protein. Per rotamer the chromophore centre, the transition
 * dipole and the weight are kept -- and every atom, so nothing downstream has
 * to re-derive them -- which is what a FRET *rate* distribution over a pair of
 * labels needs (\f$R_{ij}\f$, \f$\kappa^2_{ij}\f$, \f$w_i w_j\f$), and not
 * just a mean position.
 *
 * A **sibling** of #IMP::bff::ProbeAccessibleVolume, not a subclass: both are
 * #IMP::bff::States, so every distance and FRET helper takes either, while an
 * ensemble is spared the grid an AV has and it does not.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#include <IMP/bff/States.h>

#include <IMP/bff/HierarchyFrame.h>

#include <map>

IMPBFF_BEGIN_NAMESPACE

//! fps.json `simulation_type` of a screened 1:1 rotamer library.
IMPBFFEXPORT extern const char* const SIMULATION_TYPE_R1;

// --------------------------------------------------------------------------
// The ensemble
// --------------------------------------------------------------------------

//! How a library is screened at a site: the scoring parameters.
/*! FRETpredict's, with its defaults. They are a value rather than eight
    arguments because they travel together -- a site, a pair and a whole fps
    file are all screened with one setting of them. */
struct IMPBFFEXPORT ProbeRotamerSiteOptions {
    //! K.
    double temperature;
    //! Add the Debye-Hueckel term from the dye's charge selectors.
    bool electrostatic;
    //! `"lj"` or `"gauss"`.
    std::string potential;
    //! Mask hydrogens out of the steric term.
    bool ignore_h;
    //! Probe-probe scaling of the LJ parameters.
    double sigma_scaling, epsilon_scaling;

    ProbeRotamerSiteOptions(double temperature = 298.15,
                       bool electrostatic = false,
                       const std::string& potential = "lj",
                       bool ignore_h = true, double sigma_scaling = 0.5,
                       double epsilon_scaling = 1.0)
        : temperature(temperature), electrostatic(electrostatic),
          potential(potential), ignore_h(ignore_h),
          sigma_scaling(sigma_scaling), epsilon_scaling(epsilon_scaling) {}

    IMP_SHOWABLE_INLINE(ProbeRotamerSiteOptions,
                        out << "RotamerSiteOptions(" << potential << ", T="
                            << temperature << ")");
};
IMP_VALUES(ProbeRotamerSiteOptions, ProbeRotamerSiteOptionsList);

//! A rotamer library placed and screened at one labelling site.
/*! From #IMP::bff::States: `points` (the chromophore centre and the weight of
    each rotamer), `attachment_point` (the site's CA), `orientations` (the
    per-rotamer transition dipoles), `position_name` and `params`. Added here:
    every atom of every rotamer, the atom and residue names they are in the
    order of, the interaction energies the weights came from, and the
    partition function. */
class IMPBFFEXPORT ProbeRotamerEnsemble : public States {
    std::vector<double> atoms_;
    std::vector<std::string> atom_names_, resnames_;
    std::vector<double> energies_;
    double partition_;
    std::string library_, chain_;
    int residue_;

public:
    //! \param[in] points flat (x, y, z, w) per rotamer -- centre and weight
    /*! \param[in] attachment_point the site's CA
        \param[in] orientations flat transition dipoles, three per rotamer
        \param[in] position_name the site's name, e.g. "A344"
        \param[in] params how the ensemble was produced
        \param[in] atoms flat, `n_rotamers * n_atoms * 3`, in the protein frame
        \param[in] atom_names,resnames one per atom of a rotamer
        \param[in] energies one per rotamer, the interaction energy scored
        \param[in] partition the Boltzmann partition function Z
        \param[in] library the library's *name* -- never a rendering of it
        \param[in] chain,residue the site */
    ProbeRotamerEnsemble(
            const std::vector<double>& points = std::vector<double>(),
            const std::vector<double>& attachment_point =
                    std::vector<double>(),
            const std::vector<double>& orientations = std::vector<double>(),
            const std::string& position_name = "",
            const std::map<std::string, std::string>& params =
                    std::map<std::string, std::string>(),
            const std::vector<double>& atoms = std::vector<double>(),
            const std::vector<std::string>& atom_names =
                    std::vector<std::string>(),
            const std::vector<std::string>& resnames =
                    std::vector<std::string>(),
            const std::vector<double>& energies = std::vector<double>(),
            double partition = 0.0, const std::string& library = "",
            const std::string& chain = "", int residue = 0);

    //! Every atom of every rotamer, flat `n_rotamers * n_atoms * 3`.
    void get_atoms(double** out_view, int* n_out_view) const;
    //! The interaction energy of each rotamer, kT.
    void get_energies(double** out_view, int* n_out_view) const;
    //! The chromophore centres, flat `n_rotamers * 3` -- `points` without the
    //! weight column.
    void get_centres(double** out_view, int* n_out_view) const;
    //! The normalised weights, one per rotamer -- `points`' fourth column.
    void get_weights(double** out_view, int* n_out_view) const;

    std::vector<std::string> get_atom_names() const { return atom_names_; }
    std::vector<std::string> get_resnames() const { return resnames_; }
    double get_partition() const { return partition_; }
    std::string get_library() const { return library_; }
    std::string get_chain() const { return chain_; }
    int get_residue() const { return residue_; }
    int get_n_rotamers() const { return get_n_points(); }
    //! Atoms per rotamer; zero when the ensemble carries no atoms.
    int get_n_atoms() const;

    //! How many conformers this site's ensemble is really carrying.
    /*!
        Kish's effective sample size, \f$(\sum w)^2 / \sum w^2\f$. The library
        is one fixed sample of the *free* dye; a site reweights that sample and
        never resamples it, so what the site gets is only as good as the
        overlap between the sample and the conformations the site allows. The
        rotamer count does not show this -- a site can keep twenty rotamers and
        put three quarters of the weight on one of them.

        Measured on the bundled pins (2026-08-24): at Hsp90 residue 637 the
        cutoff-30 library carries **1.7**, and the answer it gives differs from
        the cutoff-10 library's by 0.068 in E and 0.32 in
        \f$\langle\kappa^2\rangle\f$. At the pp11 sites, which are open, the
        same comparison moves E by 0.008. Below roughly ten, treat what this
        ensemble says as unsupported and load a deeper cutoff.
    */
    double get_effective_sample_size() const;

    //! \f$R_{ij}\f$, \f$\kappa^2_{ij}\f$ and \f$w_i w_j\f$ against \p other.
    /*! \param[in] other any states: an AV cloud has no dipoles, and the
                   isotropic \f$\kappa^2 = 2/3\f$ is used for it
        \param[in] use_dipoles false ignores both ensembles' dipoles and puts
                   \f$\kappa^2 = 2/3\f$ on every pair -- the fps.json and
                   #IMP::bff::ProbeNetworkRestraint convention, which is what
                   makes a rotamer distance comparable with an AV's */
    FRETPairGeometry pair_geometry(const States& other,
                                   bool use_dipoles = true) const;

    //! The pair's FRET rate distribution and its three averages.
    /*! \param[in] other the acceptor states
        \param[in] forster_radius \f$R_0\f$ at \f$\kappa^2 = 2/3\f$, A
        \param[in] tau0 donor lifetime, ns; negative leaves `k_fret` empty */
    FRETPairEfficiencies pair_distribution(const States& other,
                                           double forster_radius,
                                           double tau0 = -1.0) const;

    //! The same, with \f$R_0\f$ computed from the two dyes' spectra.
    /*! FRETpredict's convention: the pair's own \f$\langle\kappa^2\rangle\f$
        enters \f$R_0\f$, which is then applied with the isotropic formula.
        \param[in] donor,acceptor dye names in the bundled spectra library
        \throw ValueException when either dye is unknown */
    FRETPairEfficiencies pair_distribution_from_probes(
            const States& other, const std::string& donor,
            const std::string& acceptor, double tau0 = -1.0) const;

    //! Place \p library on `(chain, residue)` of \p frame and screen it.
    /*!
        The scoring is FRETpredict's (Lennard-Jones or Gaussian, optionally
        Debye-Hueckel): the labelled residue's own atoms and, by default,
        hydrogens are not obstacles.

        Which atom is the chromophore centre, and which two span the
        transition dipole, come from the library's registry metadata (`r` and
        `mu`). A library with neither falls back to its first atom for the
        centre and to the vector from its first atom to its second for the
        dipole -- a direction, not the dye's, so an off-registry library
        should carry selectors.

        \param[in] frame the protein
        \param[in] chain,residue the labelled site; an empty chain matches any
        \param[in] library the rotamer library, with its metadata
        \param[in] options the screening parameters
        \param[in] position_name the ensemble's name; empty makes it
                   `<chain><residue>`
        \throw ValueException when the site has no CA, N or C, or when a
               selector matches no atom
    */
    static ProbeRotamerEnsemble from_frame(
            const ProteinFrame& frame, const std::string& chain, int residue,
            const ProbeRotamerLibrary& library,
            const ProbeRotamerSiteOptions& options = ProbeRotamerSiteOptions(),
            const std::string& position_name = "");

    //! The same, reading the structure and the library from their files.
    /*! A separate name and not an overload: SWIG turns keyword arguments off
        for an overloaded function, and `from_site(pdb, "A", 452, dye,
        position_name="A452")` is how every caller of this writes it. */
    /*!
        \param[in] structure a PDB (or multi-MODEL PDB) path
        \param[in] chain,residue the labelled site
        \param[in] library a registry name, a locator or a path
        \param[in] options the screening parameters
        \param[in] position_name the ensemble's name
        \param[in] frame_index which model of a multi-MODEL PDB
        \throw ValueException when \p frame_index is past the last model
    */
    static ProbeRotamerEnsemble from_site(
            const std::string& structure, const std::string& chain,
            int residue, const std::string& library,
            const ProbeRotamerSiteOptions& options = ProbeRotamerSiteOptions(),
            const std::string& position_name = "", int frame_index = 0);

    IMP_SHOWABLE_INLINE(ProbeRotamerEnsemble,
                        out << "RotamerEnsemble(" << get_n_rotamers()
                            << " rotamers, \"" << library_ << "\", Z="
                            << partition_ << ")");
};
IMP_VALUES(ProbeRotamerEnsemble, ProbeRotamerEnsembles);

// --------------------------------------------------------------------------
// Weighted averaging over an ensemble
// --------------------------------------------------------------------------

//! Per-frame weights from a pair of partition functions.
/** Named `rotamer_frame_weights` for its first caller; nothing in it is
    particular to a rotamer. Two partition functions per frame -- from
    rotamer libraries, from an umbrella sampling, from anything that can say
    how much of its ensemble each frame leaves -- weight the frames by their
    product. */
/** `z` is `(n_frames, 2)` flattened -- the donor's and acceptor's partition
    function in each frame. The weight of a frame is the product, normalised
    over frames. A single `(2,)` row means one frame and weight 1.

    When every product is zero the weights are uniform rather than undefined:
    a frame in which neither dye has any accessible conformer says nothing
    about the others. */
IMPBFFEXPORT void frame_weights_from_partitions(
        double* z_values, int n_frames, int n_pair,
        double** out_view, int* n_out_view);

//! Weighted mean, standard deviation and standard error of `values`.
/** Non-finite values and their weights are dropped first -- a frame where the
    dye could not be placed contributes nothing rather than poisoning the mean.
    The weights are renormalised **after** that removal, so the surviving
    frames still sum to one.

    Returns `(mean, sd, se)`; all three are NaN when nothing is finite. The
    standard error divides by the *count* of surviving frames, not by the
    effective count. */
IMPBFFEXPORT std::vector<double> weighted_average_sd_se(
        double* values, int n_values,
        double* weights, int n_weights);

//! The effective number of contributing frames, `exp` of the entropy.
/** Zero weights are dropped, then `exp(-sum w log(w / uniform))`. Equals the
    frame count for uniform weights and falls toward 1 as the weight
    concentrates on one frame. Zero when nothing contributes. */
IMPBFFEXPORT double effective_frame_fraction(double* weights, int n_weights);

IMPBFF_END_NAMESPACE


#endif
