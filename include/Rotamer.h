#ifndef IMPBFF_ROTAMER_H
#define IMPBFF_ROTAMER_H

/**
 *  \file IMP/bff/Rotamer.h
 *  \brief A rotamer library placed on a residue, screened, and what two such
 *         libraries predict.
 *
 * Five former headers, in the order a rotamer analysis runs:
 *
 * 1. **Energy** (formerly `RotamerEnergy.h`) -- the interaction energy of
 *    each rotamer with its protein.
 * 2. **Site** (formerly `RotamerSite.h`) -- placing a library on a residue:
 *    the backbone frame, the atom selectors, the library registry.
 * 3. **Ensemble** (formerly `RotamerEnsemble.h`) -- a library placed and
 *    screened at one labelling site, as a #IMP::bff::States.
 * 4. **fps** (formerly `RotamerFps.h`) -- fps.json and ensembles: a labelled
 *    pair in, the distances it predicts out.
 * 5. **FRET** (formerly `RotamerFret.h`) -- FRET over a trajectory from two
 *    screened libraries.
 *
 * `RotamerLibrary.h` stays separate: it is the library file itself (the
 * `.drot.pto` reader/writer), which everything here consumes.
 *
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

// -------- from RotamerEnergy.h --------
/**
 *  (formerly IMP/bff/RotamerEnergy.h, now a section of this file)
 *  \brief The interaction energy of each rotamer with its protein.
 *
 * A rotamer library is a list of candidate dye conformers; scoring turns it
 * into a weighted ensemble by asking how hard each one collides with the
 * structure it is attached to. That is an all-pairs problem — every dye atom
 * against every protein atom, for every conformer — and it is the inner loop
 * of `RotamerFRET`.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#include <IMP/bff/bff_config.h>

#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! Which functional form the steric term takes.
enum RotamerPotential {
    ROTAMER_POTENTIAL_LJ = 0,     //!< 12-6 Lennard-Jones
    ROTAMER_POTENTIAL_GAUSS = 1   //!< a soft Gaussian well, for a coarse dye
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
    \param[in] potential a #RotamerPotential
    \param[in] steric_cutoff pairs beyond this contribute nothing, Angstrom
    \param[in] coulomb_cutoff the same for the screened Coulomb term
    \param[in] debye_length screening length of the Debye-Huckel term, Angstrom
    \param[in] coulomb_prefactor the `q_i q_j / r` coefficient
    \return two values per rotamer: steric energy, then electrostatic energy.
            Separate rather than summed because the Boltzmann factor treats them
            differently -- the steric term is divided by RT and the screened
            Coulomb term is already in units of it.
*/
IMPBFFEXPORT std::vector<double> rotamer_interaction_energies(
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

 //IMPBFF_ROTAMERENERGY_H

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
#include <IMP/bff/RotamerLibrary.h>

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

// `selector_resnames(selectors)` lives in Scoring.h beside the rest of the
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
IMPBFFEXPORT std::string rotamer_library_registry();

//! One registry entry of `libraries.json`, plus the resolved spelling.
/*! The returned JSON object is the registry's own entry with `name` (the
    key), `library_name` (the name as asked) and `cutoff` added. The registry
    is read once and cached. An unknown library raises. */
IMPBFFEXPORT std::string rotamer_library_metadata(const std::string& name);

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
IMPBFFEXPORT std::string rotamer_library_metadata_for_path(
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
IMPBFFEXPORT std::string resolve_rotamer_library_path(
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
IMPBFFEXPORT RotamerLibrary load_rotamer_library(const std::string& name,
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
IMPBFFEXPORT std::vector<std::string> infer_rotamer_resnames(
        const std::vector<std::string>& atom_names,
        const std::string& metadata_json);

IMPBFF_END_NAMESPACE

 //IMPBFF_ROTAMERSITE_H

// -------- from RotamerEnsemble.h --------
/**
 *  (formerly IMP/bff/RotamerEnsemble.h, now a section of this file)
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
 * A **sibling** of #IMP::bff::AccessibleVolume, not a subclass: both are
 * #IMP::bff::States, so every distance and FRET helper takes either, while an
 * ensemble is spared the grid an AV has and it does not.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#include <IMP/bff/States.h>
#include <IMP/bff/FRET.h>
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
struct IMPBFFEXPORT RotamerSiteOptions {
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

    RotamerSiteOptions(double temperature = 298.15,
                       bool electrostatic = false,
                       const std::string& potential = "lj",
                       bool ignore_h = true, double sigma_scaling = 0.5,
                       double epsilon_scaling = 1.0)
        : temperature(temperature), electrostatic(electrostatic),
          potential(potential), ignore_h(ignore_h),
          sigma_scaling(sigma_scaling), epsilon_scaling(epsilon_scaling) {}

    IMP_SHOWABLE_INLINE(RotamerSiteOptions,
                        out << "RotamerSiteOptions(" << potential << ", T="
                            << temperature << ")");
};
IMP_VALUES(RotamerSiteOptions, RotamerSiteOptionsList);

//! A rotamer library placed and screened at one labelling site.
/*! From #IMP::bff::States: `points` (the chromophore centre and the weight of
    each rotamer), `attachment_point` (the site's CA), `orientations` (the
    per-rotamer transition dipoles), `position_name` and `params`. Added here:
    every atom of every rotamer, the atom and residue names they are in the
    order of, the interaction energies the weights came from, and the
    partition function. */
class IMPBFFEXPORT RotamerEnsemble : public States {
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
    RotamerEnsemble(
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
    static RotamerEnsemble from_frame(
            const ProteinFrame& frame, const std::string& chain, int residue,
            const RotamerLibrary& library,
            const RotamerSiteOptions& options = RotamerSiteOptions(),
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
    static RotamerEnsemble from_site(
            const std::string& structure, const std::string& chain,
            int residue, const std::string& library,
            const RotamerSiteOptions& options = RotamerSiteOptions(),
            const std::string& position_name = "", int frame_index = 0);

    IMP_SHOWABLE_INLINE(RotamerEnsemble,
                        out << "RotamerEnsemble(" << get_n_rotamers()
                            << " rotamers, \"" << library_ << "\", Z="
                            << partition_ << ")");
};
IMP_VALUES(RotamerEnsemble, RotamerEnsembles);

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

 //IMPBFF_ROTAMERENSEMBLE_H

// -------- from RotamerFps.h --------
/**
 *  (formerly IMP/bff/RotamerFps.h, now a section of this file)
 *  \brief fps.json and rotamer ensembles: reading a labelled pair out of a
 *         file, and writing predicted distances back into one.
 *
 * The fps.json format names the same thing several ways -- a position's chain
 * is `chain_identifier`, `chain` or `segid`, a distance's donor is
 * `position1_name`, `donor_position`, `donor_position_name` or
 * `dye1_position` -- because it has been written by several programs. The
 * aliases are the format's, so reading them is this module's job and not
 * every caller's.
 *
 * Entries cross the language boundary as JSON text, which is what the rest of
 * the fps layer does (#IMP::bff::read_fps_json, #IMP::bff::write_fps_json):
 * an entry carries whatever keys its writer put there, and a typed struct
 * would either lose them or have to grow a field per program. What *is* typed
 * is the part this module reasons about -- which chain, which residue, which
 * library, which dye.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */

#include <utility>

IMPBFF_BEGIN_NAMESPACE

//! A labelling position of an fps.json file.
struct IMPBFFEXPORT RotamerPosition {
    std::string name;
    //! Empty when the file names none, which matches any chain.
    std::string chain;
    int residue;
    //! The anchor atom; `CA` unless the file says otherwise.
    std::string atom_name;
    //! The dye, and the rotamer library; either may be empty.
    std::string dye, library;
    //! `donor`, `acceptor`, ... as the file spells it; may be empty.
    std::string role;

    RotamerPosition() : residue(0), atom_name("CA") {}

    IMP_SHOWABLE_INLINE(RotamerPosition,
                        out << "RotamerPosition(" << name << ", " << chain
                            << residue << ")");
};
IMP_VALUES(RotamerPosition, RotamerPositions);

//! A FRET distance between two labelling positions of an fps.json file.
struct IMPBFFEXPORT RotamerDistance {
    std::string name;
    std::string donor_position, acceptor_position;
    //! The dyes and their libraries, taken from the distance entry when it
    //! names them and from the two positions when it does not.
    std::string donor, acceptor, libname_1, libname_2;

    RotamerDistance() {}

    IMP_SHOWABLE_INLINE(RotamerDistance,
                        out << "RotamerDistance(" << name << ": "
                            << donor_position << " -> " << acceptor_position
                            << ")");
};
IMP_VALUES(RotamerDistance, RotamerDistances);

//! One distance of an fps.json file, with both its positions and the document.
struct IMPBFFEXPORT RotamerFpsSelection {
    RotamerPosition donor, acceptor;
    RotamerDistance distance;
    //! The whole document, as JSON text: every position, every distance, and
    //! the score sets and extra sections merged.
    std::string positions, distances, extra;

    RotamerFpsSelection() : positions("{}"), distances("{}"), extra("{}") {}

    IMP_SHOWABLE_INLINE(RotamerFpsSelection,
                        out << "RotamerFpsSelection(" << distance.name << ")");
};
IMP_VALUES(RotamerFpsSelection, RotamerFpsSelections);

//! Parse one fps.json position entry.
/*! \param[in] name the entry's key
    \param[in] payload_json the entry */
IMPBFFEXPORT RotamerPosition rotamer_position_from_payload(
        const std::string& name, const std::string& payload_json);

//! Parse one fps.json distance entry, filling from its positions.
/*! A distance that names no dye or library takes the ones its positions name,
    which is how an fps file written position-first still describes a pair.

    \param[in] name the entry's key
    \param[in] payload_json the entry
    \param[in] positions_json every position of the document */
IMPBFFEXPORT RotamerDistance rotamer_distance_from_payload(
        const std::string& name, const std::string& payload_json,
        const std::string& positions_json);

//! Read one FRET distance and its two positions out of an fps.json file.
/*!
    \param[in] path the fps.json file
    \param[in] distance_name which distance; empty takes the first
    \throw ValueException when the file has no distances, when \p
           distance_name is not one of them, or when a distance names a
           position the file does not have
*/
IMPBFFEXPORT RotamerFpsSelection read_rotamer_fps(
        const std::string& path, const std::string& distance_name = "");

//! The fps.json entry of a rotamer position (`simulation_type` `R1`).
/*!
    \param[in] chain,residue,library what the position is
    \param[in] atom_name the anchor atom
    \param[in] dye the dye's name; empty leaves the key out
    \param[in] temperature K; NaN leaves the key out
    \param[in] electrostatic 1 true, 0 false, -1 leaves the key out
    \param[in] potential `lj` or `gauss`; empty leaves the key out
    \return the entry, as JSON text

    A key that is left out is not the same as a key written empty: the reader
    falls back on its own default for the first and takes the second at its
    word.
*/
IMPBFFEXPORT std::string rotamer_position_payload(
        const std::string& chain, int residue, const std::string& library,
        const std::string& atom_name = "CA", const std::string& dye = "",
        double temperature = std::numeric_limits<double>::quiet_NaN(),
        int electrostatic = -1, const std::string& potential = "");

//! The `R1` entries describing existing ensembles, keyed as they are.
/*! Each ensemble's own `params` -- the record of how it was screened -- fills
    the temperature, the potential and the electrostatic flag, so a file
    written from ensembles says how to reproduce them. */
IMPBFFEXPORT std::string rotamer_positions_payload(
        const std::map<std::string, RotamerEnsemble>& ensembles,
        const std::string& atom_name = "CA");

//! Predicted fps.json distance entries between rotamer ensembles.
/*!
    \param[in] ensembles the screened ensembles, keyed by position name
    \param[in] pairs which of them to measure, donor first
    \param[in] forster_radius \f$R_0\f$ at \f$\kappa^2 = 2/3\f$, A
    \param[in] distance_type `RDAMean` (\f$\langle R_{DA}\rangle\f$),
               `RDAMeanE` (the FRET-averaged one) or `Rmp` (between mean
               positions) -- from the full pair matrix, with no sampling
    \param[in] error the error bar to write; negative uses \p error_fraction
    \param[in] error_fraction of the distance, when no \p error is given
    \param[in] kappa2 `isotropic` puts \f$\kappa^2 = 2/3\f$ on every pair,
               which is the fps.json and #IMP::bff::ProbeNetworkRestraint
               convention and is what makes the number comparable with an
               AV's; `dipoles` uses the ensembles' own per-pair
               \f$\kappa^2\f$, which is the orientation-resolved answer
    \return the entries, keyed `<donor>_<acceptor>`, as JSON text
    \throw ValueException for an unknown \p distance_type or \p kappa2, or a
           pair naming an ensemble that is not there
*/
IMPBFFEXPORT std::string distances_from_ensembles(
        const std::map<std::string, RotamerEnsemble>& ensembles,
        const std::vector<std::pair<std::string, std::string> >& pairs,
        double forster_radius,
        const std::string& distance_type = "RDAMeanE", double error = -1.0,
        double error_fraction = 0.05,
        const std::string& kappa2 = "isotropic");

//! One ensemble per fps.json position that names a rotamer library.
/*!
    A position's `rotamer_library` (or `library_name` / `libname` / `library`)
    field selects the library; \p library_map, keyed by position name,
    overrides it or supplies one for a position that names none -- which is
    how an AV-only fps file is screened as rotamers. Positions with neither
    are skipped rather than guessed at.

    The structure is read once and each library at most once, however many
    positions share them.

    \param[in] fps_json the fps.json path
    \param[in] structure the PDB the positions refer to
    \param[in] library_map position name -> library name
    \param[in] options the screening parameters, for every position alike
    \param[in] frame_index which model of a multi-MODEL PDB
*/
IMPBFFEXPORT std::map<std::string, RotamerEnsemble> rotamer_ensembles_from_fps(
        const std::string& fps_json, const std::string& structure,
        const std::map<std::string, std::string>& library_map =
                std::map<std::string, std::string>(),
        const RotamerSiteOptions& options = RotamerSiteOptions(),
        int frame_index = 0);

//! Write (or merge into) an fps.json with rotamer (`R1`) positions.
/*!
    \param[in] path the file to write
    \param[in] positions_json,distances_json the entries, as JSON objects
    \param[in] merge_into an existing fps.json whose positions and distances
               are kept; entries of the same name are replaced. Its score sets
               and extra sections are carried over untouched.
    \param[in] validate check against the fps.json schema before writing
    \throw ValueException when validation fails
*/
IMPBFFEXPORT void write_rotamer_fps(const std::string& path,
                                    const std::string& positions_json,
                                    const std::string& distances_json = "{}",
                                    const std::string& merge_into = "",
                                    bool validate = true);

IMPBFF_END_NAMESPACE

 //IMPBFF_ROTAMERFPS_H

// -------- from RotamerFret.h --------
/**
 *  (formerly IMP/bff/RotamerFret.h, now a section of this file)
 *  \brief FRET over a trajectory from two screened rotamer libraries.
 *
 * The driver FRETpredict is: place a rotamer library at each of two labelled
 * residues in every frame, screen both against the protein, and report the
 * pair's \f$\langle\kappa^2\rangle\f$ and its efficiency in the three
 * averaging limits, frame by frame.
 *
 * **The parameter names here are FRETpredict's**, deliberately -- `fixed_R0`,
 * `ign_H`, `libname_1`, `r0lib` -- because the parity harness
 * (`test/cgprobe/rotamer/expensive_test_fretpredict_parity.py`) hands the *same*
 * keyword dictionary to this class and to FRETpredict and compares the files
 * they write. A rename here would be a rename of the experiment.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */


IMPBFF_BEGIN_NAMESPACE

//! What one frame of a trajectory says about a labelled pair.
struct IMPBFFEXPORT FRETFrameResult {
    //! The two partition functions, donor first: how much of each library the
    //! site leaves accessible.
    double z_donor, z_acceptor;
    double kappa2_avg;
    //! The three averaging limits (see #IMP::bff::FRETPairEfficiencies).
    double static_efficiency, dynamic1, dynamic2;

    FRETFrameResult()
        : z_donor(0), z_acceptor(0), kappa2_avg(2.0 / 3.0),
          static_efficiency(0), dynamic1(0), dynamic2(0) {}

    IMP_SHOWABLE_INLINE(FRETFrameResult,
                        out << "FRETFrameResult(E=" << static_efficiency
                            << ", k2=" << kappa2_avg << ")");
};
IMP_VALUES(FRETFrameResult, FRETFrameResults);

//! FRET over a trajectory from two screened rotamer libraries.
class IMPBFFEXPORT RotamerFRET {
    std::vector<ProteinFrame> frames_;
    std::vector<int> residues_;
    std::vector<std::string> chains_;
    std::string donor_, acceptor_, libname_1_, libname_2_, r0lib_;
    std::string output_prefix_;
    RotamerSiteOptions site_;
    RotamerLibrary lib_1_, lib_2_;
    bool fixed_r0_;
    double r0_, z_cutoff_;
    std::vector<double> user_weights_, weights_;

    std::vector<double> z_values_, k2_values_, estatic_, edynamic1_,
            edynamic2_;

    void load_libraries();
    void write_summary(const std::string& prefix,
                       const std::vector<double>& k2,
                       const std::vector<double>& es,
                       const std::vector<double>& ed1,
                       const std::vector<double>& ed2,
                       const std::vector<double>& weights) const;

    //! The frames-in constructor #from_frames uses.
    RotamerFRET(const std::vector<ProteinFrame>& frames,
                const std::vector<int>& residues,
                const std::vector<std::string>& chains,
                const std::string& donor, const std::string& acceptor,
                const std::string& libname_1, const std::string& libname_2,
                const RotamerSiteOptions& site, bool fixed_R0, double r0,
                const std::string& r0lib, double z_cutoff,
                const std::string& output_prefix,
                const std::vector<double>& user_weights);

public:
    //! \param[in] protein a PDB (or multi-MODEL PDB) path
    /*! \param[in] residues the two labelled residue numbers
        \param[in] chains their chains; empty matches any
        \param[in] donor,acceptor the dye names the spectra are looked up by
        \param[in] libname_1,libname_2 the two rotamer libraries
        \param[in] temperature K
        \param[in] electrostatic add the Debye-Hueckel term
        \param[in] potential `lj` or `gauss`
        \param[in] ign_H mask hydrogens out of the steric term
        \param[in] sigma_scaling,epsilon_scaling dye-probe LJ scaling
        \param[in] fixed_R0 keep \p r0 rather than computing one per frame
        \param[in] r0 the Förster radius, **Angström** (it was nanometres
                   until 2026-08-25, when R0 became one unit across the
                   package: a fixed 5.5 nm is 55 here)
        \param[in] r0lib a dye-library CIF to take the spectra from
        \param[in] z_cutoff frames where either partition function falls below
                   this report NaN rather than a number the screening does not
                   support
        \param[in] output_prefix what #save writes its files under
        \param[in] max_frames stop after this many; negative reads all
        \param[in] user_weights one per frame, for the summary table
        \throw ValueException unless exactly two residues are given */
    RotamerFRET(const std::string& protein,
                const std::vector<int>& residues = std::vector<int>(),
                const std::vector<std::string>& chains =
                        std::vector<std::string>(),
                const std::string& donor = "AlexaFluor 488",
                const std::string& acceptor = "AlexaFluor 594",
                const std::string& libname_1 = "AlexaFluor 488 C1R cutoff30",
                const std::string& libname_2 = "AlexaFluor 594 C1R cutoff30",
                double temperature = 300.0, bool electrostatic = false,
                const std::string& potential = "lj", bool ign_H = true,
                double sigma_scaling = 0.5, double epsilon_scaling = 1.0,
                bool fixed_R0 = false, double r0 = 54.0,
                const std::string& r0lib = "", double z_cutoff = 0.05,
                const std::string& output_prefix = "res",
                int max_frames = -1,
                const std::vector<double>& user_weights =
                        std::vector<double>());

    //! The same, over frames a caller already has.
    /*! The door for the formats this module does not read itself: an RMF
        trajectory becomes `ProteinFrame` values (`protein_frames_from_rmf`)
        and comes in here.

        A named factory and not a second constructor, because SWIG turns
        keyword arguments off for anything overloaded and the constructor's
        twenty parameters are keyword arguments in every caller. */
    static RotamerFRET from_frames(
            const std::vector<ProteinFrame>& frames,
                const std::vector<int>& residues,
                const std::vector<std::string>& chains =
                        std::vector<std::string>(),
                const std::string& donor = "AlexaFluor 488",
                const std::string& acceptor = "AlexaFluor 594",
                const std::string& libname_1 = "AlexaFluor 488 C1R cutoff30",
                const std::string& libname_2 = "AlexaFluor 594 C1R cutoff30",
                double temperature = 300.0, bool electrostatic = false,
                const std::string& potential = "lj", bool ign_H = true,
                double sigma_scaling = 0.5, double epsilon_scaling = 1.0,
                bool fixed_R0 = false, double r0 = 54.0,
                const std::string& r0lib = "", double z_cutoff = 0.05,
            const std::string& output_prefix = "res",
            const std::vector<double>& user_weights = std::vector<double>());

    //! What one frame says: both partition functions, ⟨κ²⟩ and the three
    //! efficiencies.
    /*! \note With `fixed_R0` false this **sets** the object's `r0` to the one
        this frame's \f$\langle\kappa^2\rangle\f$ gives, which is what
        FRETpredict does and why `r0` is readable afterwards. */
    FRETFrameResult frame_fret(const ProteinFrame& frame);

    //! Every frame, into the value arrays.
    /*! A frame where either partition function is at or below `z_cutoff`
        keeps its partition functions and reports NaN for the rest: the
        screening left the dye nowhere to be, and a number computed from that
        would be a number about the fallback, not about the structure. */
    void trajectory_analysis();

    //! Write the FRETpredict-style output files.
    /*!
        `<prefix>-Z-<r1>-<r2>.dat` and one each for `w_s`, `k2`, `Es`, `Ed1`,
        `Ed2`, plus a labelled `-data-` table of (average, SD, SE) per
        quantity. A text table and not a pickled DataFrame: a `.pkl` cannot be
        read without the library that wrote it.

        \param[in] output_prefix overrides the constructor's
        \throw ValueException when `user_weights` is not one per frame
    */
    void save(const std::string& output_prefix = "");

    //! Recompute the summary table from the saved files with other weights.
    /*!
        \param[in] boltzmann_weights weight each frame by the product of its
                   two partition functions (read back from the `-Z-` file)
        \param[in] user_weights one per frame, multiplied into the above
        \param[in] output_prefix where the table goes; the inputs are always
                   read from the constructor's prefix
    */
    void reweight(bool boltzmann_weights = false,
                  const std::vector<double>& user_weights =
                          std::vector<double>(),
                  const std::string& output_prefix = "");

    //! #trajectory_analysis followed by #save.
    void run();

    //! `(n_frames, 2)` partition functions, donor first.
    void get_z_values(double** out_view, int* n_out_view) const;
    //! One per frame; NaN where the screening left nothing.
    void get_k2_values(double** out_view, int* n_out_view) const;
    void get_estatic_values(double** out_view, int* n_out_view) const;
    void get_edynamic1_values(double** out_view, int* n_out_view) const;
    void get_edynamic2_values(double** out_view, int* n_out_view) const;

    //! The Förster radius in Å -- the fixed one, or the last frame's.
    double get_r0() const { return r0_; }
    //! The two labelled residues, and the chains they are on.
    std::vector<int> get_residues() const { return residues_; }
    std::vector<std::string> get_chains() const { return chains_; }
    //! The two dyes, and the two libraries screened at their sites.
    std::string get_donor() const { return donor_; }
    std::string get_acceptor() const { return acceptor_; }
    std::string get_libname_1() const { return libname_1_; }
    std::string get_libname_2() const { return libname_2_; }
    int get_n_frames() const { return static_cast<int>(frames_.size()); }
    std::string get_output_prefix() const { return output_prefix_; }

    IMP_SHOWABLE_INLINE(RotamerFRET,
                        out << "RotamerFRET(" << frames_.size() << " frames, "
                            << donor_ << " -> " << acceptor_ << ")");
};

//! Build a #IMP::bff::RotamerFRET from an fps.json distance.
/*! The file names the two positions, their dyes and their libraries; anything
    it does not name is left at the constructor's default.

    \param[in] fps_path the fps.json file
    \param[in] protein the structure its positions refer to
    \param[in] distance_name which distance; empty takes the first
    \param[in] output_prefix what #IMP::bff::RotamerFRET::save writes under
    \param[in] fixed_R0,r0 as in the constructor (\p r0 in Å)
    \param[in] temperature,electrostatic as in the constructor
*/
IMPBFFEXPORT RotamerFRET rotamer_fret_from_fps(
        const std::string& fps_path, const std::string& protein,
        const std::string& distance_name = "",
        const std::string& output_prefix = "res", bool fixed_R0 = false,
        double r0 = 54.0, double temperature = 300.0,
        bool electrostatic = false);

IMPBFF_END_NAMESPACE

 //IMPBFF_ROTAMERFRET_H

#endif  // IMPBFF_ROTAMER_H
