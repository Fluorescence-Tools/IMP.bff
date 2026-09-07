/**
 *  \file IMP/bff/RotamerEnsemble.h
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
#ifndef IMPBFF_ROTAMERENSEMBLE_H
#define IMPBFF_ROTAMERENSEMBLE_H

#include <IMP/bff/bff_config.h>
#include <IMP/bff/AVModel.h>
#include <IMP/bff/FRETPair.h>
#include <IMP/bff/HierarchyFrame.h>
#include <IMP/bff/RotamerLibrary.h>

#include <map>
#include <string>
#include <vector>

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

#endif //IMPBFF_ROTAMERENSEMBLE_H
