/**
 *  \file IMP/bff/Scoring.h
 *  \brief Stage-2 scoring: the CHARMM36 table, LJ, Boltzmann, and the
 *         rotamer-protein score that turns conformers into a weighted ensemble.
 *
 * The inner kernels (all-pairs steric/electrostatic energy, the pair energy
 * matrix, the per-frame LJ over an explicit pair list) are in
 * \ref Rotamer.h (its RotamerEnergy section). What lives here is the orchestration around them: the
 * CHARMM36 parameter table and the Lorentz-Berthelot combining rules, the
 * Boltzmann weight, the axis-aligned bounding-box pre-filter, and the
 * end-to-end rotamer score that builds masks, assembles parameters, calls the
 * kernel, and returns normalised weights.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_SCORING_H
#define IMPBFF_SCORING_H

#include <IMP/bff/bff_config.h>
#include <IMP/bff/ProbeLibrary.h>

#include <IMP/bff/Base.h>

#include <IMP/Model.h>
#include <IMP/Restraint.h>
#include <IMP/core/Cosine.h>

#include <vector>
#include <string>
#include <map>
#include <utility>

IMPBFF_BEGIN_NAMESPACE

// --------------------------------------------------------------------------
// CHARMM36 Lennard-Jones parameters
// --------------------------------------------------------------------------

//! The element type of a PDB-style atom name: S, N, O, H, or C.
/*! Uppercased, and the first of the five a dye carries that the name starts
    with -- `"SD"` is S, `"NZ"` is N, anything else is carbon. */
IMPBFFEXPORT std::string atom_type(const std::string& atom_name);

//! The CHARMM36 LJ parameters for the five elements bff dyes carry.
/*!
    \param[out] out_view,n_out_view two values: rmin_half (A), then epsilon
               (kcal/mol). A view rather than a pair: it marshals as a numpy
               array, which unpacks and indexes like a tuple.
    Unknown elements fall back to carbon, as \ref atom_type arranges for.
*/
IMPBFFEXPORT void charmm36_lj(const std::string& element,
                              double** out_view, int* n_out_view);

//! Lorentz-Berthelot cross parameters for two elements.
/*!
    \c rmin = rmin_half_i + rmin_half_j (Angstrom),
    \c eps = sqrt(eps_i * eps_j) (kcal/mol, positive well depth).
    \param[out] out_view,n_out_view two values: rmin, then eps.
*/
IMPBFFEXPORT void lj_cross(const std::string& elem_i,
                           const std::string& elem_j,
                           double** out_view, int* n_out_view);

//! Per-element LJ parameter arrays.
struct LJArrays {
    std::vector<double> rmin_half;  //!< one per element, A
    std::vector<double> epsilon;    //!< one per element, kcal/mol
};

//! (rmin_half, epsilon) for a sequence of element symbols.
/*! Unknown elements fall back to carbon. */
IMPBFFEXPORT LJArrays lj_parameter_arrays(
        const std::vector<std::string>& elements);

//! The same two arrays, scaled by the dye-probe factors.
IMPBFFEXPORT LJArrays scaled_parameters(
        const std::vector<std::string>& elements, double sigma_scaling,
        double epsilon_scaling);

//! 12-6 Lennard-Jones energy eps * ((rmin/r)^12 - 2 (rmin/r)^6).
/*!
    Vectorised over flat arrays. \c r is clamped to \c r_floor to keep the
    energy finite. \c repulsive_only zeroes the attractive tail (r >= rmin),
    \c cutoff zeroes pairs beyond that distance.
*/
IMPBFFEXPORT std::vector<double> lj_energy(const std::vector<double>& r,
                                            const std::vector<double>& rmin,
                                            const std::vector<double>& eps,
                                            bool repulsive_only = false,
                                            double cutoff = 0.0,
                                            double r_floor = 0.01);

//! Repulsive-only LJ energy of one pair (scalar).
IMPBFFEXPORT double lj_score(double r, double rmin, double eps);

//! The flat cross parameters of two element lists, `n_a * n_b` entries.
IMPBFFEXPORT LJArrays cross_lj_params(
        const std::vector<std::string>& elements_a,
        const std::vector<std::string>& elements_b);

//! Sum of the repulsive-only LJ energy over every (a, b) coordinate pair.
/*! The reference the AABB-prefiltered kernel is checked against. */
IMPBFFEXPORT double lj_pairs_sum(const std::vector<double>& coords_a,
                                 const std::vector<double>& coords_b,
                                 const std::vector<double>& rmin_ij,
                                 const std::vector<double>& eps_ij,
                                 double r_cutoff = 12.0);

//! The conformer-pair energy matrix between two sets, AABB-prefiltered.
/*!
    Builds the Lorentz-Berthelot cross parameters for
    `elements_a x elements_b` and calls the kernel in Rotamer.h.

    \param[in] coords_a flat, `n_a_conf * n_a_atoms * 3`
    \param[in] coords_b flat, `n_b_conf * n_b_atoms * 3`
    \param[in] elements_a,elements_b one per atom of each side
    \param[in] r_cutoff pairs beyond this contribute nothing, A
    \param[in] aabb_pad padding on each bounding box, A
    \param[out] out_view,n_out_view flat `n_a_conf * n_b_conf`, row-major
*/
IMPBFFEXPORT void pair_energy_matrix(
        const std::vector<double>& coords_a,
        const std::vector<double>& coords_b,
        const std::vector<std::string>& elements_a,
        const std::vector<std::string>& elements_b,
        int n_a_conf, int n_b_conf,
        double** out_view, int* n_out_view,
        double r_cutoff = 12.0, double aabb_pad = 3.5);

// --------------------------------------------------------------------------
// Boltzmann weights
// --------------------------------------------------------------------------

//! Normalised Boltzmann weights from energies (log-sum-exp stabilised).
/*!
    \param energies in kcal/mol (or consistent units with kT)
    \param temperature in Kelvin
    \return weights summing to 1.0
*/
IMPBFFEXPORT std::vector<double> boltzmann_weights(
        const std::vector<double>& energies, double temperature = 298.15);

//! Aggregate frame weights into cluster weights.
/*!
    \param assignments cluster index per frame (n_frames,)
    \param frame_weights normalised weight per frame (n_frames,)
    \param n_clusters total number of clusters
    \return normalised cluster weights (n_clusters,)
*/
IMPBFFEXPORT std::vector<double> cluster_weights(
        const std::vector<int>& assignments,
        const std::vector<double>& frame_weights,
        int n_clusters);

// --------------------------------------------------------------------------
// Axis-aligned bounding-box pre-filter
// --------------------------------------------------------------------------

//! Build padded AABBs for a batch of coordinate frames.
/*!
    \param coords flat, n_frames * n_atoms * 3
    \param n_frames, n_atoms, pad shapes and padding in Angstrom
    \return flat n_frames * 6: [xmin, ymin, zmin, xmax, ymax, zmax] per frame
*/
IMPBFFEXPORT std::vector<double> aabb_build(
        const std::vector<double>& coords,
        int n_frames, int n_atoms, double pad);

//! Build a single padded AABB for a static coordinate set.
/*!
    \param coords flat, n_atoms * 3
    \param n_atoms, pad
    \return 6 values: [xmin, ymin, zmin, xmax, ymax, zmax]
*/
IMPBFFEXPORT std::vector<double> aabb_build_single(
        const std::vector<double>& coords, int n_atoms, double pad);

//! Vectorised AABB intersection test of many frames vs one reference.
/*!
    \param rot_boxes n_frames * 6
    \param ref_box 6
    \return bool mask, n_frames: true if frame could clash
*/
IMPBFFEXPORT std::vector<int> aabb_intersects_reference(
        const std::vector<double>& rot_boxes,
        const std::vector<double>& ref_box);

//! What #BoundingBoxFilter::filter_frames returns.
struct AABBFilterResult {
    std::vector<double> coords;  //!< the surviving frames, flat
    std::vector<int> mask;       //!< one per input frame: 1 survived
};

//! One energy per frame beside the mask that decided which were scored.
struct EnergyMaskResult {
    std::vector<double> energies;  //!< one per frame; 0 where masked out
    std::vector<int> mask;         //!< one per frame: 1 was scored
};

//! An AABB clash pre-filter for rotamer coordinate batches.
class IMPBFFEXPORT BoundingBoxFilter {
    double pad_;

public:
    explicit BoundingBoxFilter(double pad = 3.5) : pad_(pad) {}

    //! Padded AABBs, one per frame: flat `n_frames * 6`.
    std::vector<double> build(const std::vector<double>& coords, int n_frames,
                              int n_atoms) const;
    //! One padded AABB for a static coordinate set: 6 values.
    std::vector<double> build_single(const std::vector<double>& coords,
                                     int n_atoms) const;
    //! Whether two 6-value boxes overlap.
    static bool intersects(const std::vector<double>& box_a,
                           const std::vector<double>& box_b);
    //! The overlap mask of many frames' boxes against one reference box.
    std::vector<int> intersects_reference(
            const std::vector<double>& rot_boxes,
            const std::vector<double>& ref_box) const;
    //! The frames whose box overlaps the reference's, and the survival mask.
    AABBFilterResult filter_frames(const std::vector<double>& coords,
                                   int n_frames, int n_atoms,
                                   const std::vector<double>& reference_coords,
                                   int n_ref_atoms) const;
};

// --------------------------------------------------------------------------
// The FRETpredict selector mini-language, and the masks built from it
// --------------------------------------------------------------------------

//! Does one `NAME [and resname RES]` selector match one atom?
/*! An empty \p resname matches only a selector without a resname clause. */
IMPBFFEXPORT bool selector_matches(const std::string& selector,
                                   const std::string& atom_name,
                                   const std::string& resname = "");

//! Every resname clause of a selector list, uppercased.
IMPBFFEXPORT std::vector<std::string> selector_resnames(
        const std::vector<std::string>& selectors);

//! The atom names a selector list selects.
/*!
    Without per-atom \p resnames, this is the atom-name clause of each selector
    uppercased. With them, it is every atom any selector matches. A list whose
    resname clauses name more than one distinct residue selects nothing -- the
    selection is ambiguous, and guessing one residue would mis-mask silently.
*/
IMPBFFEXPORT std::vector<std::string> selector_atom_names(
        const std::vector<std::string>& selectors,
        const std::vector<std::string>& atom_names = std::vector<std::string>(),
        const std::vector<std::string>& resnames = std::vector<std::string>());

//! One per atom: 1 where the atom is a hydrogen.
IMPBFFEXPORT std::vector<int> hydrogen_mask(
        const std::vector<std::string>& atom_names);

//! One per atom: 1 where the atom is masked out of a steric term.
/*!
    The labelling \p site_residue (with \p site_chain, matched case-blind;
    \p -1 / empty mean no site) is always masked, so the dye is not sterically
    blocked by the residue it is attached to. \p mask_backbone additionally
    masks the backbone names (CA, C, N, O) -- what the rotamer side wants and
    the protein side does not.
*/
IMPBFFEXPORT std::vector<int> site_mask(
        const std::vector<std::string>& atom_names,
        const std::vector<int>& residue_indices,
        int site_residue, const std::string& site_chain,
        const std::vector<std::string>& chain_ids,
        bool mask_backbone = true);

//! Point charges of a protein's atoms, by residue name.
/*! ARG CZ and LYS NZ +1, ASP CG and GLU CD -1, HIS ND1/NE2 +0.25; 0 elsewhere. */
IMPBFFEXPORT std::vector<double> protein_charge_mask(
        const std::vector<std::string>& atom_names,
        const std::vector<std::string>& resnames);

//! Point charges of a dye's atoms, from its positive/negative selectors.
/*! A selected negative atom is -1.0 and a selected positive one +0.5; a list
    whose resname clauses name more than one distinct residue selects nothing. */
IMPBFFEXPORT std::vector<double> rotamer_charge_mask(
        const std::vector<std::string>& atom_names,
        const std::vector<std::string>& positive,
        const std::vector<std::string>& negative,
        const std::vector<std::string>& resnames = std::vector<std::string>());

// --------------------------------------------------------------------------
// Rotamer-protein score (the end-to-end orchestration)
// --------------------------------------------------------------------------

//! The result of rotamer-protein scoring.
struct RotamerScoreResult {
    std::vector<double> weights;   //!< normalised Boltzmann weights
    double partition;              //!< partition function before normalisation
    std::vector<double> energies;   //!< raw potential + electrostatic energies

    const std::vector<double>& get_weights() const { return weights; }
    const std::vector<double>& get_energies() const { return energies; }
    double get_partition() const { return partition; }
    void show(std::ostream& out) const {
        out << "RotamerScoreResult(partition=" << partition
            << ", n_weights=" << weights.size() << ")";
    }
};

//! A list of rotamer score results (the SWIG plural type).
typedef std::vector<RotamerScoreResult> RotamerScoreResults;

//! Compute Boltzmann weights for a rotamer ensemble against a protein.
/*!
    The end-to-end scorer: builds the hydrogen/backbone/site masks from the
    atom and residue names, assembles the per-pair Lorentz-Berthelot
    parameters, calls the all-pairs steric/electrostatic kernel, and returns
    normalised Boltzmann weights multiplied by the library's own.

    \param[in] rotamer_coords flat, `n_rotamers * n_probe_atoms * 3`
    \param[in] protein_coords flat, `n_protein_atoms * 3`
    \param[in] protein_atom_names,protein_resnames one per protein atom
    \param[in] rotamer_atom_names one per dye atom
    \param[in] positive,negative the dye's charge selectors (its metadata's
               `positive`/`negative`), used only when \p electrostatic
    \param[in] rotamer_resnames one per dye atom, or empty
    \param[in] protein_residue_indices,protein_chain_ids one per protein atom,
               or empty -- what the \p site is resolved against
    \param[in] site_residue the labelled residue number, or -1 for none
    \param[in] site_chain the labelled residue's chain; empty matches any
    \param[in] rotamer_weights the library's own per-rotamer populations, or
               empty for uniform
    \param[in] temperature in K
    \param[in] ignore_h mask hydrogens out of the steric term
    \param[in] electrostatic add the Debye-Hueckel term from the charges
    \param[in] potential "lj" or "gauss"
    \param[in] sigma_scaling,epsilon_scaling dye-probe scaling of the LJ
               parameters
    \return RotamerScoreResult -- uniform weights with partition 1.0 when a
            mask leaves nothing to score
    \throw ValueException for any other \p potential
*/
IMPBFFEXPORT RotamerScoreResult get_rotamer_score(
        const std::vector<double>& rotamer_coords,
        const std::vector<double>& protein_coords,
        const std::vector<std::string>& protein_atom_names,
        const std::vector<std::string>& protein_resnames,
        const std::vector<std::string>& rotamer_atom_names,
        const std::vector<std::string>& positive = std::vector<std::string>(),
        const std::vector<std::string>& negative = std::vector<std::string>(),
        const std::vector<std::string>& rotamer_resnames =
                std::vector<std::string>(),
        const std::vector<int>& protein_residue_indices = std::vector<int>(),
        const std::vector<std::string>& protein_chain_ids =
                std::vector<std::string>(),
        int site_residue = -1, const std::string& site_chain = "",
        const std::vector<double>& rotamer_weights = std::vector<double>(),
        double temperature = 300.0, bool ignore_h = true,
        bool electrostatic = false, const std::string& potential = "lj",
        double sigma_scaling = 0.5, double epsilon_scaling = 1.0);

//! Iterative mean-field update of one dye's rotamer cluster weights.
/*! Eq. 39-42 of IMP's RotamerCalculator: the dye-protein interaction energies
    are computed once and the iteration is a matrix-vector product. \p K = 0
    returns the prior. */
IMPBFFEXPORT std::vector<double> rotamer_mean_field_weights(
        const std::vector<double>& rotamer_coords,
        const std::vector<double>& initial_weights,
        const std::vector<double>& protein_coords,
        const std::vector<std::string>& probe_elements,
        const std::vector<std::string>& protein_elements,
        double K = 1.0, int n_iter = 10, double aabb_pad = 3.5,
        double r_cutoff = 12.0);

//! Mean-field weights for several dyes at once, cross-interactions included.
/*!
    The single-dye #rotamer_mean_field_weights sees the protein only. With two
    or more labels on one structure each dye is also part of the others'
    environment, so the update carries a conformer-pair energy matrix per dye
    pair and iterates all of them together:
    \f$\log q_d \leftarrow \log q_d - K(E^{bb}_d + \sum_{e \neq d} E^{sc}_{de} q_e)\f$,
    renormalised by log-sum-exp each round.

    \param[in] rotamer_coords_list one flat `n_conf * n_atoms * 3` per dye
    \param[in] initial_weights_list one weight vector per dye, same order
    \param[in] protein_coords flat `n_atoms * 3`
    \param[in] probe_elements_list one element per atom, per dye
    \param[in] protein_elements one per protein atom
    \param[in] K inverse temperature of the update
    \param[in] n_iter fixed-point iterations
    \param[in] aabb_pad,r_cutoff as #pair_energy_matrix
    \return one normalised weight vector per dye
    \throw ValueException when the per-dye lists disagree in length
*/
IMPBFFEXPORT std::vector<std::vector<double> >
rotamer_mean_field_weights_multi_probe(
        const std::vector<std::vector<double> >& rotamer_coords_list,
        const std::vector<std::vector<double> >& initial_weights_list,
        const std::vector<double>& protein_coords,
        const std::vector<std::vector<std::string> >& probe_elements_list,
        const std::vector<std::string>& protein_elements,
        double K = 1.0, int n_iter = 10, double aabb_pad = 3.5,
        double r_cutoff = 12.0);

// --------------------------------------------------------------------------
// Walkers over a typed force-field system
// --------------------------------------------------------------------------

//! `{site_id: element}` from a system's sites.
IMPBFFEXPORT std::map<std::string, std::string> site_element_map(
        const ProbeForceFieldSystem& system);

//! One non-excluded site pair with its Lorentz-Berthelot parameters.
struct LJSitePair {
    std::string site_a, site_b;
    double rmin, eps;

    IMP_SHOWABLE_INLINE(LJSitePair, out << "LJSitePair(" << site_a << "-"
                                        << site_b << ")");
};

//! Every non-excluded site pair of a system, with LJ cross parameters.
/*! The exclusions are the system's own (bonds, angle and dihedral end pairs). */
IMPBFFEXPORT std::vector<LJSitePair> get_lj_pair_sites(
        const ProbeForceFieldSystem& system);

//! The \c IMP::core::Cosine a torsion type spells, converted from CHARMM.
/*!
    cgprobe stores torsions in the CHARMM convention
    \f$V = k(1 + \cos(n\phi - \delta))\f$, and \c IMP::core::Cosine scores
    \f$k(1 - \cos(n\phi - \delta'))\f$ -- so the phase needs
    \f$\delta' = \delta + \pi\f$. Getting that wrong put every conjugated
    torsion's minimum at 90 degrees (PRD-108); it is one conversion, in one
    place, for that reason.
*/
IMPBFFEXPORT IMP::core::Cosine* torsion_cosine(const FFTorsionType& type);

//! Bonded and element-aware nonbonded restraints for a typed dye system.
/*!
    Bonds and angles become harmonics, torsions the CHARMM cosine above,
    impropers a harmonic about the *current* geometry's dihedral, and every
    non-excluded site pair a lower-bound harmonic at its Lennard-Jones
    \f$R_{min}\f$. A site with no particle is skipped rather than an error:
    a system may describe more than the caller decorated.

    \param[in] model the model the restraints score in
    \param[in] system the typed system (bonds, angles, torsions, impropers)
    \param[in] site_ids,particles the decorated sites, parallel; a site id
               that appears twice takes its last particle
    \param[in] nonbonded add the repulsion (#create_steric_restraint). False
               leaves it to the caller, which is what a run wants when it
               scores rigid moves against the repulsion *alone*.
    \throw ValueException when the two sequences disagree in length

    \note A bond whose equilibrium length is not set (zero, which no real bond
    has) is restrained about the geometry as it stands, and an angle likewise.
    That is what a system built without a template carries, and refusing it
    would refuse the systems this module builds itself.
*/
IMPBFFEXPORT IMP::Restraints create_probe_restraints(
        IMP::Model* model, const ProbeForceFieldSystem& system,
        const std::vector<std::string>& site_ids,
        const IMP::ParticleIndexes& particles, bool nonbonded = true);

//! The steric term alone: one soft-sphere restraint over the whole system.
/*!
    Every site pair the system does not exclude, in a single
    #IMP::container::PairsRestraint over a #IMP::core::SoftSpherePairScore.

    **This is the repulsion, for molecular dynamics and Monte Carlo alike.**
    It is differentiable, so a dynamics run can use it; it is one restraint
    over a container rather than one per pair, so a long run can afford it;
    and it depends only on where the spheres are, so a rigid move can be
    scored against it and nothing else -- a rigid move cannot change a bond,
    an angle or a torsion, and scoring those during one is work whose answer
    never changes. One repulsion serves both: soft spheres for the Monte-Carlo
    step *and* for dynamics, rather than these plus per-pair Lennard-Jones
    lower bounds elsewhere, which would be two implementations of one piece of
    physics.

    The Lennard-Jones parameters remain what
    #IMP::bff::IntramolecularEnergy *evaluates*: an energy of a
    conformation is a different question from keeping two atoms apart, and a
    sphere overlap answers the second.

    The exclusions are the system's own (#ProbeForceFieldSystem::get_exclusions),
    so this and the bonded terms cannot disagree about which pairs are 1-2,
    1-3 or 1-4.

    \param[in] model,system,site_ids,particles as for #create_probe_restraints
    \param[in] k the soft-sphere force constant; negative takes the system's
               own `nonbonded.k`
    \return the restraint; null when the system's non-bonded term is off or
            it has no non-excluded pair
*/
IMPBFFEXPORT IMP::Restraint* create_steric_restraint(
        IMP::Model* model, const ProbeForceFieldSystem& system,
        const std::vector<std::string>& site_ids,
        const IMP::ParticleIndexes& particles, double k = -1.0);

//! Place a guest rigidly around a host and keep the best-scoring pose.
/*!
    A random search, not an optimisation: each trial rotates the guest about a
    random axis by a random angle and puts its centre a fixed distance from
    the host's, in a random direction; the pose that scores lowest is the one
    left in the model. That is enough to start a simulation somewhere
    plausible rather than wherever the input files happened to put the two
    components -- which for a dye and a protein read from separate files is
    usually on top of each other.

    \param[in] scoring_function what a trial is judged by
    \param[in] model the model both sets of particles belong to
    \param[in] host,guest the two particle sets; only \p guest moves
    \param[in] distance the separation of the two centres, A
    \param[in] n_trials how many poses to try
    \param[in] seed the random seed, so a run repeats
    \return the best score, or NaN when there is nothing to place

    The guest keeps its own shape throughout: every trial applies one rotation
    and one translation to the same input coordinates, so this cannot distort
    what it is placing.
*/
IMPBFFEXPORT double place_guest_by_score(
        IMP::ScoringFunction* scoring_function, IMP::Model* model,
        const IMP::ParticleIndexes& host, const IMP::ParticleIndexes& guest,
        double distance = 62.5, int n_trials = 200, int seed = 42);

//! Native-contact (Go) restraints that hold a component in its own shape.
/*!
    Every heavy-atom pair of a component that is closer than \p cutoff and is
    not already within two bonds gets a harmonic at the distance it currently
    has. That is what keeps a component's fold while its parts are free to
    move: the bonded terms fix the local geometry, these fix the tertiary one.

    Hydrogens are left out (a name beginning with `H`): they add pairs without
    adding shape, and they are the majority of the atoms.

    \param[in] model,system the model and the typed system
    \param[in] site_ids,particles the decorated sites, parallel
    \param[in] site_atom_names the atom name of each site id, for the hydrogen
               test; a site missing from it counts as heavy
    \param[in] component which component to restrain
    \param[in] only_sites when non-empty, a pair is kept only if at least one
               of its two sites is in this set -- how a mostly-rigid component
               releases a few atoms and restrains only what they touch
    \param[in] k the force constant
    \param[in] cutoff the contact distance, A
    \return one restraint per contact, named `go_<component>_<a>_<b>`

    \note The equilibrium is `max(d, 1 A)`: two atoms that a structure has
    placed on top of each other would otherwise get a harmonic at zero, which
    is a singularity a minimiser walks straight into.
*/
IMPBFFEXPORT IMP::Restraints create_go_restraints(
        IMP::Model* model, const ProbeForceFieldSystem& system,
        const std::vector<std::string>& site_ids,
        const IMP::ParticleIndexes& particles,
        const std::map<std::string, std::string>& site_atom_names,
        const std::string& component,
        const std::vector<std::string>& only_sites = std::vector<std::string>(),
        double k = 3.0, double cutoff = 6.0);

//! `_ff_lj_type` entries (`LJ_<elem>`) for a set of elements.
IMPBFFEXPORT std::map<std::string, FFLJType> get_lj_type_table(
        const std::vector<std::string>& elements);

//! Evaluates the internal LJ energy of dye conformations.
/*!
    The pair list is every non-excluded site pair of \p system (see
    #get_lj_pair_sites), so bonded neighbours do not sterically block each
    other. \c evaluate takes one frame; \c evaluate_batch one per frame;
    \c evaluate_batch_filtered pre-filters with a #BoundingBoxFilter and gives
    non-overlapping frames zero energy.
*/
class IMPBFFEXPORT IntramolecularEnergy {
    std::vector<LJSitePair> pairs_;
    std::vector<int> idx_a_, idx_b_;
    std::vector<double> rmin_, eps_;

public:
    IntramolecularEnergy(const ProbeForceFieldSystem& system);

    std::vector<LJSitePair> get_pairs() const { return pairs_; }

    double evaluate(const std::vector<double>& coords, int n_atoms) const;
    std::vector<double> evaluate_batch(const std::vector<double>& coords,
                                       int n_frames, int n_atoms) const;
    EnergyMaskResult evaluate_batch_filtered(
            const std::vector<double>& coords, int n_frames, int n_atoms,
            const std::vector<double>& reference_coords, int n_ref_atoms,
            double pad = 3.5) const;
};

IMPBFF_END_NAMESPACE

#endif // IMPBFF_SCORING_H
