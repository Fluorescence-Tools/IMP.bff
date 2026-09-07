#ifndef IMPBFF_QUENCHING_H
#define IMPBFF_QUENCHING_H

/**
 *  \file IMP/bff/Quenching.h
 *  \brief Collisional quenching of a diffusing dye: the tables, the fields, the decay.
 *
 * One header for the physics below #IMP::bff::QuenchedDonorDecay, in the
 * order it is used:
 *
 * 1. **Tables** (formerly `PETQuenching.h`) -- which moieties quench
 *    (#IMP::bff::Quencher), how hard for a given dye (#IMP::bff::PETParameters),
 *    and what a diffusing dye sees at a residue (#IMP::bff::ResidueQuenching).
 * 2. **Stamping** (formerly `QuenchingGrid.h`) -- spheres of influence onto an
 *    accessible-volume grid: `stamp_spheres` and the combines built on it.
 * 3. **Fields** (formerly `QuenchingMap.h`) -- mobility, quenching-rate and
 *    FRET-rate maps on that grid, each a named composition of a kernel above.
 * 4. **The race** (formerly `QuenchedDecay.h`) -- walk, read the rate along
 *    the walk, race photons against it.
 *
 * `QuenchingModel.h` stays a header of its own: it composes these with the
 * accessible volume, the walk and the grid solver, and it reaches
 * `InteractionTerms.h`, which reaches the tables here -- so the model cannot
 * live in the same file as the tables without including itself.
 *
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

// -------- from PETQuenching.h --------
/**
 *  (formerly IMP/bff/PETQuenching.h, now a section of this file)
 *  \brief Photoinduced electron transfer: which moieties quench, and how hard.
 *
 * The chemistry, as three separable things:
 *
 * - **Identity.** #IMP::bff::Quencher — which residue, which atoms. Deliberately
 *   *not* CB: electron transfer happens at the indole ring, the phenol, the
 *   thioether or the thiol, and stamping a rate on CB puts it up to 4 A from
 *   the chemistry.
 * - **Rate.** #IMP::bff::PETParameters — a **pair** property, because \f$k_Q\f$
 *   depends on the redox potentials of dye *and* quencher. The bundled table
 *   was measured for one xanthene dye; asking for another returns the same
 *   numbers with `measured_for` set to that one, so a transfer is visible as an
 *   assumption rather than passing for a measurement.
 * - **Field.** #IMP::bff::ResidueQuenching — what a diffusing dye sees per
 *   residue: a stickiness factor, a rate, a contact radius and the atoms that
 *   define the centre.
 *
 * The published rates are *starting points to be calibrated against measured
 * lifetimes*, which is what `kQ_scale` is for; PRD-111 found it recoverable to
 * a few percent from six sites jointly and not at all from one.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#include <IMP/bff/bff_config.h>

#include <IMP/bff/Base.h>

#include <limits>
#include <map>
#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! The dye the bundled PET table was measured for -- a xanthene, Alexa488-like.
/*! Peulen et al., J. Phys. Chem. B 2017, 121, 8211. */
IMPBFFEXPORT extern const char* const REFERENCE_DYE;

//! Radius of the probe sphere, Angstrom, when a caller does not say.
IMPBFFEXPORT extern const double DEFAULT_PROBE_RADIUS;

//! A quenching moiety: which residue, which atoms. **No rate.**
/*!
    The rate lives in #IMP::bff::PETParameters, because it takes two partners to
    define one.
*/
struct IMPBFFEXPORT Quencher {
    //! Residue name, e.g. `"TRP"`. flrCIF `_flr_poly_probe_position.comp_id`.
    std::string comp_id;
    //! The redox-active atoms. flrCIF `_flr_poly_probe_position.atom_id`.
    std::vector<std::string> atom_ids;
    //! Set when this is a *specific* residue in a structure rather than a type.
    std::string asym_id;
    //! Set with `asym_id`; ignored when this describes a type.
    int seq_id;

    Quencher(std::string comp_id = "",
             const std::vector<std::string>& atom_ids =
                     std::vector<std::string>(),
             std::string asym_id = "", int seq_id = -1);

    //! True when this describes a residue *type* rather than one residue.
    bool get_is_typed() const { return asym_id.empty() && seq_id < 0; }

    //! The same moiety, located at one residue of a structure.
    Quencher at(std::string asym_id, int seq_id) const;

    IMP_SHOWABLE(Quencher);
};

//! PET between one dye and one quencher type -- a **pair** property.
struct IMPBFFEXPORT PETParameters {
    //! Chromophore name the parameters apply to.
    std::string dye;
    //! Quencher residue type.
    std::string comp_id;
    //! \f$k_Q\f$ in 1/ns at contact.
    double rate_constant;
    //! From the *dye surface* to the quenching centre, Angstrom -- roughly van
    //! der Waals contact plus the offset from the moiety centroid to its outer
    //! atoms. Surface-relative on purpose, so the geometry transfers between
    //! dyes of different size even when \f$k_Q\f$ does not.
    double contact_distance;
    //! \f$r_C\f$, the exponential decay length of the through-space rate. NaN
    //! for a hard contact-sphere model.
    double attenuation_length;
    //! The dye the values were actually measured with. When it differs from
    //! `dye` the parameters were **transferred**, which is an assumption about
    //! redox chemistry rather than a measurement.
    std::string measured_for;

    PETParameters(std::string dye = "", std::string comp_id = "",
                  double rate_constant = 0.0, double contact_distance = 0.0,
                  double attenuation_length =
                          std::numeric_limits<double>::quiet_NaN(),
                  std::string measured_for = "");

    //! True when these values were measured with a *different* dye.
    bool get_is_transferred() const { return measured_for != dye; }

    //! The same pair with \f$k_Q\f$ multiplied -- what a calibration turns.
    PETParameters scaled(double rate_scale) const;

    IMP_SHOWABLE(PETParameters);
};

//! What a diffusing dye sees near one residue type.
struct IMPBFFEXPORT ResidueQuenching {
    //! Diffusion scaling near the residue (unspecific stickiness), in [0, 1].
    double slow_factor;
    //! \f$k_Q\f$, 1/ns. Zero for a residue that does not quench.
    double kQ;
    //! Contact radius, Angstrom. **NaN means "inherit the model-wide critical
    //! distance"**, which is how a project sets one radius for everything and
    //! overrides it per residue type where it matters.
    double quench_radius;
    //! The atoms whose centroid is the quenching centre.
    std::vector<std::string> quench_atoms;

    ResidueQuenching(double slow_factor = 1.0, double kQ = 0.0,
                     double quench_radius =
                             std::numeric_limits<double>::quiet_NaN(),
                     const std::vector<std::string>& quench_atoms =
                             std::vector<std::string>());

    IMP_SHOWABLE(ResidueQuenching);
};

IMP_VALUES(Quencher, Quenchers);
IMP_VALUES(PETParameters, PETParametersList);
IMP_VALUES(ResidueQuenching, ResidueQuenchings);

//! The reference \f$k_Q\f$ and dye-surface contact distance per quencher.
struct IMPBFFEXPORT PETReference {
    double kQ;
    double contact_distance;
    PETReference(double kQ = 0.0, double contact_distance = 0.0)
        : kQ(kQ), contact_distance(contact_distance) {}
    IMP_SHOWABLE_INLINE(PETReference, out << "PETReference(kQ=" << kQ << ")");
};
IMP_VALUES(PETReference, PETReferences);

//! The twenty standard amino-acid residue names.
IMPBFFEXPORT std::vector<std::string> standard_amino_acid_residues();

//! The redox-active atoms of each residue type. Not CB — see #IMP::bff::Quencher.
IMPBFFEXPORT std::map<std::string, std::vector<std::string> > quencher_atoms();

//! The published PET chemistry: \f$k_Q\f$ and contact distance per quencher.
IMPBFFEXPORT std::map<std::string, PETReference> pet_quenching_reference();

//! The redox-active moieties, read from the one table that defines them.
/*! Identity only — see #IMP::bff::reference_pet_parameters for the rates, which
    depend on the dye. */
IMPBFFEXPORT std::map<std::string, Quencher> reference_quenchers();

//! The published PET chemistry for one dye against each quencher type.
/*!
    \param[in] dye chromophore the parameters are wanted for
    \param[in] rate_scale multiplies every \f$k_Q\f$
    \param[in] attenuation_length \f$r_C\f$ for the exponential through-space
               form; NaN keeps the hard contact-sphere model
*/
IMPBFFEXPORT std::map<std::string, PETParameters> reference_pet_parameters(
        std::string dye = "AlexaFluor488", double rate_scale = 1.0,
        double attenuation_length = std::numeric_limits<double>::quiet_NaN());

//! The full per-residue interaction table, with defaults filled in.
/*!
    Accepts a partial table keyed by residue name. Unknown residue names are
    kept, so a non-standard residue can be given a rate. Slow factors are
    clamped to [0, 1], rates to non-negative, and a non-positive or non-finite
    radius becomes NaN — "inherit the model-wide critical distance".
*/
IMPBFFEXPORT std::map<std::string, ResidueQuenching>
normalize_amino_acid_quenching(
        const std::map<std::string, ResidueQuenching>& table =
                std::map<std::string, ResidueQuenching>());

//! A full interaction table built from #IMP::bff::pet_quenching_reference.
/*!
    \param[in] kQ_scale dye-specific multiplier on every reference \f$k_Q\f$.
               Probes that are harder to reduce or oxidise use a value below one.
    \param[in] slow_factor diffusion scaling applied near every residue
    \param[in] probe_radius Angstrom. The trajectory tracks the dye *centre*, so
               the reference surface contact distances are offset by this radius
               to give centre-to-centre quench radii.
*/
IMPBFFEXPORT std::map<std::string, ResidueQuenching>
amino_acid_quenching_defaults(double kQ_scale = 1.0, double slow_factor = 1.0,
                              double probe_radius = 3.5);

//! The stickiness factor of each residue, in \p residue_names order.
IMPBFFEXPORT std::vector<double> slow_factors_for_residues(
        const std::vector<std::string>& residue_names,
        const std::map<std::string, ResidueQuenching>& table);

//! The quenching rate (1/ns) of each residue.
IMPBFFEXPORT std::vector<double> quenching_rates_for_residues(
        const std::vector<std::string>& residue_names,
        const std::map<std::string, ResidueQuenching>& table);

//! The contact radius of each residue, inheriting \p critical_distance.
/*! A NaN `quench_radius` in the table means "use the model-wide critical
    distance". */
IMPBFFEXPORT std::vector<double> quench_radii_for_residues(
        const std::vector<std::string>& residue_names,
        const std::map<std::string, ResidueQuenching>& table,
        double critical_distance = 0.0);

//! One slow centre and one quench centre per residue, plus its type.
class IMPBFFEXPORT ResidueSites {
    std::vector<double> slow_centers_;    //!< flat, three per residue
    std::vector<double> quench_centers_;  //!< flat, three per residue
    std::vector<std::string> residue_names_;

public:
    //! Append one residue: its slow centre, its quench centre and its type.
    void add(const double* slow, const double* quench,
             const std::string& residue_name);

    unsigned int size() const {
        return static_cast<unsigned int>(residue_names_.size());
    }
    std::vector<std::string> get_residue_names() const { return residue_names_; }
    //! The centres as numpy views, `(n, 3)` once reshaped.
    void get_slow_centers(double** out_view, int* n_out_view) const;
    void get_quench_centers(double** out_view, int* n_out_view) const;

    IMP_SHOWABLE_INLINE(ResidueSites,
                        out << "ResidueSites(" << residue_names_.size()
                            << " residues)");
};
IMP_VALUES(ResidueSites, ResidueSitesList);

//! Group atoms by residue and locate its slow and quench centres.
/*!
    The slow centre is CB, or CA when there is no CB, or the residue's first
    atom. The quench centre is the centroid of the atoms the table names as
    redox-active, falling back to the slow centre when the residue has none.

    Residues are keyed by `(chain, res_id, res_name)`. **Keying on `res_id`
    alone is wrong** and was a real defect in QuEst: residue numbers restart per
    chain, so in a homodimer every number occurs twice and two residues' atoms
    were folded into one centre.

    \param[in] chains,res_ids,res_names,atom_names one per atom
    \param[in] coords,n_atoms,n_dim the atoms, `(N, 3)`
    \param[in] table whose `quench_atoms` decide the quench centres
*/
IMPBFFEXPORT ResidueSites residue_sites(
        const std::vector<std::string>& chains,
        const std::vector<int>& res_ids,
        const std::vector<std::string>& res_names,
        const std::vector<std::string>& atom_names,
        double* coords, int n_atoms, int n_dim,
        const std::map<std::string, ResidueQuenching>& table =
                std::map<std::string, ResidueQuenching>());

//! Per-atom quenching rate \f$k_Q\f$ and characteristic distance \f$r_C\f$.
/*!
    Atoms not named in \p quencher do not quench, and come back as zeros.

    \param[in] res_names,atom_names one per atom
    \param[in] parameters `{comp_id: PETParameters}` for one dye; an atom
               quenches when its residue is in here *and* its name is one of
               that residue's redox-active atoms
    \param[out] out_kQ,n_out_kQ,out_rC,n_out_rC one value per atom
*/
IMPBFFEXPORT void atomic_quenching_parameters(
        const std::vector<std::string>& res_names,
        const std::vector<std::string>& atom_names,
        const std::map<std::string, PETParameters>& parameters,
        double** out_kQ, int* n_out_kQ, double** out_rC, int* n_out_rC);

IMPBFF_END_NAMESPACE

// -------- from QuenchingGrid.h --------
/**
 *  (formerly IMP/bff/QuenchingGrid.h, now a section of this file)
 *  \brief Stamping spheres of influence onto an accessible-volume grid.
 *
 * A quencher slows a dye near it and quenches it on contact. Both are stamped
 * the same way -- a sphere of a given radius around each centre, combined into
 * whatever is already there -- and differ only in how they combine: stickiness
 * *multiplies* (identity 1) and rates *add* (identity 0), because rates of
 * parallel channels add.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */


IMPBFF_BEGIN_NAMESPACE

//! How a stamped sphere combines with what is already on the grid.
enum GridCombine {
    GRID_COMBINE_MULTIPLY = 0,  //!< stickiness: factors multiply, identity 1
    GRID_COMBINE_ADD = 1        //!< rates: parallel channels add, identity 0
};

//! Voxel indices and integer radii of a set of sphere centres.
/*!
    Uses `floor`, **not** truncation. `int()` truncates toward zero, so a centre
    on the negative side of the origin would round *up* while every other map
    here rounds down -- the walk's own occupancy test and the trajectory sampler
    both floor. The two disagreed by one voxel per axis for every centre with a
    negative offset, which is half the grid.

    \param[in] rs centre coordinates, flat, three per centre
    \param[in] r0 grid origin (the attachment point), three doubles
    \param[in] dg voxel edge
    \param[in] ng voxels per axis
    \param[in] radius per-centre radius in Angstrom
    \param[out] ix0,iy0,iz0 voxel indices of each centre
    \param[out] radius_idx radius of each centre in voxels
*/
IMPBFFEXPORT void center_grid_indices(
        const std::vector<double>& rs,
        const std::vector<double>& r0,
        double dg,
        int ng,
        const std::vector<double>& radius,
        std::vector<int>& ix0,
        std::vector<int>& iy0,
        std::vector<int>& iz0,
        std::vector<int>& radius_idx
);

//! Stamp per-centre values into spheres on the accessible volume.
/*!
    Visits only the voxels inside each centre's bounding box rather than testing
    every voxel against every centre, and keeps the outer loop over x-slabs so
    the work is free of write races.

    Sphere membership is strict, \f$d^2 < r^2\f$, so a voxel exactly on the
    boundary is outside.

    Voxels where \p density is zero keep the identity: the dye cannot be there,
    so nothing is stamped.

    \param[in] density accessible-volume density, flat, ng^3
    \param[in] ng voxels per axis
    \param[in] radius per-centre radius
    \param[in] rs centre coordinates, flat
    \param[in] r0 grid origin
    \param[in] dg voxel edge
    \param[in] values per-centre value -- a factor, or a rate
    \param[in] combine GRID_COMBINE_MULTIPLY or GRID_COMBINE_ADD
*/
IMPBFFEXPORT void stamp_spheres(
        const std::vector<double>& density,
        int ng,
        const std::vector<double>& radius,
        const std::vector<double>& rs,
        const std::vector<double>& r0,
        double dg,
        const std::vector<double>& values,
        int combine
,
        double** out_view, int* n_out_view);

//! Per-voxel diffusion scaling from overlapping sticky spheres.
/*!
    Stickiness **multiplies** where spheres overlap. Voxels outside the
    accessible volume keep 1.0, so the factor is only meaningful where the walk
    can go.

    \param[in] density binary occupancy of the accessible volume, flat ng^3
    \param[in] ng voxels per axis
    \param[in] dg voxel edge, A
    \param[in] slow_radius radius of each sticky sphere
    \param[in] rs sphere centres, flat
    \param[in] r0 the grid anchor
    \param[in] slow_fact factor in [0, 1] per centre
*/
IMPBFFEXPORT void slow_factor_grid(
        const std::vector<double>& density,
        int ng,
        double dg,
        const std::vector<double>& slow_radius,
        const std::vector<double>& rs,
        const std::vector<double>& r0,
        const std::vector<double>& slow_fact,
        double** out_view, int* n_out_view);

//! Per-voxel quenching rate (1/ns) from overlapping quencher spheres.
/*!
    Rates **add** where contact spheres overlap, which is the composition rule
    for independent PET channels. Same geometry and the same indexing as
    slow_factor_grid(); only the accumulator differs.

    \param[in] values per-centre rate constant, 1/ns
*/
IMPBFFEXPORT void quenching_rate_grid(
        const std::vector<double>& density,
        int ng,
        double dg,
        const std::vector<double>& radius,
        const std::vector<double>& rs,
        const std::vector<double>& r0,
        const std::vector<double>& values,
        double** out_view, int* n_out_view);

//! The contact ("slow") part of an accessible volume, as a binary mask.
/*!
    A thin composition over split_contact_volume(), which labels every voxel;
    for a binary density the contact part is what the PET model wants.

    \param[out] out_view_i,n_out_view_i flat ng^3, 1 in contact and 0 elsewhere
*/
IMPBFFEXPORT void av_contact_mask(
        const std::vector<double>& density,
        int ng,
        double dg,
        const std::vector<double>& slow_radius,
        const std::vector<double>& rs,
        const std::vector<double>& r0,
        int** out_view_i, int* n_out_view_i);

IMPBFF_END_NAMESPACE

// -------- from QuenchingMap.h --------
/**
 *  (formerly IMP/bff/QuenchingMap.h, now a section of this file)
 *  \brief Mobility, quenching-rate and FRET-rate fields on an AV grid.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */


IMPBFF_BEGIN_NAMESPACE

//! Slow the mobility field wherever the dye is in contact with an atom.
/*!
    The factor is applied **once per contacting atom**, so it compounds: a voxel
    touching 45 atoms at a factor of 0.985 ends at 0.5, and at 0.9 it ends at
    9e-3. Only values within a whisker of 1.0 are meaningful, which is a
    property of the model rather than a taste.

    \param[in] d_map the unmodified mobility field, flat ng^3
    \param[in] density accessible-volume density; zero voxels stay zero
    \param[in] axis voxel offsets from the grid centre, length ng
    \param[in] r0 grid origin
    \param[in] atoms_xyz obstacle coordinates, flat, three per atom
    \param[in] min_distance_sq squared contact distance
    \param[in] factor per-contact slowing factor
*/
IMPBFFEXPORT void slow_near_atoms(
        const std::vector<double>& d_map,
        const std::vector<double>& density,
        const std::vector<double>& axis,
        const std::vector<double>& r0,
        const std::vector<double>& atoms_xyz,
        double min_distance_sq,
        double factor,
        double** out_view, int* n_out_view
);

//! Quenching-rate field: \f$k(r) = 1/\tau_0 + \sum_a k_a e^{-(|r-r_a| - r_{dye})/r_{C,a}}\f$.
/*!
    Voxels outside the accessible volume stay at **zero**, not at \f$1/\tau_0\f$:
    the dye cannot be there, so it has no decay rate there.

    \param[in] density accessible-volume density
    \param[in] axis voxel offsets, length ng
    \param[in] r0 grid origin
    \param[in] atoms_xyz obstacle coordinates, flat
    \param[in] kQ per-atom rate constant; zero means the atom does not quench
    \param[in] rC per-atom attenuation length
    \param[in] probe_radius subtracted from the centre-to-centre distance
    \param[in] inv_tau0 the radiative floor
*/
IMPBFFEXPORT void quenching_map(
        const std::vector<double>& density,
        const std::vector<double>& axis,
        const std::vector<double>& r0,
        const std::vector<double>& atoms_xyz,
        const std::vector<double>& kQ,
        const std::vector<double>& rC,
        double probe_radius,
        double inv_tau0,
        double** out_view, int* n_out_view
);

//! FRET-rate field of a donor volume against a whole acceptor volume.
/*!
    **This is the harmonic mean** -- the acceptor-weighted mean transfer *time*
    is accumulated and inverted, not the mean rate. That is the static limit, in
    which the acceptor does not move within the donor's excited-state lifetime.
    The arithmetic mean of rates is the fast-exchange limit and is what
    ``fret_rate_trace`` computes along a trajectory. The two are different
    physics and differ measurably; neither is a variant of the other.

    \param[in] density_d,density_a the two densities, flat
    \param[in] axis_d,axis_a voxel offsets of each grid
    \param[in] r0_d,r0_a the two grid origins
    \param[in] r0_6 the Forster radius to the sixth power
    \param[in] kf radiative rate
    \param[in] step stride over the acceptor grid; 1 visits every voxel
*/
IMPBFFEXPORT void fret_map(
        const std::vector<double>& density_d,
        const std::vector<double>& density_a,
        const std::vector<double>& axis_d,
        const std::vector<double>& axis_a,
        const std::vector<double>& r0_d,
        const std::vector<double>& r0_a,
        double r0_6,
        double kf,
        int step,
        double** out_view, int* n_out_view
);

//! Voxel-centre offsets from the grid anchor, in Angstrom.
/*!
    \f$(i - (ng-1)/2)\,dg\f$ with the **integer** centre offset, the same one
    grid_center_index() uses, so a map built here indexes the way the Brownian
    walk samples it.
*/
IMPBFFEXPORT std::vector<double> grid_axis(int ng, double dg);

//! A diffusion map from a radial profile, evaluated on an integer-radius axis.
/*!
    \f$f\f$ is a profile of the distance from the grid anchor, sampled once per
    integer Angstrom: \p radial[k] is the value at \f$r = k\f$. Every voxel gets
    \f$\p radial[\operatorname{round}(|r|)]\f$, clamped to the table -- the old
    Python took a *callable* and there is no C++ spelling of one that does not
    call back into the interpreter per voxel, so the caller lands \f$f\f$ on a
    radius grid first. This is what a turnover lamp (``stretched linker'')
    mobility looks like before #slow_near_atoms adds the local crowding.

    \param[in] density binary occupancy, flat ng^3, whose grid side gives ng
    \param[in] dg voxel edge, A
    \param[in] radial \f$f\f$ at integer radii in Angstrom; length 0 is the
               identity (returns 1.0 everywhere)
    \param[out] out_view,n_out_view flat ng^3
*/
IMPBFFEXPORT void radial_diffusion_map(
        const std::vector<double>& density, double dg,
        const std::vector<double>& radial, double** out_view, int* n_out_view);

//! The mobility field: a base coefficient, slowed by nearby atoms.
/*!
    \param[in] density binary occupancy of the accessible volume, flat ng^3
    \param[in] r0 the grid anchor
    \param[in] dg voxel edge, A
    \param[in] atoms_xyz obstacle coordinates, flat
    \param[in] free_diffusion the unhindered dye diffusion coefficient, A^2/ns,
               used wherever \p base is not given
    \param[in] min_distance contact distance, roughly dye radius + 2 vdW
    \param[in] slow_factor factor applied once per contacting atom, in [0, 1]
    \param[in] base an optional per-voxel base coefficient replacing the
               constant \p free_diffusion -- a radial profile, say. Length 0
               takes the constant.
*/
IMPBFFEXPORT void diffusion_coefficient_map(
        const std::vector<double>& density,
        const std::vector<double>& r0,
        double dg,
        const std::vector<double>& atoms_xyz,
        double free_diffusion,
        double min_distance,
        double slow_factor,
        const std::vector<double>& base,
        double** out_view, int* n_out_view);

//! Total decay rate per voxel: intrinsic plus exponential PET.
/*!
    \f$k(r) = 1/\tau_0 + \sum_a k_{Q,a}\,e^{-(|r - r_a| - r_{dye})/r_{C,a}}\f$

    The distance is measured from the dye **surface**. Unlike the contact-sphere
    law in quenching_rate_grid() there is no cut-off: a distant atom contributes
    exponentially little rather than nothing. Voxels outside the accessible
    volume stay at zero -- **not** at \f$1/\tau_0\f$; the solver reads this as
    the rate field on a domain masked by the same occupancy, so an unreachable
    voxel carries no rate at all.

    \param[in] tau0 unquenched lifetime, ns; non-positive drops the floor
*/
IMPBFFEXPORT void quenching_rate_map(
        const std::vector<double>& density,
        const std::vector<double>& r0,
        double dg,
        const std::vector<double>& atoms_xyz,
        const std::vector<double>& kQ,
        const std::vector<double>& rC,
        double tau0,
        double probe_radius,
        double** out_view, int* n_out_view);

//! An effective FRET rate for every donor voxel, from the acceptor cloud.
/*!
    A fixed donor position does not have *a* FRET rate: it has a distribution of
    them, one per accessible acceptor position. This approximates that sum by a
    single exponential whose rate is the reciprocal of the **mean transfer
    time** -- the *harmonic* mean, dominated by the slow, distant acceptor
    positions. fret_rate_trace() takes the arithmetic mean instead, because it
    models a fast-exchanging acceptor where the *rates* average.

    \param[in] kf the donor's radiative rate, \f$1/\tau_0\f$
    \param[in] acceptor_step stride over the acceptor grid; the cost is the
               product of the two grids' sizes
    \throw ValueException for a non-positive \p kf
*/
IMPBFFEXPORT void fret_rate_map(
        const std::vector<double>& density_donor,
        const std::vector<double>& density_acceptor,
        const std::vector<double>& r0_donor,
        const std::vector<double>& r0_acceptor,
        double dg_donor,
        double dg_acceptor,
        double forster_radius,
        double kf,
        int acceptor_step,
        double** out_view, int* n_out_view);

IMPBFF_END_NAMESPACE

// -------- from QuenchedDecay.h --------
/**
 *  (formerly IMP/bff/QuenchedDecay.h, now a section of this file)
 *  \brief Walk, read the quenching rate along the walk, and race photons —
 *         without the trajectory ever leaving C++.
 *
 * A PET-quenching prediction is three steps -- walk, read the rate map along
 * the walk, race photons against that trace -- and driving them as three calls
 * puts an array between each. The trajectory is the largest object in the
 * chain and **nobody wants it**: a
 * 5 000 000-step walk is 20 million doubles crossing the boundary twice to
 * produce a few thousand photons.
 *
 * This is the same computation with the intermediates kept where they are made.
 * It is not an approximation and not a different model: given the same seeds it
 * reproduces the three-call path exactly, which is what
 * `test/quenching/test_fused_decay.py` asserts.
 *
 * The separate calls remain, and should: a trajectory is what a correlation
 * function or a visualisation needs, and the rate trace is worth looking at on
 * its own. This is the path for when only the decay is wanted, which is the
 * common case and the one inside every fitting loop.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */


IMPBFF_BEGIN_NAMESPACE

//! Photon trace of a dye diffusing in its accessible volume, fused.
/*!
    Runs \p walk_seeds.size() independent walks, concatenates the quenching rate
    each one saw, and races \p n_photons excitations against the result.
    Concatenated rather than averaged, because each photon draws its own start
    frame — a longer record is exactly what that wants.

    The rate is read at the voxel the walk is standing on, with no coordinate
    round trip. It is the same voxel that
    `floor((xyz - x0)/dg + (ng-1)//2)` names: the attachment point cancels, and
    the walk's own occupancy test has already resolved the index.

    \param[in] occupancy,n_occupancy flat ng^3, nonzero where the dye may be
    \param[in] mobility,n_mobility flat ng^3 step-variance scaling, or length 0
    \param[in] rate_map,n_rate_map flat ng^3 quenching rate, 1/ns

    The three grids are taken as raw buffers so numpy's arrays pass straight
    through; converting them into `std::vector`s costs ~34 ns per element, and
    there are three of them.
    \param[in] ng,dg grid size and voxel edge
    \param[in] t_max,t_step walk duration and step, ns
    \param[in] diffusion_coefficient A^2/ns, per-component variance 2 D dt
    \param[in] walk_seeds one per independent trajectory
    \param[in] tau0 intrinsic lifetime, ns
    \param[in] n_photons excitation events
    \param[in] photon_seed reproducible when non-negative
    \param[out] stats five values: frames, accepted, rejected, mean rate,
                collision fraction — everything the caller would have derived
                from the trajectory it no longer receives
    \return two values per photon: delay time in ns (0 if quenched), and 1 or 0
            for emitted. Empty if no walk found a starting voxel.
*/
IMPBFFEXPORT void quenched_donor_photons(
        int* occupancy, int n_occupancy,
        double* mobility, int n_mobility,
        double* rate_map, int n_rate_map,
        int ng, double dg, double t_max, double t_step,
        double diffusion_coefficient,
        const std::vector<int>& walk_seeds,
        double tau0, int n_photons, int photon_seed,
        std::vector<double>& stats,
        double** out_view, int* n_out_view);

IMPBFF_END_NAMESPACE


#endif  // IMPBFF_QUENCHING_H
