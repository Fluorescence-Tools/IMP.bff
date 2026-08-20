/**
 *  \file IMP/bff/PETQuenching.h
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
#ifndef IMPBFF_PETQUENCHING_H
#define IMPBFF_PETQUENCHING_H

#include <IMP/bff/bff_config.h>

#include <IMP/value_macros.h>
#include <IMP/showable_macros.h>

#include <limits>
#include <map>
#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! The dye the bundled PET table was measured for -- a xanthene, Alexa488-like.
/*! Peulen et al., J. Phys. Chem. B 2017, 121, 8211. */
IMPBFFEXPORT extern const char* const REFERENCE_DYE;

//! Radius of the dye sphere, Angstrom, when a caller does not say.
IMPBFFEXPORT extern const double DEFAULT_DYE_RADIUS;

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
               Dyes that are harder to reduce or oxidise use a value below one.
    \param[in] slow_factor diffusion scaling applied near every residue
    \param[in] dye_radius Angstrom. The trajectory tracks the dye *centre*, so
               the reference surface contact distances are offset by this radius
               to give centre-to-centre quench radii.
*/
IMPBFFEXPORT std::map<std::string, ResidueQuenching>
amino_acid_quenching_defaults(double kQ_scale = 1.0, double slow_factor = 1.0,
                              double dye_radius = 3.5);

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

#endif //IMPBFF_PETQUENCHING_H
