/**
 *  \file IMP/bff/InteractionTerms.h
 *  \brief The channels that deactivate an excited dye, one object each.
 *
 * A term is **one deactivation channel**, and it answers one question: the rate
 * constant of every state of the dye, in 1/ns. Parallel channels add, which is
 * what #IMP::bff::total_rate does and why these are objects rather than numbers
 * computed inside an observable — the field solver already relies on summing a
 * quenching map and a FRET map, and burying either inside its own function is
 * what left `fret_rate_trace` and `fret_rate_map` as two incompatible calls.
 *
 * Terms are **representation-agnostic**: they consume #IMP::bff::States, so one
 * implementation serves an accessible volume, a rotamer library, a
 * coarse-grained model and an MD trajectory alike. A term that needs
 * orientations says so, and a representation that cannot supply them (an
 * accessible volume) is told rather than silently averaged.
 *
 * Arity distinguishes state-dependent from emergent:
 *
 * - **1-body** — #IMP::bff::RadiativeTerm. The dye and the solvent only.
 * - **2-body** — #IMP::bff::PETTerm (a dye against quencher atoms),
 *   #IMP::bff::FRETTerm (a dye against a dye).
 * - **N-body** — homo-FRET and multi-chromophore transfer, where the rate is
 *   not a sum over pairs. Not implemented; the interface is shaped to take it.
 *
 * Every term takes the same two participants and reports through `arity`
 * whether it reads the second, so #IMP::bff::total_rate sums them without
 * knowing which kind each one is.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_INTERACTIONTERMS_H
#define IMPBFF_INTERACTIONTERMS_H

#include <IMP/bff/bff_config.h>
#include <IMP/bff/States.h>
#include <IMP/bff/ProbeLibrary.h>
#include <IMP/bff/Quenching.h>

#include <IMP/Object.h>
#include <IMP/bff/Base.h>

#include <limits>
#include <map>
#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! One deactivation channel.
class IMPBFFEXPORT InteractionTerm : public IMP::Object {
public:
    InteractionTerm(std::string name) : IMP::Object(name) {}

    //! How many participants the term is defined over.
    /*! 1 for a process of the dye alone, 2 for a pair, more for a genuinely
        many-body process. */
    virtual int get_arity() const = 0;

    //! True when the term needs transition dipoles.
    /*! A caller can then tell whether a representation without them (an
        accessible volume) forces an isotropic assumption instead of resolving
        the orientation. */
    virtual bool get_needs_orientations() const { return false; }

    //! Rate constants in 1/ns, one per state of \p first.
    /*! \param[in] first the dye whose states are being rated
        \param[in] second the second participant; ignored by a 1-body term */
    virtual std::vector<double> rate_constants(
            const States& first,
            const States& second = IMP::bff::States()) const = 0;

    IMP_OBJECT_METHODS(InteractionTerm);
};
IMP_OBJECTS(InteractionTerm, InteractionTerms);

//! 1-body: the dye's own decay, \f$1/\tau_0\f$.
/*!
    The floor every other channel adds to. Depends on the dye and its medium,
    not on where it is or what is near it — which is exactly what "1-body"
    means here.
*/
class IMPBFFEXPORT RadiativeTerm : public InteractionTerm {
    double lifetime_;

public:
    //! \param[in] lifetime the unquenched fluorescence lifetime, ns
    RadiativeTerm(double lifetime);

    double get_lifetime() const { return lifetime_; }
    int get_arity() const override { return 1; }
    std::vector<double> rate_constants(
            const States& first,
            const States& second = IMP::bff::States()) const override;
};

//! 2-body: photoinduced electron transfer between a dye and quencher atoms.
/*!
    \f$k(r) = \sum_a k_{Q,a} \exp(-(|r - r_a| - r_{dye}) / r_{C,a})\f$

    The quenching atoms are part of the term, not a per-call argument: resolving
    which of a structure's atoms quench, and how hard, is a function of the
    residue and atom names and the dye's #IMP::bff::PETParameters, and it does
    not change from one set of dye states to the next. Doing it once at
    construction is also what lets the rate be a plain loop.
*/
class IMPBFFEXPORT PETTerm : public InteractionTerm {
    std::vector<double> coords_;    //!< flat, three per quenching atom
    std::vector<double> kQ_;        //!< 1/ns, per atom; zero where inactive
    std::vector<double> rC_;        //!< Angstrom, per atom
    double probe_radius_;

public:
    //! \param[in] parameters `{comp_id: PETParameters}` for **one** dye
    /*! \param[in] res_names,atom_names one per atom, naming it
        \param[in] coords,n_atoms,n_dim the atoms, `(N, 3)`
        \param[in] probe_radius subtracted from the centre-to-centre distance,
                   because the tabulated contact distances are measured from the
                   dye *surface* */
    PETTerm(const std::map<std::string, PETParameters>& parameters,
            const std::vector<std::string>& res_names,
            const std::vector<std::string>& atom_names,
            double* coords, int n_atoms, int n_dim,
            double probe_radius = 3.5);

    double get_probe_radius() const { return probe_radius_; }
    //! How many atoms carry a non-zero rate.
    unsigned int get_n_active() const;
    //! Per-atom \f$k_Q\f$, as a numpy view.
    void get_kQ(double** out_view, int* n_out_view) const;
    //! Per-atom \f$r_C\f$, as a numpy view.
    void get_rC(double** out_view, int* n_out_view) const;

    int get_arity() const override { return 2; }
    std::vector<double> rate_constants(
            const States& first,
            const States& second = IMP::bff::States()) const override;
};

//! 2-body: Förster transfer between two dyes.
/*!
    \f$k(r) = (1/\tau_0) (R_0/r)^6 \kappa^2 / (2/3)\f$

    **R0 is derived, not supplied** — from the two dyes' spectra, the donor's
    quantum yield, the medium's refractive index and \f$\kappa^2\f$. Passing it
    in was how `forster_radius=52.0` came to be a default in a dozen signatures.
*/
class IMPBFFEXPORT FRETTerm : public InteractionTerm {
    Probe donor_, acceptor_;
    double refractive_index_;
    double kappa2_;             //!< NaN resolves it from the participants
    double r_min_;

public:
    //! \param[in] donor,acceptor the two species
    /*! \param[in] refractive_index of the medium between them
        \param[in] kappa2 orientation factor; NaN falls back to the isotropic
                   2/3, reported through `get_used_isotropic_kappa2`
        \param[in] r_min closest approach of the two dye centres, Angstrom */
    FRETTerm(const Probe& donor, const Probe& acceptor,
             double refractive_index = 1.4,
             double kappa2 = std::numeric_limits<double>::quiet_NaN(),
             double r_min = 7.0);

    //! \f$R_0\f$ in Angstrom, derived from the pair and the medium.
    double get_forster_radius() const;
    //! True when no \f$\kappa^2\f$ was given, so the isotropic 2/3 is assumed.
    bool get_used_isotropic_kappa2() const { return kappa2_ != kappa2_; }
    double get_kappa2() const { return kappa2_; }

    int get_arity() const override { return 2; }
    bool get_needs_orientations() const override { return true; }
    std::vector<double> rate_constants(
            const States& first,
            const States& second = IMP::bff::States()) const override;
};

//! Sum the channels.
/*!
    Parallel deactivation channels add, which is the property that makes a
    quenching map and a FRET map summable on the same grid. Doing it here rather
    than inside each observable is what keeps that true when a channel is added.

    \throw ValueException when \p terms is empty
*/
IMPBFFEXPORT std::vector<double> total_rate(
        const InteractionTerms& terms, const States& first,
        const States& second = IMP::bff::States());

IMPBFF_END_NAMESPACE

#endif //IMPBFF_INTERACTIONTERMS_H
