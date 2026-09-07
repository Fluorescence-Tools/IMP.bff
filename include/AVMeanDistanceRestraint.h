/**
 *  \file IMP/bff/AVMeanDistanceRestraint.h
 *  \brief A FRET restraint on the distance between two mean dye positions.
 *
 * The cheap half of #IMP::bff::ProbeNetworkRestraint. That one rebuilds both
 * accessible volumes on every evaluation, which is the whole cost of an
 * AV-restrained optimisation; this one scores the separation of the two volumes'
 * *mean positions*, converts it to the modelled observable through a cached
 * #IMP::bff::FRETDistanceConverter, and scores that with the measurement's
 * asymmetric chi-squared.
 *
 * That is an approximation, and its content is that the shape of the volume does
 * not change when the structure moves — true while the volumes are carried as
 * rigid-body members and false as soon as the linker's environment changes.
 *
 * **Derivatives are computed.** The score depends on the coordinates only
 * through \f$d_{mp}\f$, so \f$dS/dd_{mp}\f$ comes from a 1-D central difference
 * through the cached converter and maps onto the two particles analytically as
 * \f$\partial d_{mp}/\partial x = \pm\hat r\f$. That is what lets an IMP
 * optimiser dock by minimisation rather than only by Monte-Carlo.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_AVMEANDISTANCERESTRAINT_H
#define IMPBFF_AVMEANDISTANCERESTRAINT_H

#include <IMP/bff/bff_config.h>
#include <IMP/bff/AV.h>
#include <IMP/bff/ProbeNetworkRestraint.h>
#include <IMP/bff/States.h>

#include <IMP/OptimizerState.h>
#include <IMP/core/Harmonic.h>
#include <IMP/Restraint.h>
#include <IMP/RestraintSet.h>
#include <IMP/particle_index.h>

IMPBFF_BEGIN_NAMESPACE

//! Score the separation of two mean dye positions against a measurement.
class IMPBFFEXPORT AVMeanDistanceRestraint : public IMP::Restraint {
    IMP::ParticleIndex p1_, p2_;
    AVPairDistanceMeasurement measurement_;
    FRETDistanceConverter converter_;
    double weight_;
    double max_force_;

    //! The chi-squared at a given mean-position separation.
    double score_at(double d_mp) const;

public:
    //! \param[in] m the model the two particles belong to
    /*! \param[in] p1,p2 the two dye particles; their `XYZ` coordinates are the
               mean positions
        \param[in] measurement the experimental distance, its asymmetric error,
               its Förster radius and which distance convention it is in
        \param[in] sigma the per-component width of the separation vector, for
               the \f$R_{mp} \to \langle R_{DA}\rangle\f$ conversion
        \param[in] weight multiplies the score
        \param[in] max_force FPS's `MaxForce`: past
               \f$\Delta r_{max} = F_{max}\sigma_{exp}^2/2\f$ the restraint
               becomes **linear** rather than growing as a parabola
               (#IMP::bff::chi2_score_capped). **0, the default, is no cap**
               and is the plain asymmetric \f$\chi^2\f$ this restraint has
               always scored; FPS's docking default is 400. */
    /*! \note The two particles are \c ParticleIndexAdaptor, which is how
        `IMP::core`'s own restraints take theirs: a `Particle`, an `AV`
        decorator or a bare index all convert. */
    AVMeanDistanceRestraint(IMP::Model* m, IMP::ParticleIndexAdaptor p1,
                            IMP::ParticleIndexAdaptor p2,
                            const AVPairDistanceMeasurement& measurement,
                            double sigma = 6.0, double weight = 1.0,
                            double max_force = 0.0);

    double get_weight_factor() const { return weight_; }
    void set_weight_factor(double w) { weight_ = w; }
    //! FPS's `MaxForce`; 0 is uncapped.
    double get_max_force() const { return max_force_; }
    void set_max_force(double f) { max_force_ = f; }
    const FRETDistanceConverter& get_converter() const { return converter_; }
    //! The measurement this restraint scores against.
    const AVPairDistanceMeasurement& get_measurement() const {
        return measurement_;
    }
    //! The score at a given mean-position separation, without moving anything.
    /*! What proves the linear tail is linear: evaluate at two separations past
        the knee and divide. */
    double get_score_at(double d_mp) const { return weight_ * score_at(d_mp); }

    virtual double unprotected_evaluate(
            IMP::DerivativeAccumulator* accum) const override;
    virtual IMP::ModelObjectsTemp do_get_inputs() const override;

    IMP_OBJECT_METHODS(AVMeanDistanceRestraint);
};

//! A flat-bottom distance restraint between two labelling positions.
/*!
    The restraint form an MD engine wants, and the one this module was missing.
    #IMP::bff::AVMeanDistanceRestraint scores a chi-squared, which pulls at
    every separation and never stops; a flat-bottom well is **zero** inside the
    experimental error bars and only pushes back outside them, so a trajectory
    can move freely through everything the data permits and is steered only
    where it does not.

    The functional form is AMBER's `&rst`, so a restraint written here and one
    written for an MD engine are the same restraint:

    \f[
    E(d) = \begin{cases}
      2k_2(r_1-r_2)(d-r_1) + k_2(r_1-r_2)^2 & d < r_1 \\
      k_2 (d-r_2)^2                         & r_1 \le d < r_2 \\
      0                                     & r_2 \le d \le r_3 \\
      k_3 (d-r_3)^2                         & r_3 < d \le r_4 \\
      2k_3(r_4-r_3)(d-r_4) + k_3(r_4-r_3)^2 & d > r_4
    \end{cases}
    \f]

    Outside \f$[r_1, r_4]\f$ the well continues **linearly**, tangent to the
    parabola: the force is capped there rather than growing without bound,
    which is what keeps a badly-placed starting structure from blowing up the
    first integration step.

    \note The four bounds are distances between the two particles. They are in
          whatever convention the caller converted to -- for an experimental
          FRET distance that conversion is
          #IMP::bff::rmp_flat_bottom_bounds(), which is not a constant offset.
*/
class IMPBFFEXPORT AVFlatBottomRestraint : public IMP::Restraint {
    IMP::ParticleIndex p1_, p2_;
    double r1_, r2_, r3_, r4_, k2_, k3_;

public:
    //! \param[in] m the model the two particles belong to
    /*! \param[in] p1,p2 the two particles the distance is measured between
        \param[in] r1,r2,r3,r4 the well: harmonic in `[r1,r2]` and `[r3,r4]`,
                   flat in `[r2,r3]`, linear outside. Must be non-decreasing.
        \param[in] k2,k3 the force constants of the lower and upper walls
        \param[in] name the restraint's name
        \throw ValueException when the bounds are not non-decreasing */
    AVFlatBottomRestraint(IMP::Model* m, IMP::ParticleIndexAdaptor p1,
                          IMP::ParticleIndexAdaptor p2, double r1, double r2,
                          double r3, double r4, double k2, double k3,
                          std::string name = "AVFlatBottomRestraint%1%");

    //! The bounds, as `(r1, r2, r3, r4)`.
    std::vector<double> get_bounds() const;
    //! Move the well. Used when the volumes are rebuilt during a run.
    /*! \throw ValueException when the bounds are not non-decreasing */
    void set_bounds(double r1, double r2, double r3, double r4);
    //! The force constants, as `(k2, k3)`.
    std::vector<double> get_force_constants() const;
    //! The two particles the distance is measured between.
    IMP::ParticleIndexes get_particle_indexes() const;
    //! The distance between the two particles right now.
    double get_distance() const;

    virtual double unprotected_evaluate(
            IMP::DerivativeAccumulator* accum) const override;
    virtual IMP::ModelObjectsTemp do_get_inputs() const override;

    IMP_OBJECT_METHODS(AVFlatBottomRestraint);
};
IMP_OBJECTS(AVFlatBottomRestraint, AVFlatBottomRestraints);

//! One labelling position as a particle an MD engine can move.
/*! What Olga writes into an AMBER restart file as a dummy atom, here as an
    IMP particle: an `XYZR` with a mass, tethered to its attachment atom by a
    stiff harmonic bond at the offset the volume itself has, so it follows the
    site as the structure moves. The FRET restraints are then ordinary
    two-particle restraints between these. */
struct IMPBFFEXPORT ProbeParticle {
    IMP::ParticleIndex particle;    //!< the probe
    IMP::ParticleIndex attachment;  //!< the atom it is tethered to
    std::string position_name;      //!< the fps.json position it stands for
    double tether_k;                //!< tether force constant, kcal/mol/A^2

    //! \note There is deliberately no `offset` here. The tether length is a
    //!       *live* quantity -- IMP::bff::AVRebuildOptimizerState moves it when
    //!       it rebuilds the volumes -- and it is held in one place, the
    //!       tether itself. Ask IMP::bff::MDRestraintSystem::get_tether_length().
    IMP_SHOWABLE_INLINE(ProbeParticle,
                        out << "ProbeParticle(" << position_name << ")");
};
typedef std::vector<ProbeParticle> ProbeParticles;

//! Everything needed to run an fps.json network on an MD engine.
/*! \see IMP::bff::md_flat_bottom_restraints() */
class IMPBFFEXPORT MDRestraintSystem {
    IMP::Pointer<IMP::RestraintSet> restraints_;
    IMP::Pointer<ProbeNetworkRestraint> network_;
    AVFlatBottomRestraints wells_;
    IMP::Vector<IMP::Pointer<IMP::core::Harmonic> > tethers_;
    ProbeParticles probes_;
    std::vector<std::string> pair_names_;
    std::vector<double> bounds_;

public:
    MDRestraintSystem() {}
    MDRestraintSystem(IMP::RestraintSet* rs, ProbeNetworkRestraint* network,
                      const AVFlatBottomRestraints& wells,
                      const IMP::Vector<IMP::Pointer<IMP::core::Harmonic> >& tethers,
                      const ProbeParticles& probes,
                      const std::vector<std::string>& pair_names,
                      const std::vector<double>& bounds)
            : restraints_(rs), network_(network), wells_(wells),
              tethers_(tethers),
              probes_(probes), pair_names_(pair_names), bounds_(bounds) {}

    //! The FRET wells and the probe tethers, as one set to add to a model.
    IMP::RestraintSet* get_restraints() const { return restraints_; }

    //! The network that owns the volumes, for rebuilding them mid-run.
    /*! \see IMP::bff::AVRebuildOptimizerState */
    ProbeNetworkRestraint* get_network() const { return network_; }

    //! Just the FRET wells, in get_pair_names() order.
    /*! A restraint read back out of a #IMP::RestraintSet arrives as a base
        #IMP::Restraint, so a caller that wants the distance or the bounds of
        one would have to downcast it. They are kept here instead. */
    const AVFlatBottomRestraints& get_wells() const { return wells_; }
    //! The harmonic holding each probe to its site, in probe order.
    /*! Exposed so that a rebuild can move a tether when the volume's offset
        from its attachment changes. */
    const IMP::Vector<IMP::Pointer<IMP::core::Harmonic> >& get_tethers() const {
        return tethers_;
    }

    //! The current tether length of probe \p i, A.
    /*! The one place it lives: a rebuild moves the tether, and reading it from
        anywhere else would read a value that was true when the system was
        built. */
    double get_tether_length(unsigned int i) const {
        return i < tethers_.size() ? tethers_[i]->get_mean() : 0.0;
    }

    //! One probe per labelling position, in fps.json name order.
    const ProbeParticles& get_probes() const { return probes_; }
    //! The restrained pairs, in the order get_bounds() reports them.
    const std::vector<std::string>& get_pair_names() const {
        return pair_names_;
    }
    //! The wells, flat: `r1, r2, r3, r4` per pair.
    const std::vector<double>& get_bounds() const { return bounds_; }

    IMP_SHOWABLE_INLINE(MDRestraintSystem,
                        out << "MDRestraintSystem(" << probes_.size()
                            << " probes, " << pair_names_.size() << " pairs)");
};

//! Build flat-bottom FRET restraints an IMP molecular-dynamics run can use.
/*!
    The one call between an fps.json and a restrained MD trajectory. For every
    distance in \p score_set it

    1. builds both accessible volumes on the starting structure,
    2. converts the measured distance and its asymmetric errors to
       mean-position bounds with #IMP::bff::rmp_flat_bottom_bounds(),
    3. creates a probe particle per labelling position, tethered to its
       attachment atom at the volume's own offset, and
    4. adds an #IMP::bff::AVFlatBottomRestraint between the two probes.

    The force constants follow the error bars, as the reference implementation
    does: \f$k = F_{max} / (2\sigma)\f$ on each side, so a model one error bar
    outside the well feels \p f_max and no restraint can pull harder than that
    however wrong it is.

    The volumes are computed **once**, on the starting structure. That is the
    approximation this buys its speed with, and it is the same one
    #IMP::bff::AVMeanDistanceRestraint makes: the shape of a volume is taken
    not to change as the structure moves. For a trajectory that stays near its
    start -- which a restrained refinement does -- the error is small; for one
    that does not, rebuild and call again.

    \param[in] hier the structure, already at its starting coordinates
    \param[in] fps_json the labelling and distance file
    \param[in] score_set which distances to restrain; empty uses every one
    \param[in] f_max the force at one error bar, kcal/mol/A
    \param[in] tether_k the force constant tethering a probe to its atom
    \param[in] probe_mass the mass given to a probe particle, Da
    \param[in] tether_atom the atom a probe is tethered to, by name. Empty
               tethers to the volume's own attachment atom, which is what an
               all-atom run wants; `"CA"` tethers to the alpha carbon of the
               same residue, which is what a run that moves only alpha carbons
               needs -- a probe tethered to an atom nothing optimises cannot
               follow the structure.
    \return the restraint set, the probes and the bounds that were derived
    \throw IOException when \p fps_json cannot be read
    \throw ValueException when a position it names is not in \p hier
*/
IMPBFFEXPORT MDRestraintSystem md_flat_bottom_restraints(
        const IMP::core::Hierarchy& hier, const std::string& fps_json,
        std::string score_set = "", double f_max = 5.0,
        double tether_k = 10.0, double probe_mass = 100.0,
        std::string tether_atom = "");

// --------------------------------------------------------------------------
// out to OpenMM
// --------------------------------------------------------------------------

//! The OpenMM `CustomBondForce` energy expression for the flat-bottom well.
/*!
    The same piecewise function #IMP::bff::AVFlatBottomRestraint evaluates,
    written in OpenMM's expression language with `r1, r2, r3, r4, k2, k3` as
    per-bond parameters. One string, so the restraint that runs here and the
    restraint that runs there cannot drift apart:

    ```python
    force = openmm.CustomBondForce(IMP.bff.openmm_flat_bottom_energy())
    for name in ("r1", "r2", "r3", "r4", "k2", "k3"):
        force.addPerBondParameter(name)
    ```

    \note The expression is in **OpenMM's units** -- `r` in nm and the energy in
          kJ/mol -- so the parameters written by
          #IMP::bff::write_openmm_restraints() are converted, not the ones this
          module scores with.
*/
IMPBFFEXPORT std::string openmm_flat_bottom_energy();

//! Write an MD restraint system as an OpenMM setup document.
/*!
    A JSON document describing everything an OpenMM script has to add to a
    `System` built from the same structure: the probe particles, their tethers,
    and one flat-bottom bond per measured pair.

    **Positions and force constants are converted to OpenMM's units** -- nm and
    kJ/mol -- and the document says so in a `units` block, because a restraint
    exported in the wrong unit is a restraint that runs and is wrong.

    Attachment atoms are identified by **chain, residue number and atom name**,
    not by index. An index depends on how the reader built its topology --
    whether it kept hydrogens, waters, altlocs -- and a spec that names atoms by
    position in someone else's file is a spec that silently restrains the wrong
    atoms.

    `examples/labels/plot_fret_restrained_md.py` writes one and shows the dozen
    lines of OpenMM that consume it.

    \param[in] system what #IMP::bff::md_flat_bottom_restraints() returned
    \param[in] path where to write the document
    \param[in] hier the structure the attachment atoms are named from
    \throw IOException when \p path cannot be opened
*/
IMPBFFEXPORT void write_openmm_restraints(const MDRestraintSystem& system,
                                          const std::string& path,
                                          const IMP::core::Hierarchy& hier);

//! Write a runnable OpenMM script that sets up and runs these restraints.
/*!
    The document #IMP::bff::write_openmm_restraints() writes still needs someone
    to write the OpenMM around it. This writes that too: a **self-contained
    Python script** with the restraint table embedded, which loads \p pdb_path,
    builds a `System` with \p force_field, adds the probes, their tethers and
    the flat-bottom wells, minimises, and runs.

    It is a starting point that runs, not a production protocol -- the
    integrator, the solvent model and the run length are the first things a
    caller will change, and they are at the top of the generated file for that
    reason.

    \param[in] system what #IMP::bff::md_flat_bottom_restraints() returned
    \param[in] path where to write the `.py`
    \param[in] hier the structure the attachment atoms are named from
    \param[in] pdb_path the structure the generated script loads
    \param[in] force_field OpenMM force-field XMLs, space separated
    \param[in] n_steps how many steps the generated script runs
    \throw IOException when \p path cannot be opened
*/
IMPBFFEXPORT void write_openmm_script(
        const MDRestraintSystem& system, const std::string& path,
        const IMP::core::Hierarchy& hier, const std::string& pdb_path,
        const std::string& force_field = "amber14-all.xml amber14/tip3pfb.xml",
        int n_steps = 100000);

//! Rebuild the accessible volumes during a run and re-derive the wells.
/*!
    `md_flat_bottom_restraints` computes both volumes once, on the starting
    structure, and every well's bounds follow from those shapes. That is the
    approximation the whole scheme buys its speed with, and it decays: as the
    structure moves, a volume's shape and its offset from the attachment change,
    and the \f$R_{mp}\f$ that reproduces a measured \f$\langle R_{DA}\rangle\f$
    changes with them.

    This is an #IMP::OptimizerState that resamples the volumes every \p period
    steps and updates each well's bounds and each probe's tether length from the
    geometry the trajectory has actually reached. A run with it is an
    approximation refreshed at a stated interval; a run without it is the same
    approximation held fixed for the whole trajectory.

    \note Resampling every volume is the expensive part of an AV calculation, so
          \p period trades accuracy for speed directly. It is a *periodic*
          correction, not a per-step one, deliberately.
*/
class IMPBFFEXPORT AVRebuildOptimizerState : public IMP::OptimizerState {
    IMP::Pointer<ProbeNetworkRestraint> network_;
    AVFlatBottomRestraints wells_;
    ProbeParticles probes_;
    IMP::Vector<IMP::Pointer<IMP::core::Harmonic> > tethers_;
    std::vector<std::string> pair_names_;
    int n_updates_;

public:
    //! \param[in] m the model
    /*! \param[in] system what md_flat_bottom_restraints() returned; it carries
               the network whose volumes are resampled
        \param[in] period how many optimiser steps between rebuilds */
    AVRebuildOptimizerState(IMP::Model* m, const MDRestraintSystem& system,
                            unsigned int period = 1000);

    //! How many times the volumes have been rebuilt.
    int get_number_of_updates() const { return n_updates_; }

    //! Resample now and re-derive every well, whatever the period says.
    void update_now();

    virtual void do_update(unsigned int call) override;

    IMP_OBJECT_METHODS(AVRebuildOptimizerState);
};

//! Make each volume a member of its attachment atom's rigid body.
/*! The coordinates of an accessible volume are the mean of its density, so a
    volume is resampled before it is added: the position it is added *at* is
    the one the rigid body will carry it by. Volumes whose source atom is not
    in a rigid body are left alone. */
IMPBFFEXPORT void add_avs_to_rigid_bodies(const AVs& avs);

//! Give each volume a radius and a mass, so it occupies space.
/*! The radius is the largest of the volume's three; the mass is twice that.
    Without this a volume is a point, and nothing keeps two of them -- or a
    volume and the structure -- from sharing the same coordinates. */
IMPBFFEXPORT void set_av_xyzr_mass(const AVs& avs);

//! Every restraint an fps.json network asks for, as one set.
/*!
    Everything a PMI `RestraintBase` wrapper would do beyond PMI's own
    bookkeeping is here. IMP.pmi is not one of this module's modules, so a
    caller that wants a `RestraintBase` subclasses one around this set, in the
    program that already imports PMI.

    \param[in] hier the structure the labelled residues are in
    \param[in] fps_json the network: positions, distances and score sets
    \param[in] name names the restraint, and is the score set's key
    \param[in] score_set which set of distances to score; empty means all
    \param[in] mean_position_restraint score the separation of the volumes'
               *mean positions* (one #IMP::bff::AVMeanDistanceRestraint per
               distance) instead of rebuilding both volumes per evaluation.
               That is an approximation whose content is that the shape of a
               volume does not change when the structure moves, so the volumes
               are made rigid-body members first.
    \param[in] sigma_DA the per-component width for the mean-position branch
    \param[in] occupy_volume give the volumes a radius and a mass
    \param[in] weight multiplies the set's score
    \return the set, owned by the caller (an IMP::Object; keep a reference)
    \throw IOException when \p fps_json cannot be read
*/
IMPBFFEXPORT IMP::RestraintSet* probe_network_restraint_set(
        const IMP::core::Hierarchy& hier, const std::string& fps_json,
        std::string name = "ProbeNetworkRestraint",
        std::string score_set = "", bool mean_position_restraint = false,
        double sigma_DA = 6.0, bool occupy_volume = true, double weight = 1.0);

IMPBFF_END_NAMESPACE

#endif //IMPBFF_AVMEANDISTANCERESTRAINT_H
