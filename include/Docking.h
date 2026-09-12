/**
 *  \file IMP/bff/Docking.h
 *  \brief What a docking run is told, and what it reports.
 *
 * The values a rigid-body FRET docking passes around -- its control
 * parameters, one experiment-against-model distance, and the outcome -- plus
 * the two things that read a scored assembly: the per-pair table and the CSV
 * it is written as.
 *
 * The *engine* (build the assembly, run the Monte-Carlo or the minimisation,
 * write the RMF) is a program and lives in `bin/`; what is here is what a
 * program and a library user both need.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_DOCKING_H
#define IMPBFF_DOCKING_H

#include <IMP/bff/bff_config.h>
#include <IMP/bff/ProbeNetworkRestraint.h>
#include <IMP/bff/VdwRadii.h>

#include <IMP/Pointer.h>
#include <IMP/RestraintSet.h>
#include <IMP/algebra/Transformation3D.h>
#include <IMP/algebra/Vector3D.h>
#include <IMP/atom/Hierarchy.h>
#include <IMP/core/rigid_bodies.h>
#include <IMP/core/RestraintsScoringFunction.h>
#include <IMP/OptimizerState.h>
#include <IMP/container/ListSingletonContainer.h>

#include <map>
#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! One experimental distance against its model value (FPS-style diagnostics).
struct IMPBFFEXPORT PairDistance {
    std::string name, position1, position2;
    double distance_exp, distance_model;
    //! `RDAMean`, `RDAMeanE` or `Rmp`, as the fps.json spells it.
    std::string distance_type;
    double forster_radius;
    double error_neg, error_pos;
    //! FPS's **bond**: both ends are plain atoms, so this is not a FRET pair.
    /*! A crosslink, an EPR anchor, a covalent tie between subunits. It is
        scored like any other distance -- what it changes is that its two
        anchor atoms stop clashing and that its energy is reported separately
        as `E_bond` (#IMP::bff::ProbeNetworkRestraint::get_is_bond). */
    bool is_bond;

    PairDistance()
        : distance_exp(0), distance_model(0), distance_type("RDAMean"),
          forster_radius(52.0), error_neg(0), error_pos(0), is_bond(false) {}

    //! Model minus data, Å: a model that is too large is positive.
    double get_residual() const { return distance_model - distance_exp; }
    //! The asymmetric \f$\chi^2\f$ of this pair (#IMP::bff::chi2_score).
    double get_chi2() const;
    double get_efficiency_model() const;
    double get_efficiency_exp() const;
    //! The pair as a JSON object, with the four derived numbers included.
    std::string get_json() const;

    IMP_SHOWABLE_INLINE(PairDistance,
                        out << "PairDistance(" << name << ": model "
                            << distance_model << " vs " << distance_exp
                            << ")");
};
IMP_VALUES(PairDistance, PairDistances);

//! One atom of the frame a fixed labelling position was measured in.
/*!
    A labelling position may be an `XYZ` -- a mean dye position measured once,
    on one structure, rather than a volume to be simulated. Such a coordinate
    is only meaningful in the frame it was measured in, so FPS lets the
    position carry the atoms that *define* that frame (`FilterEngine.cs`,
    `ReferenceAtom` / `ReadRefAtoms`). Screening then superimposes those atoms
    onto each library structure and carries the fixed coordinate along.

    \p x, \p y, \p z are where the atom sat in the measuring frame. The
    remaining fields are how the *same* atom is found again in another
    structure: chain / residue / atom name, or -- as the legacy FPS files
    spell it -- the PDB serial number, which is what `ATOM` lines trailing an
    LP line carry.
*/
struct IMPBFFEXPORT ReferenceAtom {
    std::string chain_identifier;
    int residue_seq_number;
    std::string atom_name;
    //! PDB serial in the file the frame was measured in; non-positive is unset.
    /*! A serial identifies an atom only within the file it was written in.
        It is honoured because the legacy files have nothing else, and it is
        tried *after* chain/residue/atom name for the same reason. */
    int atom_serial;
    double x, y, z;

    ReferenceAtom()
        : residue_seq_number(0), atom_serial(0), x(0), y(0), z(0) {}

    IMP::algebra::Vector3D get_coordinates() const {
        return IMP::algebra::Vector3D(x, y, z);
    }

    IMP_SHOWABLE_INLINE(ReferenceAtom,
                        out << "ReferenceAtom(" << chain_identifier << ":"
                            << residue_seq_number << ":" << atom_name << ")");
};
IMP_VALUES(ReferenceAtom, ReferenceAtoms);

//! Where a fixed position lands on one structure, and how well its frame fits.
struct IMPBFFEXPORT ReferenceFit {
    //! The fps.json position this frame belongs to.
    std::string position;
    //! How many of the declared reference atoms this structure contains.
    int n_atoms;
    //! Sum of squared deviations after the fit, Å².
    /*! FPS pools these over *all* positions before taking the root, so the
        sum -- not the root -- is what a caller adds up (#reference_rmsd). */
    double squared_deviation;
    //! \f$\sqrt{\text{squared\_deviation} / n}\f$: this frame's own fit, Å.
    double rmsd;
    //! The fixed coordinate carried onto this structure.
    IMP::algebra::Vector3D coordinates;
    //! False when no frame could be fitted; \p coordinates is then as given.
    bool fitted;

    ReferenceFit()
        : n_atoms(0), squared_deviation(0), rmsd(0), coordinates(0, 0, 0),
          fitted(false) {}

    IMP_SHOWABLE_INLINE(ReferenceFit,
                        out << "ReferenceFit(" << position << ", " << n_atoms
                            << " atoms, rmsd " << rmsd << ")");
};
IMP_VALUES(ReferenceFit, ReferenceFits);

//! The `reference_atoms` of one fps.json position object, given as JSON text.
/*! An absent or empty array gives an empty list, which is the ordinary case:
    only fixed (`XYZ`) positions normally carry a frame.
    \throw ValueException when \p position_json is not JSON */
IMPBFFEXPORT ReferenceAtoms fps_reference_atoms(
        const std::string& position_json);

//! Superimpose a frame onto a structure and carry a point through it.
/*!
    The Kabsch fit of `FilterEngine.CalculateChi2`: centre both sets on their
    centroids, take the SVD of the cross-covariance, and guard the reflection
    by flipping the third singular direction when the determinant is negative.
    IMP's #IMP::algebra::get_transformation_aligning_first_to_second is that
    algorithm with that guard, so it is what is used rather than a second copy
    of it.

    \param[in] atoms the frame; its coordinates are the *source* side
    \param[in] structure where the same atoms are looked up
    \param[in] point the coordinate to carry along -- the position's fixed
               mean dye position
    \return the fit, with \p coordinates the transported point

    Fewer than three matched atoms leaves the rotation about the remaining
    axis unconstrained, so the transported point would be arbitrary rather
    than merely uncertain; the fit is then refused (`fitted` false, the point
    returned as given). FPS does not check this -- it has no case where it
    happens, having read its frames from files a human wrote.
*/
IMPBFFEXPORT ReferenceFit fit_reference_atoms(
        const ReferenceAtoms& atoms, IMP::atom::Hierarchy structure,
        const IMP::algebra::Vector3D& point = IMP::algebra::Vector3D(0, 0, 0));

//! Fit every frame an fps.json declares onto one structure.
/*!
    \param[in] positions_json the file's `Positions` section, as JSON text
    \param[in] structure the library structure
    \return one entry per position carrying a non-empty `reference_atoms`, in
            the order the section lists them; positions without a frame are
            not reported, because there is nothing to say about them
*/
IMPBFFEXPORT ReferenceFits fit_reference_positions(
        const std::string& positions_json, IMP::atom::Hierarchy structure);

//! FPS's pooled `RefRMSD`: \f$\sqrt{\sum d^2 / \sum n}\f$ over the fits.
/*! Pooled over positions rather than averaged over them -- a frame of twenty
    atoms and one of four do not weigh the same. No fits at all is 0, which is
    what FPS reports for a file with no reference atoms. */
IMPBFFEXPORT double reference_rmsd(const ReferenceFits& fits);

//! Control parameters for rigid-body FRET docking.
struct IMPBFFEXPORT DockingParameters {
    //! Which sampler #IMP::bff::dock runs.
    /*! One interface, four backends -- the names #IMP::bff::MCMCSampler takes:

        - `"metropolis"` (aliases `"mc"`, `"walk"`): a single-chain random
          walk over the pose vector. The plain one, and the default.
        - `"stretch"` (aliases `"emcee"`, `"ensemble"`, `"affine"`): the
          affine-invariant ensemble move. Scale-free, so it does not need the
          mover amplitudes tuned.
        - `"slice"` (aliases `"zeus"`, `"ensemble_slice"`): the ensemble slice
          sampler. Every move is accepted, so there is no acceptance rate to
          tune either.
        - `"de"`: differential-evolution MCMC.

        The pose is `(tx, ty, tz, rx, ry, rz)` per mobile body -- a
        translation and a rotation vector -- so an ensemble sampler over six
        numbers per body is an ordinary small problem, and the same objective
        serves all four. */
    std::string sampler;
    //! Ensemble size for `"stretch"`, `"slice"` and `"de"`; 0 picks 4 * ndim.
    int n_walkers;
    //! Outer sampler iterations, and steps per iteration.
    int n_frames, mc_steps;
    //! kT, and the per-step rigid-body mover amplitudes (Å, radian).
    double mc_temperature, max_translation, max_rotation;
    //! The simulated-annealing schedule and its temperature bounds.
    bool simulated_annealing;
    double sa_tmin, sa_tmax;
    //! How many best-scoring models to write.
    int n_best;
    //! Weight of the excluded-volume (clash) term.
    double ev_weight;
    //! Initial random shuffle amplitude, Å; zero does not shuffle.
    double shuffle_max_translation;
    //! Score the separation of the AV *mean positions* rather than rebuilding
    //! both volumes on every evaluation.
    /*! Recommended for sampling. The approximation's content is that the shape
        of a volume does not change when the structure moves -- true while the
        volumes are carried as rigid-body members, false as soon as the
        linker's environment changes. */
    bool mean_position_restraint;
    //! Width of the mean-position transfer function, Å.
    double sigma_da;
    //! A named score set of the fps.json `chi2` section; empty means all.
    std::string score_set;
    //! The `body_id` held fixed; the others move.
    int fixed_body;
    //! Clash detection on one coarse bead per residue rather than all atoms.
    /*! Minimisation only. About three times cheaper per step, and the step is
        the dominant cost. */
    bool coarse_clash;
    //! FPS-style refinement cycles after docking: re-sample the AVs in the
    //! docked context (inter-body occlusion moves the mean positions) and
    //! minimise again. Zero is off, and off is the default because a cycle
    //! re-runs the one-time AV calculation.
    int refine_av_cycles;
    //! Export the full \f$p(R_{DA})\f$ distributions after docking.
    bool save_distributions;
    //! Record the run's path as an RMF trajectory.
    /*! On by default, and that is not an arbitrary default: a sampling run's
        trajectory *is* what distinguishes it from a minimisation, and a walk
        whose path was thrown away cannot be judged for mixing, for the basin
        it settled in, or for whether it moved at all. #IMP::bff::dock writes a
        frame each time the walk improves -- a frame per *proposal* would be a
        file the size of the run. Ignored where the build has no RMF. */
    bool save_trajectory;
    //! FPS's `OptimizeSelected`: `"Selected"`, `"All"` or `"SelectedThenAll"`.
    /*!
        **"Selected" means *distances*, not molecules** -- the commonest
        misreading of FPS's protocol (`SpringEngine.cs:428` gates on
        `Distance.IsSelected`; `Molecule.Selected` is an unrelated flag about
        which body is randomised). Here the selection is the fps.json score
        set: `"Selected"` scores only #score_set's distances, `"All"` scores
        every distance in the file whatever the score set says, and
        `"SelectedThenAll"` runs two phases of #n_frames / 2 each -- the score
        set first, then everything (`SpringEngine.cs:255`).

        **Clashes are never gated.** FPS evaluates the excluded volume over all
        molecules in every phase (`SpringEngine.cs:326`), and so does this.

        Unrecognised text is `"Selected"`, which is the default and reproduces
        what this function did before the flag existed.
    */
    std::string optimize_selected;
    //! FPS's `MaxForce`: past \f$F_{max}\sigma^2/2\f$ a restraint goes linear.
    /*! The Huber knee of #IMP::bff::chi2_score_capped. FPS ships **400** for
        Dock/Sample/Screening and 10000 for Refine/Error estimation, i.e. the
        robust tail is on for docking and effectively off for refinement --
        the single largest behavioural difference between its modes. 0 is a
        pure parabola. Note which restraints it touches: the knee is at
        \f$F_{max}\sigma^2/2\f$, so it robustifies **tight** restraints
        (\f$\sigma\f$ = 0.1 Å, knee 2 Å) and leaves loose ones alone
        (\f$\sigma\f$ = 3 Å, knee 1800 Å). */
    double max_force;
    //! FPS's `ClashTolerance`, Å: the overlap that costs one \f$\chi^2\f$ unit.
    /*! It becomes a spring constant exactly as an error bar does:
        \f$k_{clash} = 2/\text{ClashTolerance}^2\f$ (`SpringEngine.cs:147`).
        FPS ships 1.0 Å for Dock/Sample/Screening and 0.5 Å for
        Refine/Error estimation -- refinement penalises a clash four times
        harder per angstrom. **The clash term is not force-capped**:
        `drmaxclash` is computed at `SpringEngine.cs:148` and never used, and
        the native kernel takes `kclash` alone, so in the shipping FPS path
        distance restraints are capped and clashes are not. Non-positive keeps
        IMP's historical soft-sphere constant \f$k = 1\f$. */
    double clash_tolerance;
    //! Which van der Waals radii the excluded-volume term measures overlap by.
    /*! `"imp"` (the default) is the radius each particle carries -- after
        `IMP::atom::read_pdb` the CHARMM-derived **united-atom** set, which
        carries implicit hydrogens; `"olga"` is the name-keyed Bondi-scale
        table (#IMP::bff::olga_vdw_radius), which is the set FPS's
        #clash_tolerance was calibrated against. Same spellings and same
        machinery as #IMP::bff::AV::set_radii_source, deliberately: "how big is
        an atom" must have one answer per run even though the volume and the
        clash term ask it separately.

        **Why it exists.** `k_clash = 2/ClashTolerance^2` is FPS's constant and
        FPS applies it to Bondi radii. Applied to IMP's larger ones the same
        interface reads as a far worse clash. Measured statically at the input
        pose over HIV-RT's protein--DNA interface at k = 8: `"imp"` gives 268
        overlapping atom pairs, 90.49 A of total overlap and energy 198.40;
        `"olga"` gives 30, 7.59 A and 11.11 -- a factor of 17.9 in energy.

        \see #clash_radii_scale, and #coarse_clash, which this refuses to
        combine with. */
    std::string clash_radii_source;
    //! Multiplies whatever radius #clash_radii_source gives. Default 1.0.
    /*! Source and scale compose: the source picks the table, the scale
        multiplies it. A caller who only wants smaller spheres uses the scale;
        one who wants FPS's model picks the source.

        **One global scale cannot stand in for the source**, and the numbers
        say so rather than the reasoning. Bisected over the same HIV-RT
        interface, the IMP-radii scale that reproduces Olga's **30 pairs** is
        **0.8206**, the one that reproduces its **7.595 A** of overlap is
        **0.8386**, and the one that reproduces its **energy 11.114** is
        **0.8476** -- three different answers, and at the first of them the
        energy is 5.447, off by a factor of two. The ratio of IMP's radii to
        Bondi's is itself per-element (C 0.808, N 0.838, P 0.841, O 0.889,
        S 0.900 on this structure) and IMP's carbon alone spans 1.70--2.275 A
        where Bondi's is a single 1.70, so one factor can match any *one*
        statistic of the contact distribution and none of the others. It is a
        knob for loosening a clash term, not a radii set.

        \throw ValueException on a non-positive value. */
    double clash_radii_scale;

    DockingParameters()
        : sampler("metropolis"), n_walkers(0), n_frames(500), mc_steps(10),
          mc_temperature(1.0),
          max_translation(4.0), max_rotation(0.1),
          simulated_annealing(false), sa_tmin(1.0), sa_tmax(2.5), n_best(20),
          ev_weight(1.0), shuffle_max_translation(10.0),
          mean_position_restraint(true), sigma_da(6.0), fixed_body(0),
          coarse_clash(true), refine_av_cycles(0), save_distributions(false),
          save_trajectory(true), optimize_selected("Selected"),
          max_force(400.0), clash_tolerance(1.0),
          clash_radii_source("imp"), clash_radii_scale(1.0) {}

    IMP_SHOWABLE_INLINE(DockingParameters,
                        out << "DockingParameters(" << n_frames << " x "
                            << mc_steps << " steps)");
};
IMP_VALUES(DockingParameters, DockingParametersList);

//! The outcome of a docking, scoring or refinement run.
struct IMPBFFEXPORT DockingResult {
    double score;
    int n_avs, n_distances;
    //! FPS's `Ebond`: the part of #score that the bond distances contribute.
    /*!
        **A subset of #score, never an addition to it.** FPS accumulates it in
        the same loop that accumulates the total (`SpringEngine.cs:445`) and
        writes it as the `chi2_bond` column; adding the two together
        double-counts every crosslink. Zero when the file declares no bonds,
        which is the ordinary case.

        In this module's units, i.e. half a \f$\chi^2\f$ per distance, so
        `e_bond / score` is directly the fraction of the restraint energy that
        is not FRET.
    */
    double e_bond;
    //! How many of #pairs are bonds.
    int n_bonds;
    //! FPS's `Eclash`: the excluded-volume term, **not** part of #score's
    //! restraint sum but part of the number FPS stores in `E`.
    /*! FPS keeps the two apart in `SimulationResult` (`SpringEngine.cs:465`
        writes `Eclash`, `:292` writes `E`) and then writes them as two columns
        `chi2` and `chi2_clash`. Here #score is what the scoring function
        returned -- restraints **and** clash, because that is the objective the
        minimiser descended -- so `score - e_clash` is FPS's `E`. Recorded
        rather than derived because a table with a clash column and no way to
        fill it is the reason the column gets dropped. */
    double e_clash;
    //! Did the descent stop moving before the budget ran out?
    /*! **Not FPS's definition, deliberately.** FPS reports
        `niter < MaxIterations` (`SpringEngine.cs:295`), which is
        unconditionally `true` under `SelectedThenAll` (its D5) and says
        nothing about the gradient. Here it is measured: the score at the end
        of two consecutive optimisation chunks differing by less than
        \f$10^{-6}\max(1, |E|)\f$. A run that was stopped early is never
        converged. */
    bool converged;
    /* `PairDistances`, not `std::vector<PairDistance>`: SWIG returns a
       bare std::vector **member** of a temporary as a pointer into freed
       memory, and a freed vector reports size 0 -- so
       `score_structures(...).pairs[0]` raised IndexError, reading as
       "this structure has no distances", while binding the result first
       worked. IMP_VALUES wraps this type properly. */
    PairDistances pairs;
    std::string output_dir, rmf_file, score_csv;
    std::vector<std::string> best_pdbs;
    //! Whatever the run wants to report beyond the fixed fields, as JSON.
    std::string extra;
    //! The docked state, as JSON: one `{"body_id", "t", "q"}` per rigid body.
    /*! Translation and rotation quaternion. Re-applied against the same input
        PDBs they reconstruct the pose exactly, which is what lets a run be
        stored and continued. */
    std::string poses;

    //! \param[in] score the total, and what it was over
    /*! \param[in] n_avs,n_distances how many volumes and distances were used
        \param[in] pairs the experiment-against-model table
        \param[in] output_dir,rmf_file,score_csv,best_pdbs what was written
        \param[in] extra whatever the run reports beyond these, as JSON
        \param[in] poses the docked state, as a JSON array */
    DockingResult(double score = 0.0, int n_avs = 0, int n_distances = 0,
                  const std::vector<PairDistance>& pairs =
                          std::vector<PairDistance>(),
                  const std::string& output_dir = "",
                  const std::string& rmf_file = "",
                  const std::string& score_csv = "",
                  const std::vector<std::string>& best_pdbs =
                          std::vector<std::string>(),
                  const std::string& extra = "{}",
                  const std::string& poses = "[]", double e_bond = 0.0,
                  int n_bonds = 0, double e_clash = 0.0, bool converged = false)
        : score(score), n_avs(n_avs), n_distances(n_distances), e_bond(e_bond),
          n_bonds(n_bonds), e_clash(e_clash), converged(converged),
          pairs(pairs), output_dir(output_dir),
          rmf_file(rmf_file), score_csv(score_csv), best_pdbs(best_pdbs),
          extra(extra), poses(poses) {}

    //! The whole result as JSON, pairs and their derived numbers included.
    std::string get_json() const;

    IMP_SHOWABLE_INLINE(DockingResult,
                        out << "DockingResult(score=" << score << ", "
                            << n_distances << " distances)");
};
IMP_VALUES(DockingResult, DockingResults);

//! A cooperative stop for a run that may take minutes.
/*! Subclass it (in C++ or, through the director, in Python) and answer true
    to end a run early, keeping the pose reached. It is asked between chunks
    of optimisation and not inside the optimiser's own callback, because
    raising out of an `IMP::OptimizerState` mid-step is not reliable. */
class IMPBFFEXPORT DockingStop : public IMP::Object {
public:
    DockingStop(std::string name = "DockingStop%1%") : IMP::Object(name) {}
    //! True ends the run at the next chunk boundary.
    virtual bool should_stop() { return false; }
    IMP_OBJECT_METHODS(DockingStop);
};

//! Append `frame,score` to a CSV while an optimiser runs.
/*! What a progress plot reads. The file is written as the run goes, not at
    the end, so a caller watching it sees the descent live. */
class IMPBFFEXPORT ScoreTrace : public IMP::OptimizerState {
    IMP::Pointer<IMP::ScoringFunction> scoring_function_;
    std::string path_;
    mutable int step_;

public:
    ScoreTrace(IMP::Model* m, IMP::ScoringFunction* scoring_function,
               const std::string& path);
    virtual void do_update(unsigned int call_number) override;
    IMP_OBJECT_METHODS(ScoreTrace);
};

//! The container the excluded-volume term searches for close pairs.
/*!
    \param[in] m,root the model and the structure
    \param[in] coarse one enlarged backbone bead per residue instead of every
               atom. The atomistic close-pair search dominates the cost of a
               minimisation step, and IMP's excluded volume is rigid-body
               aware, so only inter-body pairs are scored either way.
    \param[in] bead_radius the radius those beads are given, Å
    \param[in] radii_source which van der Waals radii the returned spheres
               carry. #IMP::bff::AV_RADII_IMP (the default) is each particle's
               own -- today's behaviour, and the radius `IMP::core::XYZR`
               reports. #IMP::bff::AV_RADII_OLGA is the Bondi-scale table FPS's
               `ClashTolerance` was calibrated against.
    \param[in] radii_scale multiplies whatever the source gives; 1.0 is the
               default and changes nothing

    \note **Nothing is written on the structure's own radii.** Under the
    default (`AV_RADII_IMP`, scale 1.0) this returns the atoms themselves, as
    it always has. Under any other choice each atom gets a **shadow sphere**:
    a new particle at the same coordinate, added to the same rigid body, whose
    only content is the radius the clash term should measure by. The atoms keep
    theirs. That matters because #IMP::bff::AV under `AV_RADII_IMP` inflates
    its obstacles by exactly those radii and #IMP::bff::dock_minimize resamples
    volumes between minimisations (`refine_av_cycles`): rewriting the radii
    here would silently resize every volume computed afterwards. It is also
    FPS's defect D12 -- `SpringEngine` mutates session-global `Molecule`
    objects and never reverts -- and not reproducing it was already a decision
    for #IMP::bff::set_bond_anchor_radii.

    \note The shadows are not hierarchy leaves, so nothing that walks the
    structure (RMSD, pose capture, PDB output) sees them. They *are* rigid-body
    members, which is what carries them when a body moves and what propagates
    the clash gradient back to it.

    \throw ValueException when \p coarse is combined with a non-default
    \p radii_source or \p radii_scale. A coarse bead is one sphere standing for
    a whole residue, so a per-atom van der Waals table has nothing to say about
    it and scaling it is only `bead_radius` spelled twice. Turn the coarse
    representation off (`DockingParameters::coarse_clash`, `--no-coarse-clash`)
    or set \p bead_radius directly.
    \throw ValueException on a non-positive \p radii_scale.
*/
IMPBFFEXPORT IMP::container::ListSingletonContainer* clash_container(
        IMP::Model* m, IMP::atom::Hierarchy root, bool coarse = true,
        double bead_radius = 2.5,
        AVRadiiSource radii_source = AV_RADII_IMP, double radii_scale = 1.0);

//! One structure of a screened library, with what it scored and how it failed.
/*!
    A rank alone says which structure is least bad; it does not say *which*
    restraint the others miss, and that is the question a screen is run to
    answer. The diagnostics are FPS's (`FilterEngine.CalculateChi2` filling a
    `FilteringResult`): the reduced \f$\chi^2\f$, the σ-violation counts, the
    number of distances with no model value, the reference-frame fit, and the
    per-pair table itself.
*/
struct IMPBFFEXPORT ScreenedStructure {
    std::string path;
    //! The assembly's total score, and what the ranking is by.
    /*! Half the \f$\chi^2\f$ sum plus the excluded-volume term: a Gaussian
        restraint is worth \f$-\log L = (d/\sigma)^2/2\f$, which is what
        `AVPairDistanceMeasurement::score_model` returns, so that FRET mixes
        with IMP's other terms at the right weight. Infinite as soon as one
        distance has no model value. */
    double score;
    //! \f$\chi^2\f$ per *scored* distance -- FPS's `E`, and what it ranks by.
    /*! Distances with no model value are left out of both the sum and the
        count, exactly as FPS leaves them out; #score, by contrast, goes
        infinite on the first of them. NaN when none could be scored.
        Carries no clash term and no factor of a half, so on a single body
        with no clashes `2 * score == chi2_r * (n - invalid_r)` (measured on
        T4L/3GUN, `chi2_C2_33p`: **22.5732 and 0.684037** over 33 pairs --
        re-measured 2026-09-01, when the accessible contact volume the file
        asks for stopped being ignored and the obstacle radii became Olga's;
        it read 22.5277 / 1.36532 before either). */
    double chi2_r;
    //! Distances deviating by more than 1, 2 and 3 times their error.
    /*! The error on the side the deviation falls: a model that is too long is
        judged against `error_pos`, one that is too short against `error_neg`.
        The counts nest -- a 3σ outlier is counted in all three. */
    int sigma1, sigma2, sigma3;
    //! Distances with no model value at all (`InvalidR`).
    /*! An accessible volume that came out empty at one end -- a buried site,
        or a linker that cannot leave its attachment atom. Such a distance is
        neither satisfied nor violated, so it is counted rather than scored. */
    int invalid_r;
    //! Pooled Kabsch fit of every declared reference frame on this structure, Å.
    /*! 0 when the file declares no `reference_atoms`, which is FPS's own
        answer for that case; NaN when frames are declared and none could be
        fitted here. A large value means the fixed positions were transported
        onto a structure they do not belong to, and every distance touching
        one of them is then worth as much as that fit. */
    double ref_rmsd;
    //! The experiment-against-model table of this structure.
    /*! What makes a screen readable: which restraint this structure misses,
        by how much, and in which direction. */
    /* `PairDistances`, not `std::vector<PairDistance>`: SWIG returns a
       bare std::vector **member** of a temporary as a pointer into freed
       memory, and a freed vector reports size 0 -- so
       `score_structures(...).pairs[0]` raised IndexError, reading as
       "this structure has no distances", while binding the result first
       worked. IMP_VALUES wraps this type properly. */
    PairDistances pairs;

    ScreenedStructure()
        : score(0), chi2_r(0), sigma1(0), sigma2(0), sigma3(0), invalid_r(0),
          ref_rmsd(0) {}
    ScreenedStructure(const std::string& path, double score)
        : path(path), score(score), chi2_r(score), sigma1(0), sigma2(0),
          sigma3(0), invalid_r(0), ref_rmsd(0) {}

    IMP_SHOWABLE_INLINE(ScreenedStructure,
                        out << "ScreenedStructure(" << path << ", " << score
                            << ")");
};
IMP_VALUES(ScreenedStructure, ScreenedStructures);

//! An assembled model: rigid bodies, the FRET network, and the clash term.
/*! What `score`, `dock` and `refine` all start from. The value holds the
    model alive -- an `IMP::Model` is reference counted and everything else
    here indexes into it -- so a caller keeps the assembly, not its parts. */
class IMPBFFEXPORT DockingAssembly {
    IMP::Pointer<IMP::Model> model_;
    IMP::atom::Hierarchy root_;
    IMP::Pointer<ProbeNetworkRestraint> network_;
    IMP::Pointer<IMP::RestraintSet> restraints_;
    IMP::Pointer<IMP::core::RestraintsScoringFunction> scoring_function_;
    IMP::core::RigidBodies bodies_;
    std::vector<int> body_of_pdb_;
    bool mean_position_;
    double sigma_da_;
    std::string fps_json_path_;
    double max_force_ = 0.0;

public:
    DockingAssembly() : mean_position_(true), sigma_da_(6.0) {}
    DockingAssembly(IMP::Model* model, IMP::atom::Hierarchy root,
                    ProbeNetworkRestraint* network, IMP::RestraintSet* restraints,
                    IMP::core::RestraintsScoringFunction* scoring_function,
                    const IMP::core::RigidBodies& bodies,
                    const std::vector<int>& body_of_pdb, bool mean_position,
                    double sigma_da, const std::string& fps_json_path = "",
                    double max_force = 0.0);

    IMP::Model* get_model() const { return model_; }
    IMP::atom::Hierarchy get_root() const { return root_; }
    ProbeNetworkRestraint* get_network() const { return network_; }
    IMP::RestraintSet* get_restraints() const { return restraints_; }
    IMP::core::RestraintsScoringFunction* get_scoring_function() const {
        return scoring_function_;
    }
    IMP::core::RigidBodies get_rigid_bodies() const { return bodies_; }
    //! The `body_id` each input PDB became, in the order they were given.
    std::vector<int> get_body_of_pdb() const { return body_of_pdb_; }
    bool get_mean_position() const { return mean_position_; }
    double get_sigma_da() const { return sigma_da_; }
    //! The fps.json the network actually read.
    /*! Not necessarily what the caller passed: a legacy C# labelling file is
        converted beside itself first, and a run that wants to read the score
        sets back (`OptimizeSelected = All`) needs the file that was used. */
    std::string get_fps_json_path() const { return fps_json_path_; }
    //! The `MaxForce` its mean-position restraints were built with; 0 = uncapped.
    double get_max_force() const { return max_force_; }

    //! The total score of the assembly as it stands.
    double evaluate() const { return scoring_function_->evaluate(false); }

    IMP_SHOWABLE_INLINE(DockingAssembly,
                        out << "DockingAssembly(" << bodies_.size()
                            << " bodies)");
};

//! Assemble a model from PDBs and an fps.json: bodies, FRET, excluded volume.
/*!
    Each PDB becomes one rigid body. Waters and hydrogens are dropped, as they
    are for every AV and clash calculation; nothing else is pruned, because
    IMP.bff resolves labelling positions by body, so a chain id repeated
    across bodies is not ambiguous and a complete molecule (both strands of a
    duplex, say) stays whole.

    \param[in] pdb_paths one or more structures, in `body_id` order
    \param[in] fps_json_path the labelling and distance file
    \param[in] score_set a named score set; empty uses every distance
    \param[in] mean_position_restraint score the separation of the AV mean
               positions rather than rebuilding both volumes every evaluation.
               The volumes are resampled once and attached to their rigid
               bodies, so they move with the structure.
    \param[in] ev_weight the excluded-volume weight
    \param[in] sigma_da the mean-position transfer function's width, Å
    \param[in] clash_tolerance FPS's `ClashTolerance`, Å: the soft-sphere
               constant becomes \f$k = 2/\text{ClashTolerance}^2\f$.
               **Non-positive (the default) keeps IMP's historical \f$k = 1\f$**
               -- this door is the *scoring* door and its numbers are pinned
               (HIV-RT `resolved` at 59.2547, of which 24.8005 is the clash
               term -- the score moved with Olga's radii on 2026-09-01, the
               clash term did not, because `clash_container` reads the
               particles' own radii and not the volume's),
               so FPS's docking constant is opted into rather than imposed. The
               docking protocol asks for it through
               #IMP::bff::DockingParameters::clash_tolerance.
    \param[in] max_force FPS's `MaxForce` for the mean-position restraints;
               0 (the default) is a pure parabola. See
               #IMP::bff::chi2_score_capped.
    \param[in] clash_radii_source `"imp"` (the default, the particles' own
               united-atom radii) or `"olga"` (the Bondi-scale table FPS's
               clash constant was calibrated against). The clash term here is
               all-atom, so this always applies -- there is no coarse bead on
               the scoring door.
    \param[in] clash_radii_scale multiplies whatever that source gives; 1.0 is
               the default and moves no pinned number
    \throw IOException when a file is missing
    \throw ValueException when no PDB is given, or on an unknown radii source

    \note A position of `simulation_type` `XYZ` or `ATOM` becomes a **point**
    rather than a volume, and every distance touching one is scored as
    \f$R_{mp}\f$ whatever its `distance_type` says -- FPS's rule
    (`FilterEngine.cs:305-307`). An `XYZ` point is attached to the rigid body
    its `body_id` names, so it moves with that body; an `ATOM` point *is* an
    atom of the structure and moves for free.

    \note \p fps_json_path may be the original FPS / C# labelling file rather
    than an fps.json; it is converted beside itself once and the converted
    copy is what the network reads.

    \note This does **not** need `IMP.pmi`. The PMI wrapper it replaces was
    the only reason the scoring path pulled a sampler's dependency, and the
    lazy name it was reached by is what made the whole engine unrunnable
    (`okf/log.md`, 2026-08-25).
*/
IMPBFFEXPORT DockingAssembly create_docking_assembly(
        const std::vector<std::string>& pdb_paths,
        const std::string& fps_json_path, const std::string& score_set = "",
        bool mean_position_restraint = true, double ev_weight = 1.0,
        double sigma_da = 6.0, double clash_tolerance = 0.0,
        double max_force = 0.0,
        const std::string& clash_radii_source = "imp",
        double clash_radii_scale = 1.0);

//! FPS's bond rule: the anchor atoms of atom-to-atom restraints stop clashing.
/*!
    `SpringEngine.cs:119-134`. A distance whose two ends are both plain atoms
    is a **bond** -- a crosslink, a covalent tie -- and a bonded pair is *meant*
    to be pulled to a covalent distance, so its two atoms must not repel each
    other with a full van der Waals radius. FPS drops both to
    `AtomData.vdWRNoClash = 0.4 Å` (`StaticData.cs:33`).

    \param[in] network the network whose used distances are tested
    \param[in] radius the radius the anchors are given, Å
    \return how many atoms had their radius changed

    \note **FPS's own defect (D12) is not reproduced.** There, `Molecule` is a
    reference type and the engine's setter copies the *list*, not the
    molecules, so the 0.4 Å is written into the session-global structure and
    leaks into every later clone, mode and run, never reverted -- deselect the
    bond, change mode, re-dock, and the atom still cannot clash until the file
    is reloaded from disk. Here the radius is written on the particles of one
    `IMP::Model`, which #IMP::bff::create_docking_assembly creates per call from
    its own `read_pdb`, so the mutation cannot outlive the assembly.

    \note It has to be applied **after** anything that rewrites radii.
    #IMP::bff::clash_container sets a bead radius on every backbone atom it
    selects, which would undo this for an anchor that happens to be a `CA`.
*/
IMPBFFEXPORT int set_bond_anchor_radii(ProbeNetworkRestraint* network,
                                       double radius = 0.4);

//! Score one assembled model against its FRET restraints.
/*!
    \param[in] assembly what to score
    \param[in] output_csv where to write the per-pair table; empty writes none
    \return the score, the table, and how much of the file was used
*/
IMPBFFEXPORT DockingResult score_assembly(const DockingAssembly& assembly,
                                          const std::string& output_csv = "");

//! The rigid-body reference frames, as a compact JSON array.
/*! One `{"body_id", "t": [x, y, z], "q": [w, x, y, z]}` per body: a
    translation and a rotation quaternion. Re-applied against the same input
    PDBs they reconstruct the pose exactly, which is what lets a docked state
    be stored in a project file and continued. */
IMPBFFEXPORT std::string capture_poses(const DockingAssembly& assembly);

//! Put such a state back onto an assembly built from the same inputs.
/*! Bodies the JSON does not name are left where they are.
    \throw ValueException when \p poses_json is not an array of poses */
IMPBFFEXPORT void apply_poses(const DockingAssembly& assembly,
                              const std::string& poses_json);

//! The used distances **as they are actually scored**, with FPS's point rule.
/*!
    `ProbeNetworkRestraint::get_used_distances` reports what the file says. This
    reports what is scored, and the two differ in exactly one way: a distance
    with an `XYZ` or `ATOM` end comes back with `distance_type` set to `Rmp`,
    because that end is a point and there is no distribution at it for
    \f$\langle R_{DA}\rangle\f$ to average over (`FilterEngine.cs:305-307`).

    Applying the rule to the *measurement* rather than at each of the four
    places a model distance is formed is deliberate: the transfer function, the
    restraint, the CSV's `distance_type` column and the table then all agree by
    construction, and the column says what was done rather than what was asked
    for.

    \param[in] network the network; null gives an empty map
*/
IMPBFFEXPORT std::map<std::string, AVPairDistanceMeasurement>
        scored_measurements(ProbeNetworkRestraint* network);

//! The experiment-against-model table of a scored AV network.
/*!
    \param[in] restraint the network whose used distances are read
    \param[in] mean_position take the model distance from the separation of the
               two AV *mean positions*, converted through the transfer
               function, rather than from the volumes themselves
    \param[in] sigma_da the transfer function's width, Å (mean-position only)
    \return one entry per used distance, in the network's own order

    A pair whose model distance cannot be formed -- a position the network has
    no AV for -- comes back with a NaN model distance rather than being
    dropped: a table with a row missing looks like a table of a smaller
    experiment.
*/
IMPBFFEXPORT std::vector<PairDistance> collect_pair_distances(
        ProbeNetworkRestraint* restraint, bool mean_position = false,
        double sigma_da = 6.0);

//! The same table, at coordinates a caller holds rather than the network's.
/*!
    The minimisation path carries each AV's mean position on a *proxy*
    particle that is a member of the rigid body, so the volumes themselves do
    not move with the structure and reading them would report the pose the run
    started from. This takes the particles that did move.

    \param[in] distances the used distances, keyed as the network keys them
    \param[in] m the model \p particles belong to
    \param[in] position_names,particles one particle per fps.json position,
               in the same order
    \param[in] sigma_da the transfer function's width, Å
    \return one entry per distance; a distance whose positions are not in
            \p position_names gets a NaN model distance
*/
IMPBFFEXPORT std::vector<PairDistance> pair_distances_at_positions(
        const std::map<std::string, AVPairDistanceMeasurement>& distances,
        IMP::Model* m, const std::vector<std::string>& position_names,
        const IMP::ParticleIndexes& particles, double sigma_da = 6.0);

//! Score structures against an fps.json, in one call.
/*! Builds the assembly and scores it -- what a caller wants when it has files
    rather than a model in hand.

    \param[in] pdb_paths one structure per rigid body
    \param[in] fps_json_path the labelling and distance file
    \param[in] score_set a named score set; empty uses every distance
    \param[in] mean_position_restraint score the volumes' mean positions
               rather than rebuilding them (faster, an approximation)
    \param[in] sigma_da the mean-position transfer width, Å
    \param[in] output_csv where to write the per-pair table; empty writes none
    \param[in] clash_radii_source `"imp"` or `"olga"`; see
               #IMP::bff::create_docking_assembly. The default leaves the
               pinned score (HIV-RT `resolved` 59.0404) where it is.
    \param[in] clash_radii_scale multiplies that source's radii; 1.0 default
*/
IMPBFFEXPORT DockingResult score_structures(
        const std::vector<std::string>& pdb_paths,
        const std::string& fps_json_path, const std::string& score_set = "",
        bool mean_position_restraint = false, double sigma_da = 6.0,
        const std::string& output_csv = "",
        const std::string& clash_radii_source = "imp",
        double clash_radii_scale = 1.0);

//! Dock by FRET-restrained energy minimisation (conjugate gradients).
/*!
    The deterministic alternative to the Monte-Carlo sampler (`imp_bff dock`),
    and the FPS approach: each volume's mean position rides on its rigid body
    as a plain point-member *proxy*, the proxies are pulled toward the
    experimental distances by #IMP::bff::AVMeanDistanceRestraint, and a coarse
    excluded-volume term keeps the bodies apart.

    \note The proxies exist for a reason worth keeping written down. Attaching
    the dye particles themselves as body members puts IMP into
    adjoint-derivative mode, which drops the gradients a minimiser needs; a
    plain XYZ point member propagates them the classic way.

    \param[in] pdb_paths one structure per rigid body
    \param[in] fps_json_path the labelling and distance file
    \param[in] output_dir where `docked.pdb`, `scores.csv` and
               `convergence.csv` are written
    \param[in] params the run's control parameters; `n_frames` is the
               iteration budget
    \param[in] stop asked between chunks; null never stops
    \param[in] initial_poses a docked state to resume from, as JSON; empty
               starts from the input pose (and shuffles, if asked to)
    \return the score, the table, the pose and what was written
*/
IMPBFFEXPORT DockingResult dock_minimize(
        const std::vector<std::string>& pdb_paths,
        const std::string& fps_json_path, const std::string& output_dir,
        const DockingParameters& params = DockingParameters(),
        DockingStop* stop = NULL, const std::string& initial_poses = "");

//! Dock by FRET-restrained Monte Carlo.
/*!
    The stochastic alternative to #IMP::bff::dock_minimize, and what
    `imp_bff dock` runs. Each mobile rigid body gets an
    `IMP::core::RigidBodyMover`; the sampler proposes, the network scores, and
    Metropolis at \c mc_temperature decides. Unlike the minimiser this needs no
    point-member proxies -- nothing differentiates the score, so the volumes'
    own mean positions can be scored directly, which is why the assembly is
    built with \c mean_position_restraint.

    **The best pose is what comes back**, not the last one. A Monte-Carlo walk
    ends wherever it happens to be, and at a temperature that samples properly
    that is routinely worse than the best state it passed through; the
    optimiser is told to restore the best, and the returned pose, PDB and table
    are all read from it.

    \param[in] pdb_paths one structure per rigid body
    \param[in] fps_json_path the labelling and distance file
    \param[in] output_dir where `docked.pdb`, `best_*.pdb`, `scores.csv` and
               `convergence.csv` are written
    \param[in] params `n_frames` outer iterations of `mc_steps` each,
               `mc_temperature` the Metropolis kT, `max_translation` and
               `max_rotation` the mover amplitudes, `shuffle_max_translation`
               the initial randomisation (0 does not shuffle), `n_best` how
               many best-scoring models to write, `fixed_body` the body held
               still as the reference frame
    \param[in] stop asked between frames; null never stops
    \param[in] initial_poses a docked state to resume from, as JSON; empty
               starts from the input pose (and shuffles, if asked to)
    \return the best score, its table, its pose and what was written
*/
IMPBFFEXPORT DockingResult dock(const std::vector<std::string>& pdb_paths,
                                const std::string& fps_json_path,
                                const std::string& output_dir,
                                const DockingParameters& params =
                                        DockingParameters(),
                                DockingStop* stop = NULL,
                                const std::string& initial_poses = "");

//! Repeat a docking run from independent random starts and report the spread.
/*!
    The precision an experiment's distances actually pin down, as opposed to
    the score one run happened to reach: each trial shuffles and docks again,
    and what the trials disagree about is what the data does not determine.
    This is the *sampling* spread; #IMP::bff::fps_bootstrap is the other
    question, how much the distance uncertainties move the answer.

    Trials run in sequence, each in `output_dir/trial_000` and so on. The
    returned #IMP::bff::DockingResult is the best trial's, with the spread in
    its `extra` JSON: `n_trials`, `score_mean`, `score_std`, `best_trial`, and
    `trial_scores` in trial order.

    \param[in] pdb_paths,fps_json_path,output_dir as for #dock
    \param[in] params the per-trial parameters; `n_frames` is one trial's
               budget, not the total
    \param[in] n_trials how many independent starts
    \param[in] minimize run #dock_minimize per trial rather than #dock. The
               default: independent random starts followed by a deterministic
               descent is what the spread of a *fit* means.
    \param[in] stop asked between trials; null never stops
    \throw ValueException when \p n_trials is below one
*/
IMPBFFEXPORT DockingResult estimate_docking_errors(
        const std::vector<std::string>& pdb_paths,
        const std::string& fps_json_path, const std::string& output_dir,
        const DockingParameters& params = DockingParameters(), int n_trials = 10,
        bool minimize = true, DockingStop* stop = NULL);

//! Refine a pose in place: minimise, write the structure and the table.
/*! No shuffle and no proxies -- the volumes' own mean positions are scored,
    so this is the local polish after a pose is roughly right.

    \param[in] pdb_paths,fps_json_path,output_dir as for #dock_minimize
    \param[in] score_set a named score set; empty uses every distance
    \param[in] steps the iteration budget
    \param[in] ev_weight the excluded-volume weight */
IMPBFFEXPORT DockingResult refine_docking(const std::vector<std::string>& pdb_paths,
                                  const std::string& fps_json_path,
                                  const std::string& output_dir,
                                  const std::string& score_set = "",
                                  int steps = 500, double ev_weight = 1.0);

//! Score a library of structures and rank them, best first.
/*!
    \param[in] pdb_paths the structures; a directory is expanded to the
               `*.pdb` files directly in it
    \param[in] fps_json_path the labelling and distance file
    \param[in] score_set a named score set; empty uses every distance
    \param[in] output_csv where to write the ranked table; empty writes none
    \param[in] mean_position_restraint score mean positions rather than
               rebuilding both volumes for every structure
    \return one entry per structure, ascending by score, unscorable ones
            (NaN) last

    Each entry carries FPS's screening diagnostics as well as the score: the
    reduced \f$\chi^2\f$, the 1σ/2σ/3σ violation counts, the number of
    distances with no model value, the reference-frame fit and the per-pair
    table (#IMP::bff::ScreenedStructure).

    A structure that cannot be scored does not stop the library: it is
    recorded with a NaN score and the run goes on. Read the NaNs -- a table
    that is *all* NaN is not a library of bad models, it is a broken run.
*/
IMPBFFEXPORT std::vector<ScreenedStructure> screen_structures(
        const std::vector<std::string>& pdb_paths,
        const std::string& fps_json_path, const std::string& score_set = "",
        const std::string& output_csv = "",
        bool mean_position_restraint = false);

//! Write the score and the per-pair table as CSV.
/*! One `total_score` row, then a header and one row per pair. */
IMPBFFEXPORT void write_score_csv(const std::string& path, double score,
                                  const std::vector<PairDistance>& pairs);

//! Every screened structure's pair table, stacked into one CSV.
/*! One row per (structure, distance), the structure's path in the first
    column: what says *which* restraint the runners-up of a screen miss, in a
    shape a spreadsheet or a dataframe can group by.
    \throw IOException when \p path cannot be written */
IMPBFFEXPORT void write_screening_pairs_csv(
        const std::string& path,
        const std::vector<ScreenedStructure>& structures);

/* ------------------------------------------------------------------------
 * Error estimation -- the parametric bootstrap (PRD-121 G2)
 * ------------------------------------------------------------------------ */

//! How a replica's synthetic distance is drawn around the truth.
/*!
    Both models put the **mode** at the target and give the two sides the
    asymmetric errors `error_neg` / `error_pos`. They differ in how much
    probability each side gets, and that is not a detail: it decides whether
    the density is continuous and how far the ensemble drifts.

    For \f$\sigma_+ = 10\f$, \f$\sigma_- = 3\f$ (an ordinary FRET error bar):

    | | #BFF_SPLIT_NORMAL | #FPS_SIGN_SPLIT_NORMAL |
    |---|---|---|
    | \f$P(X > 0)\f$ | \f$\sigma_+/(\sigma_++\sigma_-)\f$ = 0.7692 | 0.5 |
    | density at \f$0^\pm\f$ | equal -- **continuous** | ratio \f$\sigma_-/\sigma_+\f$ = 0.3 -- **jumps** |
    | mean | \f$2(\sigma_+-\sigma_-)/\sqrt{2\pi}\f$ = **+5.585** | \f$(\sigma_+-\sigma_-)/\sqrt{2\pi}\f$ = **+2.792** |
    | sd | \f$\sqrt{\sigma_+^2-\sigma_+\sigma_-+\sigma_-^2 - \mu^2}\f$ = 6.914 | 6.836 |

    Read the mean row before choosing. Fixing FPS's discontinuity does **not**
    remove the outward drift -- it doubles it, because giving the wide side its
    proper share of the mass moves more draws onto the wide side. Neither model
    is unbiased; a perturbation that is unbiased *and* respects two different
    error bars does not exist, because the mode and the mean of an asymmetric
    density are different points. What #BFF_SPLIT_NORMAL buys is that it is a
    real density with the stated quantiles, so "1 sigma" means the same thing
    on both sides of it.
*/
enum FPSPerturbationModel {
    //! The two-piece (split) normal: continuous, half-masses \f$\sigma_\pm/(\sigma_++\sigma_-)\f$.
    /*! \f$f(x) \propto \phi(x/\sigma_+)\f$ for \f$x>0\f$ and
        \f$\phi(x/\sigma_-)\f$ for \f$x<0\f$, normalised once over both halves.
        The default, because it is the object the two error bars describe. */
    BFF_SPLIT_NORMAL = 0,
    //! FPS's sign-split normal (`ErrorEstimation.cs:57`), each half mass 1/2.
    /*! `z ~ N(0,1)`, scaled by `ErrPlus` when `z > 0` and by `ErrMinus` when
        `z < 0`. Not a normalised density with those two scales -- it is two
        half-normals glued with equal weight, so the density jumps by
        \f$\sigma_-/\sigma_+\f$ at the target. Select it to reproduce an FPS
        error estimate; do not select it to make one. */
    FPS_SIGN_SPLIT_NORMAL = 1
};

//! Draw \p n offsets \f$R' - R\f$ from a perturbation model.
/*! The engine's own draw, exposed so its distribution can be measured rather
    than asserted. \p seed selects the stream; the same seed gives the same
    draws on every platform (`std::mt19937`).
    \param[in] error_neg,error_pos the two error bars, Å; non-positive is 0
    \param[in] n how many draws
    \param[in] seed the `std::mt19937` seed
    \param[in] model which of the two densities */
IMPBFFEXPORT std::vector<double> sample_distance_perturbations(
        double error_neg, double error_pos, int n, unsigned int seed,
        FPSPerturbationModel model = BFF_SPLIT_NORMAL);

//! The closed-form mean of a perturbation offset, Å. See #FPSPerturbationModel.
IMPBFFEXPORT double perturbation_mean(double error_neg, double error_pos,
                                      FPSPerturbationModel model);
//! The closed-form standard deviation of a perturbation offset, Å.
IMPBFFEXPORT double perturbation_sd(double error_neg, double error_pos,
                                    FPSPerturbationModel model);
//! The closed-form \f$P(R' > R)\f$ of a perturbation offset.
IMPBFFEXPORT double perturbation_upper_mass(double error_neg, double error_pos,
                                            FPSPerturbationModel model);

//! FPS's shipped Error-estimation parameters (`ProjectData.cs:101-115`).
/*! `MaxForce = 10000` (so the Huber tail is effectively off), `ClashTolerance
    = 0.5` Å (four times harder per angstrom than docking) and
    `OptimizeSelected = All`. The iteration budget is *not* FPS's 100 000: that
    counts damped-Verlet steps of an integrator this module does not have, and
    #IMP::bff::DockingParameters::n_frames counts conjugate-gradient steps. It
    is left at the shipped default and is the one number a caller should set.

    Two settings are **not** FPS's.
    #IMP::bff::DockingParameters::shuffle_max_translation is 0 because this mode
    re-optimises from the parent pose and never restarts randomly
    (`SimulationJobManager.cs:96-104`), and
    #IMP::bff::DockingParameters::coarse_clash is **off** because FPS's clash is
    all-atom over real van der Waals radii; the coarse one-bead-per-residue term
    is a docking-speed approximation, and a mode whose pose moves by fractions
    of an angstrom has nothing to buy with it.

    \warning **Read #BootstrapResult::parent_e_clash before believing the
    spread.** With `ev_weight = 1` these are *clash-dominated* parameters:
    \f$k_{clash} = 2/0.5^2 = 8\f$, and IMP's united-atom radii are larger than
    the Bondi set FPS uses, so a genuine protein--DNA interface registers as a
    large overlap. Measured on HIV-RT `resolved`, parent = the docked pose, 10
    replicas, 300 steps:

    | `ClashTolerance` | `ev_weight` | parent score (clash) | RMSD vs parent |
    |---|---|---|---|
    | 0.5 | 1.0 | 134.80 (102.92) | **0.000 ± 0.000 Å** |
    | 1.0 | 1.0 | 59.72 (27.89) | **0.000 ± 0.000 Å** |
    | 1.0 | 0.1 | 34.58 (2.97) | 0.056 ± 0.127 Å |
    | 1.0 | 0.0 | 14.83 (0) | **1.935 ± 1.133 Å** |

    When the clash term carries most of the score the pose is held by geometry,
    not by the FRET data, and the bootstrap correctly reports that the *data*
    constrain nothing further -- but a reader who takes 0.000 Å as "the FRET
    network locates this subunit to a thousandth of an angstrom" has read it
    backwards. The number that says which regime a run is in is the clash
    fraction of the parent score. */
IMPBFFEXPORT DockingParameters fps_error_estimation_parameters();

//! What an error-estimation run is told beyond the docking parameters.
struct IMPBFFEXPORT BootstrapParameters {
    //! Replicas per parent. FPS's repetitions spinner defaults to 10.
    int n_replicas;
    //! The `std::mt19937` seed; the run is reproducible from it.
    unsigned int seed;
    //! Which density the synthetic distances are drawn from.
    FPSPerturbationModel perturbation;
    //! Give every scored distance noise, not only the score set's.
    /*!
        **FPS's defect D2, as a switch.** `ErrorEstimation.SetState` rewrites
        the target of *all* distances with the parent's model distance
        (`:23-29`) but `PerturbDistances` skips the deselected ones (`:54`) --
        and error estimation ships `OptimizeSelected = All`
        (`ProjectData.cs:113`). A deselected distance therefore enters the
        re-optimisation as a **zero-noise restraint already satisfied exactly**
        at the parent pose: a pin, holding the replica where it started and
        biasing the reported uncertainty **downward**.

        `true` (the default) perturbs every distance the run will score, which
        is what a parametric bootstrap means: every datum that constrains the
        fit must be resampled. `false` reproduces FPS, and is only useful for
        comparing against one of its numbers.

        Invisible when the score set is everything -- there are then no
        deselected distances to pin with.
    */
    bool perturb_deselected;
    //! Keep each replica's output directory instead of only its numbers.
    bool keep_replica_dirs;

    BootstrapParameters()
        : n_replicas(10), seed(1), perturbation(BFF_SPLIT_NORMAL),
          perturb_deselected(true), keep_replica_dirs(false) {}

    IMP_SHOWABLE_INLINE(BootstrapParameters,
                        out << "BootstrapParameters(" << n_replicas
                            << " replicas, seed " << seed << ")");
};
IMP_VALUES(BootstrapParameters, BootstrapParametersList);

//! One error-estimation replica: a synthetic dataset, re-fitted.
struct IMPBFFEXPORT BootstrapReplica {
    int replica;
    //! The score the replica reached **against its own perturbed targets**.
    double score;
    //! The bond subset of #score, and the excluded-volume term beside it.
    double e_bond, e_clash;
    bool converged;
    //! \f$\chi^2\f$ per distance of the re-fitted pose against the *truth*.
    /*! The unperturbed parent distances, i.e. what the replica gave up by
        chasing noise. Roughly 1 per distance when the error model is honest
        and the fit is not over-determined. */
    double chi2_r_truth;
    //! RMSD of the replica against the parent, Å -- **the reported statistic**.
    /*! All atoms of all bodies, unweighted, no superposition: the bodies are
        already in a common frame, and superposing would remove exactly the
        rigid-body displacement the estimate is about. */
    double rmsd_to_parent;
    //! The same number under FPS's sign error (see #IMP::bff::pose_rmsd).
    double rmsd_to_parent_fps;
    //! The replica's pose, as #IMP::bff::capture_poses writes it.
    std::string poses;
    //! The replica's experiment-against-model table; `distance_exp` is the
    //! *perturbed* target it was fitted to, not the experiment.
    PairDistances pairs;
    //! Its own directory, when #BootstrapParameters::keep_replica_dirs is set.
    std::string output_dir;

    BootstrapReplica()
        : replica(0), score(0), e_bond(0), e_clash(0), converged(false),
          chi2_r_truth(0), rmsd_to_parent(0), rmsd_to_parent_fps(0) {}

    IMP_SHOWABLE_INLINE(BootstrapReplica,
                        out << "BootstrapReplica(" << replica << ", rmsd "
                            << rmsd_to_parent << ")");
};
IMP_VALUES(BootstrapReplica, BootstrapReplicas);

//! What a parametric bootstrap reports.
struct IMPBFFEXPORT BootstrapResult {
    //! The parent's pose -- the "truth" every replica was drawn from -- and
    //! its score against the **experimental** distances.
    /*! Stated the same way however the parent arrived (docked here, or handed
        in): the mean-position restraint sum over every distance in the file
        plus the excluded-volume term, in this module's half-\f$\chi^2\f$
        units. That makes it comparable with #BootstrapReplica::score, which is
        the same objective against perturbed targets. */
    double parent_score;
    std::string parent_poses;
    //! The parent's own model distances, i.e. the targets with residuals zeroed.
    /*! `distance_exp` is the experiment as the file states it and
        `distance_model` is what the parent pose produces; the second column is
        what the replicas were drawn around. */
    PairDistances truth;
    BootstrapReplicas replicas;
    //! The excluded-volume part of #parent_score.
    /*! **The first number to read.** When it is most of #parent_score the
        replicas are held by the clash term rather than by the data, and the
        spread below collapses toward zero however good the FRET network is --
        measured on HIV-RT at 103 of 135, spread 0.000 Å, against 0 of 15 and
        1.935 Å with the excluded volume switched off. */
    double parent_e_clash;
    //! Mean, sample standard deviation and maximum of #BootstrapReplica::rmsd_to_parent.
    double rmsd_mean, rmsd_sd, rmsd_max;
    //! How many distances carried noise, and how many were zero-noise pins.
    /*! #n_pinned is non-zero only under FPS's D2
        (#BootstrapParameters::perturb_deselected false with a score set that
        leaves distances out), and every one of them drags the estimate down. */
    int n_perturbed, n_pinned;
    //! The settings and the closed-form perturbation moments, as JSON.
    std::string extra;

    BootstrapResult()
        : parent_score(0), parent_e_clash(0), rmsd_mean(0), rmsd_sd(0),
          rmsd_max(0),
          n_perturbed(0), n_pinned(0), extra("{}") {}

    //! The whole result as JSON, replicas and their numbers included.
    std::string get_json() const;

    IMP_SHOWABLE_INLINE(BootstrapResult,
                        out << "BootstrapResult(" << replicas.size()
                            << " replicas, rmsd " << rmsd_mean << " +- "
                            << rmsd_sd << " A)");
};
IMP_VALUES(BootstrapResult, BootstrapResults);

//! FPS's error estimation: perturb the model's own distances and re-fit.
/*!
    A **parametric bootstrap**, and not the same quantity as
    `imp_bff dock-errors`. That command re-docks from independent random
    starts and measures how reproducible the *optimiser* is; this one holds the
    optimiser fixed, resamples the *data* from the stated error model, and
    measures how much the answer moves. One is a property of the search, the
    other of the measurement.

    Per replica, following `ErrorEstimation.cs` and
    `SimulationJobManager.cs:96-104`:

    1. the parent's own model distances become the targets -- **residuals
       zeroed**, so the experiment is discarded for the run and the "truth" is
       the docked model itself;
    2. each target is perturbed once by #sample_distance_perturbations;
    3. the structure is re-minimised **from the parent pose**, not from a fresh
       random start -- there is no shuffle in this mode;
    4. the targets are reset before the next replica, so the draws are i.i.d.
       and do not compound.

    What a caller reads off is the spread of
    #BootstrapReplica::rmsd_to_parent. It is a **precision** figure conditional
    on the restraint set and the error model, with no model-misspecification
    term in it -- because step 1 threw the residuals away.

    \param[in] pdb_paths one structure per rigid body
    \param[in] fps_json_path the labelling and distance file
    \param[in] output_dir where the parent and the replicas are written
    \param[in] parent_poses the docked state to bootstrap around, as
               #IMP::bff::capture_poses writes it; **empty docks first**, into
               `<output_dir>/parent`, and uses that
    \param[in] params the docking parameters each replica is re-minimised
               under -- #fps_error_estimation_parameters is FPS's own set
    \param[in] bootstrap the replica count, the seed and the two decisions
               above
    \throw IOException when an input is missing

    \note Each replica is a full #IMP::bff::dock_minimize against a perturbed
    copy of the file, written beside it. That rebuilds the assembly and
    resamples the volumes per replica, which is the expensive way; it is done
    that way so a replica goes through exactly the code path a docking run
    does, and cannot drift from it.
*/
IMPBFFEXPORT BootstrapResult fps_bootstrap(
        const std::vector<std::string>& pdb_paths,
        const std::string& fps_json_path, const std::string& output_dir,
        const std::string& parent_poses = "",
        const DockingParameters& params = DockingParameters(),
        const BootstrapParameters& bootstrap = BootstrapParameters());

//! RMSD between two poses of one assembly, over every atom, unweighted.
/*!
    No superposition: two poses of the same rigid bodies are already in one
    frame, and the displacement between them is what an error estimate is
    about. With \f$U = R_a - R_b\f$ (a difference of rotation matrices, not a
    rotation) and \f$t = t_a - t_b\f$, an atom at body-local \f$r\f$ moves by
    \f$Ur + t\f$, so
    \f$\mathrm{RMSD} = \sqrt{\frac{1}{N}\sum_i |U r_i + t|^2}\f$.

    \param[in] assembly the bodies and their atoms
    \param[in] poses_a,poses_b the two states; an empty string is the
               assembly's current pose
    \param[in] fps_sign_convention accumulate \f$|U r_i - t|^2\f$ instead

    \note **FPS's sign error, as a switch** (`SimulationResult.cs:40-72`). FPS
    accumulates `SquareNormDiff(U*r, t)` = \f$|Ur - t|^2\f$ where the
    displacement it derived is \f$Ur + t\f$, so every RMSD it has ever printed
    has the translation entering with the wrong sign. The two differ by
    \f$\frac{4}{N}\sum_i (U r_i)\cdot t\f$, which vanishes exactly when the
    body-local coordinates sum to zero -- so the size of the error is set by
    how far each body's frame origin is from the unweighted centroid of its
    atoms. FPS's origin is the **mass-weighted** centre of mass, so the sum
    does not vanish and the error is real but small; a port whose origin
    happened to be the plain centroid would find the two identical and
    conclude, wrongly, that there was nothing there.
*/
IMPBFFEXPORT double pose_rmsd(const DockingAssembly& assembly,
                              const std::string& poses_a,
                              const std::string& poses_b,
                              bool fps_sign_convention = false);

//! The rigid transform putting one pose of an assembly onto another.
/*!
    A Kabsch superposition over every atom, unweighted, **computed here rather
    than remembered**. FPS's `SaveForm` reads `SimulationResult.BestFitRotation`
    and no `Save*` method ever computes it: it is a side effect of the GUI's
    RMSD column (`SimulationResult.cs:47` called from
    `MainForm.DisplayNewStructures`), so an exported overlay is a chain of
    *pairwise* fits under the default "RMSD vs previous", and an all-zero
    matrix -- which `AngleAndAxis` turns into a spurious 180 degree rotation --
    if the best-fit box was unticked when the grid last refreshed. That is not
    reproduced. #IMP::bff::write_fps_overlay and friends call this against the
    reference they were given.

    \param[in] assembly the bodies and their atoms
    \param[in] poses the state to move
    \param[in] reference_poses the state to move it onto
    \return the transformation taking the first onto the second
*/
IMPBFFEXPORT IMP::algebra::Transformation3D pose_superposition(
        const DockingAssembly& assembly, const std::string& poses,
        const std::string& reference_poses);

IMPBFF_END_NAMESPACE

#endif //IMPBFF_DOCKING_H
