/**
 *  \file IMP/bff/QuenchingModel.h
 *  \brief A labelled site's donor decay, two ways.
 *
 * The model QuEst was built on: a dye sphere diffuses on an accessible-volume
 * grid, is slowed near sticky residues, and is quenched by photo-induced
 * electron transfer whenever it comes within contact distance of a redox-active
 * side chain. The observable is the donor's fluorescence decay, and comparing
 * it with a measured lifetime is what calibrates the accessible *contact*
 * volume.
 *
 * Both classes here take an accessible volume and the structure around it and
 * answer what the donor emits. They differ in *how the dye moves*:
 *
 * - #IMP::bff::DynamicAccessibleVolume is the **field** picture. It stamps a
 *   mobility field, a PET field and (with an acceptor) a FRET field onto the
 *   volume's grid, then integrates the excited-state population on it. No shot
 *   noise, and the whole distribution at once.
 * - #IMP::bff::QuenchedDonorDecay is the **particle** picture. It runs a
 *   Brownian walk through the volume, reads the quenching rate along it, and
 *   races photons against that rate. It resolves the *time dependence* of the
 *   rate, which a field average cannot.
 *
 * The two are not alternatives to pick by taste: the field picture's stationary
 * state and the walk's occupancy are the same distribution, and where they
 * disagree the walk is resolving something the field has averaged away.
 *
 * Not to be confused with the explicit all-atom dye under a force field. That
 * is a different model at a different scale; this one is a sphere on a grid,
 * cheap enough to scan every position in a protein.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_QUENCHINGMODEL_H
#define IMPBFF_QUENCHINGMODEL_H

#include <IMP/bff/bff_config.h>
#include <IMP/bff/AVModel.h>
#include <IMP/bff/ProbeDiffusion.h>
#include <IMP/bff/GridDiffusionSolver.h>
#include <IMP/bff/LifetimeSpectrum.h>
#include <IMP/bff/PETQuenching.h>

#include <IMP/bff/Base.h>

#include <limits>
#include <map>
#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! The obstacles a quenching model needs of a structure, named.
/*!
    Parallel arrays rather than a structured dtype: a structured dtype is a
    numpy idea, and the fields it carries are what this wants handed to it.
*/
struct IMPBFFEXPORT ObstacleAtoms {
    std::vector<std::string> chains, res_names, atom_names;
    std::vector<int> res_ids;
    //! Flat, three per atom.
    std::vector<double> coords;

    unsigned int size() const {
        return static_cast<unsigned int>(res_names.size());
    }
    void get_coords(double** out_view, int* n_out_view) const;

    IMP_SHOWABLE_INLINE(ObstacleAtoms,
                        out << "ObstacleAtoms(" << size() << " atoms)");
};
IMP_VALUES(ObstacleAtoms, ObstacleAtomsList);

//! An accessible volume with diffusion, quenching and FRET rate fields.
/*!
    The **field** picture of a quenched dye. Each field is built on first use and
    kept: they are grid-sized, and none of them changes unless a parameter does.
*/
class IMPBFFEXPORT DynamicAccessibleVolume {
    AccessibleVolume av_;
    ObstacleAtoms atoms_;
    double tau0_, probe_radius_, free_diffusion_, contact_distance_, slow_factor_;
    std::string flux_form_;

    mutable std::vector<double> diffusion_map_, quenching_rate_map_;
    mutable std::vector<double> fret_rate_map_, occupancy_;

    //! The domain mask for the solver: where the dye may be.
    std::vector<double> bounds() const;
    //! The step the solver runs at: half the explicit limit when unset.
    double resolve_t_step(double t_step) const;

public:
    //! \param[in] av the volume to decorate
    /*! \param[in] atoms the obstacles
        \param[in] tau0 unquenched donor lifetime, ns
        \param[in] probe_radius A
        \param[in] free_diffusion unhindered diffusion coefficient, A^2/ns
        \param[in] contact_distance dye-to-atom distance counted as contact, for
                   the mobility field
        \param[in] slow_factor mobility multiplier applied per contacting atom
        \param[in] flux_form `"smoluchowski"` or `"ito"` — see
                   #IMP::bff::equilibrium_occupancy
        \throw ValueException for any other \p flux_form */
    DynamicAccessibleVolume(const AccessibleVolume& av = AccessibleVolume(),
                            const ObstacleAtoms& atoms = ObstacleAtoms(),
                            double tau0 = 4.0, double probe_radius = 3.5,
                            double free_diffusion = 8.0,
                            double contact_distance = 6.5,
                            double slow_factor = 0.985,
                            const std::string& flux_form = "smoluchowski");

    int get_ng() const { return av_.get_ng(); }
    double get_dg() const { return av_.get_grid_step(); }
    double get_tau0() const { return tau0_; }
    std::string get_flux_form() const { return flux_form_; }

    //! Binary occupancy of the accessible volume, flat.
    void get_density(double** out_view, int* n_out_view) const;
    //! The grid anchor — the attachment point.
    void get_x0(double** out_view, int* n_out_view) const;
    //! Where the dye may be: the domain mask for the solver.
    void get_bounds(double** out_view, int* n_out_view) const;

    //! The mobility field: free diffusion, slowed near atoms. Built on demand.
    void get_diffusion_map(double** out_view, int* n_out_view) const;
    //! The PET field. \throw ValueException before update_quenching_map()
    void get_quenching_rate_map(double** out_view, int* n_out_view) const;
    //! The FRET field; a zero-length view until an acceptor is set.
    void get_fret_rate_map(double** out_view, int* n_out_view) const;
    //! Quenching, plus FRET where there is an acceptor.
    void get_rate_map(double** out_view, int* n_out_view) const;
    //! The equilibrium distribution of the dye, normalised. Built on demand.
    /*! Under the default `"smoluchowski"` this **is** the AV density, uniform
        over the accessible voxels, and does not depend on the mobility. */
    void get_occupancy(double** out_view, int* n_out_view) const;

    //! Build the mobility field.
    /*! \param[in] base an optional per-voxel base coefficient replacing the
               constant `free_diffusion` — a radial profile, say. Empty takes
               the constant. */
    void update_diffusion_map(
            const std::vector<double>& base = std::vector<double>());

    //! Build the PET field from the dye's pair parameters.
    /*! \param[in] quencher `{comp_id: PETParameters}` for one dye
        \param[in] rC overrides every per-atom attenuation length with one
                   electron-transfer length, which is how ChiSurf drove it;
                   NaN keeps the per-atom values */
    void update_quenching_map(
            const std::map<std::string, PETParameters>& quencher,
            double rC = std::numeric_limits<double>::quiet_NaN());

    //! Build the FRET field against an acceptor's volume.
    /*! The donor's radiative rate is \f$1/\tau_0\f$. ChiSurf read it from the
        `foerster_radius` keyword by a copy-paste slip, so passing a Förster
        radius also set the radiative rate to it. */
    void update_fret_map(const DynamicAccessibleVolume& acceptor,
                         double forster_radius = 52.0, int acceptor_step = 2);

    //! The equilibrium occupancy, in closed form.
    /*! Propagating to it instead is possible but slow and, in the `"ito"` form
        on a real site where the compounding slow factor makes `D` span orders
        of magnitude, may not converge at all: on T4L site 132 it was still
        drifting after 40 000 iterations. */
    void update_occupancy();

    //! The equilibrium occupancy the long way: propagate until it stops moving.
    void update_occupancy_by_iteration(double t_step = -1.0,
                                       int n_steps = 20000,
                                       double tolerance = 1e-8,
                                       int n_check = 100);

    //! Integrate the donor decay on the grid, from the equilibrium start.
    /*!
        \param[in] t_step negative takes half the explicit stability limit —
                   stable with room for the map to change
        \return a #GridDiffusionResult whose `time` axis is built from the
                *resolved* step actually integrated, so no caller has to parse a
                step out of a flat buffer.
    */
    GridDiffusionResult donor_decay(double t_max, double t_step,
                                    int n_out) const;

    IMP_SHOWABLE_INLINE(DynamicAccessibleVolume,
                        out << "DynamicAccessibleVolume(ng=" << av_.get_ng()
                            << ", tau0=" << tau0_ << ")");
};
IMP_VALUES(DynamicAccessibleVolume, DynamicAccessibleVolumes);

//! A labelled site's donor decay, from the structure and the PET chemistry.
/*!
    The **particle** picture: a Brownian walk through the accessible volume, the
    quenching rate read along it, and photons raced against that rate.
*/
class IMPBFFEXPORT QuenchedDonorDecay {
    AccessibleVolume av_;
    ObstacleAtoms atoms_;
    double tau0_, critical_distance_, slow_radius_, probe_radius_;
    double diffusion_coefficient_, slow_fact_, t_step_, t_max_;
    int n_photons_, n_trajectories_, random_seed_;
    std::map<std::string, ResidueQuenching> table_;

    mutable ResidueSites sites_;
    mutable std::vector<double> quenching_rate_map_, slow_factor_map_;
    mutable ProbeDiffusionSimulation walk_;
    mutable std::vector<double> delays_;
    mutable std::vector<int> emitted_;
    mutable bool has_sites_, has_grids_, has_walk_, has_photons_;
    //! Frames, mean rate and collision fraction, from whichever path ran.
    mutable double n_frames_, mean_k_quench_, collision_fraction_;

    std::vector<int> density() const;
    std::vector<double> x0() const;
    void ensure_grids() const;
    void ensure_walk() const;
    void ensure_photons() const;

public:
    //! \param[in] av the donor's accessible volume
    /*! \param[in] atoms the obstacles
        \param[in] tau0 unquenched donor lifetime, ns
        \param[in] quenching_table per-residue interaction table; the defaults of
                   #IMP::bff::amino_acid_quenching_defaults when empty
        \param[in] critical_distance contact radius inherited by residues whose
                   table entry leaves `quench_radius` unset
        \param[in] slow_radius radius of each residue's sticky sphere
        \param[in] n_photons excitation events for the photon Monte-Carlo
        \param[in] n_trajectories negative takes
                   #IMP::bff::default_trajectory_count
        \param[in] random_seed negative draws freely */
    QuenchedDonorDecay(
            const AccessibleVolume& av = AccessibleVolume(),
            const ObstacleAtoms& atoms = ObstacleAtoms(), double tau0 = 4.0,
            const std::map<std::string, ResidueQuenching>& quenching_table =
                    std::map<std::string, ResidueQuenching>(),
            double critical_distance = 7.0, double slow_radius = 10.0,
            double probe_radius = DEFAULT_PROBE_RADIUS,
            double diffusion_coefficient = 40.0, double slow_fact = 0.05,
            double t_step = 0.004, double t_max = 10000.0,
            int n_photons = 100000, int n_trajectories = -1,
            int random_seed = -1);

    double get_tau0() const { return tau0_; }
    double get_dg() const { return av_.get_grid_step(); }
    int get_n_photons() const { return n_photons_; }
    double get_t_step() const { return t_step_; }
    std::map<std::string, ResidueQuenching> get_quenching_table() const {
        return table_;
    }
    void get_density(double** out_view, int* n_out_view) const;
    void get_x0(double** out_view, int* n_out_view) const;

    //! The slow and quench centres of every residue in the structure.
    const ResidueSites& get_sites() const;

    //! Stamp the quenching-rate and stickiness fields onto the AV grid.
    void update_grids();
    void get_quenching_rate_map(double** out_view, int* n_out_view) const;
    void get_slow_factor_map(double** out_view, int* n_out_view) const;
    //! Replace the stickiness field directly (what update_grids() stamps).
    void set_slow_factor_map(const std::vector<double>& m);

    //! The photon seed, distinct from the walk's so the draws are not a replay.
    int get_photon_seed() const;

    //! Whether a walk has been run -- by either path.
    /*! The fused path builds the walk record without a trajectory, so this is
        true after it too: what it reports is that the simulation has happened,
        not that coordinates exist. */
    bool get_has_walk() const { return has_walk_; }

    //! Run the Brownian walk. \return whether a trajectory was produced
    bool simulate_diffusion();
    //! The walk, run on first use.
    const ProbeDiffusionSimulation& get_diffusion() const;
    //! The quenching rate the dye sees, frame by frame, 1/ns.
    void get_k_quench(double** out_view, int* n_out_view) const;

    //! Race photons against the quenching rate along the trajectory.
    /*! \param[in] k_quench the rate trace to race against; empty takes the
               walk's own */
    void simulate_photons(
            const std::vector<double>& k_quench = std::vector<double>());

    //! The walk, the rate along it, and the photon race — in one call.
    /*!
        Identical to simulate_diffusion() followed by simulate_photons(), photon
        for photon. What it skips is the **trajectory**: the split path builds an
        `(n_frames, 3)` array, reads a rate map along it and keeps the trace, and
        at the default `t_max` that is 20 million doubles to produce a few
        thousand photons. Nothing downstream of the decay wants the coordinates.

        Use the split calls when the trajectory *is* the point — a correlation
        function, a visualisation, an inspection of the rate trace. Use this one
        inside a fitting loop.
    */
    void photons_fused();

    //! Sampled delay times, ns. Zero where the excitation was quenched.
    void get_delays(double** out_view, int* n_out_view) const;
    //! Per event: 1 when a photon was emitted, 0 when the dye was quenched.
    void get_emitted(int** out_view_i, int* n_out_view_i) const;

    //! Frames the walk produced.
    /*! Reported by the fused path too, which never materialises them. */
    int get_n_frames() const;
    //! Mean quenching rate along the walk, 1/ns.
    double get_mean_k_quench() const;
    //! Fraction of frames inside a quenching sphere.
    double get_collision_fraction() const;

    //! Emitted photons over excitation events.
    double get_quantum_yield() const;
    //! The species-averaged lifetime of the emitted photons, ns.
    double get_fluorescence_lifetime() const;

    //! The decay as (amplitude, rate) pairs — the neutral output.
    /*!
        Prefer this to decay_histogram(). A spectrum carries no bin width and no
        time range, so whatever owns the instrument can convolve, bin and add
        noise on its own terms.

        Built from the *per-frame* total rate along the trajectory: each frame is
        a state the dye occupies, weighted equally because the walk visits them
        in proportion to their occupancy. That makes this the **static
        approximation** — exact only if the dye held each position for a whole
        excited-state lifetime, which is precisely what a diffusion simulation
        exists to deny. The spectrum is marked `exact=false` for that reason, and
        decay_histogram() is the one that resolves the averaging.

        \param[in] n_species coarse-grain to at most this many species; the
                   trajectory has one state per frame — often 10^5 — and almost
                   none are distinguishable. Zero keeps every frame.
        \throw ValueException when the walk produced no rate trace
    */
    LifetimeSpectrum lifetime_spectrum(int n_species = 128) const;

    //! Histogram of the photons that were actually emitted.
    /*!
        A binned curve is **not** the neutral output — it has chosen a bin width
        and a time range, both of which are instrument settings. Use
        lifetime_spectrum() unless the binning is the point.

        **Only emitted photons.** A quenched excitation comes back with
        `dt = 0`, so histogramming the whole trace piles every non-emitted photon
        into the first bin — a spike of photons that never existed, and a curve
        whose total is the excitation count rather than `QY * N`.

        \param[out] out_view,n_out_view `n_bins + 1` edges, then `n_bins` counts
    */
    void decay_histogram(int n_bins, double t_min, double t_max,
                         double** out_view, int* n_out_view) const;

    //! Per-frame FRET rate against an acceptor whose walk has also been run.
    /*! Both dyes are resolved in time and the frames are paired; the pairing
        uses the shorter of the two walks. */
    void fret_rate_trace_paired(const QuenchedDonorDecay& acceptor,
                                double forster_radius, double kappa2,
                                double r_min, double** out_view,
                                int* n_out_view) const;

    //! Per-frame FRET rate against a static acceptor cloud.
    /*! The fast-acceptor limit — see #IMP::bff::fret_rate_trace.
        \param[in] points flat, three per point */
    void fret_rate_trace_cloud(const std::vector<double>& points,
                               double forster_radius, double kappa2,
                               double r_min, double** out_view,
                               int* n_out_view) const;

    //! \f$1 - QY_{DA}/QY_D\f$, both from the photon Monte-Carlo.
    /*!
        Taking the ratio of two simulated quantum yields rather than an analytic
        formula keeps the quenching in: the donor is quenched by PET in both
        terms, so what is left is the FRET.

        \param[in] fret_rates the per-frame FRET rate, from one of the two trace
                   methods above
    */
    double fret_efficiency(const std::vector<double>& fret_rates) const;

    IMP_SHOWABLE_INLINE(QuenchedDonorDecay,
                        out << "QuenchedDonorDecay(tau0=" << tau0_
                            << ", n_photons=" << n_photons_ << ")");
};
IMP_VALUES(QuenchedDonorDecay, QuenchedDonorDecays);

IMPBFF_END_NAMESPACE

#endif //IMPBFF_QUENCHINGMODEL_H
