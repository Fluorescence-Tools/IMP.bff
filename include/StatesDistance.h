/**
 *  \file IMP/bff/StatesDistance.h
 *  \brief Distances between two labels, whatever represents them.
 *
 * Everything here takes #IMP::bff::States, so one implementation serves an
 * accessible volume, a rotamer library, a Gaussian, a coarse-grained ensemble
 * and an MD trajectory. That is the point: this was the third of **six**
 * implementations of dye-pair distances in the package, and they were checked
 * against each other before being merged — on T4L A132 x A65, 2180 x 3087
 * points, 2e5 samples, \f$\langle R_{DA}\rangle\f$ came out 51.945 against
 * 51.940 A and \f$\langle R_{DA}\rangle_E\f$ 51.696 against 51.692 A, with
 * \f$R_{mp}\f$ identical. The spread is Monte-Carlo sampling noise.
 *
 * There is **one sampler** now, #IMP::bff::random_distances. The Python this
 * replaces drew from `numpy.random.RandomState`, so every quantity built on it
 * agreed with the compiled kernel only to about \f$1/\sqrt{n}\f$ — two
 * estimators of the same integral rather than two computations of the same
 * number. Anything pinned against the old stream has to be regenerated.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_STATESDISTANCE_H
#define IMPBFF_STATESDISTANCE_H

#include <IMP/bff/bff_config.h>
#include <IMP/bff/AVModel.h>

#include <IMP/value_macros.h>
#include <IMP/showable_macros.h>

#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! Monte-Carlo samples a pair distance is estimated from, unless stated.
/*! In the header rather than the source so SWIG emits it as a module constant:
    it is a default the Python surface states too. */
const int N_DISTANCE_SAMPLES = 50000;

//! \f$\langle R_{DA}\rangle\f$ -- the weighted mean inter-state distance, A.
IMPBFFEXPORT double states_average_distance(const States& s1, const States& s2,
                                            int n_samples = 50000);

//! \f$R_E\f$ -- the FRET-averaged distance, A.
/*! The *efficiency* is averaged and then inverted, not the distance. */
IMPBFFEXPORT double states_mean_fret_distance(const States& s1,
                                              const States& s2,
                                              double forster_radius = 52.0,
                                              int n_samples = 50000);

//! \f$R_{mp}\f$ -- the distance between the two mean positions, A.
IMPBFFEXPORT double distance_between_mean_positions(const States& s1,
                                                    const States& s2);

//! The width of the pair-distance distribution, A.
IMPBFFEXPORT double standard_deviation_of_distances(const States& s1,
                                                    const States& s2,
                                                    int n_samples = 50000);

//! \f$(R_{mp}, \langle R_{DA}\rangle, R_E, \sigma_R)\f$ in one pass.
/*! An empty cloud on either side gives \f$R_{mp}\f$ for all three distances and
    zero width -- the mean position is the one thing a label without a volume
    still has. */
IMPBFFEXPORT std::vector<double> av_pair_statistics(
        const States& s1, const States& s2, double forster_radius = 52.0,
        int n_samples = 50000);

//! The model distance of a named type: `"Rmp"`, `"RDAMean"` or `"RDAMeanE"`.
/*! \throw ValueException for any other name */
IMPBFFEXPORT double model_distance(const States& s1, const States& s2,
                                   const std::string& distance_type,
                                   double forster_radius = 52.0,
                                   int n_samples = 50000);

//! \f$p(R_{DA})\f$ over \p axis (bin edges), weighted by the pair weights.
IMPBFFEXPORT void histogram_rda(const States& s1, const States& s2,
                                const std::vector<double>& axis, int n_samples,
                                bool normalize, double** out_view,
                                int* n_out_view);

//! Fit \f$R_{mp} \to \langle R_{DA}\rangle\f$ (or \f$R_E\f$) by translation.
/*!
    The two clouds are rigidly separated along the line joining their mean
    positions and the target distance is recomputed at seven offsets, then a
    polynomial of \p degree is fitted through them. That is what a docking run
    needs: a cheap \f$R_{mp}\f$ during sampling, corrected afterwards.

    An empty cloud, a zero total weight, or coincident mean positions give the
    identity polynomial — there is nothing to correct.

    \return coefficients, highest power first, as `numpy.polyfit` returns them
*/
IMPBFFEXPORT std::vector<double> fit_transfer_polynomial(
        const States& s1, const States& s2, const std::string& distance_type,
        double forster_radius = 52.0, int degree = 3, int n_samples = 10000);

//! \f$R_{mp} + \sigma^2 / R_{mp}\f$ -- the Gaussian mean-position correction.
/*!
    **\p sigma is the per-component width of the separation vector**, not of one
    cloud and not of the distance distribution. Writing the separation as
    \f$d = \Delta + \varepsilon\f$ with
    \f$\varepsilon \sim N(0, \sigma^2 I_3)\f$, the first-order term averages
    away and the **two** transverse components give
    \f$\langle|\varepsilon_\perp|^2\rangle = 2\sigma^2\f$, so the correction is
    \f$\sigma^2/R_{mp}\f$ — not half of it, which is what this returned under
    one of its two former names until 2026-08-18.

    **Zero where \p rmp is zero or negative**: the expansion is in
    \f$\sigma/R_{mp}\f$ and says nothing at coincident mean positions. Clamping
    the denominator instead returns 3.6e11 A there, a number that then
    propagates as if it meant something.
*/
IMPBFFEXPORT double gaussian_rmp_to_rda_mean(double rmp, double sigma);

//! Evaluate a transfer polynomial whose coefficients run **lowest power first**.
/*!
    \f$c_0 + c_1 x + c_2 x^2 + \dots\f$, the order a human writes by hand.
    #IMP::bff::polynomial_transfer takes the opposite order -- highest power
    first, which is what `numpy.polyfit` returns and therefore what a *fitted*
    calibration carries.

    The two are not interchangeable and the difference is silent: at
    `rmp = 45` with `[0, 1, 0.02]` this gives **85.5** and the other gives
    **45.02**. Both were called `polynomial_transfer`, in two modules, for as
    long as both existed.
*/
IMPBFFEXPORT double polynomial_transfer_ascending(
        double rmp, const std::vector<double>& coeffs);

//! \f$R_{mp}\f$ between two **point clouds**, rather than two `States`.
/*!
    Written over arrays so a rotamer library, a coarse-grained ensemble or an MD
    frame reaches it without building a `States` first.

    \param[in] points_a,points_b flat `(n, 3)` clouds, A
    \param[in] weights_a,weights_b per-point weights; uniform when empty
    \throw ValueException when either cloud is empty
*/
IMPBFFEXPORT double mean_position_distance(
        const std::vector<double>& points_a, const std::vector<double>& points_b,
        const std::vector<double>& weights_a = std::vector<double>(),
        const std::vector<double>& weights_b = std::vector<double>());

//! Where a dye is, as something that produces states on demand.
/*!
    Two of them: an accessible volume computed from a structure's obstacles, and
    an isotropic Gaussian that needs no structure at all. They exist so a caller
    can hold *a label* and ask it for states, without deciding which model
    produced them.

    The volume is computed **lazily**, on the first call to
    #IMP::bff::LabelDistribution::get_accessible_volume, because building one is
    the expensive thing a label does and a caller often holds several before
    asking any of them anything.
*/
class IMPBFFEXPORT LabelDistribution {
protected:
    std::string simulation_type_;
    std::vector<double> origin_;
    double simulation_grid_resolution_;
    std::string position_name_;
    mutable AccessibleVolume av_;
    mutable bool computed_;

    virtual void do_compute() const = 0;

public:
    LabelDistribution(const std::string& simulation_type = "AV1",
                      const std::vector<double>& origin = std::vector<double>(),
                      double simulation_grid_resolution = 0.5,
                      const std::string& position_name = "");
    virtual ~LabelDistribution() {}

    //! The states, computing them on first use.
    const AccessibleVolume& get_accessible_volume() const;

    std::string get_simulation_type() const { return simulation_type_; }
    std::string get_position_name() const { return position_name_; }
    double get_simulation_grid_resolution() const {
        return simulation_grid_resolution_;
    }
    void get_origin(double** out_view, int* n_out_view) const;
    int get_n_points() const { return get_accessible_volume().get_n_points(); }
};

//! A label whose states come from an accessible-volume search.
class IMPBFFEXPORT LabelDistributionAV : public LabelDistribution {
    std::vector<double> atoms_xyzr_;   //!< flat, four per obstacle
    std::vector<double> source_xyz_;
    double linker_length_, linker_width_, r1_, r2_, r3_;

    void do_compute() const override;

public:
    //! \param[in] atoms_xyzr,n_atoms,n_cols obstacles, `(N, 4)`
    /*! \param[in] source_xyz where the linker is tied
        \param[in] linker_length,linker_width the linker
        \param[in] r1,r2,r3 dye radii; a zero `r2` means the AV1 model
        \param[in] simulation_grid_resolution voxel spacing, A
        \param[in] position_name a human-readable label */
    LabelDistributionAV(double* atoms_xyzr, int n_atoms, int n_cols,
                        const std::vector<double>& source_xyz,
                        double linker_length = 20.0, double linker_width = 0.5,
                        double r1 = 3.5, double r2 = 0.0, double r3 = 0.0,
                        double simulation_grid_resolution = 1.5,
                        const std::string& position_name = "");
};

//! A label modelled as an isotropic 3-D Gaussian. No structure needed.
class IMPBFFEXPORT DyeDistributionNormal : public LabelDistribution {
    double width_;
    int n_points_;
    int seed_;

    void do_compute() const override;

public:
    //! \param[in] origin the mean position
    /*! \param[in] width the per-axis standard deviation, A
        \param[in] n_points states to draw
        \param[in] seed for reproducibility -- the Python this replaces drew
                   from the global `numpy.random` state, so two labels built in
                   the same process were not reproducible at all
        \param[in] position_name a human-readable label */
    DyeDistributionNormal(const std::vector<double>& origin, double width = 6.0,
                          int n_points = 50000, int seed = 0,
                          const std::string& position_name = "");

    double get_width() const { return width_; }
};

//! \f$R_{mp} \to \langle R_{DA}\rangle\f$ under a Gaussian pair model.
/*!
    A lookup table over centre-to-centre distance: at each one the pair distance
    is distributed as the sum of two mirrored normals of width \p sigma, and the
    three conventions -- mean distance, FRET-averaged distance, mean efficiency
    -- are integrated over it once and interpolated afterwards.

    This is what a docking run needs when it has an \f$R_{mp}\f$ and an
    experiment that measured something else.
*/
class IMPBFFEXPORT FRETDistanceConverter {
    double forster_radius_, sigma_;
    std::vector<double> distances_, efficiencies_;
    std::vector<double> d_mean_, d_mean_fret_, e_mean_;

    void update_efficiencies();
    void update_lookup();

public:
    FRETDistanceConverter(double forster_radius = 52.0, double sigma = 6.0,
                          double distance_min = 1.0, double distance_max = 100.0,
                          int n_distances = 128);

    double get_forster_radius() const { return forster_radius_; }
    void set_forster_radius(double v);
    double get_sigma() const { return sigma_; }
    void set_sigma(double v);

    //! \f$\langle R_{DA}\rangle\f$ at a centre-to-centre distance.
    double get_distance_mean(double dist_center_center) const;
    //! \f$R_E\f$ at a centre-to-centre distance.
    double get_distance_mean_fret(double dist_center_center) const;
    //! \f$\langle E\rangle\f$ at a centre-to-centre distance.
    double get_fret_efficiency_mean(double dist_center_center) const;

    //! The distance of a named convention, by the `DYE_PAIR_DISTANCE_*` codes.
    double get_effective_distance(double value, int distance_type) const;

    IMP_SHOWABLE_INLINE(FRETDistanceConverter,
                        out << "FRETDistanceConverter(R0=" << forster_radius_
                            << ", sigma=" << sigma_ << ")");
};
IMP_VALUES(FRETDistanceConverter, FRETDistanceConverters);

IMPBFF_END_NAMESPACE

#endif //IMPBFF_STATESDISTANCE_H
