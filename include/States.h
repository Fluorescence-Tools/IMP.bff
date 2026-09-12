#ifndef IMPBFF_STATES_H
#define IMPBFF_STATES_H

/**
 *  \file IMP/bff/States.h
 *  \brief A label's states, whatever represents them: the cloud, the kernels
 *         over it, and the distances between two.
 *
 * Representation-neutral. A label's configuration space is a weighted point
 * cloud (#IMP::bff::States); an accessible volume (ProbeAccessibleVolume.h) and a rotamer
 * ensemble (Rotamer.h) are two ways of enumerating one, and everything here
 * treats them identically. Three sections:
 *
 * 1. **Kernels** (formerly `AVDistance.h`) -- distances and reductions over
 *    point clouds, and the fps.json distance vocabulary.
 * 2. **States** (formerly the base class in `ProbeAccessibleVolume.h`) -- the cloud as a
 *    C++ value, with the three distances that are functions of it alone.
 * 3. **Distances between two labels** (formerly `StatesDistance.h`, less the
 *    accessible-volume-backed #IMP::bff::LabelDistribution, which is in
 *    ProbeAccessibleVolume.h with the representation it wraps) -- the pair statistics,
 *    the FRET distance converter, the transfer polynomials.
 *
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

// -------- from AVDistance.h --------
/**
 *  (formerly IMP/bff/AVDistance.h, now a section of this file)
 *  \brief Distances and reductions over accessible-volume point clouds.
 *
 * These take **point arrays** rather than the ``AV`` decorators the
 * `av_distance` family in ProbeAccessibleVolumeDecorator.h takes, so they serve a rotamer library, a
 * coarse-grained ensemble or an MD trajectory as readily as an accessible
 * volume -- anything that supplies states.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#include <IMP/bff/bff_config.h>

#include <vector>

IMPBFF_BEGIN_NAMESPACE

/// Different types of distances between two accessible volumes
typedef enum{
    PROBE_PAIR_DISTANCE_E,            /// Mean FRET averaged distance R_E
    PROBE_PAIR_DISTANCE_MEAN,         /// Mean distance <R_DA>
    PROBE_PAIR_DISTANCE_MP,           /// Distance between AV mean positions
    PROBE_PAIR_EFFICIENCY,            /// Mean FRET efficiency
    PROBE_PAIR_DISTANCE_DISTRIBUTION, /// Distance distribution
    PROBE_PAIR_XYZ_DISTANCE,          /// Distance between XYZ of dye particles
    PROBE_PAIR_DISTANCE_MIN           /// Closest approach of the two clouds
} ProbePairMeasures;


//! The fps.json name of a probe-pair distance type -- `RDAMean`, `RDAMeanE`,
//! `Rmp`, `Efficiency`, `pRDA`, `Rmin`.
/*! The vocabulary is the file format's, so it belongs beside the enum and not
    inside a reader: a table that reports what it scored needs the same names
    the file used. An unrecognised type gives the empty string. */
IMPBFFEXPORT std::string probe_pair_distance_type_name(int distance_type);

//! The type an fps.json distance name denotes; -1 when the name is unknown.
IMPBFFEXPORT int probe_pair_distance_type(const std::string& name);


//! What split_contact_volume() found in a voxel.
enum AVVoxel {
    AV_VOXEL_EMPTY = 0,       //!< no density: outside the accessible volume
    AV_VOXEL_CONTACT = 1,     //!< inside a slow centre's sphere
    AV_VOXEL_FREE = 2         //!< accessible, and away from every slow centre
};

//! Weighted mean position of a point cloud.
/*!
    \param[in] points flat, four per point: x, y, z, weight
    \return three coordinates; the origin for an empty or zero-weight cloud
*/
IMPBFFEXPORT std::vector<double> points_weighted_mean(
        const std::vector<double>& points);

//! Random distance and weight-product samples between two clouds.
/*!
    Draws one point from each cloud per sample. **The two clouds are
    independent**, so any pairing samples the joint distribution.

    The clouds arrive straight from numpy, one `(n, 4)` array per cloud via
    IN_ARRAY2. The result is published as an `(n_samples, 2)` managed numpy
    view: column 0 is the distance, column 1 the weight product. Both shapes
    are part of the contract, so a caller never reshapes a flat buffer to
    recover the picture the data already has.

    \param[in] p1,p2 flat clouds, four per point. Via IN_ARRAY2, carrying the
                (n, 4) shape; the values are read as `4 * n` in row order.
    \param[in] n_samples samples to draw
    \param[in] seed for reproducibility
    \param[out] output,n_output1,n_output2 `(n_samples, 2)` view. The buffer
                is malloc'd and numpy adopts it.
*/
IMPBFFEXPORT void random_distances(
        double* p1, int n_p1, int n_p1c,
        double* p2, int n_p2, int n_p2c,
        int n_samples,
        int seed,
        double** output, int* n_output1, int* n_output2
);

//! Weighted mean inter-point distance \f$\langle R_{DA}\rangle\f$.
IMPBFFEXPORT double average_distance(
        const std::vector<double>& p1,
        const std::vector<double>& p2,
        int n_samples = 50000,
        int seed = 0
);

//! FRET-averaged distance \f$\langle R_{DA}\rangle_E\f$.
/*!
    Averages the *efficiency* over the sampled pairs and converts back, which is
    not the same as averaging the distance -- \f$1/r^6\f$ weights close pairs far
    more heavily, so this is always the shorter of the two.

    \return 0 when the mean efficiency saturates at 1, infinity when it reaches 0
*/
IMPBFFEXPORT double mean_fret_distance(
        const std::vector<double>& p1,
        const std::vector<double>& p2,
        double forster_radius = 52.0,
        int n_samples = 50000,
        int seed = 0
);

//! The four distance statistics of a weighted distance sample.
/*!
    One reduction rather than four passes, and one place where the conventions
    live: \f$\langle R_{DA}\rangle\f$, \f$R_E\f$, \f$\langle E\rangle\f$ and
    the width of the distance distribution were computed in three modules that
    disagreed at the limits.

    Note what is **not** here: \f$R_{mp}\f$, the distance between the clouds'
    mean positions. It is not a function of the distribution of pair distances
    -- \f$|\langle a\rangle - \langle b\rangle|\f$ cannot be recovered from
    \f$|a - b|\f$ -- so it takes the clouds, not a sample of them. A slot in
    this tuple returned the mean distance under that name until 2026-07-28;
    on 148l E15/E90 the two were 8 % apart.

    \param[in] distances sampled pair distances
    \param[in] weights per-sample weight
    \param[in] forster_radius \f$R_0\f$
    \return four values: mean distance, \f$R_E\f$, mean efficiency, standard
            deviation. \f$R_E\f$ is 0 when the mean efficiency saturates at 1
            and infinity when it reaches 0; all four are 0 for zero total weight.
*/
IMPBFFEXPORT std::vector<double> distance_sample_statistics(
        const std::vector<double>& distances,
        const std::vector<double>& weights,
        double forster_radius = 52.0
);

//! Smallest distance between any point of one cloud and any of the other.
/*!
    A *bound*, where every other distance type here is an average: two labels
    can approach no closer than this, whatever the linkers do. That is the
    quantity a crosslink-style upper limit is written against, and the reason
    it is not a mean.

    Points of zero weight are skipped -- a zero-weight point is not in the
    volume. The double loop is exact and \f$O(n_1 n_2)\f$; at the few thousand
    points an AV carries that is milliseconds, and there is no sampled version
    because a minimum estimated from a sample is biased high.

    \param[in] p1,p2 flat clouds, four per point: x, y, z, weight
    \return the minimum distance; infinity when either cloud has no weight
*/
IMPBFFEXPORT double minimum_distance(const std::vector<double>& p1,
                                     const std::vector<double>& p2);

//! Fraction of a cloud's weight lying within \p radius of any reference point.
/*!
    How much of a volume touches something -- a patch of surface, a set of
    quenching residues, a binding partner. The answer is a weight fraction in
    [0, 1], so it does not depend on how finely the volume was rastered.

    \param[in] points the cloud, flat, four per point
    \param[in] reference_coords the reference positions, flat, **three** per
               point (they carry no weight)
    \param[in] radius the contact radius, A
    \return the fraction of \p points weight within \p radius of at least one
            reference position; 0 for an empty cloud or no references
    \throw ValueException when \p reference_coords is not a multiple of three
*/
IMPBFFEXPORT double cloud_overlap(const std::vector<double>& points,
                                  const std::vector<double>& reference_coords,
                                  double radius);

//! The model distance between two clouds, in a named convention.
/*!
    One dispatch for the point-cloud forms of the #IMP::bff::ProbePairMeasures
    conventions, so a caller that has two clouds rather than two
    #IMP::bff::ProbeAccessibleVolumeDecorator decorators does not have to pick the right reduction itself.

    \param[in] p1,p2 flat clouds, four per point
    \param[in] distance_type a #IMP::bff::ProbePairMeasures value
    \param[in] forster_radius \f$R_0\f$, for the two FRET conventions
    \param[in] n_samples,seed for the sampled conventions; `Rmp` and `Rmin` are
               exact and ignore both
*/
IMPBFFEXPORT double cloud_model_distance(
        const std::vector<double>& p1, const std::vector<double>& p2,
        int distance_type = PROBE_PAIR_DISTANCE_MEAN,
        double forster_radius = 52.0, int n_samples = 50000, int seed = 0);

//! The distance distribution between two clouds: the `pRDA` convention.
/*!
    `pRDA` is the one entry in #IMP::bff::ProbePairMeasures that is **not a
    distance**, and #IMP::bff::cloud_model_distance() refuses it for that
    reason. This is what it asks for instead: the weighted histogram of pair
    distances, \f$p(R_{DA})\f$, which is what a fluorescence decay is actually
    a function of.

    \param[in] p1,p2 flat clouds, four per point
    \param[in] axis the histogram's **bin edges**, ascending; a distance
               outside them is dropped rather than piled into an end bin
    \param[in] n_samples,seed the pair sample
    \param[in] normalize divide by the total weight, so the result sums to 1
               over the bins that were in range
    \return one value per bin -- `axis.size() - 1` of them
    \throw ValueException when \p axis has fewer than two edges or is not
           ascending
*/
IMPBFFEXPORT std::vector<double> cloud_distance_distribution(
        const std::vector<double>& p1, const std::vector<double>& p2,
        const std::vector<double>& axis, int n_samples = 50000, int seed = 0,
        bool normalize = true);

//! The mean-position separation that reproduces a given model distance.
/*!
    The inverse of #IMP::bff::cloud_model_distance in its first argument: slide
    the two clouds apart, rigidly and along the line joining their mean
    positions, until the modelled observable equals \p target_distance, and
    report the \f$R_{mp}\f$ at which that happens.

    This is the number that turns a measurement into a restraint an ordinary
    two-particle potential can carry. An experiment reports
    \f$\langle R_{DA}\rangle\f$; a distance restraint between two mean
    positions needs \f$R_{mp}\f$; the two differ by several angstrom and the
    difference depends on the shapes of both volumes, so it cannot be a
    constant. #IMP::bff::effective_distance is the *forward* approximation of
    the same relation, parametrised rather than measured; this is the exact
    answer for the two volumes in hand.

    The pairs are sampled **once** and then translated, so the function being
    bisected is smooth: resampling per iteration would add noise of order
    \f$1/\sqrt{n}\f$ to a root that is being located to \p accuracy. Bracket by
    expansion, then bisect -- which is stricter than the fixed-point iteration
    of the reference implementation, whose step assumes
    \f$dR_{model}/dR_{mp} = 1\f$ and which stops when it stops improving rather
    than when it has converged.

    \param[in] p1,p2 flat clouds, four per point
    \param[in] target_distance the measured value to reproduce
    \param[in] distance_type which convention \p target_distance is in
    \param[in] forster_radius \f$R_0\f$
    \param[in] accuracy stop when the model distance is this close, A
    \param[in] n_samples,seed the fixed pair sample
    \return the \f$R_{mp}\f$ that reproduces \p target_distance; 0 when either
            cloud is empty
    \throw ValueException when \p target_distance is unreachable by any
           separation -- a mean distance shorter than the volumes allow, or an
           efficiency outside (0, 1)
*/
IMPBFFEXPORT double rmp_from_model_distance(
        const std::vector<double>& p1, const std::vector<double>& p2,
        double target_distance,
        int distance_type = PROBE_PAIR_DISTANCE_MEAN,
        double forster_radius = 52.0, double accuracy = 0.01,
        int n_samples = 50000, int seed = 0);

//! Asymmetric chi-squared contribution of one distance restraint.
/*!
    \f$\chi^2 = (d_m - d_e)^2 / \sigma^2\f$. **The residual is model minus
    data**, which is the definition the whole package uses: a model distance
    that is too large is a positive deviation. \f$\sigma\f$ is the
    experimental error on that side -- a positive residual is judged against
    \p error_pos, a negative one against \p error_neg. Squaring makes the sign
    invisible in the value; what the convention decides is which error bar
    applies, and with asymmetric errors that is the whole answer.

    A non-positive error contributes nothing rather than dividing by zero.
*/
IMPBFFEXPORT double chi2_score(double model_distance,
                               double experimental_distance,
                               double error_neg, double error_pos);

//! The same \f$\chi^2\f$ with FPS's force cap: harmonic near zero, linear far.
/*!
    FPS's `SpringEngine.ForceAndTorque` (`SpringEngine.cs:432-443`) does not
    score a plain parabola. With \f$k = 2/\sigma^2\f$ and
    \f$\Delta r_{max} = F_{max}/k = F_{max}\sigma^2/2\f$
    (`SpringEngine.cs:141-142`) it scores

    \f[
      E(\Delta r) = \begin{cases}
        (\Delta r/\sigma)^2                      & |\Delta r| \le \Delta r_{max}\\
        F_{max}\,(|\Delta r| - \Delta r_{max}/2) & |\Delta r| >  \Delta r_{max}
      \end{cases}
    \f]

    -- a Huber-like robustification, joined \f$C^1\f$: past the knee the force
    is exactly \p max_force and the cost grows only linearly, so one grossly
    violated restraint cannot dominate a docking run. **A pure harmonic does
    not reproduce FPS under docking's `MaxForce = 400`**, which is the
    situation every random restart starts in.

    \p max_force is stated in FPS's units -- it caps \f$dE/d|\Delta r|\f$ of
    the \f$\chi^2\f$ this returns, so the knee sits at
    \f$F_{max}\sigma^2/2\f$ exactly where FPS puts it. A restraint that scores
    half a \f$\chi^2\f$ (which is what a Gaussian \f$-\log L\f$ is, and what
    #IMP::bff::AVPairDistanceMeasurement::score_model returns) therefore feels
    half that slope; the *knee position*, which is what changes a docked
    answer, is unaffected by the factor.

    The cap engages on tight restraints and leaves loose ones alone: with
    \f$\sigma = 0.1\f$ Å and `max_force = 400` the knee is at 2 Å, with
    \f$\sigma = 3\f$ Å it is at 1800 Å and the cap is inert.

    \param[in] max_force the force cap; **non-positive means no cap**, and
               #IMP::bff::chi2_score is then returned unchanged
    \see IMP::bff::chi2_score for the sign convention (model minus data, and
         which error bar that selects)
*/
IMPBFFEXPORT double chi2_score_capped(double model_distance,
                                      double experimental_distance,
                                      double error_neg, double error_pos,
                                      double max_force);

//! Single-pair FRET efficiency \f$1/(1 + (r/R_0)^6)\f$.
IMPBFFEXPORT double fret_efficiency(double distance,
                                    double forster_radius = 52.0);

//! The distance a single-pair FRET efficiency implies.
IMPBFFEXPORT double distance_from_fret_efficiency(
        double efficiency, double forster_radius = 52.0);

//! Occupied voxels of a density as a weighted point cloud.
/*!
    Published as an `(n, 4)` managed numpy view: x, y, z, density. The shape is
    part of the contract, stated once here rather than by a Python `.reshape()`.

    \param[in] density `(nx, ny, nz)` numpy grid, via IN_ARRAY3
    \param[in] nx,ny,nz grid shape
    \param[in] dg voxel edge
    \param[in] r0 grid anchor
    \param[in] threshold keep voxels strictly above this
    \param[out] output,n_output1,n_output2 `(n, 4)` view; the buffer is
                malloc'd and numpy adopts it.
*/
IMPBFFEXPORT void density_to_points(
        double* density, int nx, int ny, int nz,
        double dg,
        const std::vector<double>& r0,
        double threshold = 0.0,
        double** output = 0, int* n_output1 = 0, int* n_output2 = 0
);
//! Label each voxel of an accessible volume contact, free, or empty.
/*!
    A voxel is *contact* when it lies inside any slow centre's sphere. One
    label array rather than a pair of masks: the two are complementary within
    the occupied volume, so returning both invites them to disagree.

    \p r0 is **the grid anchor**: the position of voxel \f$(ng-1)/2\f$ on each
    axis, not the corner and not the attachment atom. The inverse map uses
    `floor` and the *integer* offset, both deliberately: the float corner
    differs on every even \p ng -- the normal case -- putting the stamped
    spheres half a voxel from the density they mask, and truncation would round
    a negative offset up while every other index map rounds down.

    \param[in] density flat ng^3
    \param[in] ng voxels per axis
    \param[in] dg voxel edge
    \param[in] rad per-centre slow radius
    \param[in] rs centre coordinates, flat
    \param[in] r0 grid anchor
    \param[out] output_i,dim1,dim2,dim3 an `(ng, ng, ng)` #AVVoxel view; the
                buffer is malloc'd and numpy adopts it.
*/
IMPBFFEXPORT void split_contact_volume(
        const std::vector<double>& density,
        int ng,
        double dg,
        const std::vector<double>& rad,
        const std::vector<double>& rs,
        const std::vector<double>& r0,
        int** output_i, int* dim1, int* dim2, int* dim3
);

//! The contact and free regions of an accessible volume, as two uint8 masks.
/*!
    A voxel is *contact* when it lies inside any slow centre's sphere. The two
    masks are complementary within the occupied volume; returning both invites
    them to disagree, which is why this is one call publishing the pair.

    `uint8`, not `float64`: at `ng = 92` a float64 pair costs 12.5 MB per
    labelling site against 1.6 MB. A residue scan builds one per site.

    \param[in] density `(ng, ng, ng)`, straight from numpy via IN_ARRAY3
    \param[in] radius a scalar slow radius, broadcast over every centre, or one
                value per centre
    \param[out] contact,free each `(ng, ng, ng)` uint8 views; both buffers are
                malloc'd and numpy adopts them
*/
IMPBFFEXPORT void split_contact_volume_masks(
        double* density, int ng, int ng2, int ng3,
        double dg,
        const std::vector<double>& radius,
        double* rs, int n_rs, int n_rsc,
        const std::vector<double>& r0,
        unsigned char** contact, int* contact_dim1, int* contact_dim2,
        int* contact_dim3,
        unsigned char** free, int* free_dim1, int* free_dim2, int* free_dim3
);

IMPBFF_END_NAMESPACE

// -------- from ProbeAccessibleVolume.h (the States base class) --------

#include <IMP/bff/IMPCompatibility.h>

#include <map>
#include <string>

IMPBFF_BEGIN_NAMESPACE

//! A weighted set of states of one label.
/*!
    The states are `(x, y, z, weight)` per point, flat. Everything a
    representation has to supply is here and nothing else is: the grid an
    accessible volume was enumerated on belongs to
    #IMP::bff::ProbeAccessibleVolume, and the conformers a rotamer library carries
    belong to the rotamer ensemble. Neither is part of what a distance or a
    rate needs.

    The three distances live here rather than one level down for the same
    reason: they are functions of the cloud, so a rotamer ensemble and an
    accessible volume answer them identically.
*/
class IMPBFFEXPORT States {
protected:
    std::vector<double> points_;             //!< flat, four per state
    std::vector<double> attachment_point_;   //!< three, or empty when unknown
    std::vector<double> orientations_;       //!< flat, three per state; empty
                                             //!< when the representation has none
    std::string position_name_;
    //! How these states were produced — representation parameters, not dye or
    //! site properties. Provenance: written by whatever built the states, and
    //! stringified because it is read back as a record rather than as numbers.
    std::map<std::string, std::string> params_;

public:
    //! \param[in] points flat (x, y, z, w) per state; may be empty
    /*! \param[in] attachment_point where the label is tied to the structure
        \param[in] orientations flat transition dipoles, three per state
        \param[in] position_name a human-readable label, e.g. "donor_72"
        \param[in] params how these states were produced */
    States(const std::vector<double>& points = std::vector<double>(),
           const std::vector<double>& attachment_point = std::vector<double>(),
           const std::vector<double>& orientations = std::vector<double>(),
           const std::string& position_name = "",
           const std::map<std::string, std::string>& params =
                   std::map<std::string, std::string>());

    //! The cloud, as a numpy view over this object's own buffer.
    void get_points(double** out_view, int* n_out_view) const;
    //! Where the label is tied; a zero-length view when it is not known.
    void get_attachment_point(double** out_view, int* n_out_view) const;
    //! The transition dipoles; a zero-length view when there are none.
    void get_orientations(double** out_view, int* n_out_view) const;
    //! Weight-averaged position, three values.
    /*! Falls back to the attachment point when the cloud is empty or its
        weights sum to zero, because that is the one position a label always
        has. With no attachment point either, the origin. */
    void get_mean_position(double** out_view, int* n_out_view) const;

    //! The cloud itself, for a C++ caller. Python gets the numpy view above.
    const std::vector<double>& get_points_vector() const { return points_; }

    int get_n_points() const { return static_cast<int>(points_.size() / 4); }
    bool get_has_volume() const { return !points_.empty(); }
    bool get_has_orientations() const { return !orientations_.empty(); }
    std::string get_position_name() const { return position_name_; }
    void set_position_name(const std::string& n) { position_name_ = n; }
    std::map<std::string, std::string> get_params() const { return params_; }
    void set_params(const std::map<std::string, std::string>& p) { params_ = p; }

    //! Replace the cloud.
    void set_points(const std::vector<double>& points);
    //! Replace the transition dipoles; one per state, or empty for none.
    void set_orientations(const std::vector<double>& orientations);
    void set_attachment_point(const std::vector<double>& xyz);

    //! Distance between the two mean positions, A.
    /*! The cheapest of the three and the least meaningful: it is the distance
        between two *averages*, which is not the average of the distance and can
        sit several Angstrom from either of the others. */
    double dRmp(const States& other) const;

    //! \f$\langle R_{DA}\rangle\f$ -- the mean over sampled point pairs, A.
    double dRDA(const States& other, int n_samples = 50000) const;

    //! \f$R_E\f$ -- the FRET-averaged distance, A.
    /*! The *efficiency* is averaged and then inverted, not the distance:
        \f$1/R^6\f$ weights close pairs far more heavily, so this is always the
        shorter of the two averages and it is the one a measured efficiency
        corresponds to. */
    double dRDAE(const States& other, double forster_radius,
                 int n_samples = 50000) const;

    //! \f$p(R_{DA})\f$ on \p axis (bin edges), normalised to sum 1.
    void pRDA(const States& other, const std::vector<double>& axis,
              int n_samples, double** out_view, int* n_out_view) const;

    IMP_SHOWABLE_INLINE(States,
                        out << "States(" << get_n_points() << " states"
                            << (position_name_.empty()
                                        ? std::string()
                                        : ", \"" + position_name_ + "\"")
                            << ")");
};
IMP_VALUES(States, StatesList);

IMPBFF_END_NAMESPACE

// -------- from StatesDistance.h --------
/**
 *  (formerly IMP/bff/StatesDistance.h, now a section of this file)
 *  \brief Distances between two labels, whatever represents them.
 *
 * Everything here takes #IMP::bff::States, so one implementation serves an
 * accessible volume, a rotamer library, a Gaussian, a coarse-grained ensemble
 * and an MD trajectory. That is the point: this was the third of **six**
 * implementations of probe-pair distances in the package, and they were checked
 * against each other before being merged — on T4L A132 x A65, 2180 x 3087
 * points, 2e5 samples, \f$\langle R_{DA}\rangle\f$ came out 51.945 against
 * 51.940 A and \f$\langle R_{DA}\rangle_E\f$ 51.696 against 51.692 A, with
 * \f$R_{mp}\f$ identical. The spread is Monte-Carlo sampling noise.
 *
 * There is **one sampler**, #IMP::bff::random_distances, so two callers
 * asking for the same quantity compute the same number rather than two
 * estimates of one integral agreeing to \f$1/\sqrt{n}\f$.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */



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
/*! An empty \p axis (the default) spans the sample: bin edges on
    \f$[0, R_{max}]\f$ of \p n_samples pairs, \p n_samples + 1 edges (one more
    than the default sample count) -- the shape a caller wanting "just a
    histogram" gets. \p normalize scales the counts to sum to 1. */
IMPBFFEXPORT void histogram_rda(const States& s1, const States& s2,
                                const std::vector<double>& axis =
                                        std::vector<double>(),
                                int n_samples = 50000, bool normalize = true,
                                double** out_view = NULL,
                                int* n_out_view = NULL);

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

//! Apply a named transfer function to a mean-position distance.
/*!
    The dispatch an fps.json calibration names, in one place. `"None"` and the
    `Rmp` convention return \p rmp untouched; `"Gaussian"` applies
    gaussian_rmp_to_rda_mean(); `"Polynomial"` evaluates \p coeffs through
    polynomial_transfer_ascending() and **falls back to the Gaussian** when no
    coefficients were given -- a calibration that names a polynomial and carries
    none is a calibration that has not been fitted yet, and \p sigma_rda is the
    parametric stand-in for it.

    \param[in] rmp distance between the two mean positions, A
    \param[in] transfer_function_type `"None"`, `"Gaussian"` or `"Polynomial"`;
               an unrecognised name returns \p rmp, because a transfer function
               nobody implements is no transfer function
    \param[in] sigma_rda per-component width of the separation vector; a
               non-positive value disables the Gaussian correction
    \param[in] coeffs polynomial coefficients, **lowest power first**
*/
IMPBFFEXPORT double effective_distance(
        double rmp, const std::string& transfer_function_type,
        double sigma_rda = 0.0,
        const std::vector<double>& coeffs = std::vector<double>());

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

    //! The distance of a named convention, by the `PROBE_PAIR_DISTANCE_*` codes.
    double get_effective_distance(double value, int distance_type) const;

    IMP_SHOWABLE_INLINE(FRETDistanceConverter,
                        out << "FRETDistanceConverter(R0=" << forster_radius_
                            << ", sigma=" << sigma_ << ")");
};
IMP_VALUES(FRETDistanceConverter, FRETDistanceConverters);

// --------------------------------------------------------------------------
// Empirical corrections from a computed distance to a measured one
// --------------------------------------------------------------------------

//! Evaluate a polynomial by Horner's method.
/*!
    \param[in] x abscissa
    \param[in] coefficients **highest power first** -- what `np.polyfit`
               returns, and therefore what a fitted calibration carries. The
               opposite order was in use elsewhere and evaluated the same
               calibration to a different number (85.5 against 45.02 at
               `x = 45`, `c = [0, 1, 0.02]`); a caller holding ascending
               coefficients must reverse them.
    \return the polynomial's value; 0 for no coefficients
*/
IMPBFFEXPORT double polynomial_transfer(
        double x, const std::vector<double>& coefficients);

//! Evaluate a polynomial at many abscissae, as a managed view.
/*! The direct vector form; a caller with many abscissae uses this rather than
    dispatching a scalar across them. */
IMPBFFEXPORT void polynomial_transfer_vector(
        const std::vector<double>& x, const std::vector<double>& coefficients,
        double** out_view, int* n_out_view);

IMPBFF_END_NAMESPACE

#endif  // IMPBFF_STATES_H
