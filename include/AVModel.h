/**
 *  \file IMP/bff/AVModel.h
 *  \brief What a label's configuration space is, whatever enumerates it.
 *
 * #IMP::bff::States is the shared answer: **positions, weights and
 * orientations**. An accessible volume, a rotamer library, a coarse-grained
 * conformer set and an MD trajectory are four ways of producing the same three
 * things, so distances, kappa^2 and the interaction terms are written once and
 * serve all of them. A representation that resolves transition dipoles says so;
 * one that cannot — an accessible volume — is *asked*, rather than silently
 * averaged over an isotropic assumption.
 *
 * #IMP::bff::AV is the *solver* — it grows a path map from an attachment point
 * and produces a density. #IMP::bff::AccessibleVolume is the *result*: those
 * states plus the grid they were enumerated on. They are deliberately separate,
 * because a cloud does not have to come from a solve. It can be read from a
 * file, sampled from a rotamer library, or taken from a trajectory frame, and
 * everything downstream should not care which.
 *
 * #IMP::bff::ACV adds the accessible **contact** volume: the part of the cloud
 * lying within reach of a surface, weighted up by a trapped fraction because a
 * dye that touches the protein stays there longer than free diffusion would put
 * it there. That is one number applied to a geometric split, not a model of
 * sticking — see the class documentation.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_AVMODEL_H
#define IMPBFF_AVMODEL_H

#include <IMP/bff/bff_config.h>

#include <IMP/value_macros.h>
#include <IMP/showable_macros.h>

#include <map>
#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! A weighted set of states of one label.
/*!
    The states are `(x, y, z, weight)` per point, flat. Everything a
    representation has to supply is here and nothing else is: the grid an
    accessible volume was enumerated on belongs to
    #IMP::bff::AccessibleVolume, and the conformers a rotamer library carries
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

//! States enumerated as a voxel grid, plus the grid itself.
/*!
    A density grid may be carried alongside the cloud — some callers have one,
    some only have points — and the cloud is derived from the grid on demand
    when only the grid was given.

    The grid is cubic: #IMP::bff::PathMapHeader sizes it from the linker length
    in all three axes, so `ng` is one number rather than three.

    There is no padding: the kernel returns exactly the points it found, so
    the count is the cloud's length and a separate `n_points` cannot disagree
    with it.
*/
class IMPBFFEXPORT AccessibleVolume : public States {
protected:
    std::vector<double> density_;       //!< flat ng^3, empty if not carried
    std::vector<double> grid_origin_;   //!< the first voxel's centre
    double grid_step_;
    int ng_;

public:
    //! \param[in] points flat (x, y, z, w) per point; may be empty
    /*! \param[in] density flat ng^3; empty when only a cloud is known
        \param[in] grid_origin the first voxel's centre, three values
        \param[in] grid_step voxel spacing, A
        \param[in] position_name a human-readable label, e.g. "donor_72"
        \param[in] attachment_point where the label is tied to the structure
        \param[in] orientations flat transition dipoles, three per point
        \param[in] params how this volume was produced */
    AccessibleVolume(
            const std::vector<double>& points = std::vector<double>(),
            const std::vector<double>& density = std::vector<double>(),
            const std::vector<double>& grid_origin = std::vector<double>(),
            double grid_step = 1.5,
            const std::string& position_name = "",
            const std::vector<double>& attachment_point = std::vector<double>(),
            const std::vector<double>& orientations = std::vector<double>(),
            const std::map<std::string, std::string>& params =
                    std::map<std::string, std::string>());

    //! The density grid; a zero-length view when none is carried.
    void get_density(double** out_view, int* n_out_view) const;
    void get_grid_origin(double** out_view, int* n_out_view) const;

    //! The grid itself, for a C++ caller. Python gets the numpy views above.
    const std::vector<double>& get_density_vector() const { return density_; }
    const std::vector<double>& get_grid_origin_vector() const {
        return grid_origin_;
    }

    int get_ng() const { return ng_; }
    double get_grid_step() const { return grid_step_; }

    //! Rebuild the cloud from the density grid.
    /*! \throws IMP::ValueException when no grid is carried */
    void update_points();

    IMP_SHOWABLE_INLINE(AccessibleVolume,
                        out << "AccessibleVolume(" << get_n_points()
                            << " points"
                            << (position_name_.empty()
                                        ? std::string()
                                        : ", \"" + position_name_ + "\"")
                            << ")");
};
IMP_VALUES(AccessibleVolume, AccessibleVolumes);

//! Write a volume to \p path; the format follows the file extension.
/*!
    The overload for the value type #IMP::bff::compute_av() returns, so the
    array front door and the decorator front door are written the same way.

    | extension | what is written |
    |---|---|
    | `.xyz` | the point cloud, weight in a fifth column |
    | `.pqr` | the point cloud, weight in the charge column |
    | `.dx` | the density grid, as OpenDX |

    \param[in] av the volume
    \param[in] path the output path, whose extension picks the format
    \throw ValueException on an unknown extension, or on a grid format when
           the volume carries only a cloud
    \throw IOException when \p path cannot be opened
*/
IMPBFFEXPORT void write_av(const AccessibleVolume& av, const std::string& path);

//! Accessible contact volume: the cloud split by proximity to a surface.
/*!
    Points within \p slow_radius of any slow centre are the *contact* volume and
    are weighted up so that they carry \p trapped_fraction of the total. The
    rest share what is left.

    This is a geometric split with one number on it, not a model of sticking.
    The trapped fraction is a free parameter fitted per construct; it does not
    come from the structure, and a contact volume built with the default 0.8 is
    a statement about the shape of the answer, not a prediction of it.
*/
class IMPBFFEXPORT ACV : public AccessibleVolume {
    std::vector<double> slow_centers_;   //!< flat, three per centre
    std::vector<double> slow_radius_;    //!< one per centre, or one for all
    double trapped_fraction_;
    std::vector<double> contact_density_;

    void update_contact_density();

public:
    ACV(const std::vector<double>& points = std::vector<double>(),
        const std::vector<double>& density = std::vector<double>(),
        const std::vector<double>& grid_origin = std::vector<double>(),
        double grid_step = 1.5,
        const std::vector<double>& slow_centers = std::vector<double>(),
        const std::vector<double>& slow_radius = std::vector<double>(),
        double trapped_fraction = 0.8,
        const std::string& position_name = "");

    //! Build an ACV from a solved volume, which must carry a density grid.
    static ACV from_accessible_volume(const AccessibleVolume& av,
                                      const std::vector<double>& slow_centers,
                                      const std::vector<double>& slow_radius,
                                      double trapped_fraction = 0.8);

    void get_slow_centers(double** out_view, int* n_out_view) const;
    void get_slow_radius(double** out_view, int* n_out_view) const;
    //! The re-weighted grid; a zero-length view before it is built.
    void get_contact_density(double** out_view, int* n_out_view) const;
    double get_trapped_fraction() const { return trapped_fraction_; }

    IMP_SHOWABLE_INLINE(ACV,
                        out << "ACV(" << get_n_points() << " points, trapped "
                            << trapped_fraction_ << ")");
};
IMP_VALUES(ACV, ACVs);

IMPBFF_END_NAMESPACE

#endif //IMPBFF_AVMODEL_H
