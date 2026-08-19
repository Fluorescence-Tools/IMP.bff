/**
 *  \file IMP/bff/AVModel.h
 *  \brief The accessible volume as an object: a point cloud, and what it answers.
 *
 * #IMP::bff::AV is the *solver* -- it grows a path map from an attachment point
 * and produces a density. #BasicAV is the *result*: the cloud that comes out,
 * and the questions a labelling site is asked. They are deliberately separate,
 * because a cloud does not have to come from a solve. It can be read from a
 * file, sampled from a rotamer library, or taken from a trajectory frame, and
 * everything downstream should not care which.
 *
 * #ACV adds the accessible **contact** volume: the part of the cloud lying
 * within reach of a surface, weighted up by a trapped fraction because a dye
 * that touches the protein stays there longer than free diffusion would put it
 * there. That is one number applied to a geometric split, not a model of
 * sticking -- see the class documentation.
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

#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! Accessible volume for one labelling site: a weighted point cloud.
/*!
    The cloud is `(x, y, z, weight)` per point, flat. A density grid may be
    carried alongside it -- some callers have one, some only have points -- and
    the cloud is derived from the grid on demand when only the grid was given.

    There is no padding here. The Python this replaces kept a `_padded_points`
    buffer and an `n_points` count because the numba kernel it once called
    returned a grid-sized array with unused trailing rows; the C++ kernel
    returns exactly the points it found, so the count is the cloud's length and
    the two cannot disagree.
*/
class IMPBFFEXPORT BasicAV {
protected:
    std::vector<double> points_;        //!< flat, four per point
    std::vector<double> density_;       //!< flat ng^3, empty if not carried
    std::vector<double> grid_origin_;   //!< the first voxel's centre
    double grid_step_;
    int ng_;
    std::string position_name_;

public:
    //! \param[in] points flat (x, y, z, w) per point; may be empty
    /*! \param[in] density flat ng^3; empty when only a cloud is known
        \param[in] grid_origin the first voxel's centre, three values
        \param[in] grid_step voxel spacing, A
        \param[in] position_name a human-readable label, e.g. "donor_72" */
    BasicAV(const std::vector<double>& points = std::vector<double>(),
            const std::vector<double>& density = std::vector<double>(),
            const std::vector<double>& grid_origin = std::vector<double>(),
            double grid_step = 1.5,
            const std::string& position_name = "");

    //! The cloud, as a numpy view over this object's own buffer.
    void get_points(double** out_view, int* n_out_view) const;
    //! The density grid; a zero-length view when none is carried.
    void get_density(double** out_view, int* n_out_view) const;
    void get_grid_origin(double** out_view, int* n_out_view) const;
    //! Density-weighted mean position, three values.
    void get_mean_position(double** out_view, int* n_out_view) const;

    int get_n_points() const { return static_cast<int>(points_.size() / 4); }
    int get_ng() const { return ng_; }
    double get_grid_step() const { return grid_step_; }
    std::string get_position_name() const { return position_name_; }
    void set_position_name(const std::string& n) { position_name_ = n; }

    //! Replace the cloud.
    void set_points(const std::vector<double>& points);
    //! Rebuild the cloud from the density grid.
    /*! \throws IMP::ValueException when no grid is carried */
    void update_points();

    //! Distance between the two mean positions, A.
    /*! The cheapest of the three and the least meaningful: it is the distance
        between two *averages*, which is not the average of the distance and can
        sit several Angstrom from either of the others. */
    double dRmp(const BasicAV& other) const;

    //! \f$\langle R_{DA}\rangle\f$ -- the mean over sampled point pairs, A.
    double dRDA(const BasicAV& other, int n_samples = 50000) const;

    //! \f$R_E\f$ -- the FRET-averaged distance, A.
    /*! The *efficiency* is averaged and then inverted, not the distance:
        \f$1/R^6\f$ weights close pairs far more heavily, so this is always the
        shorter of the two averages and it is the one a measured efficiency
        corresponds to. */
    double dRDAE(const BasicAV& other, double forster_radius,
                 int n_samples = 50000) const;

    //! \f$p(R_{DA})\f$ on \p axis (bin edges), normalised to sum 1.
    void pRDA(const BasicAV& other, const std::vector<double>& axis,
              int n_samples, double** out_view, int* n_out_view) const;

    IMP_SHOWABLE_INLINE(BasicAV,
                        out << "BasicAV(" << get_n_points() << " points"
                            << (position_name_.empty()
                                        ? std::string()
                                        : ", \"" + position_name_ + "\"")
                            << ")");
};
IMP_VALUES(BasicAV, BasicAVs);

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
class IMPBFFEXPORT ACV : public BasicAV {
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

    //! Build an ACV from a solved BasicAV, which must carry a density grid.
    static ACV from_basic_av(const BasicAV& av,
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
