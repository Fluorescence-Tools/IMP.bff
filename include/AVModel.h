#ifndef IMPBFF_AVMODEL_H
#define IMPBFF_AVMODEL_H

/**
 *  \file IMP/bff/AVModel.h
 *  \brief The accessible volume as a representation of a label's states,
 *         and the label distributions built on it.
 *
 *  One representation of #IMP::bff::States (which lives in States.h with
 *  everything representation-neutral): #IMP::bff::AccessibleVolume enumerates
 *  the states on a grid, #IMP::bff::ACV adds the contact volume, and
 *  #IMP::bff::LabelDistribution (formerly in `StatesDistance.h`) is a label
 *  that produces such a volume on demand. The other representation is the
 *  rotamer ensemble, in Rotamer.h.
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

#include <IMP/bff/States.h>
#include <IMP/bff/Base.h>

#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

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

// -------- from StatesDistance.h (the label distributions) --------
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
    //! The cloud as a numpy view, via the underlying states (which may compute
    //! on first access).
    void get_points(double** out_view, int* n_out_view) const;
    //! Weight-averaged position, three values (via the underlying states).
    void get_mean_position(double** out_view, int* n_out_view) const;
    //! Distance between the two labels' mean positions, A (via the states).
    double dRmp(const LabelDistribution& other) const;
    //! Mean over sampled point pairs, A (via the states).
    double dRDA(const LabelDistribution& other, int n_samples = 50000) const;
    //! FRET-averaged distance, A (via the states).
    double dRDAE(const LabelDistribution& other, double forster_radius,
                 int n_samples = 50000) const;
    int get_n_points() const { return get_accessible_volume().get_n_points(); }
};

//! A label whose states come from an accessible-volume search.
class IMPBFFEXPORT LabelDistributionAV : public LabelDistribution {
    std::vector<double> atoms_xyzr_;   //!< flat, four per obstacle
    std::vector<double> source_xyz_;
    double linker_length_, linker_width_, r1_, r2_, r3_;

    void do_compute() const override;

public:
    //! \param[in] atoms_xyz flat `(N, 3)` obstacle coordinates
    /*! \param[in] atoms_vdw flat `(N,)` obstacle van der Waals radii, A
        \param[in] source_xyz where the linker is tied
        \param[in] linker_length,linker_width the linker
        \param[in] r1,r2,r3 dye radii; a zero `r2` means the AV1 model
        \param[in] simulation_grid_resolution voxel spacing, A
        \param[in] position_name a human-readable label */
    LabelDistributionAV(const std::vector<double>& atoms_xyz,
                        const std::vector<double>& atoms_vdw,
                        const std::vector<double>& source_xyz =
                                std::vector<double>(),
                        double linker_length = 20.0, double linker_width = 0.5,
                        double r1 = 3.5, double r2 = 0.0, double r3 = 0.0,
                        double simulation_grid_resolution = 1.5,
                        const std::string& position_name = "");
};

//! A label modelled as an isotropic 3-D Gaussian. No structure needed.
class IMPBFFEXPORT ProbeDistributionNormal : public LabelDistribution {
    double width_;
    int n_points_;
    int seed_;

    void do_compute() const override;

public:
    //! \param[in] origin the mean position
    /*! \param[in] width the per-axis standard deviation, A
        \param[in] n_points states to draw
        \param[in] seed for reproducibility. A generator seeded here, rather
                   than a process-global one, is what makes two labels built in
                   one process reproducible
        \param[in] position_name a human-readable label */
    ProbeDistributionNormal(const std::vector<double>& origin, double width = 6.0,
                          int n_points = 50000, int seed = 0,
                          const std::string& position_name = "");

    double get_width() const { return width_; }
};

IMPBFF_END_NAMESPACE

#endif  // IMPBFF_AVMODEL_H
