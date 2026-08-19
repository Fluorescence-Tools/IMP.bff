/**
 * \file AVModel.cpp
 * \brief The accessible volume as an object: a point cloud, and what it answers.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/AVModel.h>
#include <IMP/bff/AVDistance.h>
#include <IMP/bff/internal/GridShape.h>
#include <IMP/bff/internal/OutputView.h>
#include <IMP/exception.h>

#include <cmath>
#include <cstring>

IMPBFF_BEGIN_NAMESPACE

BasicAV::BasicAV(const std::vector<double>& points,
                 const std::vector<double>& density,
                 const std::vector<double>& grid_origin,
                 double grid_step,
                 const std::string& position_name)
        : points_(points), density_(density), grid_origin_(grid_origin),
          grid_step_(grid_step), ng_(internal::cube_side(density.size())),
          position_name_(position_name) {
    if (!points_.empty() && points_.size() % 4 != 0) {
        IMP_THROW("a point cloud is four values per point (x, y, z, weight); "
                          << points_.size() << " is not a multiple of four",
                  IMP::ValueException);
    }
    if (!density_.empty() && ng_ == 0) {
        IMP_THROW("a density grid must be cubic; " << density_.size()
                          << " values are not a whole cube",
                  IMP::ValueException);
    }
    // Only a grid was given: derive the cloud, which is what every consumer
    // actually reads.
    if (points_.empty() && !density_.empty() && grid_origin_.size() == 3) {
        update_points();
    }
}

void BasicAV::get_points(double** out_view, int* n_out_view) const {
    internal::copy_to_view(points_, out_view, n_out_view);
}

void BasicAV::get_density(double** out_view, int* n_out_view) const {
    internal::copy_to_view(density_, out_view, n_out_view);
}

void BasicAV::get_grid_origin(double** out_view, int* n_out_view) const {
    internal::copy_to_view(grid_origin_, out_view, n_out_view);
}

void BasicAV::get_mean_position(double** out_view, int* n_out_view) const {
    internal::copy_to_view(points_weighted_mean(points_), out_view, n_out_view);
}

void BasicAV::set_points(const std::vector<double>& points) {
    if (!points.empty() && points.size() % 4 != 0) {
        IMP_THROW("a point cloud is four values per point (x, y, z, weight)",
                  IMP::ValueException);
    }
    points_ = points;
}

void BasicAV::update_points() {
    if (density_.empty() || grid_origin_.size() != 3) {
        IMP_THROW("no density grid available to convert", IMP::ValueException);
    }
    double* buf = nullptr;
    int n = 0;
    density_to_points(density_, ng_, ng_, ng_, grid_step_, grid_origin_, 0.0,
                      &buf, &n);
    points_.assign(buf, buf + n);
    std::free(buf);       // the view was never published; this side owns it
}

double BasicAV::dRmp(const BasicAV& other) const {
    const std::vector<double> a = points_weighted_mean(points_);
    const std::vector<double> b = points_weighted_mean(other.points_);
    const double dx = a[0] - b[0], dy = a[1] - b[1], dz = a[2] - b[2];
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

double BasicAV::dRDA(const BasicAV& other, int n_samples) const {
    return average_distance(points_, other.points_, n_samples, 0);
}

double BasicAV::dRDAE(const BasicAV& other, double forster_radius,
                      int n_samples) const {
    return mean_fret_distance(points_, other.points_, forster_radius,
                              n_samples, 0);
}

void BasicAV::pRDA(const BasicAV& other, const std::vector<double>& axis,
                   int n_samples, double** out_view, int* n_out_view) const {
    // `axis` is bin *edges*, so the histogram has one fewer bin than it has
    // edges -- the Python this replaces passed the same array to
    // `np.histogram(bins=axis)`, which reads it the same way.
    const std::size_t n_bins = axis.size() > 1 ? axis.size() - 1 : 0;
    double* out = internal::new_double_view(n_bins, out_view, n_out_view);
    if (out == nullptr || n_bins == 0) return;

    double* d = nullptr;
    int n = 0;
    random_distances(points_, other.points_, n_samples, 0, &d, &n);
    if (d == nullptr) return;

    double total = 0.0;
    for (int i = 0; i < n / 2; ++i) {
        const double r = d[2 * i + 0], w = d[2 * i + 1];
        if (r < axis.front() || r > axis.back()) continue;
        // Upper edge closed, as numpy's histogram has it: a sample exactly at
        // the top lands in the last bin rather than nowhere.
        std::size_t b = 0;
        while (b + 1 < n_bins && r >= axis[b + 1]) ++b;
        out[b] += w;
        total += w;
    }
    std::free(d);
    if (total > 0.0) {
        for (std::size_t b = 0; b < n_bins; ++b) out[b] /= total;
    }
}

// --------------------------------------------------------------------------
// ACV
// --------------------------------------------------------------------------

ACV::ACV(const std::vector<double>& points, const std::vector<double>& density,
         const std::vector<double>& grid_origin, double grid_step,
         const std::vector<double>& slow_centers,
         const std::vector<double>& slow_radius, double trapped_fraction,
         const std::string& position_name)
        : BasicAV(points, density, grid_origin, grid_step, position_name),
          slow_centers_(slow_centers), slow_radius_(slow_radius),
          trapped_fraction_(trapped_fraction) {
    if (!density_.empty() && !slow_centers_.empty()) {
        update_contact_density();
    }
}

ACV ACV::from_basic_av(const BasicAV& av,
                       const std::vector<double>& slow_centers,
                       const std::vector<double>& slow_radius,
                       double trapped_fraction) {
    double* d = nullptr; int nd = 0;
    double* o = nullptr; int no = 0;
    av.get_density(&d, &nd);
    av.get_grid_origin(&o, &no);
    ACV out(std::vector<double>(), std::vector<double>(d, d + nd),
            std::vector<double>(o, o + no), av.get_grid_step(),
            slow_centers, slow_radius, trapped_fraction,
            av.get_position_name());
    std::free(d);
    std::free(o);
    return out;
}

void ACV::update_contact_density() {
    if (density_.empty() || grid_origin_.size() != 3 || slow_centers_.empty()) {
        return;
    }
    const std::size_t n_centre = slow_centers_.size() / 3;
    std::vector<double> radius = slow_radius_;
    if (radius.empty()) radius.assign(n_centre, 10.0);
    if (radius.size() != n_centre) radius.assign(n_centre, radius[0]);

    int* label = nullptr; int n_label = 0;
    split_contact_volume(density_, ng_, grid_step_, radius, slow_centers_,
                         grid_origin_, &label, &n_label);
    if (label == nullptr) return;

    std::size_t n_contact = 0, n_free = 0;
    for (int i = 0; i < n_label; ++i) {
        if (label[i] == AV_VOXEL_CONTACT) ++n_contact;
        else if (label[i] == AV_VOXEL_FREE) ++n_free;
    }

    // Each side is scaled so that it carries its share of the total: the
    // contact voxels together hold `trapped_fraction`, the rest hold the
    // remainder. Dividing by the *count* is what makes it a share rather than a
    // multiplier -- two contact voxels each get half of the trapped fraction.
    contact_density_.assign(density_.size(), 0.0);
    std::vector<double> merged(density_.size(), 0.0);
    const double c_scale = n_contact > 0 ? trapped_fraction_ / n_contact : 0.0;
    const double f_scale =
            n_free > 0 ? (1.0 - trapped_fraction_) / n_free : 0.0;
    for (int i = 0; i < n_label; ++i) {
        if (label[i] == AV_VOXEL_CONTACT) {
            contact_density_[i] = density_[i] * c_scale;
            merged[i] = contact_density_[i];
        } else if (label[i] == AV_VOXEL_FREE) {
            merged[i] = density_[i] * f_scale;
        }
    }
    std::free(label);

    double total = 0.0;
    for (std::size_t i = 0; i < merged.size(); ++i) total += merged[i];
    if (total > 0.0) {
        for (std::size_t i = 0; i < merged.size(); ++i) merged[i] /= total;
    }
    density_ = merged;
    update_points();
}

void ACV::get_slow_centers(double** out_view, int* n_out_view) const {
    internal::copy_to_view(slow_centers_, out_view, n_out_view);
}

void ACV::get_slow_radius(double** out_view, int* n_out_view) const {
    internal::copy_to_view(slow_radius_, out_view, n_out_view);
}

void ACV::get_contact_density(double** out_view, int* n_out_view) const {
    internal::copy_to_view(contact_density_, out_view, n_out_view);
}

IMPBFF_END_NAMESPACE
