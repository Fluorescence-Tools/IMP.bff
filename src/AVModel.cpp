/**
 * \file AVModel.cpp
 * \brief The accessible volume: the grid-enumerated cloud, the contact volume,
 *        and the label distributions that produce one on demand.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#include <IMP/bff/AVModel.h>
#include <IMP/bff/AVBuilder.h>
#include <IMP/bff/StructureIO.h>
#include <IMP/bff/internal/GridShape.h>
#include <IMP/bff/Base.h>

#include <cstdlib>
#include <cstring>

IMPBFF_BEGIN_NAMESPACE

// --------------------------------------------------------------------------
// AccessibleVolume
// --------------------------------------------------------------------------

AccessibleVolume::AccessibleVolume(
        const std::vector<double>& points, const std::vector<double>& density,
        const std::vector<double>& grid_origin, double grid_step,
        const std::string& position_name,
        const std::vector<double>& attachment_point,
        const std::vector<double>& orientations,
        const std::map<std::string, std::string>& params)
        : States(points, attachment_point, orientations, position_name, params),
          density_(density), grid_origin_(grid_origin), grid_step_(grid_step),
          ng_(internal::cube_side(density.size())) {
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

void AccessibleVolume::get_density(double** out_view, int* n_out_view) const {
    internal::copy_to_view(density_, out_view, n_out_view);
}

void AccessibleVolume::get_grid_origin(double** out_view,
                                       int* n_out_view) const {
    internal::copy_to_view(grid_origin_, out_view, n_out_view);
}

void AccessibleVolume::update_points() {
    if (density_.empty() || grid_origin_.size() != 3) {
        IMP_THROW("no density grid available to convert", IMP::ValueException);
    }
    double* buf = nullptr;
    int n = 0, nc = 0;
    const int ng = ng_;
    density_to_points(
            const_cast<double*>(density_.data()), ng, ng, ng,
            grid_step_, grid_origin_, 0.0, &buf, &n, &nc);
    set_points(std::vector<double>(buf, buf + static_cast<std::size_t>(n) * nc));
    std::free(buf);       // the view was never published; this side owns it
}

// --------------------------------------------------------------------------
// ACV
// --------------------------------------------------------------------------

ACV::ACV(const std::vector<double>& points, const std::vector<double>& density,
         const std::vector<double>& grid_origin, double grid_step,
         const std::vector<double>& slow_centers,
         const std::vector<double>& slow_radius, double trapped_fraction,
         const std::string& position_name)
        : AccessibleVolume(points, density, grid_origin, grid_step,
                           position_name),
          slow_centers_(slow_centers), slow_radius_(slow_radius),
          trapped_fraction_(trapped_fraction) {
    // A single `slow_radius` is broadcast over every centre: a caller that has
    // one sticky-sphere radius for every residue should not have to repeat it.
    if (slow_radius_.size() == 1 && !slow_centers_.empty()) {
        slow_radius_.assign(slow_centers_.size() / 3, slow_radius_[0]);
    }
    if (!density_.empty() && !slow_centers_.empty()) {
        update_contact_density();
    }
}

ACV ACV::from_accessible_volume(const AccessibleVolume& av,
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
    // The attachment point and the provenance belong to the label, not to the
    // grid, so they survive the split.
    double* a = nullptr; int na = 0;
    av.get_attachment_point(&a, &na);
    out.set_attachment_point(std::vector<double>(a, a + na));
    std::free(a);
    out.set_params(av.get_params());
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

    int* label = nullptr; int n_label = 0, d2 = 0, d3 = 0;
    split_contact_volume(density_, ng_, grid_step_, radius, slow_centers_,
                         grid_origin_, &label, &n_label, &d2, &d3);
    if (label == nullptr) return;
    const std::size_t n_labels =
            static_cast<std::size_t>(n_label) * d2 * d3;

    std::size_t n_contact = 0, n_free = 0;
    for (std::size_t i = 0; i < n_labels; ++i) {
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

void write_av(const AccessibleVolume& av, const std::string& path) {
    const std::size_t dot = path.rfind('.');
    std::string ext = dot == std::string::npos ? "" : path.substr(dot + 1);
    for (auto& c : ext) c = (char) std::tolower((unsigned char) c);

    if (ext == "xyz") {
        write_points_xyz(path, av.get_points_vector(), "He",
                         "accessible volume, column 5 is the weight");
        return;
    }
    if (ext == "pqr") {
        write_points_pqr(path, av.get_points_vector(), 1.0);
        return;
    }
    if (ext == "dx") {
        const std::vector<double> d = av.get_density_vector();
        if (d.empty()) {
            IMP_THROW("this volume carries only a cloud, so it has no grid to "
                      "write as OpenDX -- use .xyz or .pqr", ValueException);
        }
        write_opendx(path, d, av.get_ng(), av.get_ng(), av.get_ng(),
                     av.get_grid_origin_vector(), av.get_grid_step());
        return;
    }
    IMP_THROW("cannot write " << path << ": the extension must be one of "
                              << "xyz, pqr, dx", ValueException);
}

// -------- from StatesDistance.cpp (the label distributions) --------
// --------------------------------------------------------------------------
// label distributions
// --------------------------------------------------------------------------

LabelDistribution::LabelDistribution(const std::string& simulation_type,
                                     const std::vector<double>& origin,
                                     double simulation_grid_resolution,
                                     const std::string& position_name)
    : simulation_type_(simulation_type), origin_(origin),
      simulation_grid_resolution_(simulation_grid_resolution),
      position_name_(position_name), computed_(false) {}

const AccessibleVolume& LabelDistribution::get_accessible_volume() const {
    if (!computed_) {
        do_compute();
        computed_ = true;
    }
    return av_;
}

void LabelDistribution::get_origin(double** out_view, int* n_out_view) const {
    internal::copy_to_view(origin_, out_view, n_out_view);
}

void LabelDistribution::get_points(double** out_view,
                                   int* n_out_view) const {
    get_accessible_volume().get_points(out_view, n_out_view);
}

void LabelDistribution::get_mean_position(double** out_view,
                                          int* n_out_view) const {
    get_accessible_volume().get_mean_position(out_view, n_out_view);
}

double LabelDistribution::dRmp(const LabelDistribution& other) const {
    return get_accessible_volume().dRmp(other.get_accessible_volume());
}

double LabelDistribution::dRDA(const LabelDistribution& other,
                               int n_samples) const {
    return get_accessible_volume().dRDA(other.get_accessible_volume(),
                                        n_samples);
}

double LabelDistribution::dRDAE(const LabelDistribution& other,
                                double forster_radius, int n_samples) const {
    return get_accessible_volume().dRDAE(other.get_accessible_volume(),
                                         forster_radius, n_samples);
}

LabelDistributionAV::LabelDistributionAV(
        const std::vector<double>& atoms_xyz,
        const std::vector<double>& atoms_vdw,
        const std::vector<double>& source_xyz, double linker_length,
        double linker_width, double r1, double r2, double r3,
        double simulation_grid_resolution, const std::string& position_name)
    : LabelDistribution(r2 == 0.0 ? "AV1" : "AV3", source_xyz,
                        simulation_grid_resolution, position_name),
      source_xyz_(source_xyz), linker_length_(linker_length),
      linker_width_(linker_width), r1_(r1), r2_(r2), r3_(r3) {
    const std::size_t n = atoms_xyz.size() / 3;
    if (atoms_xyz.size() % 3 != 0 || atoms_vdw.size() != n) {
        IMP_THROW("obstacle coords must be (N, 3) and vdw (N,), got "
                          << atoms_xyz.size() << " coords and "
                          << atoms_vdw.size() << " radii",
                  IMP::ValueException);
    }
    // Interleave (x, y, z, ball-radius) so compute_av sees its (N, 4) layout.
    atoms_xyzr_.resize(4 * n);
    for (std::size_t i = 0; i < n; ++i) {
        atoms_xyzr_[4 * i + 0] = atoms_xyz[3 * i + 0];
        atoms_xyzr_[4 * i + 1] = atoms_xyz[3 * i + 1];
        atoms_xyzr_[4 * i + 2] = atoms_xyz[3 * i + 2];
        atoms_xyzr_[4 * i + 3] = atoms_vdw[i];
    }
    // No explicit source means the first obstacle.
    if (source_xyz_.empty() && !atoms_xyzr_.empty()) {
        source_xyz_.assign(atoms_xyzr_.begin(), atoms_xyzr_.begin() + 3);
        origin_.assign(source_xyz_.begin(), source_xyz_.end());
    }
}

void LabelDistributionAV::do_compute() const {
    av_ = compute_av(const_cast<double*>(atoms_xyzr_.data()),
                     static_cast<int>(atoms_xyzr_.size() / 4), 4, source_xyz_,
                     linker_length_, linker_width_, r1_, r2_, r3_,
                     simulation_grid_resolution_, DEFAULT_ALLOWED_SPHERE_RADIUS,
                     0);
    av_.set_position_name(position_name_);
}

ProbeDistributionNormal::ProbeDistributionNormal(const std::vector<double>& origin,
                                             double width, int n_points,
                                             int seed,
                                             const std::string& position_name)
    : LabelDistribution("Normal", origin, 1.0, position_name), width_(width),
      n_points_(n_points), seed_(seed) {
    if (origin.size() != 3) {
        IMP_THROW("the origin is three coordinates, not " << origin.size(),
                  ValueException);
    }
}

void ProbeDistributionNormal::do_compute() const {
    // Box-Muller off a seeded engine, so two labels built in one process are
    // reproducible and so is one across runs.
    std::mt19937 engine(static_cast<unsigned int>(seed_));
    std::normal_distribution<double> normal(0.0, 1.0);

    std::vector<double> points(static_cast<std::size_t>(n_points_) * 4);
    double total = 0.0;
    for (int i = 0; i < n_points_; ++i) {
        double r2 = 0.0;
        for (int c = 0; c < 3; ++c) {
            const double d = normal(engine) * width_;
            points[i * 4 + c] = origin_[c] + d;
            r2 += d * d;
        }
        // The weight is the Gaussian height at the drawn point, which makes the
        // cloud an importance-weighted sample of the same distribution it was
        // drawn from rather than a uniform one.
        const double w = std::exp(-0.5 * r2 / (width_ * width_));
        points[i * 4 + 3] = w;
        total += w;
    }
    if (total > 0.0) {
        for (int i = 0; i < n_points_; ++i) points[i * 4 + 3] /= total;
    }
    av_ = AccessibleVolume(points, std::vector<double>(), std::vector<double>(),
                           1.0, position_name_, origin_);
}

IMPBFF_END_NAMESPACE
