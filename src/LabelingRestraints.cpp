/**
 * \file LabelingRestraints.cpp
 * \brief Chi-squared scoring of labelling data, with and without volumes.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#include <IMP/bff/LabelingRestraints.h>

#include <IMP/bff/AVDistance.h>
#include <IMP/exception.h>

#include <cmath>

IMPBFF_BEGIN_NAMESPACE

AVMeasurement::AVMeasurement(std::string av1_name, std::string av2_name,
                             double distance, double error_neg,
                             double error_pos, double forster_radius,
                             std::string distance_type)
    : av1_name(av1_name), av2_name(av2_name), distance(distance),
      error_neg(error_neg), error_pos(error_pos),
      forster_radius(forster_radius), distance_type(distance_type) {}

LabelingSite::LabelingSite(int residue_seq_number, std::string atom_name,
                           double distance, double error_neg, double error_pos,
                           double forster_radius)
    : residue_seq_number(residue_seq_number), atom_name(atom_name),
      distance(distance), error_neg(error_neg), error_pos(error_pos),
      forster_radius(forster_radius) {}

double SimpleAVNetworkRestraint::model_distance(const AVMeasurement& m) const {
    std::map<std::string, BasicAV>::const_iterator a1 = avs_.find(m.av1_name);
    std::map<std::string, BasicAV>::const_iterator a2 = avs_.find(m.av2_name);
    if (a1 == avs_.end() || a2 == avs_.end()) {
        IMP_THROW("no accessible volume registered as '"
                      << (a1 == avs_.end() ? m.av1_name : m.av2_name) << "'",
                  ValueException);
    }
    if (m.distance_type == "RDAMean") return a1->second.dRDA(a2->second);
    if (m.distance_type == "RDAMeanE")
        return a1->second.dRDAE(a2->second, m.forster_radius);
    if (m.distance_type == "Rmp") return a1->second.dRmp(a2->second);
    IMP_THROW("unknown distance_type '" << m.distance_type << "'",
              ValueException);
}

double SimpleAVNetworkRestraint::evaluate() const {
    double total = 0.0;
    for (std::size_t i = 0; i < measurements_.size(); ++i) {
        const AVMeasurement& m = measurements_[i];
        total += chi2_score(model_distance(m), m.distance, m.error_neg,
                            m.error_pos);
    }
    return total * weight_;
}

std::map<std::string, double> SimpleAVNetworkRestraint::get_model_distances()
        const {
    std::map<std::string, double> out;
    for (std::size_t i = 0; i < measurements_.size(); ++i) {
        const AVMeasurement& m = measurements_[i];
        out[m.av1_name + "_" + m.av2_name] = model_distance(m);
    }
    return out;
}

DirectLabelingRestraint::DirectLabelingRestraint(double* xyz, int n_atoms,
                                                 int n_dim, double weight)
    : n_atoms_(0), weight_(weight) {
    if (xyz == 0 || n_atoms <= 0) return;
    if (n_dim != 3) {
        IMP_THROW("coordinates must be (N, 3), not (" << n_atoms << ", "
                                                      << n_dim << ")",
                  ValueException);
    }
    n_atoms_ = n_atoms;
    xyz_.assign(xyz, xyz + static_cast<std::size_t>(n_atoms) * 3);
}

double DirectLabelingRestraint::evaluate() const {
    if (n_atoms_ == 0 || sites_.size() < 2) return 0.0;

    double total = 0.0;
    for (std::size_t i = 0; i < sites_.size(); ++i) {
        const LabelingSite& si = sites_[i];
        for (std::size_t j = i + 1; j < sites_.size(); ++j) {
            const LabelingSite& sj = sites_[j];
            if (si.residue_seq_number < 0 || sj.residue_seq_number < 0 ||
                si.residue_seq_number >= n_atoms_ ||
                sj.residue_seq_number >= n_atoms_) {
                IMP_THROW("site residue_seq_number outside the "
                              << n_atoms_ << " coordinates given",
                          ValueException);
            }
            const double* a = &xyz_[static_cast<std::size_t>(
                                        si.residue_seq_number) * 3];
            const double* b = &xyz_[static_cast<std::size_t>(
                                        sj.residue_seq_number) * 3];
            const double dx = a[0] - b[0], dy = a[1] - b[1], dz = a[2] - b[2];
            const double d = std::sqrt(dx * dx + dy * dy + dz * dz);
            // The experimental value is the *first* site's: a pair is encoded
            // as two sites carrying the same measurement, so which one is read
            // does not matter, and reading both would double-count the error.
            total += chi2_score(d, si.distance, si.error_neg, si.error_pos);
        }
    }
    return total * weight_;
}

IMPBFF_END_NAMESPACE
