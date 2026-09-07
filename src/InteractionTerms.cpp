/**
 * \file InteractionTerms.cpp
 * \brief The channels that deactivate an excited dye, one object each.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#include <IMP/bff/InteractionTerms.h>

#include <IMP/bff/FRETRateTrace.h>
#include <IMP/bff/internal/OutputView.h>
#include <IMP/bff/Base.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>

IMPBFF_BEGIN_NAMESPACE

// --------------------------------------------------------------------------
// radiative
// --------------------------------------------------------------------------

RadiativeTerm::RadiativeTerm(double lifetime)
    : InteractionTerm("RadiativeTerm%1%"), lifetime_(lifetime) {}

std::vector<double> RadiativeTerm::rate_constants(const States& first,
                                                  const States&) const {
    if (!(lifetime_ > 0.0)) {
        IMP_THROW("lifetime must be > 0, not " << lifetime_, ValueException);
    }
    return std::vector<double>(static_cast<std::size_t>(first.get_n_points()),
                               1.0 / lifetime_);
}

// --------------------------------------------------------------------------
// PET
// --------------------------------------------------------------------------

PETTerm::PETTerm(const std::map<std::string, PETParameters>& parameters,
                 const std::vector<std::string>& res_names,
                 const std::vector<std::string>& atom_names, double* coords,
                 int n_atoms, int n_dim, double probe_radius)
    : InteractionTerm("PETTerm%1%"), probe_radius_(probe_radius) {
    if (n_dim != 3) {
        IMP_THROW("quencher atoms must be (N, 3), not (" << n_atoms << ", "
                                                         << n_dim << ")",
                  ValueException);
    }
    if (res_names.size() != static_cast<std::size_t>(n_atoms) ||
        atom_names.size() != static_cast<std::size_t>(n_atoms)) {
        IMP_THROW("one residue name and one atom name per atom: "
                          << res_names.size() << ", " << atom_names.size()
                          << " against " << n_atoms << " atoms",
                  ValueException);
    }

    // Which atoms of which residues are redox-active is the one table that
    // defines it; the rate for the pair comes from `parameters`. Resolving both
    // here means the loop below never looks at a name.
    const std::map<std::string, std::vector<std::string> > active =
            quencher_atoms();

    coords_.assign(coords, coords + static_cast<std::size_t>(n_atoms) * 3);
    kQ_.assign(static_cast<std::size_t>(n_atoms), 0.0);
    rC_.assign(static_cast<std::size_t>(n_atoms), 0.0);

    for (int i = 0; i < n_atoms; ++i) {
        std::map<std::string, PETParameters>::const_iterator p =
                parameters.find(res_names[i]);
        if (p == parameters.end()) continue;
        std::map<std::string, std::vector<std::string> >::const_iterator a =
                active.find(res_names[i]);
        if (a == active.end()) continue;
        bool is_active = false;
        for (std::size_t k = 0; k < a->second.size(); ++k) {
            if (a->second[k] == atom_names[i]) { is_active = true; break; }
        }
        if (!is_active) continue;
        kQ_[i] = p->second.rate_constant;
        // An absent or zero attenuation length is a hard contact sphere, which
        // this exponential form spells as a 1 A decay -- the same fallback the
        // Python had, kept because the tabulated values are surface-relative
        // and a zero would divide.
        const double rc = p->second.attenuation_length;
        rC_[i] = (rc == rc && rc != 0.0) ? rc : 1.0;
    }
}

unsigned int PETTerm::get_n_active() const {
    unsigned int n = 0;
    for (std::size_t i = 0; i < kQ_.size(); ++i) {
        if (kQ_[i] > 0.0 && rC_[i] > 0.0) ++n;
    }
    return n;
}

void PETTerm::get_kQ(double** out_view, int* n_out_view) const {
    internal::copy_to_view(kQ_, out_view, n_out_view);
}

void PETTerm::get_rC(double** out_view, int* n_out_view) const {
    internal::copy_to_view(rC_, out_view, n_out_view);
}

std::vector<double> PETTerm::rate_constants(const States& first,
                                            const States&) const {
    const std::size_t n_states = static_cast<std::size_t>(first.get_n_points());
    std::vector<double> out(n_states, 0.0);
    if (n_states == 0) return out;

    double* points = nullptr;
    int n_points = 0;
    first.get_points(&points, &n_points);
    if (points == nullptr) return out;

    const std::size_t n_atoms = kQ_.size();
    for (std::size_t i = 0; i < n_states; ++i) {
        const double x = points[i * 4 + 0];
        const double y = points[i * 4 + 1];
        const double z = points[i * 4 + 2];
        double total = 0.0;
        for (std::size_t a = 0; a < n_atoms; ++a) {
            if (!(kQ_[a] > 0.0) || !(rC_[a] > 0.0)) continue;
            const double dx = coords_[a * 3 + 0] - x;
            const double dy = coords_[a * 3 + 1] - y;
            const double dz = coords_[a * 3 + 2] - z;
            const double d = std::sqrt(dx * dx + dy * dy + dz * dz) - probe_radius_;
            total += kQ_[a] * std::exp(-d / rC_[a]);
        }
        out[i] = total;
    }
    std::free(points);
    return out;
}

// --------------------------------------------------------------------------
// FRET
// --------------------------------------------------------------------------

FRETTerm::FRETTerm(const Probe& donor, const Probe& acceptor,
                   double refractive_index, double kappa2, double r_min)
    : InteractionTerm("FRETTerm%1%"), donor_(donor), acceptor_(acceptor),
      refractive_index_(refractive_index), kappa2_(kappa2), r_min_(r_min) {}

double FRETTerm::get_forster_radius() const {
    // The isotropic value when none was given, so R0 is a number even where
    // kappa^2 is not resolved; `get_used_isotropic_kappa2` reports which.
    const double k2 = (kappa2_ == kappa2_) ? kappa2_ : 2.0 / 3.0;
    return forster_radius(donor_, acceptor_, k2, refractive_index_);  // Angstrom
}

std::vector<double> FRETTerm::rate_constants(const States& first,
                                             const States& second) const {
    const double tau0 = donor_.lifetime;
    if (!(tau0 > 0.0)) {
        IMP_THROW(donor_.name << " has no lifetime, so a FRET rate cannot be "
                                 "expressed as 1/tau0 * (R0/r)^6",
                  ValueException);
    }
    if (second.get_n_points() == 0) {
        IMP_THROW("a FRET term needs the acceptor's states", ValueException);
    }

    double* donor_points = nullptr;
    double* acceptor_points = nullptr;
    int n_donor = 0, n_acceptor = 0;
    first.get_points(&donor_points, &n_donor);
    second.get_points(&acceptor_points, &n_acceptor);

    // The kernel takes (x, y, z) per point; the clouds carry a fourth column.
    std::vector<double> donor_xyz, acceptor_xyz;
    donor_xyz.reserve(static_cast<std::size_t>(n_donor) / 4 * 3);
    for (int i = 0; i + 3 < n_donor; i += 4) {
        donor_xyz.push_back(donor_points[i]);
        donor_xyz.push_back(donor_points[i + 1]);
        donor_xyz.push_back(donor_points[i + 2]);
    }
    acceptor_xyz.reserve(static_cast<std::size_t>(n_acceptor) / 4 * 3);
    for (int i = 0; i + 3 < n_acceptor; i += 4) {
        acceptor_xyz.push_back(acceptor_points[i]);
        acceptor_xyz.push_back(acceptor_points[i + 1]);
        acceptor_xyz.push_back(acceptor_points[i + 2]);
    }
    std::free(donor_points);
    std::free(acceptor_points);

    const double k2 = (kappa2_ == kappa2_) ? kappa2_ : 2.0 / 3.0;
    const double r_min = std::max(r_min_, 1e-6);

    double* rates = nullptr;
    int n_rates = 0;
    fret_rate_trace_kernel(donor_xyz, acceptor_xyz, get_forster_radius(), tau0,
                           r_min * r_min, k2 / (2.0 / 3.0), &rates, &n_rates);
    std::vector<double> out(rates, rates + n_rates);
    std::free(rates);
    return out;
}

// --------------------------------------------------------------------------
// the sum
// --------------------------------------------------------------------------

std::vector<double> total_rate(const InteractionTerms& terms,
                               const States& first, const States& second) {
    if (terms.empty()) {
        IMP_THROW("no interaction terms to sum", ValueException);
    }
    std::vector<double> out;
    for (unsigned int t = 0; t < terms.size(); ++t) {
        const std::vector<double> r = terms[t]->rate_constants(first, second);
        if (out.empty()) {
            out = r;
        } else if (out.size() != r.size()) {
            IMP_THROW("terms disagree about the number of states: " << out.size()
                              << " against " << r.size(),
                      ValueException);
        } else {
            for (std::size_t i = 0; i < out.size(); ++i) out[i] += r[i];
        }
    }
    return out;
}

IMPBFF_END_NAMESPACE
