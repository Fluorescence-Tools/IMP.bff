/**
 * \file Potentials.cpp
 * \brief The coarse-grained protein potentials, as IMP scores and restraints.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#include <IMP/bff/Potentials.h>
#include <IMP/bff/PotentialTables.h>
#include <IMP/bff/SolventAccessibleSurface.h>
#include <IMP/bff/ZMatrix.h>
#include <IMP/bff/internal/OutputView.h>
#include <IMP/bff/internal/Text.h>

#include <IMP/atom/Charged.h>
#include <IMP/atom/Selection.h>
#include <IMP/atom/StereochemistryPairFilter.h>
#include <IMP/atom/bond_graph.h>
#include <IMP/atom/hierarchy_tools.h>
#include <IMP/container/ClosePairContainer.h>
#include <IMP/container/ListSingletonContainer.h>
#include <IMP/container/PairsRestraint.h>
#include <IMP/core/AngleRestraint.h>
#include <IMP/core/DihedralRestraint.h>
#include <IMP/core/DistanceRestraint.h>
#include <IMP/core/Harmonic.h>
#include <IMP/core/SphereDistancePairScore.h>
#include <IMP/core/XYZR.h>
#include <IMP/bff/Base.h>

// for the restraint factories that came from Scoring.cpp
#include <IMP/bff/Mol2IO.h>
#include <IMP/bff/Rotamer.h>
#include <IMP/algebra/vector_generators.h>
#include <IMP/constants.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <IMP/core/HarmonicLowerBound.h>
#include <IMP/container/ListPairContainer.h>
#include <IMP/core/XYZ.h>
#include <IMP/core/internal/dihedral_helpers.h>
#include <cctype>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <algorithm>
#include <cmath>
#include <sstream>

IMPBFF_BEGIN_NAMESPACE

namespace {

inline double distance_of(const std::vector<double>& xyz, std::size_t a,
                          std::size_t b) {
    const double dx = xyz[3 * a + 0] - xyz[3 * b + 0];
    const double dy = xyz[3 * a + 1] - xyz[3 * b + 1];
    const double dz = xyz[3 * a + 2] - xyz[3 * b + 2];
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

inline double squared_distance_of(const std::vector<double>& xyz,
                                  std::size_t a, std::size_t b) {
    const double dx = xyz[3 * a + 0] - xyz[3 * b + 0];
    const double dy = xyz[3 * a + 1] - xyz[3 * b + 1];
    const double dz = xyz[3 * a + 2] - xyz[3 * b + 2];
    return dx * dx + dy * dy + dz * dz;
}

//! The atom of a residue with this type, or an unset handle.
IMP::atom::Atom residue_atom(IMP::atom::Hierarchy residue,
                             IMP::atom::AtomType type) {
    const IMP::atom::Hierarchies children = residue.get_children();
    for (unsigned int i = 0; i < children.size(); ++i) {
        if (!IMP::atom::Atom::get_is_setup(children[i])) continue;
        IMP::atom::Atom a(children[i]);
        if (a.get_atom_type() == type) return a;
    }
    return IMP::atom::Atom();
}

//! Every residue of a hierarchy, in order.
IMP::atom::Hierarchies residues_of(IMP::atom::Hierarchy hierarchy) {
    IMP::atom::Hierarchies out;
    const IMP::atom::Hierarchies all =
            IMP::atom::get_by_type(hierarchy, IMP::atom::RESIDUE_TYPE);
    for (unsigned int i = 0; i < all.size(); ++i) out.push_back(all[i]);
    return out;
}

//! `index` into `pis`, or -1 when the atom is not there.
int push_atom(IMP::atom::Hierarchy residue, IMP::atom::AtomType type,
              IMP::ParticleIndexes& pis) {
    IMP::atom::Atom a = residue_atom(residue, type);
    if (!a) return -1;
    pis.push_back(a.get_particle_index());
    return static_cast<int>(pis.size()) - 1;
}

//! The coordinates behind an index into `pis`.
inline IMP::algebra::Vector3D coordinates_at(IMP::Model* m,
                                             const IMP::ParticleIndexes& pis,
                                             int slot) {
    return IMP::core::XYZ(m, pis[slot]).get_coordinates();
}

}  // namespace

// ---------------------------------------------------------------------------
// Typing the particles
// ---------------------------------------------------------------------------

IMP::IntKey get_residue_type_key() {
    static const IMP::IntKey key("bff residue type");
    return key;
}

void load_residue_contact_types(std::string table, std::string path) {
    const PotentialTable t = read_potential_table(table, path);
    if (t.text.empty()) {
        IMP_THROW("`" << table << "` is not a table with residue names in it",
                  ValueException);
    }
    // The names are the first two fields of every line but the header, and
    // naming a Key is what registers it -- the same thing the PMF reader does
    // as it reads the file.
    std::istringstream in(t.text);
    std::string line;
    bool header = true;
    while (std::getline(in, line)) {
        if (header) {
            header = false;
            continue;
        }
        std::istringstream fields(line);
        std::string a, b;
        if (!(fields >> a >> b)) continue;
        // Naming a Key registers it. Parenthesised, because
        // `ResidueContactType(a);` declares a variable called `a`.
        (void)ResidueContactType(a).get_index();
        (void)ResidueContactType(b).get_index();
    }
}

IMP::ParticlesTemp add_residue_type_score_data(IMP::atom::Hierarchy hierarchy,
                                               IMP::atom::AtomType site,
                                               std::string table) {
    if (!table.empty()) load_residue_contact_types(table);
    const IMP::IntKey key = get_residue_type_key();
    const IMP::atom::Hierarchies residues = residues_of(hierarchy);
    IMP::ParticlesTemp out;
    int n_untyped = 0;
    for (unsigned int i = 0; i < residues.size(); ++i) {
        IMP::atom::Atom a = residue_atom(residues[i], site);
        if (!a) continue;   // a glycine has no C-beta
        const std::string name = IMP::atom::Residue(residues[i])
                                         .get_residue_type()
                                         .get_string();
        // A residue the table says nothing about takes no part: typing it
        // would put an index outside the table's square, which reads as
        // another residue's row or as nothing at all.
        if (!ResidueContactType::get_key_exists(name)) {
            ++n_untyped;
            continue;
        }
        const int type = ResidueContactType(name).get_index();
        IMP::Particle* p = a.get_particle();
        if (p->has_attribute(key)) {
            p->set_value(key, type);
        } else {
            p->add_attribute(key, type);
        }
        out.push_back(p);
    }
    if (out.empty() && n_untyped > 0) {
        IMP_THROW("none of the " << n_untyped << " residues has a type the "
                                 << "contact tables know; load a table first "
                                 << "with load_residue_contact_types()",
                  ValueException);
    }
    return out;
}

// ---------------------------------------------------------------------------
// Typed contact potentials
// ---------------------------------------------------------------------------

namespace {

//! The shipped table of this name, as something IMP's PMF reader can read.
/*! `IMP::TextInput` over a stream *refers* to it, so the stream has to outlive
    the read. A temporary of this type in a member initializer lives to the end
    of that initializer -- which is after the base constructor has read the
    table -- so the conversion below is safe and a bare
    `TextInput(std::istringstream(...))` would not be. */
struct ShippedPmf {
    std::string name, text;
    std::istringstream stream;

    explicit ShippedPmf(const std::string& name)
        : name(name), text(read_potential_table(name).text), stream(text) {
        if (text.empty()) {
            IMP_THROW("the shipped container has no text for `"
                              << name << "`; pass a table of your own",
                      IOException);
        }
    }

    operator IMP::TextInput() { return IMP::TextInput(stream, name); }
};

}  // namespace

MiyazawaJerniganPairScore::MiyazawaJerniganPairScore(double threshold)
    : P(get_residue_type_key(), threshold, ShippedPmf("mj")) {}

MiyazawaJerniganPairScore::MiyazawaJerniganPairScore(double threshold,
                                                     IMP::TextInput data_file)
    : P(get_residue_type_key(), threshold, data_file) {}

UNRESCentroidPairScore::UNRESCentroidPairScore(double threshold)
    : P(get_residue_type_key(), threshold, ShippedPmf("unres")) {}

UNRESCentroidPairScore::UNRESCentroidPairScore(double threshold,
                                               IMP::TextInput data_file)
    : P(get_residue_type_key(), threshold, data_file) {}

// ---------------------------------------------------------------------------
// LennardJonesBeadPairScore
// ---------------------------------------------------------------------------

LennardJonesBeadPairScore::LennardJonesBeadPairScore(double rm)
    : IMP::PairScore("LennardJonesBeadPairScore%1%"), rm_(rm) {}

double LennardJonesBeadPairScore::evaluate_index(
        IMP::Model* m, const IMP::ParticleIndexPair& p,
        IMP::DerivativeAccumulator* da) const {
    IMP::core::XYZ a(m, p[0]), b(m, p[1]);
    const IMP::algebra::Vector3D delta =
            a.get_coordinates() - b.get_coordinates();
    const double r2 = delta.get_squared_magnitude();
    if (r2 == 0.0) return 0.0;
    const double s = std::pow(rm_ * rm_ / r2, 3);
    const double score = s * (s - 2.0);
    if (da != nullptr) {
        // dE/dr = -12 s (s - 1) / r, along the separation.
        const double r = std::sqrt(r2);
        const IMP::algebra::Vector3D grad =
                delta * (-12.0 * s * (s - 1.0) / (r * r));
        a.add_to_derivatives(grad, *da);
        b.add_to_derivatives(-grad, *da);
    }
    return score;
}

IMP::ModelObjectsTemp LennardJonesBeadPairScore::do_get_inputs(
        IMP::Model* m, const IMP::ParticleIndexes& pis) const {
    return IMP::get_particles(m, pis);
}

// ---------------------------------------------------------------------------
// The kernels
// ---------------------------------------------------------------------------

double clash_energy(const std::vector<double>& xyz,
                    const std::vector<double>& vdw, double clash_tolerance,
                    double covalent_radius) {
    if (clash_tolerance == 0.0) {
        IMP_THROW("clash_tolerance divides the overlap and cannot be zero",
                  ValueException);
    }
    const std::size_t n = std::min(xyz.size() / 3, vdw.size());
    double e = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = i + 1; j < n; ++j) {
            const double rij = distance_of(xyz, i, j);
            const double sigma = vdw[i] + vdw[j];
            if (rij > covalent_radius && rij < sigma) {
                const double overlap = (sigma - rij) / clash_tolerance;
                e += overlap * overlap;
            }
        }
    }
    return e;
}

void GoContacts::get_energies(double** out_view, int* n_out_view) const {
    internal::copy_to_view(energies, out_view, n_out_view);
}

void GoContacts::get_distances(double** out_view, int* n_out_view) const {
    internal::copy_to_view(distances, out_view, n_out_view);
}

GoContacts go_native_contacts(const std::vector<double>& xyz, double epsilon,
                              double nn_e_factor, double cutoff) {
    GoContacts out;
    const std::size_t n = xyz.size() / 3;
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = i + 1; j < n; ++j) {
            const double rm = distance_of(xyz, i, j);
            // The truncation the Go energy assumes: the well is cut at
            // 2.5 sigma, and 17/55 + 1 is rm/sigma for that cut.
            const double int_e = (17.0 / 55.0 + 1.0) * rm;
            const bool native = rm < cutoff;
            if (native) ++out.n_native;
            out.residue_i.push_back(static_cast<int>(i));
            out.residue_j.push_back(static_cast<int>(j));
            out.distances.push_back(rm);
            out.energies.push_back(
                    int_e == 0.0 ? 0.0
                                 : (native ? epsilon : nn_e_factor * epsilon) /
                                           int_e);
        }
    }
    return out;
}

namespace {

//! The truncated-LJ energy of one contact at one distance.
inline double go_pair_energy(double epsilon, double rm, double r) {
    if (r > rm && r < rm * 2.5 && r > 0.0) {
        const double sr = 2.0 * std::pow(rm / r, 6);
        // The 0.00818 shift puts the truncated well at zero at 2.5 rm.
        return epsilon * (-sr + sr * sr / 4.0 + 0.00818);
    }
    return -epsilon;
}

}  // namespace

double go_energy(const std::vector<double>& xyz, const GoContacts& contacts) {
    double total = 0.0;
    for (std::size_t k = 0; k < contacts.energies.size(); ++k) {
        const double r = distance_of(
                xyz, static_cast<std::size_t>(contacts.residue_i[k]),
                static_cast<std::size_t>(contacts.residue_j[k]));
        total += go_pair_energy(contacts.energies[k], contacts.distances[k], r);
    }
    return total;
}

double lennard_jones_bead_energy(const std::vector<double>& xyz, double rm) {
    const std::size_t n = xyz.size() / 3;
    const double rm2 = rm * rm;
    double e = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = i + 1; j < n; ++j) {
            const double r2 = squared_distance_of(xyz, i, j);
            if (r2 == 0.0) continue;
            const double s = std::pow(rm2 / r2, 3);
            e += s * (s - 2.0);
        }
    }
    return e;
}

double generalized_born_energy(const std::vector<double>& xyz,
                               const std::vector<double>& radii,
                               const std::vector<double>& charges,
                               double epsilon, double epsilon0,
                               double cutoff) {
    const std::size_t n =
            std::min(xyz.size() / 3, std::min(radii.size(), charges.size()));
    const double cutoff2 = cutoff * cutoff;
    const double pre = 1.0 / (8.0 * M_PI) * (1.0 / epsilon0 - 1.0 / epsilon);
    double energy = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const double qi = charges[i];
        if (qi == 0.0) continue;
        for (std::size_t j = i; j < n; ++j) {
            const double qj = charges[j];
            if (qj == 0.0) continue;
            const double rij2 = squared_distance_of(xyz, i, j);
            if (rij2 > cutoff2) continue;
            const double aij2 = radii[i] * radii[j];
            if (aij2 == 0.0) continue;
            const double qij = qi * qj;
            const double d = rij2 / (4.0 * aij2);
            const double fgbr = std::sqrt(rij2 + aij2 * std::exp(-d));
            // `cutoff`, not `cutoff * cutoff`: see the note in the header.
            const double fgbc = std::sqrt(cutoff + aij2 * std::exp(-d));
            energy += qij / fgbr - qij / fgbc;
        }
    }
    return energy * pre;
}

double residue_asa(const std::vector<double>& xyz, int n_sphere, double probe,
                   double radius) {
    const std::size_t n = xyz.size() / 3;
    if (n == 0) return 0.0;
    std::vector<int> indices(n);
    for (std::size_t i = 0; i < n; ++i) indices[i] = static_cast<int>(i);
    const std::vector<double> radii(n, radius);
    const std::vector<double> points = sphere_points(n_sphere);
    const std::vector<double> areas = solvent_accessible_surface_area(
            xyz, radii, indices, points, probe, radius);
    double total = 0.0;
    for (std::size_t i = 0; i < areas.size(); ++i) total += areas[i];
    return total;
}

double ramachandran_energy(double phi, double psi,
                           const std::vector<double>& grid, int n_bins,
                           double empty_penalty) {
    if (n_bins <= 0 ||
        grid.size() < static_cast<std::size_t>(n_bins) * n_bins) {
        IMP_THROW("the Ramachandran grid needs "
                          << n_bins << " squared numbers and has "
                          << grid.size(),
                  ValueException);
    }
    if (phi != phi || psi != psi) return 0.0;   // an end of the chain

    const double resolution = 360.0 / static_cast<double>(n_bins);
    int ix = static_cast<int>((phi * 180.0 / M_PI + 180.0) / resolution);
    int iy = static_cast<int>((psi * 180.0 / M_PI + 180.0) / resolution);
    ix = std::max(0, std::min(ix, n_bins - 1));
    iy = std::max(0, std::min(iy, n_bins - 1));
    const double value = grid[static_cast<std::size_t>(ix) * n_bins + iy];

    double maximum = grid[0];
    for (std::size_t i = 1; i < grid.size(); ++i) {
        maximum = std::max(maximum, grid[i]);
    }
    // A map that never goes positive is holding pseudo-energies already.
    if (maximum <= 0.0) return -value;
    if (value <= 0.0) return empty_penalty;
    return -std::log(value / maximum);
}

// ---------------------------------------------------------------------------
// GoRestraint
// ---------------------------------------------------------------------------

GoRestraint::GoRestraint(IMP::Model* m, const IMP::ParticleIndexes& pis,
                         double epsilon, double cutoff, double nn_e_factor)
    : IMP::Restraint(m, "GoRestraint%1%"), pis_(pis), epsilon_(epsilon),
      cutoff_(cutoff), nn_e_factor_(nn_e_factor) {
    set_native_contacts();
}

namespace {

//! The beads' coordinates, flat.
std::vector<double> bead_coords(IMP::Model* m,
                                const IMP::ParticleIndexes& pis) {
    std::vector<double> xyz(pis.size() * 3, 0.0);
    for (unsigned int i = 0; i < pis.size(); ++i) {
        const IMP::algebra::Vector3D v =
                IMP::core::XYZ(m, pis[i]).get_coordinates();
        xyz[3 * i + 0] = v[0];
        xyz[3 * i + 1] = v[1];
        xyz[3 * i + 2] = v[2];
    }
    return xyz;
}

}  // namespace

void GoRestraint::set_native_contacts() {
    contacts_ = go_native_contacts(bead_coords(get_model(), pis_), epsilon_,
                                   nn_e_factor_, cutoff_);
}

int GoRestraint::get_n_non_native() const {
    return contacts_.get_n_contacts() - contacts_.n_native;
}

double GoRestraint::unprotected_evaluate(
        IMP::DerivativeAccumulator* accum) const {
    // No derivatives, like IMP::atom::CAAngleRestraint: the term is a
    // knowledge-based well over a fixed contact list, and every caller of it
    // so far samples rather than minimises.
    IMP_UNUSED(accum);
    return go_energy(bead_coords(get_model(), pis_), contacts_);
}

IMP::ModelObjectsTemp GoRestraint::do_get_inputs() const {
    return IMP::get_particles(get_model(), pis_);
}

// ---------------------------------------------------------------------------
// HydrogenBondRestraint
// ---------------------------------------------------------------------------

HydrogenBondRestraint::HydrogenBondRestraint(IMP::Model* m,
                                             IMP::atom::Hierarchy hierarchy,
                                             const std::vector<double>& table,
                                             int n_bins, double cutoff_ca,
                                             double cutoff_h,
                                             double bin_width)
    : IMP::Restraint(m, "HydrogenBondRestraint%1%"), table_(table),
      n_bins_(n_bins), cutoff_ca_(cutoff_ca), cutoff_h_(cutoff_h),
      bin_width_(bin_width), ch_(true), on_(true), oh_(true), cn_(true),
      n_hbonds_(0) {
    if (table_.empty()) {
        IMP_THROW("the hydrogen-bond lookup is empty", ValueException);
    }
    if (n_bins_ <= 0) n_bins_ = static_cast<int>(table_.size() / 4);
    if (table_.size() < static_cast<std::size_t>(4) * n_bins_) {
        IMP_THROW("the lookup needs 4 channels of " << n_bins_ << " bins and "
                                                    << "has " << table_.size()
                                                    << " numbers",
                  ValueException);
    }

    const IMP::atom::Hierarchies residues = residues_of(hierarchy);
    for (unsigned int i = 0; i < residues.size(); ++i) {
        // A residue takes part only if it has a C-alpha to measure the
        // prefilter with; the four bonding atoms may each be missing.
        IMP::atom::Atom ca = residue_atom(residues[i], IMP::atom::AT_CA);
        if (!ca) continue;
        pis_.push_back(ca.get_particle_index());
        n_.push_back(push_atom(residues[i], IMP::atom::AT_N, pis_));
        c_.push_back(push_atom(residues[i], IMP::atom::AT_C, pis_));
        o_.push_back(push_atom(residues[i], IMP::atom::AT_O, pis_));
        h_.push_back(push_atom(residues[i], IMP::atom::AT_H, pis_));
    }
}

void HydrogenBondRestraint::set_channels(bool ch, bool on, bool oh, bool cn) {
    ch_ = ch;
    on_ = on;
    oh_ = oh;
    cn_ = cn;
}

namespace {

//! One channel of the lookup at a distance, or zero past its end.
inline double hbond_bin(const std::vector<double>& table, int n_bins,
                        int channel, double dist, double bin_width) {
    const int b = static_cast<int>(dist / bin_width);
    if (b < 0 || b >= n_bins) return 0.0;
    return table[static_cast<std::size_t>(channel) * n_bins + b];
}

}  // namespace

double HydrogenBondRestraint::unprotected_evaluate(
        IMP::DerivativeAccumulator* accum) const {
    IMP_UNUSED(accum);   // a binned lookup has no useful gradient
    IMP::Model* m = get_model();
    const int n = static_cast<int>(n_.size());
    const double cutoff_h2 = cutoff_h_ * cutoff_h_;
    double e = 0.0;
    n_hbonds_ = 0;

    // The C-alpha of residue k is the entry of pis_ that precedes its four
    // bonding atoms; it was pushed first, so its slot is recoverable from the
    // count of slots taken before it. Rather than recompute that, the walk
    // below keeps the coordinates it needs.
    std::vector<IMP::algebra::Vector3D> ca(n);
    {
        int slot = 0;
        for (int k = 0; k < n; ++k) {
            ca[k] = IMP::core::XYZ(m, pis_[slot]).get_coordinates();
            // one C-alpha plus however many of N, C, O, H were present
            slot += 1 + (n_[k] >= 0) + (c_[k] >= 0) + (o_[k] >= 0) +
                    (h_[k] >= 0);
        }
    }

    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            // Applied, unlike the kernel this came from: see the header.
            if (IMP::algebra::get_distance(ca[i], ca[j]) > cutoff_ca_) {
                continue;
            }
            for (int direction = 0; direction < 2; ++direction) {
                // donor -> acceptor: the donor lends its amide H to the
                // acceptor's carbonyl O.
                const int donor = direction == 0 ? j : i;
                const int acceptor = direction == 0 ? i : j;
                if (h_[donor] < 0 || n_[donor] < 0 || o_[acceptor] < 0 ||
                    c_[acceptor] < 0) {
                    continue;
                }
                const IMP::algebra::Vector3D h =
                        coordinates_at(m, pis_, h_[donor]);
                const IMP::algebra::Vector3D o =
                        coordinates_at(m, pis_, o_[acceptor]);
                if ((h - o).get_squared_magnitude() >= cutoff_h2) continue;
                ++n_hbonds_;
                const IMP::algebra::Vector3D nn =
                        coordinates_at(m, pis_, n_[donor]);
                const IMP::algebra::Vector3D cc =
                        coordinates_at(m, pis_, c_[acceptor]);
                if (oh_) {
                    e += hbond_bin(table_, n_bins_, 2,
                                   IMP::algebra::get_distance(h, o),
                                   bin_width_);
                }
                if (on_) {
                    e += hbond_bin(table_, n_bins_, 1,
                                   IMP::algebra::get_distance(nn, o),
                                   bin_width_);
                }
                if (ch_) {
                    e += hbond_bin(table_, n_bins_, 0,
                                   IMP::algebra::get_distance(cc, h),
                                   bin_width_);
                }
                if (cn_) {
                    e += hbond_bin(table_, n_bins_, 3,
                                   IMP::algebra::get_distance(nn, cc),
                                   bin_width_);
                }
            }
        }
    }
    return e;
}

IMP::ModelObjectsTemp HydrogenBondRestraint::do_get_inputs() const {
    return IMP::get_particles(get_model(), pis_);
}

// ---------------------------------------------------------------------------
// GeneralizedBornRestraint
// ---------------------------------------------------------------------------

GeneralizedBornRestraint::GeneralizedBornRestraint(
        IMP::Model* m, const IMP::ParticleIndexes& pis, double epsilon,
        double epsilon0, double cutoff)
    : IMP::Restraint(m, "GeneralizedBornRestraint%1%"), pis_(pis),
      epsilon_(epsilon), epsilon0_(epsilon0), cutoff_(cutoff) {}

double GeneralizedBornRestraint::unprotected_evaluate(
        IMP::DerivativeAccumulator* accum) const {
    IMP_UNUSED(accum);
    IMP::Model* m = get_model();
    std::vector<double> xyz(pis_.size() * 3, 0.0);
    std::vector<double> radii(pis_.size(), 0.0), charges(pis_.size(), 0.0);
    for (unsigned int i = 0; i < pis_.size(); ++i) {
        const IMP::algebra::Vector3D v =
                IMP::core::XYZ(m, pis_[i]).get_coordinates();
        xyz[3 * i + 0] = v[0];
        xyz[3 * i + 1] = v[1];
        xyz[3 * i + 2] = v[2];
        radii[i] = IMP::core::XYZR::get_is_setup(m, pis_[i])
                           ? IMP::core::XYZR(m, pis_[i]).get_radius()
                           : 0.0;
        charges[i] = IMP::atom::Charged::get_is_setup(m, pis_[i])
                             ? IMP::atom::Charged(m, pis_[i]).get_charge()
                             : 0.0;
    }
    return generalized_born_energy(xyz, radii, charges, epsilon_, epsilon0_,
                                   cutoff_);
}

IMP::ModelObjectsTemp GeneralizedBornRestraint::do_get_inputs() const {
    return IMP::get_particles(get_model(), pis_);
}

// ---------------------------------------------------------------------------
// RamachandranRestraint
// ---------------------------------------------------------------------------

RamachandranRestraint::RamachandranRestraint(
        IMP::Model* m, IMP::atom::Hierarchy hierarchy,
        const std::vector<double>& grids, int n_channels, int n_bins,
        double empty_penalty)
    : IMP::Restraint(m, "RamachandranRestraint%1%"), grids_(grids),
      n_channels_(n_channels), n_bins_(n_bins), empty_penalty_(empty_penalty) {
    if (grids_.empty()) IMP_THROW("the Ramachandran map is empty", ValueException);
    if (n_channels_ <= 0) n_channels_ = 1;
    if (n_bins_ <= 0) {
        n_bins_ = static_cast<int>(
                std::sqrt(static_cast<double>(grids_.size() / n_channels_)) +
                0.5);
    }

    const IMP::atom::Hierarchies residues = residues_of(hierarchy);
    const IMP::atom::ResidueType pro = IMP::atom::PRO;
    const IMP::atom::ResidueType gly = IMP::atom::GLY;
    for (unsigned int i = 0; i < residues.size(); ++i) {
        n_.push_back(push_atom(residues[i], IMP::atom::AT_N, pis_));
        ca_.push_back(push_atom(residues[i], IMP::atom::AT_CA, pis_));
        c_.push_back(push_atom(residues[i], IMP::atom::AT_C, pis_));
        const IMP::atom::ResidueType t =
                IMP::atom::Residue(residues[i]).get_residue_type();
        // The shipped map is general, proline, glycine.
        int channel = 0;
        if (n_channels_ >= 3) {
            if (t == pro) channel = 1;
            else if (t == gly) channel = 2;
        }
        channel_.push_back(channel);
    }
    // phi needs the previous residue's C and psi the next residue's N.
    const int n = static_cast<int>(n_.size());
    c_prev_.assign(static_cast<std::size_t>(std::max(0, n)), -1);
    n_next_.assign(static_cast<std::size_t>(std::max(0, n)), -1);
    for (int i = 0; i < n; ++i) {
        if (i > 0) c_prev_[i] = c_[i - 1];
        if (i + 1 < n) n_next_[i] = n_[i + 1];
    }
}

namespace {

//! The two backbone dihedrals of every residue, radians; NaN at the ends.
void backbone_dihedrals_of(IMP::Model* m, const IMP::ParticleIndexes& pis,
                           const std::vector<int>& c_prev,
                           const std::vector<int>& n, const std::vector<int>& ca,
                           const std::vector<int>& c,
                           const std::vector<int>& n_next,
                           std::vector<double>& phi, std::vector<double>& psi) {
    const std::size_t count = n.size();
    phi.assign(count, internal::nan_value());
    psi.assign(count, internal::nan_value());
    const double to_rad = M_PI / 180.0;
    for (std::size_t i = 0; i < count; ++i) {
        if (n[i] < 0 || ca[i] < 0 || c[i] < 0) continue;
        const IMP::algebra::Vector3D nn = coordinates_at(m, pis, n[i]);
        const IMP::algebra::Vector3D a = coordinates_at(m, pis, ca[i]);
        const IMP::algebra::Vector3D cc = coordinates_at(m, pis, c[i]);
        if (c_prev[i] >= 0) {
            phi[i] = dihedral_deg(coordinates_at(m, pis, c_prev[i]), nn, a, cc) *
                     to_rad;
        }
        if (n_next[i] >= 0) {
            psi[i] = dihedral_deg(nn, a, cc, coordinates_at(m, pis, n_next[i])) *
                     to_rad;
        }
    }
}

}  // namespace

void RamachandranRestraint::get_phi(double** out_view,
                                    int* n_out_view) const {
    std::vector<double> phi, psi;
    backbone_dihedrals_of(get_model(), pis_, c_prev_, n_, ca_, c_, n_next_,
                          phi, psi);
    internal::copy_to_view(phi, out_view, n_out_view);
}

void RamachandranRestraint::get_psi(double** out_view,
                                    int* n_out_view) const {
    std::vector<double> phi, psi;
    backbone_dihedrals_of(get_model(), pis_, c_prev_, n_, ca_, c_, n_next_,
                          phi, psi);
    internal::copy_to_view(psi, out_view, n_out_view);
}

double RamachandranRestraint::unprotected_evaluate(
        IMP::DerivativeAccumulator* accum) const {
    IMP_UNUSED(accum);
    std::vector<double> phi, psi;
    backbone_dihedrals_of(get_model(), pis_, c_prev_, n_, ca_, c_, n_next_,
                          phi, psi);
    const std::size_t stride = static_cast<std::size_t>(n_bins_) * n_bins_;
    double e = 0.0;
    for (std::size_t i = 0; i < phi.size(); ++i) {
        if (phi[i] != phi[i] || psi[i] != psi[i]) continue;
        const std::size_t base = static_cast<std::size_t>(channel_[i]) * stride;
        if (base + stride > grids_.size()) continue;
        const std::vector<double> grid(
                grids_.begin() + static_cast<std::ptrdiff_t>(base),
                grids_.begin() + static_cast<std::ptrdiff_t>(base + stride));
        e += ramachandran_energy(phi[i], psi[i], grid, n_bins_, empty_penalty_);
    }
    return e;
}

IMP::ModelObjectsTemp RamachandranRestraint::do_get_inputs() const {
    return IMP::get_particles(get_model(), pis_);
}

// ---------------------------------------------------------------------------
// Built out of IMP's own scores
// ---------------------------------------------------------------------------

IMP::Restraint* create_clash_restraint(IMP::atom::Hierarchy hierarchy,
                                      double clash_tolerance, double slack) {
    if (clash_tolerance == 0.0) {
        IMP_THROW("clash_tolerance divides the overlap and cannot be zero",
                  ValueException);
    }
    IMP::Model* m = hierarchy.get_model();
    const IMP::atom::Hierarchies leaves = IMP::atom::get_leaves(hierarchy);
    IMP::ParticleIndexes pis;
    for (unsigned int i = 0; i < leaves.size(); ++i) {
        if (IMP::core::XYZR::get_is_setup(leaves[i])) {
            pis.push_back(leaves[i].get_particle_index());
        }
    }
    IMP_NEW(IMP::container::ListSingletonContainer, lsc, (m, pis));
    IMP_NEW(IMP::container::ClosePairContainer, cpc, (lsc, 0.0, slack));

    // "Do not count pairs closer than a covalent radius" is a distance guess at
    // "do not count bonded pairs"; IMP has the bond graph, so it can say so.
    const IMP::atom::Bonds bonds = IMP::atom::get_internal_bonds(hierarchy);
    if (!bonds.empty()) {
        IMP_NEW(IMP::atom::StereochemistryPairFilter, filter, ());
        filter->set_bonds(bonds);
        cpc->add_pair_filter(filter);
    }

    // SoftSphere is 0.5 k (sigma - r)^2 and the ported term is
    // ((sigma - r)/t)^2, so k = 2/t^2.
    const double k = 2.0 / (clash_tolerance * clash_tolerance);
    IMP_NEW(IMP::core::SoftSpherePairScore, score, (k));
    IMP_NEW(IMP::container::PairsRestraint, r, (score, cpc, "clash"));
    return r.release();
}

IMP::Restraints create_ca_internal_restraints(
        IMP::Model* m, const IMP::ParticleIndexes& cas,
        const IMP::ParticleIndexes& reference, double k_bond, double k_angle,
        double k_dihedral) {
    if (cas.size() != reference.size()) {
        IMP_THROW("the trace has " << cas.size() << " C-alphas and the "
                                   << "reference " << reference.size(),
                  ValueException);
    }
    IMP::Restraints out;
    const std::size_t n = cas.size();
    if (k_bond != 0.0) {
        for (std::size_t i = 0; i + 1 < n; ++i) {
            const double x0 = IMP::core::get_distance(
                    IMP::core::XYZ(m, reference[i]),
                    IMP::core::XYZ(m, reference[i + 1]));
            IMP_NEW(IMP::core::Harmonic, f, (x0, k_bond));
            out.push_back(new IMP::core::DistanceRestraint(m, f, cas[i],
                                                           cas[i + 1]));
        }
    }
    if (k_angle != 0.0) {
        for (std::size_t i = 0; i + 2 < n; ++i) {
            const double a0 = bond_angle_rad(
                    IMP::core::XYZ(m, reference[i]).get_coordinates(),
                    IMP::core::XYZ(m, reference[i + 1]).get_coordinates(),
                    IMP::core::XYZ(m, reference[i + 2]).get_coordinates());
            IMP_NEW(IMP::core::Harmonic, f, (a0, k_angle));
            out.push_back(new IMP::core::AngleRestraint(
                    m, f, cas[i], cas[i + 1], cas[i + 2]));
        }
    }
    if (k_dihedral != 0.0) {
        for (std::size_t i = 0; i + 3 < n; ++i) {
            const double d0 =
                    dihedral_deg(
                            IMP::core::XYZ(m, reference[i]).get_coordinates(),
                            IMP::core::XYZ(m, reference[i + 1])
                                    .get_coordinates(),
                            IMP::core::XYZ(m, reference[i + 2])
                                    .get_coordinates(),
                            IMP::core::XYZ(m, reference[i + 3])
                                    .get_coordinates()) *
                    M_PI / 180.0;
            IMP_NEW(IMP::core::Harmonic, f, (d0, k_dihedral));
            out.push_back(new IMP::core::DihedralRestraint(
                    m, f, cas[i], cas[i + 1], cas[i + 2], cas[i + 3]));
        }
    }
    return out;
}

double residue_solvent_accessible_surface(IMP::atom::Hierarchy hierarchy,
                                          IMP::atom::AtomType site,
                                          int n_sphere, double probe,
                                          double radius) {
    const IMP::atom::Hierarchies residues = residues_of(hierarchy);
    std::vector<double> centres;
    for (unsigned int i = 0; i < residues.size(); ++i) {
        IMP::atom::Atom a = residue_atom(residues[i], site);
        if (!a) continue;
        const IMP::algebra::Vector3D v =
                IMP::core::XYZ(a).get_coordinates();
        centres.push_back(v[0]);
        centres.push_back(v[1]);
        centres.push_back(v[2]);
    }
    return residue_asa(centres, n_sphere, probe, radius);
}

// -------- from Scoring.cpp (the restraint factories over a typed dye system) --------

IMP::core::Cosine* torsion_cosine(const FFTorsionType& type) {
    // CHARMM's k(1 + cos(n phi - delta)) against Cosine's k(1 - cos(...)):
    // the same curve with the phase shifted by pi.
    return new IMP::core::Cosine(type.k, type.periodicity,
                                 type.phase + IMP::algebra::PI);
}

IMP::Restraints create_probe_restraints(
        IMP::Model* model, const ProbeForceFieldSystem& system,
        const std::vector<std::string>& site_ids,
        const IMP::ParticleIndexes& particles, bool nonbonded) {
    if (site_ids.size() != particles.size()) {
        IMP_THROW("create_probe_restraints: " << site_ids.size() << " site ids "
                  << "against " << particles.size() << " particles",
                  IMP::ValueException);
    }
    std::map<std::string, IMP::ParticleIndex> site;
    for (std::size_t i = 0; i < site_ids.size(); ++i) {
        site[site_ids[i]] = particles[i];
    }
    const bool has_all = true;
    IMP::Restraints out;

    // -- bonds ------------------------------------------------------------
    const std::map<std::string, double>& bt = system.get_bond_types();
    for (std::size_t i = 0; i < system.get_bonds().size(); ++i) {
        const FFBond& b = system.get_bonds()[i];
        std::map<std::string, IMP::ParticleIndex>::const_iterator a =
                site.find(b.site_a), c = site.find(b.site_b);
        if (a == site.end() || c == site.end()) continue;
        std::map<std::string, double>::const_iterator k = bt.find(b.type_id);
        if (k == bt.end()) continue;
        // A length of zero is a length nobody set: restrain about the
        // geometry as it stands rather than pulling the sites together.
        const double length =
                b.length > 0.0
                        ? b.length
                        : IMP::core::get_distance(
                                  IMP::core::XYZ(model, a->second),
                                  IMP::core::XYZ(model, c->second));
        out.push_back(new IMP::core::DistanceRestraint(
                model, new IMP::core::Harmonic(length, k->second),
                model->get_particle(a->second),
                model->get_particle(c->second)));
    }

    // -- angles -----------------------------------------------------------
    const std::map<std::string, double>& at = system.get_angle_types();
    for (std::size_t i = 0; i < system.get_angles().size(); ++i) {
        const FFAngle& an = system.get_angles()[i];
        std::map<std::string, IMP::ParticleIndex>::const_iterator a =
                site.find(an.site_a), b = site.find(an.site_b),
                c = site.find(an.site_c);
        if (a == site.end() || b == site.end() || c == site.end()) continue;
        std::map<std::string, double>::const_iterator k = at.find(an.type_id);
        if (k == at.end()) continue;
        const double theta =
                an.theta > 0.0
                        ? an.theta
                        : bond_angle_rad(
                                  IMP::core::XYZ(model, a->second)
                                          .get_coordinates(),
                                  IMP::core::XYZ(model, b->second)
                                          .get_coordinates(),
                                  IMP::core::XYZ(model, c->second)
                                          .get_coordinates());
        out.push_back(new IMP::core::AngleRestraint(
                model, new IMP::core::Harmonic(theta, k->second),
                model->get_particle(a->second),
                model->get_particle(b->second),
                model->get_particle(c->second)));
    }

    // -- torsions ---------------------------------------------------------
    const std::map<std::string, FFTorsionType>& tt = system.get_torsion_types();
    for (std::size_t i = 0; i < system.get_dihedrals().size(); ++i) {
        const FFTorsion& t = system.get_dihedrals()[i];
        std::map<std::string, IMP::ParticleIndex>::const_iterator
                a = site.find(t.site_a), b = site.find(t.site_b),
                c = site.find(t.site_c), d = site.find(t.site_d);
        if (a == site.end() || b == site.end() || c == site.end() ||
            d == site.end()) continue;
        std::map<std::string, FFTorsionType>::const_iterator ty =
                tt.find(t.type_id);
        if (ty == tt.end()) continue;
        out.push_back(new IMP::core::DihedralRestraint(
                model, torsion_cosine(ty->second),
                model->get_particle(a->second),
                model->get_particle(b->second),
                model->get_particle(c->second),
                model->get_particle(d->second)));
    }

    // -- impropers: harmonic about the geometry as it stands now ----------
    const std::map<std::string, FFTorsionType>& it =
            system.get_improper_types();
    for (std::size_t i = 0; i < system.get_impropers().size(); ++i) {
        const FFTorsion& t = system.get_impropers()[i];
        std::map<std::string, IMP::ParticleIndex>::const_iterator
                a = site.find(t.site_a), b = site.find(t.site_b),
                c = site.find(t.site_c), d = site.find(t.site_d);
        if (a == site.end() || b == site.end() || c == site.end() ||
            d == site.end()) continue;
        std::map<std::string, FFTorsionType>::const_iterator ty =
                it.find(t.type_id);
        if (ty == it.end()) continue;
        const IMP::core::XYZ xa(model, a->second), xb(model, b->second),
                xc(model, c->second), xd(model, d->second);
        const double theta0 = IMP::core::get_dihedral(xa, xb, xc, xd);
        out.push_back(new IMP::core::DihedralRestraint(
                model, new IMP::core::Harmonic(theta0, ty->second.k),
                model->get_particle(a->second),
                model->get_particle(b->second),
                model->get_particle(c->second),
                model->get_particle(d->second)));
    }

    // -- repulsion --------------------------------------------------------
    // One soft-sphere restraint over every non-excluded pair: the same term a
    // Monte-Carlo step is scored against, and differentiable, so a dynamics
    // run uses it too. It was per-pair Lennard-Jones lower bounds here and
    // soft spheres there -- two implementations of one piece of physics, and
    // thousands of restraints where one does.
    if (nonbonded) {
        IMP::Restraint* steric =
                create_steric_restraint(model, system, site_ids, particles);
        if (steric != NULL) out.push_back(steric);
    }
    (void)has_all;
    return out;
}

namespace {

//! The mean position of a set of particles.
IMP::algebra::Vector3D centre_of(IMP::Model* model,
                                 const IMP::ParticleIndexes& ps) {
    IMP::algebra::Vector3D c(0.0, 0.0, 0.0);
    for (std::size_t i = 0; i < ps.size(); ++i) {
        c += IMP::core::XYZ(model, ps[i]).get_coordinates();
    }
    return ps.empty() ? c : c / static_cast<double>(ps.size());
}

}  // namespace

double place_guest_by_score(IMP::ScoringFunction* scoring_function,
                            IMP::Model* model,
                            const IMP::ParticleIndexes& host,
                            const IMP::ParticleIndexes& guest, double distance,
                            int n_trials, int seed) {
    if (guest.empty() || scoring_function == NULL) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    const IMP::algebra::Vector3D host_centre = centre_of(model, host);
    const IMP::algebra::Vector3D guest_centre = centre_of(model, guest);
    IMP::algebra::Vector3Ds shape;
    for (std::size_t i = 0; i < guest.size(); ++i) {
        shape.push_back(IMP::core::XYZ(model, guest[i]).get_coordinates() -
                        guest_centre);
    }

    boost::mt19937 rng(static_cast<boost::uint32_t>(seed));
    boost::uniform_real<double> unit(0.0, 1.0);

    IMP::algebra::Vector3Ds best = shape;
    IMP::algebra::Vector3D best_offset = host_centre;
    double best_score = std::numeric_limits<double>::infinity();
    for (int trial = 0; trial < std::max(1, n_trials); ++trial) {
        const IMP::algebra::Rotation3D rotation =
                IMP::algebra::get_random_rotation_3d();
        // A direction drawn uniformly on the sphere: z uniform in [-1, 1] and
        // the azimuth uniform, which is the one construction that does not
        // crowd the poles.
        const double z = 2.0 * unit(rng) - 1.0;
        const double azimuth = 2.0 * IMP::PI * unit(rng);
        const double r = std::sqrt(std::max(0.0, 1.0 - z * z));
        const IMP::algebra::Vector3D direction(r * std::cos(azimuth),
                                               r * std::sin(azimuth), z);
        const IMP::algebra::Vector3D offset = host_centre + direction * distance;

        IMP::algebra::Vector3Ds pose;
        for (std::size_t i = 0; i < shape.size(); ++i) {
            pose.push_back(offset + rotation.get_rotated(shape[i]));
        }
        for (std::size_t i = 0; i < guest.size(); ++i) {
            IMP::core::XYZ(model, guest[i]).set_coordinates(pose[i]);
        }
        const double score = scoring_function->evaluate(false);
        if (score < best_score) {
            best_score = score;
            best = pose;
        }
    }
    for (std::size_t i = 0; i < guest.size(); ++i) {
        IMP::core::XYZ(model, guest[i]).set_coordinates(best[i]);
    }
    (void)best_offset;
    return best_score;
}

IMP::Restraints create_go_restraints(
        IMP::Model* model, const ProbeForceFieldSystem& system,
        const std::vector<std::string>& site_ids,
        const IMP::ParticleIndexes& particles,
        const std::map<std::string, std::string>& site_atom_names,
        const std::string& component,
        const std::vector<std::string>& only_sites, double k, double cutoff) {
    IMP::Restraints out;
    if (site_ids.size() != particles.size()) {
        IMP_THROW("create_go_restraints: " << site_ids.size() << " site ids "
                  << "against " << particles.size() << " particles",
                  IMP::ValueException);
    }
    std::map<std::string, IMP::ParticleIndex> site;
    for (std::size_t i = 0; i < site_ids.size(); ++i) {
        site[site_ids[i]] = particles[i];
    }
    const std::set<std::string> released(only_sites.begin(), only_sites.end());

    // The component's heavy sites, sorted, so the restraints come out in the
    // site ids' order and a repeated run builds the same set.
    std::vector<std::string> heavy;
    for (std::size_t i = 0; i < system.get_sites().size(); ++i) {
        const FFSite& s = system.get_sites()[i];
        if (s.component != component) continue;
        if (site.find(s.id) == site.end()) continue;
        std::map<std::string, std::string>::const_iterator name =
                site_atom_names.find(s.id);
        if (name != site_atom_names.end() && !name->second.empty() &&
            (name->second[0] == 'H' || name->second[0] == 'h')) {
            continue;
        }
        heavy.push_back(s.id);
    }
    std::sort(heavy.begin(), heavy.end());

    for (std::size_t i = 0; i < heavy.size(); ++i) {
        for (std::size_t j = i + 1; j < heavy.size(); ++j) {
            if (!released.empty() && released.count(heavy[i]) == 0 &&
                released.count(heavy[j]) == 0) {
                continue;
            }
            // Within two bonds the bonded terms already say what the distance
            // is; a contact there would be a second opinion.
            if (system.is_within_bonds(heavy[i], heavy[j], 2)) continue;
            const IMP::core::XYZ a(model, site[heavy[i]]),
                    b(model, site[heavy[j]]);
            const double d = IMP::core::get_distance(a, b);
            if (d > cutoff) continue;
            IMP::Restraint* r = new IMP::core::DistanceRestraint(
                    model, new IMP::core::Harmonic(std::max(d, 1.0), k),
                    model->get_particle(site[heavy[i]]),
                    model->get_particle(site[heavy[j]]));
            r->set_name("go_" + component + "_" + heavy[i] + "_" + heavy[j]);
            out.push_back(r);
        }
    }
    return out;
}

IMP::Restraint* create_steric_restraint(IMP::Model* model,
                                       const ProbeForceFieldSystem& system,
                                       const std::vector<std::string>& site_ids,
                                       const IMP::ParticleIndexes& particles,
                                       double k) {
    if (site_ids.size() != particles.size()) {
        IMP_THROW("create_steric_restraint: " << site_ids.size()
                  << " site ids against " << particles.size() << " particles",
                  IMP::ValueException);
    }
    std::map<std::string, IMP::ParticleIndex> site;
    for (std::size_t i = 0; i < site_ids.size(); ++i) {
        site[site_ids[i]] = particles[i];
    }
    // A system that says its non-bonded term is off has no steric restraint,
    // rather than one nobody asked for.
    if (!system.get_nonbonded().enabled) return NULL;
    const std::set<std::pair<std::string, std::string> > excluded =
            system.get_exclusions();

    // Sorted, so the pair order is the site ids' and not a hash's: a run that
    // is repeated has to build the same container.
    std::vector<std::string> ids;
    for (std::map<std::string, IMP::ParticleIndex>::const_iterator it =
                 site.begin();
         it != site.end(); ++it) {
        ids.push_back(it->first);
    }

    // A soft sphere is a sphere: a site the caller decorated without a radius
    // would contribute nothing at all, silently. The system says what radius
    // each site has, so give it that rather than score an empty term.
    for (std::map<std::string, IMP::ParticleIndex>::const_iterator it =
                 site.begin();
         it != site.end(); ++it) {
        IMP::Particle* p = model->get_particle(it->second);
        if (IMP::core::XYZR::get_is_setup(p)) continue;
        double radius = 1.7;
        for (std::size_t i = 0; i < system.get_sites().size(); ++i) {
            if (system.get_sites()[i].id != it->first) continue;
            radius = system.get_sites()[i].radius;
            break;
        }
        IMP::core::XYZR::setup_particle(p, radius);
    }

    IMP::ParticleIndexPairs pairs;
    for (std::size_t i = 0; i < ids.size(); ++i) {
        for (std::size_t j = i + 1; j < ids.size(); ++j) {
            if (excluded.count(std::make_pair(ids[i], ids[j])) > 0) continue;
            pairs.push_back(IMP::ParticleIndexPair(site[ids[i]], site[ids[j]]));
        }
    }
    if (pairs.empty()) return NULL;

    IMP_NEW(IMP::container::ListPairContainer, container, (model, pairs));
    const double strength = k >= 0.0 ? k : system.get_nonbonded().k;
    return new IMP::container::PairsRestraint(
            new IMP::core::SoftSpherePairScore(strength), container,
            "steric");
}

IMPBFF_END_NAMESPACE
