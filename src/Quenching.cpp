/**
 * \file Quenching.cpp
 * \brief Collisional quenching: the PET tables, sphere stamping, the fields, the race.
 *
 * Sections in the order of IMP/bff/Quenching.h; each is marked with the
 * file it came from.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

// -------- from PETQuenching.cpp --------
/**
 * (formerly PETQuenching.cpp, now a section of this file)
 * \brief Photoinduced electron transfer: which moieties quench, and how hard.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#include <IMP/bff/Quenching.h>

#include <IMP/bff/internal/OutputView.h>
#include <IMP/bff/internal/Text.h>
#include <IMP/bff/Base.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <limits>

IMPBFF_BEGIN_NAMESPACE

const char* const REFERENCE_DYE = "AlexaFluor488";
const double DEFAULT_PROBE_RADIUS = 3.5;

namespace pet {
using IMP::bff::internal::nan_value;



//! Upper-cased and trimmed, which is how every residue name is compared here.
std::string residue_key(const std::string& name) {
    const std::string space = " \t\n\r";
    const std::size_t a = name.find_first_not_of(space);
    if (a == std::string::npos) return std::string();
    const std::size_t b = name.find_last_not_of(space);
    std::string out = name.substr(a, b - a + 1);
    for (std::size_t i = 0; i < out.size(); ++i) {
        out[i] = static_cast<char>(
                std::toupper(static_cast<unsigned char>(out[i])));
    }
    return out;
}

double clamp_slow_factor(double v) { return std::min(1.0, std::max(0.0, v)); }

//! A positive, finite radius, or NaN to inherit the global one.
double sane_radius(double v) {
    if (!(v > 0.0) || !std::isfinite(v)) return nan_value();
    return v;
}

//! A de-duplicated, upper-cased atom list, falling back to \p fallback.
std::vector<std::string> sane_atoms(const std::vector<std::string>& atoms,
                                    const std::vector<std::string>& fallback) {
    std::vector<std::string> out;
    for (std::size_t i = 0; i < atoms.size(); ++i) {
        const std::string name = residue_key(atoms[i]);
        if (name.empty()) continue;
        if (std::find(out.begin(), out.end(), name) == out.end()) {
            out.push_back(name);
        }
    }
    return out.empty() ? fallback : out;
}

std::vector<std::string> atoms_of(const std::string& residue) {
    const std::map<std::string, std::vector<std::string> > table =
            quencher_atoms();
    std::map<std::string, std::vector<std::string> >::const_iterator it =
            table.find(residue);
    if (it != table.end()) return it->second;
    return std::vector<std::string>(1, "CB");
}

//! `residue_names` order, one value per name, from a lookup that may miss.
const ResidueQuenching& lookup(
        const std::map<std::string, ResidueQuenching>& table,
        const std::string& residue, const ResidueQuenching& fallback) {
    std::map<std::string, ResidueQuenching>::const_iterator it =
            table.find(residue_key(residue));
    return it == table.end() ? fallback : it->second;
}

}  // namespace pet

// --------------------------------------------------------------------------
// the values
// --------------------------------------------------------------------------

Quencher::Quencher(std::string comp_id,
                   const std::vector<std::string>& atom_ids,
                   std::string asym_id, int seq_id)
    : comp_id(comp_id), atom_ids(atom_ids), asym_id(asym_id), seq_id(seq_id) {}

Quencher Quencher::at(std::string asym_id, int seq_id) const {
    return Quencher(comp_id, atom_ids, asym_id, seq_id);
}

void Quencher::show(std::ostream& out) const {
    out << "Quencher(" << comp_id;
    if (!get_is_typed()) out << " @" << asym_id << seq_id;
    out << ", atoms=[";
    for (std::size_t i = 0; i < atom_ids.size(); ++i) {
        out << (i ? ", " : "") << atom_ids[i];
    }
    out << "])";
}

PETParameters::PETParameters(std::string dye, std::string comp_id,
                             double rate_constant, double contact_distance,
                             double attenuation_length,
                             std::string measured_for)
    : dye(dye), comp_id(comp_id), rate_constant(rate_constant),
      contact_distance(contact_distance),
      attenuation_length(attenuation_length),
      // An unstated `measured_for` means "this dye" -- so `is_transferred` is
      // false by construction unless a caller says otherwise.
      measured_for(measured_for.empty() ? dye : measured_for) {
    if (rate_constant < 0.0) {
        IMP_THROW("rate_constant must be >= 0, not " << rate_constant,
                  ValueException);
    }
}

PETParameters PETParameters::scaled(double rate_scale) const {
    return PETParameters(dye, comp_id, rate_constant * rate_scale,
                         contact_distance, attenuation_length, measured_for);
}

void PETParameters::show(std::ostream& out) const {
    out << "PETParameters(" << dye << " x " << comp_id
        << ", kQ=" << rate_constant;
    if (get_is_transferred()) out << ", transferred from " << measured_for;
    out << ")";
}

ResidueQuenching::ResidueQuenching(double slow_factor, double kQ,
                                   double quench_radius,
                                   const std::vector<std::string>& quench_atoms)
    : slow_factor(slow_factor), kQ(kQ), quench_radius(quench_radius),
      quench_atoms(quench_atoms) {}

void ResidueQuenching::show(std::ostream& out) const {
    out << "ResidueQuenching(kQ=" << kQ << ", slow=" << slow_factor << ")";
}

// --------------------------------------------------------------------------
// the tables
// --------------------------------------------------------------------------

std::vector<std::string> standard_amino_acid_residues() {
    static const char* const names[] = {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"};
    return std::vector<std::string>(names, names + 20);
}

std::map<std::string, std::vector<std::string> > quencher_atoms() {
    std::map<std::string, std::vector<std::string> > out;
#define ATOMS(residue, ...)                                              \
    {                                                                    \
        const char* const a[] = {__VA_ARGS__};                           \
        out[residue] = std::vector<std::string>(a, a + sizeof(a) / sizeof(a[0])); \
    }
    ATOMS("ALA", "CB")
    ATOMS("ARG", "CZ", "NE", "NH1", "NH2")
    ATOMS("ASN", "CG", "OD1", "ND2")
    ATOMS("ASP", "CG", "OD1", "OD2")
    ATOMS("CYS", "SG")
    ATOMS("GLN", "CD", "OE1", "NE2")
    ATOMS("GLU", "CD", "OE1", "OE2")
    ATOMS("GLY", "CA")
    ATOMS("HIS", "CG", "ND1", "CD2", "CE1", "NE2")
    ATOMS("ILE", "CB")
    ATOMS("LEU", "CB")
    ATOMS("LYS", "NZ")
    ATOMS("MET", "SD")
    ATOMS("PHE", "CG", "CD1", "CD2", "CE1", "CE2", "CZ")
    ATOMS("PRO", "N", "CB", "CG", "CD")
    ATOMS("SER", "OG")
    ATOMS("THR", "OG1")
    ATOMS("TRP", "CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2", "NE1", "CG", "CD1")
    ATOMS("TYR", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH")
    ATOMS("VAL", "CB")
#undef ATOMS
    return out;
}

std::map<std::string, PETReference> pet_quenching_reference() {
    std::map<std::string, PETReference> out;
    // TRP > PRO ~ TYR > MET > HIS > CYS. Peulen et al., JPC B 2017, 121, 8211.
    out["TRP"] = PETReference(3.5, 5.0);
    out["TYR"] = PETReference(2.0, 5.0);
    out["MET"] = PETReference(1.67, 3.5);
    out["HIS"] = PETReference(1.0, 4.7);
    out["CYS"] = PETReference(0.8, 3.5);
    out["PRO"] = PETReference(2.0, 4.0);
    return out;
}

std::map<std::string, Quencher> reference_quenchers() {
    std::map<std::string, Quencher> out;
    const std::map<std::string, PETReference> reference =
            pet_quenching_reference();
    for (std::map<std::string, PETReference>::const_iterator it =
                 reference.begin(); it != reference.end(); ++it) {
        out[it->first] = Quencher(it->first, pet::atoms_of(it->first));
    }
    return out;
}

std::map<std::string, PETParameters> reference_pet_parameters(
        std::string dye, double rate_scale, double attenuation_length) {
    std::map<std::string, PETParameters> out;
    const std::map<std::string, PETReference> reference =
            pet_quenching_reference();
    for (std::map<std::string, PETReference>::const_iterator it =
                 reference.begin(); it != reference.end(); ++it) {
        // `measured_for` is the reference dye whatever `dye` is: the table was
        // measured once, and asking for another dye is a transfer, which
        // `is_transferred` then reports rather than hiding.
        out[it->first] = PETParameters(dye, it->first,
                                       it->second.kQ * rate_scale,
                                       it->second.contact_distance,
                                       attenuation_length, REFERENCE_DYE);
    }
    return out;
}

std::map<std::string, ResidueQuenching> normalize_amino_acid_quenching(
        const std::map<std::string, ResidueQuenching>& table) {
    std::map<std::string, ResidueQuenching> out;

    const std::vector<std::string> standard = standard_amino_acid_residues();
    for (std::size_t i = 0; i < standard.size(); ++i) {
        out[standard[i]] = ResidueQuenching(1.0, 0.0, pet::nan_value(),
                                            pet::atoms_of(standard[i]));
    }

    for (std::map<std::string, ResidueQuenching>::const_iterator it =
                 table.begin(); it != table.end(); ++it) {
        const std::string name = pet::residue_key(it->first);
        if (name.empty()) continue;
        const std::vector<std::string> fallback = pet::atoms_of(name);
        ResidueQuenching entry;
        entry.slow_factor = pet::clamp_slow_factor(it->second.slow_factor);
        entry.kQ = std::max(0.0, it->second.kQ);
        entry.quench_radius = pet::sane_radius(it->second.quench_radius);
        entry.quench_atoms = pet::sane_atoms(it->second.quench_atoms, fallback);
        out[name] = entry;
    }
    return out;
}

std::map<std::string, ResidueQuenching> amino_acid_quenching_defaults(
        double kQ_scale, double slow_factor, double probe_radius) {
    const double radius = std::max(0.0, probe_radius);
    std::map<std::string, ResidueQuenching> table =
            normalize_amino_acid_quenching();
    const std::map<std::string, PETReference> reference =
            pet_quenching_reference();

    for (std::map<std::string, ResidueQuenching>::iterator it = table.begin();
         it != table.end(); ++it) {
        it->second.slow_factor = pet::clamp_slow_factor(slow_factor);
        std::map<std::string, PETReference>::const_iterator ref =
                reference.find(it->first);
        if (ref == reference.end()) continue;
        it->second.kQ = std::max(0.0, ref->second.kQ * kQ_scale);
        // The table is quoted dye-surface-to-quencher; the walk tracks the dye
        // *centre*, so the radius is added on the way in.
        it->second.quench_radius =
                pet::sane_radius(radius + ref->second.contact_distance);
    }
    return table;
}

// --------------------------------------------------------------------------
// per-residue lookups
// --------------------------------------------------------------------------

std::vector<double> slow_factors_for_residues(
        const std::vector<std::string>& residue_names,
        const std::map<std::string, ResidueQuenching>& table) {
    static const ResidueQuenching fallback;
    std::vector<double> out;
    out.reserve(residue_names.size());
    for (std::size_t i = 0; i < residue_names.size(); ++i) {
        out.push_back(pet::lookup(table, residue_names[i], fallback).slow_factor);
    }
    return out;
}

std::vector<double> quenching_rates_for_residues(
        const std::vector<std::string>& residue_names,
        const std::map<std::string, ResidueQuenching>& table) {
    static const ResidueQuenching fallback;
    std::vector<double> out;
    out.reserve(residue_names.size());
    for (std::size_t i = 0; i < residue_names.size(); ++i) {
        out.push_back(pet::lookup(table, residue_names[i], fallback).kQ);
    }
    return out;
}

std::vector<double> quench_radii_for_residues(
        const std::vector<std::string>& residue_names,
        const std::map<std::string, ResidueQuenching>& table,
        double critical_distance) {
    static const ResidueQuenching fallback;
    const double global = std::isfinite(critical_distance) ? critical_distance : 0.0;
    std::vector<double> out;
    out.reserve(residue_names.size());
    for (std::size_t i = 0; i < residue_names.size(); ++i) {
        const double r = pet::lookup(table, residue_names[i], fallback).quench_radius;
        out.push_back(r == r ? r : global);
    }
    return out;
}

// --------------------------------------------------------------------------
// residues, and the atoms that quench
// --------------------------------------------------------------------------

void ResidueSites::add(const double* slow, const double* quench,
                       const std::string& residue_name) {
    for (int i = 0; i < 3; ++i) slow_centers_.push_back(slow[i]);
    for (int i = 0; i < 3; ++i) quench_centers_.push_back(quench[i]);
    residue_names_.push_back(residue_name);
}

void ResidueSites::get_slow_centers(double** out_view, int* n_out_view) const {
    internal::copy_to_view(slow_centers_, out_view, n_out_view);
}

void ResidueSites::get_quench_centers(double** out_view, int* n_out_view) const {
    internal::copy_to_view(quench_centers_, out_view, n_out_view);
}

ResidueSites residue_sites(const std::vector<std::string>& chains,
                           const std::vector<int>& res_ids,
                           const std::vector<std::string>& res_names,
                           const std::vector<std::string>& atom_names,
                           double* coords, int n_atoms, int n_dim,
                           const std::map<std::string, ResidueQuenching>& table) {
    ResidueSites out;
    if (n_atoms == 0) return out;
    if (n_dim != 3) {
        IMP_THROW("atoms must be (N, 3), not (" << n_atoms << ", " << n_dim << ")",
                  ValueException);
    }
    const std::size_t n = static_cast<std::size_t>(n_atoms);
    if (chains.size() != n || res_ids.size() != n || res_names.size() != n ||
        atom_names.size() != n) {
        IMP_THROW("one chain, residue id, residue name and atom name per atom",
                  ValueException);
    }

    const std::map<std::string, ResidueQuenching> full =
            normalize_amino_acid_quenching(table);

    // Insertion order, not sorted: the centres come back in the order the
    // residues appear in the structure, which is what every consumer indexes
    // against. A std::map would silently reorder them.
    typedef std::pair<std::pair<std::string, int>, std::string> Key;
    std::vector<Key> order;
    std::map<Key, std::vector<std::size_t> > by_residue;
    for (std::size_t i = 0; i < n; ++i) {
        const Key key(std::make_pair(pet::residue_key(chains[i]), res_ids[i]),
                      pet::residue_key(res_names[i]));
        if (by_residue.find(key) == by_residue.end()) order.push_back(key);
        by_residue[key].push_back(i);
    }

    for (std::size_t r = 0; r < order.size(); ++r) {
        const std::vector<std::size_t>& indices = by_residue[order[r]];
        const std::string residue_name = order[r].second;

        std::size_t selected = indices[0];
        for (std::size_t k = 0; k < indices.size(); ++k) {
            if (pet::residue_key(atom_names[indices[k]]) == "CB") {
                selected = indices[k];
                break;
            }
        }
        if (pet::residue_key(atom_names[selected]) != "CB") {
            for (std::size_t k = 0; k < indices.size(); ++k) {
                if (pet::residue_key(atom_names[indices[k]]) == "CA") {
                    selected = indices[k];
                    break;
                }
            }
        }

        const double* slow = coords + selected * 3;

        std::vector<std::string> wanted;
        std::map<std::string, ResidueQuenching>::const_iterator w =
                full.find(residue_name);
        wanted = (w == full.end() || w->second.quench_atoms.empty())
                         ? pet::atoms_of(residue_name)
                         : w->second.quench_atoms;

        double cx = 0.0, cy = 0.0, cz = 0.0;
        std::size_t matched = 0;
        for (std::size_t k = 0; k < indices.size(); ++k) {
            const std::string name = pet::residue_key(atom_names[indices[k]]);
            if (std::find(wanted.begin(), wanted.end(), name) == wanted.end()) {
                continue;
            }
            const double* c = coords + indices[k] * 3;
            cx += c[0]; cy += c[1]; cz += c[2];
            ++matched;
        }

        double quench[3];
        if (matched > 0) {
            quench[0] = cx / matched;
            quench[1] = cy / matched;
            quench[2] = cz / matched;
        } else {
            // A residue with no redox-active atom in the structure still has a
            // position; the slow centre is the honest one to use.
            quench[0] = slow[0]; quench[1] = slow[1]; quench[2] = slow[2];
        }
        out.add(slow, quench, residue_name);
    }
    return out;
}

void atomic_quenching_parameters(
        const std::vector<std::string>& res_names,
        const std::vector<std::string>& atom_names,
        const std::map<std::string, PETParameters>& parameters,
        double** out_kQ, int* n_out_kQ, double** out_rC, int* n_out_rC) {
    const std::size_t n = res_names.size();
    double* kQ = internal::new_double_view(n, out_kQ, n_out_kQ);
    double* rC = internal::new_double_view(n, out_rC, n_out_rC);
    if (kQ == NULL || rC == NULL) return;
    if (atom_names.size() != n) {
        IMP_THROW("one atom name per residue name: " << atom_names.size()
                          << " against " << n,
                  ValueException);
    }

    const std::map<std::string, std::vector<std::string> > active =
            quencher_atoms();
    for (std::size_t i = 0; i < n; ++i) {
        kQ[i] = 0.0;
        rC[i] = 0.0;
        const std::string residue = pet::residue_key(res_names[i]);
        std::map<std::string, PETParameters>::const_iterator p =
                parameters.find(residue);
        if (p == parameters.end()) continue;
        std::map<std::string, std::vector<std::string> >::const_iterator a =
                active.find(residue);
        if (a == active.end()) continue;
        const std::string atom = pet::residue_key(atom_names[i]);
        if (std::find(a->second.begin(), a->second.end(), atom) ==
            a->second.end()) {
            continue;
        }
        kQ[i] = p->second.rate_constant;
        // An absent or zero attenuation length is a hard contact sphere, which
        // the exponential form spells as a 1 A decay.
        const double rc = p->second.attenuation_length;
        rC[i] = (rc == rc && rc != 0.0) ? rc : 1.0;
    }
}

IMPBFF_END_NAMESPACE

// -------- from QuenchingGrid.cpp --------
/**
 * (formerly QuenchingGrid.cpp, now a section of this file)
 * \brief Stamping spheres of influence onto an accessible-volume grid.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#include <IMP/bff/States.h>

#include <cstdlib>

IMPBFF_BEGIN_NAMESPACE

void center_grid_indices(
        const std::vector<double>& rs, const std::vector<double>& r0, double dg,
        int ng, const std::vector<double>& radius, std::vector<int>& ix0,
        std::vector<int>& iy0, std::vector<int>& iz0,
        std::vector<int>& radius_idx) {
    const std::size_t n = rs.size() / 3;
    ix0.assign(n, 0); iy0.assign(n, 0); iz0.assign(n, 0); radius_idx.assign(n, 0);
    const int offset = (ng - 1) / 2;
    for (std::size_t i = 0; i < n; ++i) {
        ix0[i] = static_cast<int>(std::floor((rs[3 * i + 0] - r0[0]) / dg)) + offset;
        iy0[i] = static_cast<int>(std::floor((rs[3 * i + 1] - r0[1]) / dg)) + offset;
        iz0[i] = static_cast<int>(std::floor((rs[3 * i + 2] - r0[2]) / dg)) + offset;
        radius_idx[i] = static_cast<int>(radius[i] / dg);
    }
}

namespace {
std::vector<double> stamp_spheres_impl(
        const std::vector<double>& density, int ng,
        const std::vector<double>& radius, const std::vector<double>& rs,
        const std::vector<double>& r0, double dg,
        const std::vector<double>& values, int combine) {
    const bool multiply = (combine == GRID_COMBINE_MULTIPLY);
    const double identity = multiply ? 1.0 : 0.0;
    const std::size_t n_vox = static_cast<std::size_t>(ng) * ng * ng;
    std::vector<double> factors(n_vox, identity);

    std::vector<int> ix0, iy0, iz0, r_idx;
    center_grid_indices(rs, r0, dg, ng, radius, ix0, iy0, iz0, r_idx);
    const std::size_t n_centre = ix0.size();
    if (n_centre == 0) {
        // Nothing stamped: every voxel keeps the identity, accessible or not.
        return factors;
    }

    std::vector<double> slab(static_cast<std::size_t>(ng) * ng, identity);
    for (int ix = 0; ix < ng; ++ix) {
        std::fill(slab.begin(), slab.end(), identity);
        bool touched = false;
        for (std::size_t c = 0; c < n_centre; ++c) {
            const double value = values[c];
            if (!multiply && value == 0.0) continue;   // an added zero is a no-op
            const int ri = r_idx[c];
            const int dx = ix - ix0[c];
            const int remaining = ri * ri - dx * dx;
            if (remaining <= 0) continue;
            const int y_lo = std::max(0, iy0[c] - ri);
            const int y_hi = std::min(ng - 1, iy0[c] + ri);
            for (int iy = y_lo; iy <= y_hi; ++iy) {
                const int dy = iy - iy0[c];
                const int span2 = remaining - dy * dy;
                if (span2 <= 0) continue;
                int span = static_cast<int>(std::sqrt(static_cast<double>(span2)));
                // Membership is strict (d^2 < r^2), so drop the boundary voxel
                // when span2 is a perfect square.
                if (span * span >= span2) --span;
                if (span < 0) continue;
                const int z_lo = std::max(0, iz0[c] - span);
                const int z_hi = std::min(ng - 1, iz0[c] + span);
                for (int iz = z_lo; iz <= z_hi; ++iz) {
                    double& s = slab[static_cast<std::size_t>(iy) * ng + iz];
                    if (multiply) s *= value; else s += value;
                    touched = true;
                }
            }
        }
        if (!touched) continue;
        for (int iy = 0; iy < ng; ++iy) {
            for (int iz = 0; iz < ng; ++iz) {
                const std::size_t k =
                        (static_cast<std::size_t>(ix) * ng + iy) * ng + iz;
                if (density[k] != 0.0) {
                    factors[k] = slab[static_cast<std::size_t>(iy) * ng + iz];
                }
            }
        }
    }
    return factors;
}
}  // namespace

void stamp_spheres(const std::vector<double>& density, int ng,
        const std::vector<double>& radius, const std::vector<double>& rs,
        const std::vector<double>& r0, double dg,
        const std::vector<double>& values, int combine, double** out_view, int* n_out_view) {
    internal::copy_to_view(stamp_spheres_impl(density, ng, radius, rs, r0, dg, values, combine),
                           out_view, n_out_view);
}

void slow_factor_grid(const std::vector<double>& density, int ng, double dg,
                      const std::vector<double>& slow_radius,
                      const std::vector<double>& rs,
                      const std::vector<double>& r0,
                      const std::vector<double>& slow_fact,
                      double** out_view, int* n_out_view) {
    stamp_spheres(density, ng, slow_radius, rs, r0, dg, slow_fact,
                  GRID_COMBINE_MULTIPLY, out_view, n_out_view);
}

void quenching_rate_grid(const std::vector<double>& density, int ng, double dg,
                         const std::vector<double>& radius,
                         const std::vector<double>& rs,
                         const std::vector<double>& r0,
                         const std::vector<double>& values,
                         double** out_view, int* n_out_view) {
    stamp_spheres(density, ng, radius, rs, r0, dg, values, GRID_COMBINE_ADD,
                  out_view, n_out_view);
}

void av_contact_mask(const std::vector<double>& density, int ng, double dg,
                     const std::vector<double>& slow_radius,
                     const std::vector<double>& rs,
                     const std::vector<double>& r0,
                     int** out_view_i, int* n_out_view_i) {
    int* labels = nullptr;
    int n_labels = 0, d2 = 0, d3 = 0;
    split_contact_volume(density, ng, dg, slow_radius, rs, r0, &labels,
                         &n_labels, &d2, &d3);
    const std::size_t total = static_cast<std::size_t>(n_labels) * d2 * d3;
    int* out = internal::new_int_view(total, out_view_i, n_out_view_i);
    if (out == nullptr) { std::free(labels); return; }
    for (std::size_t i = 0; i < total; ++i) {
        out[i] = labels[i] == AV_VOXEL_CONTACT ? 1 : 0;
    }
    std::free(labels);
}

IMPBFF_END_NAMESPACE

// -------- from QuenchingMap.cpp --------
/**
 * (formerly QuenchingMap.cpp, now a section of this file)
 * \brief Mobility, quenching-rate and FRET-rate fields on an AV grid.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */



IMPBFF_BEGIN_NAMESPACE

// Named, not anonymous: IMP compiles this module as one translation unit.
namespace qmap {

//! The side of a cube of `n` voxels, rounded to the nearest integer.
/*! `cbrt` of a perfect cube can land a hair below it in double arithmetic, and
    truncating then gives `ng - 1` -- one axis short of the grid. */
int grid_side(std::size_t n) {
    return static_cast<int>(std::floor(std::cbrt(static_cast<double>(n)) + 0.5));
}

}  // namespace qmap

void slow_near_atoms(
        const std::vector<double>& d_map, const std::vector<double>& density,
        const std::vector<double>& axis, const std::vector<double>& r0,
        const std::vector<double>& atoms_xyz, double min_distance_sq,
        double factor, double** out_view, int* n_out_view) {
    const std::size_t ng = axis.size();
    const std::size_t n_atoms = atoms_xyz.size() / 3;
    double* out = internal::new_double_view(density.size(), out_view, n_out_view);
    if (out == nullptr) return;
    for (std::size_t ix = 0; ix < ng; ++ix) {
        const double x = axis[ix] + r0[0];
        for (std::size_t iy = 0; iy < ng; ++iy) {
            const double y = axis[iy] + r0[1];
            for (std::size_t iz = 0; iz < ng; ++iz) {
                const std::size_t k = (ix * ng + iy) * ng + iz;
                if (density[k] <= 0.0) continue;
                const double z = axis[iz] + r0[2];
                double slow = 1.0;
                for (std::size_t a = 0; a < n_atoms; ++a) {
                    const double dx = atoms_xyz[3 * a + 0] - x;
                    const double dy = atoms_xyz[3 * a + 1] - y;
                    const double dz = atoms_xyz[3 * a + 2] - z;
                    // Applied once per contacting atom, so it compounds.
                    if (dx * dx + dy * dy + dz * dz < min_distance_sq) slow *= factor;
                }
                out[k] = d_map[k] * slow;
            }
        }
    }
}

void quenching_map(
        const std::vector<double>& density, const std::vector<double>& axis,
        const std::vector<double>& r0, const std::vector<double>& atoms_xyz,
        const std::vector<double>& kQ, const std::vector<double>& rC,
        double probe_radius, double inv_tau0,
        double** out_view, int* n_out_view) {
    const std::size_t ng = axis.size();
    const std::size_t n_atoms = atoms_xyz.size() / 3;
    double* out = internal::new_double_view(density.size(), out_view, n_out_view);
    if (out == nullptr) return;
    for (std::size_t ix = 0; ix < ng; ++ix) {
        const double x = axis[ix] + r0[0];
        for (std::size_t iy = 0; iy < ng; ++iy) {
            const double y = axis[iy] + r0[1];
            for (std::size_t iz = 0; iz < ng; ++iz) {
                const std::size_t k = (ix * ng + iy) * ng + iz;
                // Outside the volume the rate stays zero, not 1/tau0.
                if (density[k] <= 0.0) continue;
                const double z = axis[iz] + r0[2];
                double v = inv_tau0;
                for (std::size_t a = 0; a < n_atoms; ++a) {
                    if (kQ[a] == 0.0 || rC[a] == 0.0) continue;
                    const double dx = atoms_xyz[3 * a + 0] - x;
                    const double dy = atoms_xyz[3 * a + 1] - y;
                    const double dz = atoms_xyz[3 * a + 2] - z;
                    const double d = std::sqrt(dx * dx + dy * dy + dz * dz) - probe_radius;
                    v += kQ[a] * std::exp(-d / rC[a]);
                }
                out[k] = v;
            }
        }
    }
}

void fret_map(
        const std::vector<double>& density_d, const std::vector<double>& density_a,
        const std::vector<double>& axis_d, const std::vector<double>& axis_a,
        const std::vector<double>& r0_d, const std::vector<double>& r0_a,
        double r0_6, double kf, int step,
        double** out_view, int* n_out_view) {
    const std::size_t ng_d = axis_d.size();
    const std::size_t ng_a = axis_a.size();
    double* out = internal::new_double_view(density_d.size(), out_view, n_out_view);
    if (out == nullptr) return;
    if (step < 1) step = 1;
    for (std::size_t ixd = 0; ixd < ng_d; ++ixd) {
        const double x = axis_d[ixd] + r0_d[0];
        for (std::size_t iyd = 0; iyd < ng_d; ++iyd) {
            const double y = axis_d[iyd] + r0_d[1];
            for (std::size_t izd = 0; izd < ng_d; ++izd) {
                const std::size_t kd = (ixd * ng_d + iyd) * ng_d + izd;
                if (density_d[kd] <= 0.0) continue;
                const double z = axis_d[izd] + r0_d[2];
                double t_ret = 0.0;   // accumulate transfer *time*, not rate
                double weight = 0.0;
                for (std::size_t ixa = 0; ixa < ng_a; ixa += step) {
                    const double ddx = axis_a[ixa] + r0_a[0] - x;
                    const double sx = ddx * ddx;
                    for (std::size_t iya = 0; iya < ng_a; iya += step) {
                        const double ddy = axis_a[iya] + r0_a[1] - y;
                        const double sy = ddy * ddy;
                        for (std::size_t iza = 0; iza < ng_a; iza += step) {
                            const std::size_t ka = (ixa * ng_a + iya) * ng_a + iza;
                            const double da = density_a[ka];
                            if (da <= 0.0) continue;
                            const double ddz = axis_a[iza] + r0_a[2] - z;
                            const double rda2 = sx + sy + ddz * ddz;
                            const double r2 = r0_6 / (rda2 * rda2 * rda2);
                            t_ret += da / (r2 * kf);
                            weight += da;
                        }
                    }
                }
                // Harmonic mean: the weighted mean transfer time, inverted.
                if (weight > 0.0 && t_ret > 0.0) out[kd] = weight / t_ret;
            }
        }
    }
}

std::vector<double> grid_axis(int ng, double dg) {
    std::vector<double> axis(ng > 0 ? ng : 0);
    const int centre = (ng - 1) / 2;
    for (int i = 0; i < ng; ++i) axis[i] = (i - centre) * dg;
    return axis;
}

void diffusion_coefficient_map(const std::vector<double>& density,
                               const std::vector<double>& r0, double dg,
                               const std::vector<double>& atoms_xyz,
                               double free_diffusion, double min_distance,
                               double slow_factor,
                               const std::vector<double>& base,
                               double** out_view, int* n_out_view) {
    const int ng = qmap::grid_side(density.size());
    std::vector<double> d_map;
    if (base.empty()) {
        d_map.assign(density.size(), free_diffusion);
    } else {
        if (base.size() != density.size()) {
            IMP_THROW("the base coefficient map must have as many voxels as the "
                      "density: " << base.size() << " against " << density.size(),
                      ValueException);
        }
        d_map = base;
    }
    slow_near_atoms(d_map, density, grid_axis(ng, dg), r0, atoms_xyz,
                    min_distance * min_distance, slow_factor, out_view,
                    n_out_view);
}

void quenching_rate_map(const std::vector<double>& density,
                        const std::vector<double>& r0, double dg,
                        const std::vector<double>& atoms_xyz,
                        const std::vector<double>& kQ,
                        const std::vector<double>& rC, double tau0,
                        double probe_radius, double** out_view,
                        int* n_out_view) {
    const int ng = qmap::grid_side(density.size());
    quenching_map(density, grid_axis(ng, dg), r0, atoms_xyz, kQ, rC, probe_radius,
                  tau0 > 0.0 ? 1.0 / tau0 : 0.0, out_view, n_out_view);
}

void fret_rate_map(const std::vector<double>& density_donor,
                   const std::vector<double>& density_acceptor,
                   const std::vector<double>& r0_donor,
                   const std::vector<double>& r0_acceptor, double dg_donor,
                   double dg_acceptor, double forster_radius, double kf,
                   int acceptor_step, double** out_view, int* n_out_view) {
    if (!(kf > 0.0)) {
        IMP_THROW("kf (the donor's radiative rate) must be positive, not " << kf,
                  ValueException);
    }
    const double r0_6 = std::pow(forster_radius, 6);
    fret_map(density_donor, density_acceptor,
             grid_axis(qmap::grid_side(density_donor.size()), dg_donor),
             grid_axis(qmap::grid_side(density_acceptor.size()), dg_acceptor),
             r0_donor, r0_acceptor, r0_6, kf,
             acceptor_step > 1 ? acceptor_step : 1, out_view, n_out_view);
}

void radial_diffusion_map(const std::vector<double>& density, double dg,
                           const std::vector<double>& radial,
                           double** out_view, int* n_out_view) {
    const int ng = qmap::grid_side(density.size());
    if (ng == 0) {
        IMP_THROW("the density must be a grid with an integer cube-root side",
                  ValueException);
    }
    const std::vector<double> axis = grid_axis(ng, dg);
    const int max_r = radial.empty() ? 0
                                     : static_cast<int>(radial.size()) - 1;
    double* out = internal::new_double_view(density.size(), out_view, n_out_view);
    if (out == NULL) return;
    for (int ix = 0; ix < ng; ++ix) {
        const double ax = axis[ix];
        for (int iy = 0; iy < ng; ++iy) {
            const double ay = axis[iy];
            const double s2_xy = ax * ax + ay * ay;
            for (int iz = 0; iz < ng; ++iz) {
                const std::size_t k = (static_cast<std::size_t>(ix) * ng + iy) * ng + iz;
                if (radial.empty()) {
                    out[k] = 1.0;
                    continue;
                }
                const double r2 = s2_xy + axis[iz] * axis[iz];
                int r = static_cast<int>(std::lround(std::sqrt(r2)));
                if (r > max_r) r = max_r;
                out[k] = radial[static_cast<std::size_t>(r)];
            }
        }
    }
}

IMPBFF_END_NAMESPACE

// -------- from QuenchedDecay.cpp --------
/**
 * (formerly QuenchedDecay.cpp, now a section of this file)
 * \brief The fused walk -> rate -> photons path.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#include <IMP/bff/internal/PhotonRace.h>
#include <IMP/bff/internal/RandomWalk.h>


IMPBFF_BEGIN_NAMESPACE

namespace {
std::vector<double> quenched_donor_photons_impl(
        int* occupancy, int n_occupancy,
        double* mobility, int n_mobility,
        double* rate_map, int n_rate_map,
        int ng, double dg, double t_max, double t_step,
        double diffusion_coefficient,
        const std::vector<int>& walk_seeds,
        double tau0, int n_photons, int photon_seed,
        std::vector<double>& stats) {
    stats.assign(5, 0.0);
    const int n_steps = static_cast<int>(t_max / t_step);
    const std::size_t n = static_cast<std::size_t>(ng);
    const bool has_rates = static_cast<std::size_t>(n_rate_map) == n * n * n;

    // The whole point: this vector is the only large thing allocated, it is one
    // value per step rather than four, and it never crosses into Python.
    //
    // `float`, not `double`, for two reasons. It halves the memory the photon
    // race walks -- and that race is the bottleneck once the trajectory is long,
    // making tens of millions of random reads into this array. It also matches
    // `sample_grid`, which casts the rate map to float32, so the fused path and
    // the three-call path see bit-identical rates rather than rates that agree
    // to 2.6e-8 and could in principle disagree about a photon.
    std::vector<float> k_quench;
    k_quench.reserve(static_cast<std::size_t>(std::max(0, n_steps)) *
                     std::max<std::size_t>(1, walk_seeds.size()));

    long long n_accepted = 0, n_rejected = 0;
    bool any = false;
    for (std::size_t s = 0; s < walk_seeds.size(); ++s) {
        int acc = 0, rej = 0;
        const std::size_t before = k_quench.size();
        const bool ok = internal::run_walk(
                occupancy, static_cast<std::size_t>(n_occupancy),
                mobility, static_cast<std::size_t>(n_mobility),
                ng, t_step, diffusion_coefficient, dg,
                walk_seeds[s], n_steps, acc, rej,
                [&](int /*i*/, double /*px*/, double /*py*/, double /*pz*/,
                    bool /*accepted*/, std::size_t voxel) {
                    k_quench.push_back(has_rates ? static_cast<float>(rate_map[voxel]) : 0.0f);
                });
        if (!ok) {
            k_quench.resize(before);   // a walk that never started contributes nothing
            continue;
        }
        any = true;
        n_accepted += acc;
        n_rejected += rej;
    }
    if (!any) return std::vector<double>();

    // Accumulated in float, like numpy's float32 mean, so the reported average
    // matches the three-call path exactly rather than nearly.
    float sum = 0.0f;
    long long n_contact = 0;
    for (std::size_t i = 0; i < k_quench.size(); ++i) {
        sum += k_quench[i];
        if (k_quench[i] > 0.0f) ++n_contact;
    }
    const double frames = static_cast<double>(k_quench.size());
    stats[0] = frames;
    stats[1] = static_cast<double>(n_accepted);
    stats[2] = static_cast<double>(n_rejected);
    stats[3] = frames > 0.0 ? static_cast<double>(sum / static_cast<float>(frames)) : 0.0;
    stats[4] = frames > 0.0 ? static_cast<double>(n_contact) / frames : 0.0;

    // Interleaved so one returned vector carries delay and flag, for the same
    // reason the walk packs its accept flag: a SWIG out-parameter costs ~480 ns
    // per element to convert, while a returned vector becomes a tuple and costs
    // nothing.
    const int n_ph = n_photons > 0 ? n_photons : 0;
    std::vector<double> out(static_cast<std::size_t>(n_ph) * 2, 0.0);
#pragma omp parallel for schedule(static)
    for (int i = 0; i < n_ph; ++i) {
        bool got = false;
        out[2 * i + 0] = internal::race_one_photon(
                k_quench, t_step, tau0, photon_seed, i, got);
        out[2 * i + 1] = got ? 1.0 : 0.0;
    }
    return out;
}
}  // namespace

void quenched_donor_photons(int* occupancy, int n_occupancy,
        double* mobility, int n_mobility,
        double* rate_map, int n_rate_map,
        int ng, double dg, double t_max, double t_step,
        double diffusion_coefficient,
        const std::vector<int>& walk_seeds,
        double tau0, int n_photons, int photon_seed,
        std::vector<double>& stats, double** out_view, int* n_out_view) {
    internal::copy_to_view(quenched_donor_photons_impl(occupancy, n_occupancy, mobility, n_mobility, rate_map, n_rate_map, ng, dg, t_max, t_step, diffusion_coefficient, walk_seeds, tau0, n_photons, photon_seed, stats),
                           out_view, n_out_view);
}

IMPBFF_END_NAMESPACE
