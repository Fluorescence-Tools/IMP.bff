/**
 * \file PETQuenching.cpp
 * \brief Photoinduced electron transfer: which moieties quench, and how hard.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#include <IMP/bff/PETQuenching.h>

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
