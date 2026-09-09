/**
 * \file StructureTable.cpp
 * \brief The structure record: its columns as arrays, and the residue test.
 *
 * No IMP here. Filling the table from a coordinate file is
 * src/imp/StructureReader.cpp, which is the connection layer's.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#include <IMP/bff/StructureTable.h>

#include <IMP/bff/internal/OutputView.h>

#include <algorithm>
#include <cctype>
#include <cstring>
#include <set>
#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

namespace {

//! The residue names IMP's own hierarchy calls standard, as strings.
/*! Written out rather than derived from `IMP::atom::ResidueType` because the
    question a caller asks is about the *name in the file*, and a name IMP has
    never seen has no ResidueType to ask. */
const std::set<std::string>& standard_residues() {
    static const std::set<std::string> s = {
            // the twenty, plus the two encoded by stop codons
            "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
            "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
            "TYR", "VAL", "SEC", "PYL",
            // protonation and tautomer spellings CHARMM and PDB both use
            "HSD", "HSE", "HSP", "HID", "HIE", "HIP", "CYX", "CYM", "ASH",
            "GLH", "LYN", "MSE",
            // nucleotides, DNA and RNA, in both the one- and two-letter
            // spellings deposited files carry
            "A", "C", "G", "T", "U", "I",
            "DA", "DC", "DG", "DT", "DU", "DI",
            "ADE", "CYT", "GUA", "THY", "URA",
            // the modified nucleotides a deposited RNA carries. They are
            // kept because dropping them cuts a tRNA into fragments, and a
            // structure with holes in its backbone is worse than one
            // carrying a residue whose chemistry this table does not model.
            "2DA", "2DC", "2DG", "2DT",
            "1MA", "1MG", "1MC", "1MT",
            "M2G", "OMG", "OMC", "H2U",
            "PSU", "5MC", "7MG"};
    return s;
}

}  // namespace

bool is_standard_residue(const std::string& res_name) {
    std::string name = res_name;
    while (!name.empty() && name.front() == ' ') name.erase(name.begin());
    while (!name.empty() && name.back() == ' ') name.erase(name.end() - 1);
    std::transform(name.begin(), name.end(), name.begin(), ::toupper);
    return standard_residues().count(name) > 0;
}

void StructureTable::get_xyz(double** out_view, int* n_out_view) const {
    double* out = internal::new_double_view(xyz.size(), out_view, n_out_view);
    if (out != nullptr && !xyz.empty()) {
        std::memcpy(out, xyz.data(), xyz.size() * sizeof(double));
    }
}

void StructureTable::get_radius(double** out_view, int* n_out_view) const {
    double* out =
            internal::new_double_view(radius.size(), out_view, n_out_view);
    if (out != nullptr && !radius.empty()) {
        std::memcpy(out, radius.data(), radius.size() * sizeof(double));
    }
}

void StructureTable::get_mass(double** out_view, int* n_out_view) const {
    double* out = internal::new_double_view(mass.size(), out_view, n_out_view);
    if (out != nullptr && !mass.empty()) {
        std::memcpy(out, mass.data(), mass.size() * sizeof(double));
    }
}

void StructureTable::get_bfactor(double** out_view, int* n_out_view) const {
    double* out =
            internal::new_double_view(bfactor.size(), out_view, n_out_view);
    if (out != nullptr && !bfactor.empty()) {
        std::memcpy(out, bfactor.data(), bfactor.size() * sizeof(double));
    }
}

void StructureTable::get_atom_id(int** out_view_i, int* n_out_view_i) const {
    int* out = internal::new_int_view(atom_id.size(), out_view_i, n_out_view_i);
    if (out != nullptr && !atom_id.empty()) {
        std::memcpy(out, atom_id.data(), atom_id.size() * sizeof(int));
    }
}

void StructureTable::get_res_id(int** out_view_i, int* n_out_view_i) const {
    int* out = internal::new_int_view(res_id.size(), out_view_i, n_out_view_i);
    if (out != nullptr && !res_id.empty()) {
        std::memcpy(out, res_id.data(), res_id.size() * sizeof(int));
    }
}

IMPBFF_END_NAMESPACE
