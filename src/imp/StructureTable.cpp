/**
 * \file StructureTable.cpp
 * \brief A coordinate file read into flat columns, through IMP's readers.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#include <IMP/bff/StructureTable.h>

#include <IMP/bff/internal/OutputView.h>

#include <IMP/Model.h>
#include <IMP/atom/Atom.h>
#include <IMP/atom/Chain.h>
#include <IMP/atom/Hierarchy.h>
#include <IMP/atom/Mass.h>
#include <IMP/atom/Residue.h>
#include <IMP/atom/element.h>
#include <IMP/atom/mmcif.h>
#include <IMP/atom/pdb.h>
#include <IMP/core/XYZR.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <map>
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

//! IMP prints an unrecognised atom type as `HET: C3 `; callers want `C3`.
/*! The prefix and the quoting are both IMP's, and a five-character field
    truncates `"HET: N  "` to `"HET: "` -- every ligand atom then shares one
    name and no selection can see ligands at all. Stripped here so that no
    caller has to know the spelling. */
std::string bare_atom_name(IMP::atom::Atom atom) {
    std::string name = atom.get_atom_type().get_string();
    while (!name.empty() && (name.front() == ' ' || name.back() == ' ')) {
        if (name.front() == ' ') name.erase(name.begin());
        if (!name.empty() && name.back() == ' ') name.erase(name.end() - 1);
    }
    if (name.size() >= 2 && name.front() == '"' && name.back() == '"') {
        name = name.substr(1, name.size() - 2);
    }
    if (name.size() >= 4) {
        std::string head = name.substr(0, 4);
        std::transform(head.begin(), head.end(), head.begin(), ::toupper);
        if (head == "HET:") name = name.substr(4);
    }
    while (!name.empty() && (name.front() == ' ' || name.back() == ' ')) {
        if (name.front() == ' ') name.erase(name.begin());
        if (!name.empty() && name.back() == ' ') name.erase(name.end() - 1);
    }
    return name;
}

std::string lower_extension(const std::string& path) {
    const std::size_t dot = path.find_last_of('.');
    if (dot == std::string::npos) return std::string();
    std::string ext = path.substr(dot);
    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
    return ext;
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

StructureTable read_structure_table(const std::string& path, bool keep_water,
                                    bool only_standard_residues,
                                    bool radius_no_interaction) {
    const std::string ext = lower_extension(path);
    if (ext == ".gz" || ext == ".bz2" || ext == ".zip" || ext == ".xz") {
        IMP_THROW("read_structure_table: " << path
                  << " is compressed, and IMP's readers do not decompress -- "
                     "a gzip stream read as text yields an empty structure "
                     "rather than an error. Decompress it first.",
                  IMP::ValueException);
    }
    const bool is_pdb = (ext == ".pdb" || ext == ".ent");
    const bool is_cif = (ext == ".cif" || ext == ".mmcif");
    if (!is_pdb && !is_cif) {
        IMP_THROW("read_structure_table: " << path
                  << " is not a .pdb/.ent or .cif/.mmcif file",
                  IMP::ValueException);
    }

    IMP_NEW(IMP::Model, model, ());
    // NonAlternative keeps everything but alternate locations; NonWater is
    // the same minus solvent. Alternate locations go either way, so a
    // multi-conformer file does not yield overlapping copies of a residue.
    IMP::Pointer<IMP::atom::PDBSelector> selector;
    if (keep_water) {
        selector = new IMP::atom::NonAlternativePDBSelector();
    } else {
        selector = new IMP::atom::NonWaterPDBSelector();
    }

    IMP::atom::Hierarchy hierarchy =
            is_pdb ? IMP::atom::read_pdb(path, model, selector)
                   : IMP::atom::read_mmcif(path, model, selector);

    const IMP::atom::Hierarchies atoms =
            IMP::atom::get_by_type(hierarchy, IMP::atom::ATOM_TYPE);

    StructureTable table;
    const std::size_t n = atoms.size();
    table.xyz.reserve(n * 3);
    table.radius.reserve(n);
    table.mass.reserve(n);
    table.bfactor.reserve(n);
    table.atom_id.reserve(n);
    table.res_id.reserve(n);
    table.chain.reserve(n);
    table.res_name.reserve(n);
    table.atom_name.reserve(n);
    table.element.reserve(n);

    // The radius at which the 6-12 potential is zero rather than its minimum:
    // E = eij ((Rmin/rij)^12 - 2 (Rmin/rij)^6) is zero at 2^(-1/6) Rmin.
    const double radius_scaling =
            radius_no_interaction ? std::pow(2.0, -1.0 / 6.0) : 1.0;

    // A chain decorator per *atom* rebuilds the same object for every atom of
    // a chain, and there are a handful of chains against tens of thousands of
    // atoms; the element name is likewise one of a dozen strings. Both are
    // cached by the key IMP gives them.
    const IMP::atom::ElementTable& element_table =
            IMP::atom::get_element_table();
    std::map<IMP::ParticleIndex, std::string> chain_ids;
    std::map<int, std::string> element_names;

    for (unsigned int i = 0; i < atoms.size(); ++i) {
        IMP::atom::Atom atom(atoms[i]);
        IMP::atom::Hierarchy parent = atom.get_parent();
        if (!IMP::atom::Residue::get_is_setup(parent)) continue;
        IMP::atom::Residue residue(parent);

        const std::string res_name = residue.get_residue_type().get_string();
        if (only_standard_residues && !is_standard_residue(res_name)) continue;

        IMP::atom::Hierarchy chain_particle = residue.get_parent();
        const IMP::ParticleIndex chain_key =
                chain_particle.get_particle_index();
        std::map<IMP::ParticleIndex, std::string>::const_iterator chain_it =
                chain_ids.find(chain_key);
        if (chain_it == chain_ids.end()) {
            std::string id;
            if (IMP::atom::Chain::get_is_setup(chain_particle)) {
                id = IMP::atom::Chain(chain_particle).get_id();
            }
            chain_it = chain_ids.insert(std::make_pair(chain_key, id)).first;
        }

        const IMP::atom::Element element = atom.get_element();
        std::map<int, std::string>::const_iterator element_it =
                element_names.find(static_cast<int>(element));
        if (element_it == element_names.end()) {
            element_it = element_names
                                 .insert(std::make_pair(
                                         static_cast<int>(element),
                                         element_table.get_name(element)))
                                 .first;
        }

        IMP::core::XYZR xyzr(atoms[i]);
        const IMP::algebra::Vector3D v = xyzr.get_coordinates();
        table.xyz.push_back(v[0]);
        table.xyz.push_back(v[1]);
        table.xyz.push_back(v[2]);
        table.radius.push_back(xyzr.get_radius() * radius_scaling);
        table.mass.push_back(IMP::atom::Mass(atoms[i]).get_mass());
        table.bfactor.push_back(atom.get_temperature_factor());
        table.atom_id.push_back(atom.get_input_index());
        table.res_id.push_back(residue.get_index());
        table.chain.push_back(chain_it->second);
        table.res_name.push_back(res_name);
        table.atom_name.push_back(bare_atom_name(atom));
        table.element.push_back(element_it->second);
    }

    return table;
}

IMPBFF_END_NAMESPACE
