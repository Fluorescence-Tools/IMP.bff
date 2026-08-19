/**
 * \file Mol2IO.cpp
 * \brief Reading a MOL2 component into force-field sites and bonds.
 *
 * Copyright 2007-2022 IMP Inventors. All rights reserved.
 */

#include <IMP/bff/Mol2IO.h>

#include <IMP/Model.h>
#include <IMP/atom/Atom.h>
#include <IMP/atom/Hierarchy.h>
#include <IMP/atom/Residue.h>
#include <IMP/atom/bond_decorators.h>
#include <IMP/atom/mol2.h>
#include <IMP/core/XYZ.h>
#include <IMP/log.h>

#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>

IMPBFF_BEGIN_NAMESPACE

std::string element_from_atom_name(const std::string& atom_name) {
    for (size_t i = 0; i < atom_name.size(); ++i) {
        if (std::isalpha((unsigned char) atom_name[i])) {
            return std::string(1, (char) std::toupper((unsigned char) atom_name[i]));
        }
        break;                       // only a *leading* alphabetic run counts
    }
    return "C";
}

std::map<int, std::string> read_mol2_atom_names(const std::string& path) {
    std::map<int, std::string> names;
    std::ifstream in(path.c_str());
    if (!in) return names;

    bool in_atom = false;
    std::string line;
    while (std::getline(in, line)) {
        if (line.compare(0, 13, "@<TRIPOS>ATOM") == 0) { in_atom = true; continue; }
        if (line.compare(0, 9, "@<TRIPOS>") == 0) { in_atom = false; continue; }
        if (!in_atom) continue;
        std::istringstream fields(line);
        int serial;
        std::string name;
        if (fields >> serial >> name) names[serial] = name;
    }
    return names;
}

Mol2Component read_mol2_component(const std::string& path,
                                  const std::string& component) {
    Mol2Component out;

    IMP_NEW(Model, model, ());
    const LogLevel previous = IMP::get_log_level();
    IMP::set_log_level(SILENT);
    atom::Hierarchy hierarchy;
    try {
        hierarchy = atom::read_mol2(path, model);
    } catch (...) {
        IMP::set_log_level(previous);
        throw;
    }
    IMP::set_log_level(previous);

    const std::map<int, std::string> names = read_mol2_atom_names(path);

    const atom::Hierarchies atoms = atom::get_by_type(hierarchy, atom::ATOM_TYPE);
    for (size_t i = 0; i < atoms.size(); ++i) {
        atom::Atom a(atoms[i]);
        Mol2Atom site;
        site.serial = a.get_input_index();
        site.component = component;

        std::map<int, std::string>::const_iterator it = names.find(site.serial);
        if (it != names.end() && !it->second.empty()) {
            site.atom_name = it->second;
        } else {
            // fall back to the TRIPOS type with IMP's prefixes stripped
            std::string t = a.get_atom_type().get_string();
            const char* strip[] = {"HET: ", "HET:", "ATOM: ", "ATOM:"};
            for (int s = 0; s < 4; ++s) {
                const size_t at = t.find(strip[s]);
                if (at != std::string::npos) t.erase(at, std::string(strip[s]).size());
            }
            const size_t b = t.find_first_not_of(" \t");
            const size_t e = t.find_last_not_of(" \t");
            site.atom_name = (b == std::string::npos) ? "" : t.substr(b, e - b + 1);
        }

        atom::Residue residue = atom::get_residue(a, true);
        site.resname = residue ? residue.get_residue_type().get_string() : "UNK";
        if (site.resname.empty()) site.resname = "UNK";
        site.element = element_from_atom_name(site.atom_name);

        const algebra::Vector3D c = core::XYZ(atoms[i]).get_coordinates();
        site.x = c[0]; site.y = c[1]; site.z = c[2];
        out.atoms.push_back(site);
    }

    const Particles bonds = atom::get_internal_bonds(hierarchy);
    std::vector<std::pair<int, int> > pairs;
    for (size_t i = 0; i < bonds.size(); ++i) {
        atom::Bond b(bonds[i]);
        const int s1 = atom::Atom(b.get_bonded(0).get_particle()).get_input_index();
        const int s2 = atom::Atom(b.get_bonded(1).get_particle()).get_input_index();
        pairs.push_back(s1 <= s2 ? std::make_pair(s1, s2) : std::make_pair(s2, s1));
    }
    std::sort(pairs.begin(), pairs.end());
    pairs.erase(std::unique(pairs.begin(), pairs.end()), pairs.end());
    out.bonds = pairs;

    return out;
}

IMPBFF_END_NAMESPACE
