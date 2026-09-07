/**
 * \file Mol2IO.cpp
 * \brief Reading a MOL2 component into force-field sites and bonds.
 *
 * Copyright 2007-2022 IMP Inventors. All rights reserved.
 */

#include <IMP/bff/Mol2IO.h>


#include <algorithm>
#include <cctype>
#include <fstream>
#include <map>
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

    // The Tripos file, walked the way IMP::atom::read_mol2 walks it (PRD-137
    // step 5c; test/io/test_mol2_matches_imp.py pins the answer): a MOLECULE
    // block's third header line names the residue type ("UNK" when blank);
    // an ATOM section runs to an empty line or the next '@' record, every
    // line an atom (IMP's AllMol2Selector); a BOND section likewise, a bond
    // kept only when both atoms are in the molecule. Coordinates are read as
    // doubles, as IMP's Float is.
    std::ifstream in(path.c_str());
    if (!in) IMP_THROW("Cannot open " << path, IOException);
    std::string line;
    std::string resname = "UNK";
    std::map<int, std::size_t> serial_to_atom;    // this molecule's atoms
    std::vector<std::pair<int, int> > pairs;
    while (std::getline(in, line)) {
        if (line.compare(0, 17, "@<TRIPOS>MOLECULE") == 0) {
            resname = "UNK";
            serial_to_atom.clear();
            for (int i = 0; i < 6; ++i) {
                if (in.peek() == '@') break;
                std::string header;
                if (!std::getline(in, header) || header.empty()) break;
                if (i == 2) {
                    std::istringstream fields(header);
                    std::string mol_type;
                    if (fields >> mol_type && !mol_type.empty()) resname = mol_type;
                }
            }
        } else if (line.compare(0, 13, "@<TRIPOS>ATOM") == 0) {
            while (in.peek() != '@' && std::getline(in, line) && !line.empty()) {
                std::istringstream fields(line);
                Mol2Atom site;
                std::string type_field;
                if (!(fields >> site.serial >> site.atom_name >> site.x >> site.y >> site.z)) continue;
                fields >> type_field;
                if (site.atom_name.empty()) {
                    // IMP's spelling of the Sybyl type, its prefixes stripped:
                    // ".ar"/".am" cut at the dot, the first dot erased
                    std::string n = type_field;
                    if (n.find(".ar") != std::string::npos || n.find(".am") != std::string::npos) {
                        n = n.substr(0, n.find('.'));
                    }
                    if (n.find('.') != std::string::npos) n.erase(n.find('.'), 1);
                    site.atom_name = n;
                }
                site.component = component;
                site.resname = resname.empty() ? "UNK" : resname;
                site.element = element_from_atom_name(site.atom_name);
                serial_to_atom[site.serial] = out.atoms.size();
                out.atoms.push_back(site);
            }
        } else if (line.compare(0, 13, "@<TRIPOS>BOND") == 0) {
            while (in.peek() != '@' && std::getline(in, line) && !line.empty()) {
                std::istringstream fields(line);
                int id = 0, s1 = 0, s2 = 0;
                if (!(fields >> id >> s1 >> s2)) continue;
                if (serial_to_atom.find(s1) == serial_to_atom.end() ||
                    serial_to_atom.find(s2) == serial_to_atom.end()) continue;
                pairs.push_back(s1 <= s2 ? std::make_pair(s1, s2) : std::make_pair(s2, s1));
            }
        }
    }
    std::sort(pairs.begin(), pairs.end());
    pairs.erase(std::unique(pairs.begin(), pairs.end()), pairs.end());
    out.bonds = pairs;

    return out;
}

IMPBFF_END_NAMESPACE
