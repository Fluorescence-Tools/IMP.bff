/**
 *  \file IMP/bff/CifIO.h
 *  \brief Writing compact FF topology mmCIF tables, and the template/rotamer
 *         library IO that was Python.
 *
 * The reader is \ref ForceFieldCIF.h (through the ihm C reader). This file is
 * the writer and the remaining io/cif.py surface: the FF system writer, the
 * component/dye template reader and writer, the rotamer library numpy/text
 * IO, and the string utilities (site-id splitting, range compression).
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_CIFIO_H
#define IMPBFF_CIFIO_H

#include <IMP/bff/bff_config.h>
#include <IMP/bff/DyeForceField.h>

#include <string>
#include <vector>
#include <map>
#include <utility>

IMPBFF_BEGIN_NAMESPACE

// --------------------------------------------------------------------------
// String utilities (site-id splitting, range compression)
// --------------------------------------------------------------------------

//! Split a site id like "CX4/S1" into (prefix, serial). Returns ("", -1) on failure.
IMPBFFEXPORT std::pair<std::string, int> split_site_id(const std::string& site_id);

//! Expand a group range "A1".."A3" into ["A1","A2","A3"].
IMPBFFEXPORT std::vector<std::string> expand_site_range(
        const std::string& start_id, const std::string& end_id);

//! Compress a list of integers into sorted contiguous (start, end) pairs.
IMPBFFEXPORT std::vector<std::pair<int, int>> compress_int_ranges(
        const std::vector<int>& nos);

// --------------------------------------------------------------------------
// Force-field system writer
// --------------------------------------------------------------------------

//! Write a DyeForceFieldSystem to mmCIF.
/*!
    Mirrors the reader in ForceFieldCIF.h. Writes the eighteen _ff_* categories
    plus _atom_site (from the component MOL2/PDB files) and _flr_probe_list.
    Site numbers are assigned and group members are written as integer ranges.
*/
IMPBFFEXPORT void _write_dye_forcefield_cif(const std::string& path,
                                             const DyeForceFieldSystem& system);

// --------------------------------------------------------------------------
// Component/dye template CIF reader and writer
// --------------------------------------------------------------------

//! A cgdye component feature template.
struct ComponentTemplate {
    std::string name;
    // feature_id -> {rb, md_fixed, feature_type, region_color, atoms}
    struct FeatureAtom { std::string name; int occurrence; };
    struct Feature {
        bool rb = false;
        bool md_fixed = false;
        std::string feature_type = "dof";
        std::string region_color;
        std::vector<FeatureAtom> atoms;
    };
    std::map<std::string, Feature> features;
    struct Improper { std::string center_atom; std::string type; };
    std::vector<Improper> impropers;
    // dye metadata
    std::string center_atom;
    std::string dipole_atom_1;
    std::string dipole_atom_2;
    std::vector<std::string> positive_atoms;
    std::vector<std::string> negative_atoms;

    void show(std::ostream& out) const {
        out << "ComponentTemplate(" << name << ")";
    }
};

//! A list of component templates (the SWIG plural type).
typedef std::vector<ComponentTemplate> ComponentTemplates;

//! Read a component template from mmCIF.
IMPBFFEXPORT ComponentTemplate read_component_template_cif(
        const std::string& path, bool with_dye_metadata = false);

//! Write a component template to mmCIF.
IMPBFFEXPORT void write_component_template_cif(
        const std::string& path, const ComponentTemplate& tmpl);

//! Read a dye template (component template with metadata).
IMPBFFEXPORT ComponentTemplate read_dye_template_cif(const std::string& path);

//! Write a dye template to mmCIF.
IMPBFFEXPORT void write_dye_template_cif(
        const std::string& path, const ComponentTemplate& tmpl);

//! Return {feature_id: region_color} for features with non-null region_color.
IMPBFFEXPORT std::map<std::string, std::string> region_features(
        const ComponentTemplate& tmpl);

// --------------------------------------------------------------------------
// Rotamer library IO (numpy .npy + text files)
// --------------------------------------------------------------------------

//! A rotamer library read from disk.
struct RotamerLibraryData {
    std::vector<int> id;
    std::vector<double> weight;
    std::vector<std::string> atom_names;
    // flat coords: n_rotamers * n_atoms * 3
    std::vector<double> coords;
    int n_rotamers = 0;
    int n_atoms = 0;

    void show(std::ostream& out) const {
        out << "RotamerLibraryData(" << n_rotamers << " x " << n_atoms << ")";
    }
};

//! A list of rotamer library data (the SWIG plural type).
typedef std::vector<RotamerLibraryData> RotamerLibraryDatas;

//! Read a rotamer library from numpy/text files.
IMPBFFEXPORT RotamerLibraryData read_rotamer_library(const std::string& path);

//! Write a rotamer library to numpy/text files.
IMPBFFEXPORT void write_rotamer_library(const std::string& path,
                                          const RotamerLibraryData& lib);

//! Normalize rotamer weights to sum to 1.0 in-place.
IMPBFFEXPORT void normalize_weights(RotamerLibraryData& lib);

IMPBFF_END_NAMESPACE

#endif // IMPBFF_CIFIO_H
