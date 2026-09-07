/**
 *  \file IMP/bff/ComponentTemplate.h
 *  \brief A cgprobe component's feature template, and its mmCIF form.
 *
 * The template names a component's degrees of freedom -- which atoms form each
 * feature, which features are rigid or fixed during MD, which triples are
 * impropers -- and, for a probe, the metadata that orients it (center atom,
 * dipole atoms, charged atoms). `TopologyBuild.h` turns one into topology.
 *
 * The file it is read from and written to is mmCIF, in `_cgprobe_*`
 * categories. Reading goes through IMP's parser (the vendored `ihm_format.h`
 * that `IMP::atom::read_mmcif` uses); writing goes through this module's one
 * writer, `internal/Cif.h`, because IMP has no CIF writer in C++.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#ifndef IMPBFF_COMPONENTTEMPLATE_H
#define IMPBFF_COMPONENTTEMPLATE_H

#include <IMP/bff/bff_config.h>

#include <IMP/bff/Base.h>

#include <map>
#include <ostream>
#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! A cgprobe component feature template.
struct ComponentTemplate {
    std::string name;
    // feature_id -> {rb, md_fixed, feature_type, region_color, atoms}
    struct FeatureAtom {
        std::string name;
        int occurrence;

        IMP_SHOWABLE_INLINE(FeatureAtom, out << "FeatureAtom(" << name << ", "
                                             << occurrence << ")");
    };
    struct Feature {
        bool rb = false;
        bool md_fixed = false;
        std::string feature_type = "dof";
        std::string region_color;
        std::vector<FeatureAtom> atoms;
    };
    std::map<std::string, Feature> features;
    struct Improper {
        std::string center_atom;
        std::string type;

        IMP_SHOWABLE_INLINE(Improper, out << "Improper(" << center_atom << ", "
                                          << type << ")");
    };
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
        const std::string& path, bool with_probe_metadata = false);

//! Read a dye template (component template with metadata).
IMPBFFEXPORT ComponentTemplate read_probe_template_cif(const std::string& path);

//! Return {feature_id: region_color} for features with non-null region_color.
IMPBFFEXPORT std::map<std::string, std::string> region_features(
        const ComponentTemplate& tmpl);

// --------------------------------------------------------------------------
// Rotamer library IO (numpy .npy + text files)

IMPBFF_END_NAMESPACE

#endif // IMPBFF_COMPONENTTEMPLATE_H
