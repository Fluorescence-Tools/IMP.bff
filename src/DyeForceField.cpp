/**
 * \file DyeForceField.cpp
 * \brief A coarse-grained dye/protein system: sites, bonds, and their terms.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/DyeForceField.h>

IMPBFF_BEGIN_NAMESPACE

std::vector<int> DyeForceFieldSystem::get_component_sites(
        const std::string& component) const {
    std::vector<int> out;
    for (std::size_t i = 0; i < sites_.size(); ++i) {
        if (sites_[i].component == component) {
            out.push_back(static_cast<int>(i));
        }
    }
    return out;
}

IMPBFF_END_NAMESPACE
