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

namespace {
std::vector<std::string> with_role(
        const std::map<std::string, FFComponent>& components,
        const std::string& role) {
    std::vector<std::string> out;
    for (std::map<std::string, FFComponent>::const_iterator it =
                 components.begin(); it != components.end(); ++it) {
        if (it->second.role == role) out.push_back(it->first);
    }
    return out;
}
}  // namespace

std::vector<std::string> DyeForceFieldSystem::get_fixed_components() const {
    return with_role(components_, "fixed");
}

std::vector<std::string> DyeForceFieldSystem::get_mobile_components() const {
    return with_role(components_, "mobile");
}

IMPBFF_END_NAMESPACE
