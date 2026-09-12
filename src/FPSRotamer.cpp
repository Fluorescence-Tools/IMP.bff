#include <IMP/bff/FPSRotamer.h>
#include <IMP/bff/internal/json.h>
#include <IMP/bff/FPS.h>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <sstream>


// -------- from RotamerFps.cpp --------
/**
 * (formerly RotamerFps.cpp, now a section of this file)
 * \brief fps.json and rotamer ensembles.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */




IMPBFF_BEGIN_NAMESPACE

namespace {

//! The first of \p keys the object has as a non-empty value, or "".
/*! The fps.json aliases: a caller asks for every spelling at once and takes
    whichever the writer used. A key present but empty counts as absent, which
    is what an fps file written with blank fields means. */
std::string first_text(const nlohmann::json& payload,
                       const std::vector<std::string>& keys) {
    if (!payload.is_object()) return std::string();
    for (std::size_t i = 0; i < keys.size(); ++i) {
        if (!payload.contains(keys[i]) || payload[keys[i]].is_null()) continue;
        const nlohmann::json& v = payload[keys[i]];
        if (v.is_string()) {
            const std::string s = v.get<std::string>();
            if (!s.empty()) return s;
        } else if (v.is_number()) {
            std::ostringstream out;
            out << v.get<double>();
            return out.str();
        } else if (v.is_boolean()) {
            return v.get<bool>() ? "True" : "False";
        }
    }
    return std::string();
}

std::vector<std::string> keys_of(const char* a, const char* b = NULL,
                                 const char* c = NULL, const char* d = NULL) {
    std::vector<std::string> out;
    out.push_back(a);
    if (b) out.push_back(b);
    if (c) out.push_back(c);
    if (d) out.push_back(d);
    return out;
}

nlohmann::json parse_object(const std::string& text) {
    nlohmann::json out = nlohmann::json::parse(text, NULL, false);
    if (out.is_discarded() || !out.is_object()) return nlohmann::json::object();
    return out;
}

//! One position's entry, or an empty object when the document has no such key.
nlohmann::json position_entry(const nlohmann::json& positions,
                              const std::string& name) {
    if (name.empty() || !positions.contains(name)) {
        return nlohmann::json::object();
    }
    return positions[name];
}

//! A `params` value as a number, or NaN when it is missing or not one.
double param_number(const std::map<std::string, std::string>& params,
                    const std::string& key) {
    const std::map<std::string, std::string>::const_iterator it =
            params.find(key);
    if (it == params.end() || it->second.empty()) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    char* end = NULL;
    const double v = std::strtod(it->second.c_str(), &end);
    if (end == it->second.c_str() || *end != '\0') {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return v;
}

std::string param_text(const std::map<std::string, std::string>& params,
                       const std::string& key) {
    const std::map<std::string, std::string>::const_iterator it =
            params.find(key);
    return it == params.end() ? std::string() : it->second;
}

//! A `params` flag: 1 true, 0 false, -1 when it is not recorded.
int param_flag(const std::map<std::string, std::string>& params,
               const std::string& key) {
    const std::string text = param_text(params, key);
    if (text.empty()) return -1;
    return (text == "True" || text == "true" || text == "1") ? 1 : 0;
}

}  // namespace

FPSRotamerPosition fps_rotamer_position_from_payload(const std::string& name,
                                              const std::string& payload_json) {
    const nlohmann::json payload = parse_object(payload_json);
    FPSRotamerPosition out;
    out.name = name;
    out.chain = first_text(payload, keys_of("chain_identifier", "chain",
                                            "segid"));
    const std::string residue = first_text(
            payload, keys_of("residue_seq_number", "residue", "resid"));
    out.residue = std::atoi(residue.c_str());
    const std::string atom = first_text(payload, keys_of("atom_name"));
    if (!atom.empty()) out.atom_name = atom;
    out.dye = first_text(payload, keys_of("dye_name", "dye", "chromophore"));
    out.library = first_text(payload, keys_of("rotamer_library",
                                              "library_name", "libname",
                                              "library"));
    out.role = first_text(payload, keys_of("role", "label_type"));
    return out;
}

FPSRotamerDistance fps_rotamer_distance_from_payload(
        const std::string& name, const std::string& payload_json,
        const std::string& positions_json) {
    const nlohmann::json payload = parse_object(payload_json);
    const nlohmann::json positions = parse_object(positions_json);
    FPSRotamerDistance out;
    out.name = name;
    out.donor_position = first_text(
            payload, keys_of("position1_name", "donor_position",
                             "donor_position_name", "dye1_position"));
    out.acceptor_position = first_text(
            payload, keys_of("position2_name", "acceptor_position",
                             "acceptor_position_name", "dye2_position"));
    const nlohmann::json donor_entry =
            position_entry(positions, out.donor_position);
    const nlohmann::json acceptor_entry =
            position_entry(positions, out.acceptor_position);

    out.donor = first_text(payload, keys_of("donor", "donor_dye", "dye1",
                                            "dye_1"));
    if (out.donor.empty()) out.donor = first_text(donor_entry, keys_of("dye"));
    out.acceptor = first_text(payload, keys_of("acceptor", "acceptor_dye",
                                               "dye2", "dye_2"));
    if (out.acceptor.empty()) {
        out.acceptor = first_text(acceptor_entry, keys_of("dye"));
    }

    out.libname_1 = first_text(payload, keys_of("libname_1", "donor_library",
                                                "donor_rotamer_library",
                                                "library_1"));
    if (out.libname_1.empty()) {
        out.libname_1 = first_text(donor_entry,
                                   keys_of("rotamer_library", "library"));
    }
    out.libname_2 = first_text(payload, keys_of("libname_2",
                                                "acceptor_library",
                                                "acceptor_rotamer_library",
                                                "library_2"));
    if (out.libname_2.empty()) {
        out.libname_2 = first_text(acceptor_entry,
                                   keys_of("rotamer_library", "library"));
    }
    return out;
}

FPSRotamerSelection read_fps_rotamer(const std::string& path,
                                     const std::string& distance_name) {
    const FPSDocument document = read_fps_json(path);
    const nlohmann::json positions = parse_object(document.positions);
    const nlohmann::json distances = parse_object(document.distances);
    if (distances.empty()) {
        IMP_THROW("No distances found in " << path, ValueException);
    }
    std::string name = distance_name;
    if (name.empty()) name = distances.begin().key();
    if (!distances.contains(name)) {
        IMP_THROW("Distance '" << name << "' not found in " << path,
                  ValueException);
    }

    FPSRotamerSelection out;
    out.distance = fps_rotamer_distance_from_payload(name, distances[name].dump(),
                                                 document.positions);
    if (!positions.contains(out.distance.donor_position)) {
        IMP_THROW("Donor position '" << out.distance.donor_position
                                     << "' not found in " << path,
                  ValueException);
    }
    if (!positions.contains(out.distance.acceptor_position)) {
        IMP_THROW("Acceptor position '" << out.distance.acceptor_position
                                        << "' not found in " << path,
                  ValueException);
    }
    out.donor = fps_rotamer_position_from_payload(
            out.distance.donor_position,
            positions[out.distance.donor_position].dump());
    out.acceptor = fps_rotamer_position_from_payload(
            out.distance.acceptor_position,
            positions[out.distance.acceptor_position].dump());
    out.positions = document.positions;
    out.distances = document.distances;
    // The two free-form sections, merged: what a caller wants of them is
    // "whatever else the file said", and which section a key sat in is the
    // writer's business.
    nlohmann::json extra = parse_object(document.score_sets);
    const nlohmann::json rest = parse_object(document.extra);
    for (nlohmann::json::const_iterator it = rest.begin(); it != rest.end();
         ++it) {
        extra[it.key()] = it.value();
    }
    out.extra = extra.dump();
    return out;
}

std::string fps_rotamer_position_payload(const std::string& chain, int residue,
                                     const std::string& library,
                                     const std::string& atom_name,
                                     const std::string& dye,
                                     double temperature, int electrostatic,
                                     const std::string& potential) {
    nlohmann::json payload;
    payload["chain_identifier"] = chain;
    payload["residue_seq_number"] = residue;
    payload["atom_name"] = atom_name;
    payload["simulation_type"] = SIMULATION_TYPE_R1;
    payload["rotamer_library"] = library;
    if (!dye.empty()) payload["dye_name"] = dye;
    if (!std::isnan(temperature)) payload["temperature"] = temperature;
    if (electrostatic >= 0) payload["electrostatic"] = electrostatic != 0;
    if (!potential.empty()) payload["potential"] = potential;
    return payload.dump();
}

std::string fps_rotamer_positions_payload(
        const std::map<std::string, ProbeRotamerEnsemble>& ensembles,
        const std::string& atom_name) {
    nlohmann::json out = nlohmann::json::object();
    for (std::map<std::string, ProbeRotamerEnsemble>::const_iterator it =
                 ensembles.begin();
         it != ensembles.end(); ++it) {
        const ProbeRotamerEnsemble& e = it->second;
        const std::map<std::string, std::string> params = e.get_params();
        out[it->first] = nlohmann::json::parse(fps_rotamer_position_payload(
                e.get_chain(), e.get_residue(), e.get_library(), atom_name,
                param_text(params, "dye"), param_number(params, "temperature"),
                param_flag(params, "electrostatic"),
                param_text(params, "potential")));
    }
    return out.dump();
}

std::string fps_rotamer_distances_from_ensembles(
        const std::map<std::string, ProbeRotamerEnsemble>& ensembles,
        const std::vector<std::pair<std::string, std::string> >& pairs,
        double forster_radius, const std::string& distance_type, double error,
        double error_fraction, const std::string& kappa2) {
    if (kappa2 != "isotropic" && kappa2 != "dipoles") {
        IMP_THROW("kappa2 must be 'isotropic' or 'dipoles', not '" << kappa2
                                                                   << "'",
                  ValueException);
    }
    if (distance_type != "RDAMean" && distance_type != "RDAMeanE" &&
        distance_type != "Rmp") {
        IMP_THROW("unknown distance_type '" << distance_type << "'",
                  ValueException);
    }

    nlohmann::json out = nlohmann::json::object();
    for (std::size_t p = 0; p < pairs.size(); ++p) {
        const std::string& name1 = pairs[p].first;
        const std::string& name2 = pairs[p].second;
        std::map<std::string, ProbeRotamerEnsemble>::const_iterator i1 =
                ensembles.find(name1);
        std::map<std::string, ProbeRotamerEnsemble>::const_iterator i2 =
                ensembles.find(name2);
        if (i1 == ensembles.end() || i2 == ensembles.end()) {
            IMP_THROW("no ensemble for the pair (" << name1 << ", " << name2
                                                   << ")",
                      ValueException);
        }
        const ProbeRotamerEnsemble& e1 = i1->second;
        const ProbeRotamerEnsemble& e2 = i2->second;

        // The isotropic geometry drops the dipoles, which is what puts
        // kappa2 = 2/3 on every pair.
        const FRETPairGeometry geometry =
                e1.pair_geometry(e2, kappa2 == "dipoles");
        double value = 0.0;
        if (distance_type == "Rmp") {
            value = e1.dRmp(e2);
        } else if (distance_type == "RDAMean") {
            for (std::size_t k = 0; k < geometry.R.size(); ++k) {
                value += geometry.R[k] * geometry.weight[k];
            }
        } else {
            const FRETPairEfficiencies eff =
                    fret_pair_efficiencies(geometry, forster_radius);
            const double mean_e = eff.static_efficiency;
            if (mean_e <= 0.0) {
                for (std::size_t k = 0; k < geometry.R.size(); ++k) {
                    value += geometry.R[k] * geometry.weight[k];
                }
            } else if (mean_e >= 1.0) {
                value = 0.0;
            } else {
                value = forster_radius *
                        std::pow(1.0 / mean_e - 1.0, 1.0 / 6.0);
            }
        }
        const double err = error >= 0.0 ? error : error_fraction * value;
        nlohmann::json entry;
        entry["position1_name"] = name1;
        entry["position2_name"] = name2;
        entry["distance_type"] = distance_type;
        entry["distance"] = value;
        entry["error_neg"] = err;
        entry["error_pos"] = err;
        entry["Forster_radius"] = forster_radius;
        out[name1 + "_" + name2] = entry;
    }
    return out.dump();
}

std::map<std::string, ProbeRotamerEnsemble> fps_rotamer_ensembles(
        const std::string& fps_json, const std::string& structure,
        const std::map<std::string, std::string>& library_map,
        const ProbeRotamerSiteOptions& options, int frame_index) {
    const FPSDocument document = read_fps_json(fps_json);
    const nlohmann::json positions = parse_object(document.positions);
    std::map<std::string, ProbeRotamerEnsemble> out;
    if (positions.empty()) return out;

    const std::vector<ProteinFrame> frames =
            load_protein_frames(structure, frame_index + 1);
    if (frame_index >= static_cast<int>(frames.size())) {
        IMP_THROW(structure << " has " << frames.size() << " frame(s), frame "
                            << frame_index << " requested",
                  ValueException);
    }
    const ProteinFrame& frame = frames[frame_index];

    // One read per library and not one per position: a donor and an acceptor
    // labelled with the same dye are the common case.
    std::map<std::string, ProbeRotamerLibrary> libraries;
    for (nlohmann::json::const_iterator it = positions.begin();
         it != positions.end(); ++it) {
        const FPSRotamerPosition position =
                fps_rotamer_position_from_payload(it.key(), it.value().dump());
        std::string library_name = position.library;
        const std::map<std::string, std::string>::const_iterator mapped =
                library_map.find(it.key());
        if (mapped != library_map.end()) library_name = mapped->second;
        // A position that names no library is skipped and not guessed at.
        if (library_name.empty()) continue;
        if (libraries.find(library_name) == libraries.end()) {
            libraries[library_name] = load_probe_rotamer_library(library_name);
        }
        out.insert(std::make_pair(
                it.key(),
                ProbeRotamerEnsemble::from_frame(frame, position.chain,
                                            position.residue,
                                            libraries[library_name], options,
                                            it.key())));
    }
    return out;
}

void write_fps_rotamer(const std::string& path,
                       const std::string& positions_json,
                       const std::string& distances_json,
                       const std::string& merge_into, bool validate) {
    nlohmann::json positions = nlohmann::json::object();
    nlohmann::json distances = nlohmann::json::object();
    nlohmann::json score_sets = nlohmann::json::object();
    nlohmann::json extra = nlohmann::json::object();
    if (!merge_into.empty()) {
        const FPSDocument document = read_fps_json(merge_into);
        positions = parse_object(document.positions);
        distances = parse_object(document.distances);
        score_sets = parse_object(document.score_sets);
        extra = parse_object(document.extra);
    }
    const nlohmann::json new_positions = parse_object(positions_json);
    for (nlohmann::json::const_iterator it = new_positions.begin();
         it != new_positions.end(); ++it) {
        positions[it.key()] = it.value();
    }
    const nlohmann::json new_distances = parse_object(distances_json);
    for (nlohmann::json::const_iterator it = new_distances.begin();
         it != new_distances.end(); ++it) {
        distances[it.key()] = it.value();
    }
    write_fps_json(path, positions.dump(), distances.dump(), score_sets.dump(),
                   extra.dump(), validate);
}

FRETRotamer fps_rotamer_fret(const std::string& fps_path,
                                  const std::string& protein,
                                  const std::string& distance_name,
                                  const std::string& output_prefix,
                                  bool fixed_R0, double r0, double temperature,
                                  bool electrostatic) {
    const FPSRotamerSelection selection = read_fps_rotamer(fps_path,
                                                           distance_name);
    std::vector<int> residues;
    residues.push_back(selection.donor.residue);
    residues.push_back(selection.acceptor.residue);
    std::vector<std::string> chains;
    chains.push_back(selection.donor.chain);
    chains.push_back(selection.acceptor.chain);

    const std::string donor = selection.distance.donor.empty()
                                      ? selection.donor.dye
                                      : selection.distance.donor;
    const std::string acceptor = selection.distance.acceptor.empty()
                                         ? selection.acceptor.dye
                                         : selection.distance.acceptor;
    const std::string libname_1 = selection.distance.libname_1.empty()
                                          ? selection.donor.library
                                          : selection.distance.libname_1;
    const std::string libname_2 = selection.distance.libname_2.empty()
                                          ? selection.acceptor.library
                                          : selection.distance.libname_2;
    return FRETRotamer(protein, residues, chains, donor, acceptor, libname_1,
                       libname_2, temperature, electrostatic, "lj", true, 0.5,
                       1.0, fixed_R0, r0, "", 0.05, output_prefix);
}

IMPBFF_END_NAMESPACE
