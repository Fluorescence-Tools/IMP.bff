/**
 *  \file IMP/bff/FPSReaderWriter.h
 *  \brief Simple class to read and (write) FPS.json files
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2022 IMP Inventors. All rights reserved.
 *
 */
#include <IMP/bff/internal/FPSReaderWriter.h>

IMPBFF_BEGIN_NAMESPACE


FPSReaderWriter::FPSReaderWriter(
        std::string fps_json_fn, std::string score_set):
        fn_json_(fps_json_fn), score_set_(score_set)
{
    set_filename(fps_json_fn);
}



void FPSReaderWriter::set_filename(std::string fn){
    fn_json_ = fn;
    read_json();
    update_distances();
}

void FPSReaderWriter::read_json(){
    std::ifstream in(fn_json_);
    fps_json_ = nlohmann::json::parse(in);
}

void FPSReaderWriter::update_distances(){
    std::map<std::string, IMP::bff::AVPairDistanceMeasurement> distance_map;

    // select the used distances based on the score set
    nlohmann::json distances;
    if(score_set_ == ""){
        distances = fps_json_["Distances"];
    } else{
        nlohmann::json selected_distances = fps_json_["χ²"][score_set_]["distances"];
        for(auto d : selected_distances) {
            std::string key = d;
            distances[key] = fps_json_["Distances"][key];
        }
    }

    // create AVPairDistanceMeasurement objects
    for(auto it = distances.begin(); it != distances.end(); it++){
        std::string key = it.key();
        nlohmann::json distance = distances[key];
        AVPairDistanceMeasurement d;
        d.position_1 = distance.value("position1_name", "");
        d.position_2 = distance.value("position2_name", "");
        d.forster_radius = distance.value("Forster_radius", 52.0);
        d.distance = distance.value("distance", -1.0);
        d.error_neg = distance.value("error_neg", -1.0);
        d.error_pos = distance.value("error_pos", -1.0);
        d.distance_type = DyePairMeasure_name_to_type[
                distance.value("distance_type", "RDAMeanE")
        ];
        distance_map[key] = d;
    }

    distances_ = distance_map;
}



nlohmann::json FPSReaderWriter::get_used_positions(){
    // Check which positions are used for scoring and
    // only decorate labeling sites of used positions
    std::set<std::string> used_positions_names;
    for(auto it = distances_.begin(); it != distances_.end(); it++) {
        used_positions_names.emplace(it->second.position_1);
        used_positions_names.emplace(it->second.position_2);
    }
    nlohmann::json used_positions;
    for(auto &position_name: used_positions_names){
        used_positions[position_name] = fps_json_["Positions"][position_name];
    }
    return used_positions;
}

IMPBFF_END_NAMESPACE
