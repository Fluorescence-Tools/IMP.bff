/**
 *  \file IMP/bff/FPSReaderWriter.h
 *  \brief Simple class to read and (write) FPS.json files
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2022 IMP Inventors. All rights reserved.
 *
 */

#ifndef IMPBFF_FPSREADERWRITER_H
#define IMPBFF_FPSREADERWRITER_H


#include <IMP/bff/bff_config.h>

#include <IMP/bff/AV.h>
#include <IMP/bff/internal/json.h>

#include <map>
#include <vector>
#include <algorithm>


IMPBFF_BEGIN_NAMESPACE

class FPSReaderWriter{

private:

    /// Converts distance names (also used in fps.json files) to internal distance types
    std::map<std::string, int> DyePairMeasure_name_to_type = {
            {"RDAMean",    DYE_PAIR_DISTANCE_MEAN},
            {"RDAMeanE",   DYE_PAIR_DISTANCE_E},
            {"Rmp",        DYE_PAIR_DISTANCE_MP},
            {"Efficiency", DYE_PAIR_EFFICIENCY},
            {"pRDA",    DYE_PAIR_DISTANCE_DISTRIBUTION}
    };

    std::string fn_json_;
    std::string score_set_;
    nlohmann::json fps_json_;
    std::map<std::string, AVPairDistanceMeasurement> distances_;

public:

    void read_json();

    void update_distances();

    void set_filename(std::string fn);

    std::string get_filename(){
        return fn_json_;
    }

    std::map<std::string, AVPairDistanceMeasurement> get_distances(){
        return distances_;
    }

    nlohmann::json get_used_positions();

    /// Setup the FPSReaderWriter
    /**
     * @param[in] fps_json_fn Filename of fps.json file
     */
    FPSReaderWriter(std::string fps_json_fn, std::string score_set = "");


};

IMPBFF_END_NAMESPACE

#endif //IMPBFF_FPSREADERWRITER_H
