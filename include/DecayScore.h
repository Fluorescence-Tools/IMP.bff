/**
 * \file IMP/bff/DecayScore.h
 * \brief Scoring of model fluorescence decay.
 *
 * \authors Thomas-Otavio Peulen
 * Copyright 2007-2022 IMP Inventors. All rights reserved.
 *
 */

#ifndef IMPBFF_DECAYSCORE_H
#define IMPBFF_DECAYSCORE_H

#include <IMP/bff/bff_config.h>

#include <string>
#include <iostream>

#include <IMP/bff/internal/PhotonStatistics.h>
#include <IMP/bff/DecayRange.h>
#include <IMP/bff/DecayCurve.h>


IMPBFF_BEGIN_NAMESPACE


class IMPBFFEXPORT DecayScore : public DecayRange {

private:

    DecayCurve* _model = nullptr;
    DecayCurve* _data = nullptr;
    DecayCurve* _default_model = nullptr;
    DecayCurve* _default_data = nullptr;
    std::vector<double> _weighted_residuals;
    std::string _score_type= "default";

public:

    DecayCurve* get_model();

    void set_model(DecayCurve* v);

    DecayCurve* get_data();

    void set_data(DecayCurve* v);

private:

    void update_weighted_residuals();

public:

    std::vector<double>& get_weighted_residuals();

    double get_score(int start=0, int stop=-1, const char* score_type = nullptr);

    void set_score_type(std::string v);

    std::string get_score_type();

    void set(
            DecayCurve* model, DecayCurve* data,
            std::string score_type = "poisson",
            int start=0, int stop=-1
    );


    /// Evaluate and return the score
    double score(DecayCurve* model = nullptr);

    /*!
     * score_type:
     * poisson, pearson, gauss, cnp, sswr (sum of weighted squared residuals)
     */
    DecayScore(
            DecayCurve* model= nullptr,
            DecayCurve* data= nullptr,
            std::string score_type= "poisson",
            int start=0, int stop=-1
    );

     ~DecayScore() override {
        delete _default_data;
        delete _default_model;
    }

};

IMPBFF_END_NAMESPACE


#endif //IMPBFF_DECAYSCORE_H
