/**
 * \file IMP/bff/DecayScore.h
 * \brief Scoring of model fluorescence decay.
 *
 * \authors Thomas-Otavio Peulen
 * Copyright 2007-2023 IMP Inventors. All rights reserved.
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

/**
 * \class DecayScore
 * \brief Class for scoring model fluorescence decay.
 *
 * This class provides functionality for scoring the fluorescence decay of a model
 * against experimental data. It calculates the score based on various scoring types,
 * such as Poisson, Pearson, Gauss, CNP, and SSWR (Sum of Weighted Squared Residuals).
 * The score can be calculated for a specific range of data points.
 */
class IMPBFFEXPORT DecayScore : public DecayRange {

private:
    DecayCurve* _model = nullptr; /**< Pointer to the model decay curve */
    DecayCurve* _data = nullptr; /**< Pointer to the experimental data decay curve */
    DecayCurve* _default_model = nullptr; /**< Default model decay curve */
    DecayCurve* _default_data = nullptr; /**< Default experimental data decay curve */
    std::vector<double> _weighted_residuals; /**< Vector of weighted residuals */
    std::string _score_type = "default"; /**< Type of score to calculate */
public:
    /**
     * Get the model decay curve.
     * \return Pointer to the model decay curve
     */
    DecayCurve* get_model();
    /**
     * Set the model decay curve.
     * \param v Pointer to the model decay curve
     */
    void set_model(DecayCurve* v);
    /**
     * Get the experimental data decay curve.
     * \return Pointer to the experimental data decay curve
     */
    DecayCurve* get_data();
    /**
     * Set the experimental data decay curve.
     * \param v Pointer to the experimental data decay curve
     */
    void set_data(DecayCurve* v);

private:
    /**
     * Update the vector of weighted residuals.
     */
    void update_weighted_residuals();

public:
    /**
     * Get the vector of weighted residuals.
     * \return Reference to the vector of weighted residuals
     */
    std::vector<double>& get_weighted_residuals();
    /**
     * Calculate the score based on the specified scoring type.
     * \param start Start index of the range of data points to consider (default: 0)
     * \param stop Stop index of the range of data points to consider (default: -1, i.e., all points)
     * \param score_type Type of score to calculate (default: nullptr, i.e., use the current score type)
     * \return The calculated score
     */
    double get_score(int start = 0, int stop = -1, const char* score_type = nullptr);
    /**
     * Set the score type.
     * \param v The score type to set
     */
    void set_score_type(std::string v);
    /**
     * Get the score type.
     * \return The current score type
     */
    std::string get_score_type();
    /**
     * Set the model and experimental data decay curves, score type, and range of data points.
     * \param model Pointer to the model decay curve
     * \param data Pointer to the experimental data decay curve
     * \param score_type Type of score to calculate (default: "poisson")
     * \param start Start index of the range of data points to consider (default: 0)
     * \param stop Stop index of the range of data points to consider (default: -1, i.e., all points)
     */
    void set(
        DecayCurve* model, DecayCurve* data,
        std::string score_type = "poisson",
        int start = 0, int stop = -1
    );
    /**
     * Evaluate and return the score.
     * \param model Pointer to the model decay curve (default: nullptr, i.e., use the current model)
     * \return The calculated score
     */
    double score(DecayCurve* model = nullptr);

    /**
     * @brief Constructs a DecayScore object with the specified model, experimental data, score type, and range of data points.
     * @param model Pointer to the model decay curve
     * @param data Pointer to the experimental data decay curve
     * @param score_type Type of score to calculate (default: "poisson")
     * @param start Start index of the range of data points to consider (default: 0)
     * @param stop Stop index of the range of data points to consider (default: -1, i.e., all points)
     */
    DecayScore(
        DecayCurve* model = nullptr,
        DecayCurve* data = nullptr,
        std::string score_type = "poisson",
        int start = 0,
        int stop = -1
    );

    /**
     * @brief Destroys the DecayScore object and frees any allocated memory.
     */
    ~DecayScore();

};

IMPBFF_END_NAMESPACE


#endif //IMPBFF_DECAYSCORE_H
