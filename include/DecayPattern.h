/**
 * \file IMP/bff/DecayPattern.h
 * \brief Simple Accessible Volume decorator.
 *
 * \authors Thomas-Otavio Peulen
 * Copyright 2007-2023 IMP Inventors. All rights reserved.
 */

#ifndef IMPBFF_DECAYPATTERN_H
#define IMPBFF_DECAYPATTERN_H

#include <IMP/bff/bff_config.h>
#include <iostream>
#include <algorithm>
#include <IMP/bff/DecayModifier.h>

IMPBFF_BEGIN_NAMESPACE

/**
 * @class DecayPattern
 * @brief The DecayPattern class represents a decay pattern with a constant offset and a background pattern.
 *
 * The DecayPattern class is a subclass of DecayModifier and provides functionality to add a background pattern to a decay curve.
 * It allows setting and retrieving the constant offset and the fraction of the background pattern.
 */
class IMPBFFEXPORT DecayPattern: public DecayModifier{
private:
    /// The constant offset in the model
    double constant_offset = 0.0;

    /// Fraction (area) of background pattern
    double pattern_fraction = 0.0;

public:
    /**
     * @brief Get the fraction (area) of the background pattern.
     * @return The fraction (area) of the background pattern.
     */
    double get_pattern_fraction() const;

    /**
     * @brief Set the fraction (area) of the background pattern.
     * @param v The fraction (area) of the background pattern to set.
     */
    void set_pattern_fraction(double v);

    /**
     * @brief Set the background pattern.
     * @param v The background pattern to set.
     */
    void set_pattern(DecayCurve* v);

    /**
     * @brief Get the background pattern.
     * @return The background pattern.
     */
    DecayCurve* get_pattern();

    /**
     * @brief Set the constant offset in the model.
     * @param v The constant offset to set.
     */
    void set_constant_offset(double v);

    /**
     * @brief Get the constant offset in the model.
     * @return The constant offset in the model.
     */
    double get_constant_offset() const;

    /**
     * @brief Add the background pattern to the given decay curve.
     * @param out The decay curve to add the background pattern to.
     */
    void add(DecayCurve* out) override;

    /**
     * @brief Constructor for the DecayPattern class.
     * @param constant_offset The constant offset in the model (default: 0.0).
     * @param pattern The background pattern (default: nullptr).
     * @param pattern_fraction The fraction (area) of the background pattern (default: 0.0).
     * @param start The start index of the decay curve (default: 0).
     * @param stop The stop index of the decay curve (default: -1).
     * @param active The flag indicating if the DecayPattern is active (default: true).
     */
    explicit DecayPattern(
            double constant_offset = 0.0,
            DecayCurve* pattern = nullptr,
            double pattern_fraction = 0.0,
            int start = 0,
            int stop = -1,
            bool active = true
    );
};

IMPBFF_END_NAMESPACE

#endif //IMPBFF_DECAYPATTERN_H