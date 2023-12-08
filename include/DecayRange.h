/**
 * \file IMP/bff/DecayRange.h
 * \brief Defines inspected range of fluorescence decay
 *
 * This file defines the DecayRange class, which represents an inspected range
 * of fluorescence decay. It allows setting and retrieving the start and stop
 * indices of the range, as well as setting the range using a vector of indices.
 *
 * \authors Thomas-Otavio Peulen
 * \copyright 2007-2023 IMP Inventors. All rights reserved.
 */

#ifndef IMPBFF_DECAYRANGE_H
#define IMPBFF_DECAYRANGE_H

#include <IMP/bff/bff_config.h>
#include <IMP/bff/DecayCurve.h>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

/**
 * \class DecayRange
 * \brief Represents an inspected range of fluorescence decay
 *
 * The DecayRange class represents an inspected range of fluorescence decay.
 * It allows setting and retrieving the start and stop indices of the range,
 * as well as setting the range using a vector of indices.
 */
class IMPBFFEXPORT DecayRange {
private:
    size_t _start = 0; ///< Start index of the decay range
    size_t _stop = -1; ///< Stop index of the decay range

public:
    /**
     * \brief Set the start index of the decay range
     * \param v The start index
     */
    void set_start(int v);

    /**
     * \brief Get the start index of the decay range
     * \param d The DecayCurve object (optional)
     * \return The start index
     */
    size_t get_start(DecayCurve* d = nullptr) const;

    /**
     * \brief Set the stop index of the decay range
     * \param v The stop index
     */
    void set_stop(int v);

    /**
     * \brief Get the stop index of the decay range
     * \param d The DecayCurve object (optional)
     * \return The stop index
     */
    size_t get_stop(DecayCurve* d = nullptr) const;

    /**
     * \brief Set the range using a vector of indices
     * \param v The vector of indices
     */
    void set_range(std::vector<int> v);

    /**
     * \brief Get the range as a vector of indices
     * \param d The DecayCurve object (optional)
     * \return The vector of indices
     */
    std::vector<int> get_range(DecayCurve* d = nullptr);

    /**
     * \brief Set the start and stop indices of the decay range
     * \param start The start index (default: 0)
     * \param stop The stop index (default: -1)
     */
    void set(int start = 0, int stop = -1);

    /**
     * \brief Constructor
     * \param start The start index of the decay range
     * \param stop The stop index of the decay range
     */
    DecayRange(int start, int stop);

    /**
     * \brief Destructor
     */
    virtual ~DecayRange() = default;
};

IMPBFF_END_NAMESPACE

#endif // IMPBFF_DECAYRANGE_H