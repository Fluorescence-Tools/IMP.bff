/**
 * \file IMP/bff/DecayLifetimeHandler.h
 * \brief Store and handle lifetime spectra
 *
 * This file defines the DecayLifetimeHandler class, which is responsible for storing and handling lifetime spectra.
 * A lifetime spectrum consists of a set of lifetimes and their corresponding amplitudes.
 * The class provides methods to set and retrieve the lifetime spectrum, as well as to manipulate it.
 *
 * \authors Thomas-Otavio Peulen
 * Copyright 2007-2023 IMP Inventors. All rights reserved.
 */

#ifndef IMPBFF_DECAYLIFETIMEHANDLER_H
#define IMPBFF_DECAYLIFETIMEHANDLER_H

#include <IMP/bff/bff_config.h>
#include <vector>
#include <limits>
#include <iostream>
#include <IMP/bff/DecayRoutines.h>

IMPBFF_BEGIN_NAMESPACE

class IMPBFFEXPORT DecayLifetimeHandler {
private:
    std::vector<double> _lifetime_spectrum = std::vector<double>();  //!< Lifetime spectrum / original
    std::vector<double> lt_ = std::vector<double>();  //!< Lifetime spectrum / for getter
    double amplitude_threshold = std::numeric_limits<double>::epsilon();  //!< Threshold used to discriminate lifetimes with small amplitudes
    bool use_amplitude_threshold = false;  //!< If true, lifetimes with small amplitudes are discriminated
    bool abs_lifetime_spectrum = true;  //!< If true, absolute values of lifetime spectrum are used to compute model function

public:
    /**
     * Get the amplitude threshold used to discriminate lifetimes with small amplitudes.
     * \return The amplitude threshold.
     */
    double get_amplitude_threshold();

    /**
     * Set the amplitude threshold used to discriminate lifetimes with small amplitudes.
     * \param v The amplitude threshold.
     */
    void set_amplitude_threshold(double v);

    /**
     * Check if the amplitude threshold is being used to discriminate lifetimes with small amplitudes.
     * \return True if the amplitude threshold is being used, false otherwise.
     */
    bool get_use_amplitude_threshold();

    /**
     * Set whether to use the amplitude threshold to discriminate lifetimes with small amplitudes.
     * \param v True to use the amplitude threshold, false otherwise.
     */
    void set_use_amplitude_threshold(bool v);

    /**
     * Check if the absolute values of the lifetime spectrum are being used to compute the model function.
     * \return True if the absolute values are being used, false otherwise.
     */
    bool get_abs_lifetime_spectrum() const;

    /**
     * Set whether to use the absolute values of the lifetime spectrum to compute the model function.
     * \param v True to use the absolute values, false otherwise.
     */
    void set_abs_lifetime_spectrum(bool v);

    /**
     * Set the lifetime spectrum.
     * \param v The lifetime spectrum to set.
     */
    void set_lifetime_spectrum(std::vector<double> v);

    /**
     * Add a lifetime to the lifetime spectrum.
     * \param amplitude The amplitude of the lifetime.
     * \param lifetime The lifetime value.
     */
    void add_lifetime(double amplitude, double lifetime);

    /**
     * Get a reference to the lifetime spectrum.
     * \return A reference to the lifetime spectrum.
     */
    std::vector<double>& get_lifetime_spectrum();

    /**
     * Get a view of the lifetime spectrum.
     * \param output_view A pointer to the output view of the lifetime spectrum.
     * \param n_output A pointer to the number of elements in the output view.
     */
    void get_lifetime_spectrum(double **output_view, int *n_output);

    /**
     * Construct a DecayLifetimeHandler object.
     * \param lifetime_spectrum The initial lifetime spectrum (default: empty).
     * \param use_amplitude_threshold Whether to use the amplitude threshold (default: false).
     * \param abs_lifetime_spectrum Whether to use the absolute values of the lifetime spectrum (default: false).
     * \param amplitude_threshold The amplitude threshold (default: epsilon).
     */
    DecayLifetimeHandler(
        std::vector<double> lifetime_spectrum = std::vector<double>(),
        bool use_amplitude_threshold = false,
        bool abs_lifetime_spectrum = false,
        double amplitude_threshold = std::numeric_limits<double>::epsilon()
    );
};

IMPBFF_END_NAMESPACE

#endif // IMPBFF_DECAYLIFETIMEHANDLER_H