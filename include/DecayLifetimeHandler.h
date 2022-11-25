/**
 * \file IMP/bff/DecayLifetimeHandler.h
 * \brief Store and handle lifetime spectra
 *
 * \authors Thomas-Otavio Peulen
 * Copyright 2007-2022 IMP Inventors. All rights reserved.
 *
 */

#ifndef IMPBFF_DECAYLIFETIMEHANDLER_H
#define IMPBFF_DECAYLIFETIMEHANDLER_H

#include <IMP/bff/bff_config.h>

#include <vector>
#include <limits>
#include <iostream>

#include <IMP/bff/DecayRoutines.h>

IMPBFF_BEGIN_NAMESPACE

class IMPBFFEXPORT DecayLifetimeHandler{

private:

    /// Lifetime spectrum / original
    std::vector<double> _lifetime_spectrum = std::vector<double>();

    /// Lifetime spectrum / for getter
    std::vector<double> lt_ = std::vector<double>();

    /// Threshold used to discriminate lifetimes with small amplitudes
    double amplitude_threshold = std::numeric_limits<double>::epsilon();

    /// If true lifetimes with small amplitudes are discriminated
    bool use_amplitude_threshold = false;

    /// If true absolute values of lifetime spectrum is used to compute model function
    bool abs_lifetime_spectrum = true;

public:

    double get_amplitude_threshold();

    void set_amplitude_threshold(double v);

    bool get_use_amplitude_threshold();

    void set_use_amplitude_threshold(bool v);

    bool get_abs_lifetime_spectrum() const;

    void set_abs_lifetime_spectrum(bool v);

    void set_lifetime_spectrum(std::vector<double> v);

    void add_lifetime(double amplitude, double lifetime);

    std::vector<double>& get_lifetime_spectrum();

    void get_lifetime_spectrum(double **output_view, int *n_output);

    DecayLifetimeHandler(
            std::vector<double> lifetime_spectrum = std::vector<double>(),
            bool use_amplitude_threshold = false,
            bool abs_lifetime_spectrum = false,
            double amplitude_threshold = std::numeric_limits<double>::epsilon()
            );

};

IMPBFF_END_NAMESPACE


#endif //IMPBFF_DECAYLIFETIMEHANDLER_H
