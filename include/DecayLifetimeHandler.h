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

IMPBFF_BEGIN_NAMESPACE

class DecayLifetimeHandler{

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

    double get_amplitude_threshold(){
        return std::abs(amplitude_threshold);
    }

    void set_amplitude_threshold(double v){
        amplitude_threshold = v;
    }

    bool get_use_amplitude_threshold(){
        return use_amplitude_threshold;
    }

    void set_use_amplitude_threshold(bool v){
        use_amplitude_threshold = v;
    }

    bool get_abs_lifetime_spectrum() const{
        return abs_lifetime_spectrum;
    }

    void set_abs_lifetime_spectrum(bool v){
        abs_lifetime_spectrum = v;
    }

    void set_lifetime_spectrum(std::vector<double> v) {
        _lifetime_spectrum = v;
    }

    void add_lifetime(double amplitude, double lifetime) {
        _lifetime_spectrum.emplace_back(amplitude);
        _lifetime_spectrum.emplace_back(lifetime);
    }

    std::vector<double>& get_lifetime_spectrum(){
        lt_ = _lifetime_spectrum;
        if(use_amplitude_threshold){
            discriminate_small_amplitudes(lt_.data(), lt_.size(), amplitude_threshold);
        }
        if (abs_lifetime_spectrum) {
            for(auto &v: lt_) v = std::abs(v);
        }
        return lt_;
    }

    void get_lifetime_spectrum(double **output_view, int *n_output) {
        auto lt = get_lifetime_spectrum();
        *output_view = lt.data();
        *n_output = (int) lt.size();
    }

    DecayLifetimeHandler(
            std::vector<double> lifetime_spectrum = std::vector<double>(),
            bool use_amplitude_threshold = false,
            bool abs_lifetime_spectrum = false,
            double amplitude_threshold = std::numeric_limits<double>::epsilon()
            ){
        set_use_amplitude_threshold(use_amplitude_threshold);
        set_abs_lifetime_spectrum(abs_lifetime_spectrum);
        set_amplitude_threshold(amplitude_threshold);
        set_lifetime_spectrum(lifetime_spectrum);
    }

};

IMPBFF_END_NAMESPACE


#endif //IMPBFF_DECAYLIFETIMEHANDLER_H
