/**
 * \file IMP/bff/DecayScale.h
 * \brief Simple Accessible Volume decorator.
 *
 * \authors Thomas-Otavio Peulen
 * Copyright 2007-2022 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_DECAYSCALE_H
#define IMPBFF_DECAYSCALE_H

#include <IMP/bff/bff_config.h>

#include <iostream>

#include <IMP/bff/DecayCurve.h>
#include <IMP/bff/DecayModifier.h>
#include <IMP/bff/DecayRoutines.h>  /* rescale_w_bg */



IMPBFF_BEGIN_NAMESPACE

class IMPBFFEXPORT DecayScale : public DecayModifier{

private:

    /// A constant that is subtracted from the data
    double _constant_background = 0.0;
    bool _blank_outside = true;

public:

    /// Number of photons in data between start and stop (if model is scaled to data). Otherwise
    /// user-specified number of photons
    double get_number_of_photons(){
        double re = 0.0;
        DecayCurve* d = get_data();
        for(size_t i = get_start(d); i < get_stop(d); i++)
            re += d->y[i];
        return re;
    }

    double get_constant_background() const {
        return _constant_background;
    }

    void set_constant_background(double v){
        _constant_background = v;
    }

    void set_blank_outside(double v){
        _blank_outside = v;
    }

    bool get_blank_outside(){
        return _blank_outside;
    }

    void set(
            DecayCurve* data = nullptr,
            double constant_background = 0.0,
            int start = 0, int stop = -1, bool active = true,
            bool blank_outside=true
    ){
        set_data(data);
        set_constant_background(constant_background);
        set_start(start);
        set_stop(stop);
        set_active(active);
        set_blank_outside(blank_outside);
    }

    DecayScale(
            DecayCurve* data = nullptr,
            double constant_background = 0.0,
            int start=0, int stop=-1, bool active=true,
            bool blank_outside=true
    ) : DecayModifier(data, start, stop, active){
        set_constant_background(constant_background);
    }

    void add(DecayCurve* decay){
#if IMPBFF_VERBOSE
        std::clog << "DecayScale::add" << std::endl;
#endif
        IMP_USAGE_CHECK(data != nullptr, "Experimental data not set - cannot scale model.");
        if(is_active()){
            auto d = get_data();
            double scale = 0.0;
            double* model = decay->y.data();
            double* data = d->y.data();
            double* squared_data_weights = d->ey.data();
            int n_model = decay->y.size();
            int start = (int) get_start(decay);
            int stop = (int) get_stop(decay);
            auto constant_background = get_constant_background();
            decay_rescale_w_bg(
                    model, data,
                    squared_data_weights, constant_background,
                    &scale, start, stop);
            if(_blank_outside){
                for(int i = 0; i < start; i++){
                    model[i] = 0.0;
                }
                for(int i = stop; i < n_model; i++){
                    model[i] = 0.0;
                }
            }
        }
    }

};

IMPBFF_END_NAMESPACE

#endif //IMPBFF_DECAYSCALE_H
