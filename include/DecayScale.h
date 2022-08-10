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

/// Scales the model in the specified range
//            /*!
//             *
//             *
//             * The model is either scaled to the data or to a specified number of photons
//             *
//             * @param scale_model_to_data[in] If true the model is scaled to the data
//             * @param number_of_photons[in] If scale_model_to_data is false the model is
//             * scaled to the specified number of photons
//             * @param start[in] Specifies the start index of the model range
//             * @param stop[in] Specifies the stop index of the model range
//             * @param constant_background[in] number of constant background photons in
//             * each histogram bin.
//             * @param model[in, out] Model function that is scaled. The model function is
//             * changed inplace.
//             * @param data[in] Data array to which the model function can be scaled
//             * @param squared_data_weights[in] squared weights of the data
//             */
class DecayScale : public DecayModifier{

private:

    /// A constant that is subtracted from the data
    double _constant_background = 0.0;

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

    void set(
            DecayCurve* data = nullptr,
            double constant_background = 0.0,
            int start = 0, int stop = -1, bool active = true
    ){
        set_data(data);
        set_constant_background(constant_background);
        set_start(start);
        set_stop(stop);
        set_active(active);
    }

    DecayScale(
            DecayCurve* data = nullptr,
            double constant_background = 0.0,
            int start=0, int stop=-1, bool active=true
    ) : DecayModifier(data, start, stop, active){
        set_constant_background(constant_background);
    }

    void add(DecayCurve* decay){
#if IMPBFF_VERBOSE
        std::clog << "DecayScale::add" << std::endl;
#endif
        if(is_active()){
            auto d = get_data();
            double scale = 0.0;
            double* model = decay->y.data();
            double* data = d->y.data();
            double* squared_data_weights = d->ey.data();
            int start = (int) get_start(decay);
            int stop = (int) get_stop(decay);
            auto constant_background = get_constant_background();
            rescale_w_bg(
                    model, data,
                    squared_data_weights, constant_background,
                    &scale, start, stop);

        }
    }

};

IMPBFF_END_NAMESPACE

#endif //IMPBFF_DECAYSCALE_H
