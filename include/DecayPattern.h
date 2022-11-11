/**
 * \file IMP/bff/DecayPattern.h
 * \brief Simple Accessible Volume decorator.
 *
 * \authors Thomas-Otavio Peulen
 * Copyright 2007-2022 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_DECAYPATTERN_H
#define IMPBFF_DECAYPATTERN_H

#include <IMP/bff/bff_config.h>

#include <iostream>
#include <algorithm>

# include <IMP/bff/DecayModifier.h>

IMPBFF_BEGIN_NAMESPACE


class IMPBFFEXPORT DecayPattern: public DecayModifier{

private:

    /// The constant offset in the model
    double constant_offset = 0.0;

    /// Fraction (area) of background pattern
    double pattern_fraction = 0.0;

public:

    /// Fraction (area) of background pattern
    double get_pattern_fraction() const{
        return pattern_fraction;
    }

    /// Fraction (area) of background pattern
    void set_pattern_fraction(double v){
        pattern_fraction = std::max(0.0, std::min(1.0, std::abs(v)));
    }

    void set_pattern(DecayCurve* v){
        set_data(v);
    }

    DecayCurve* get_pattern(){
        return get_data();
    }

    void set_constant_offset(double v){
        constant_offset = v;
    }

    double get_constant_offset() const{
        return constant_offset;
    }

    void add(DecayCurve* out) override{
        if(out != nullptr && is_active()) {
            // resize background pattern to size of input decay
            resize(out->size());

            int start = get_start(out);
            int stop = get_stop(out);

            double f = get_pattern_fraction();
            auto bg = get_pattern();


            // Compute pre-factors of decay and added pattern
            double ds, bs;
            if (f > 0.0) {
                ds = std::accumulate(out->y.begin() + start, out->y.begin() + stop, 0.0);
                bs = std::accumulate(bg->y.begin() + start, bg->y.begin() + stop, 0.0);
            } else {
                ds = bs = 1.0;
            }

            double fd = (1. - f);
            double fb = f * ds / bs;

            // mix constant offset, decay, and pattern
            double offset = get_constant_offset();
            for (int i = start; i < stop; i++) {
                out->y[i] = out->y[i] * fd + bg->y[i] * fb + offset;
            }
        }
    }

    DecayPattern(
            double constant_offset = 0.0,
            DecayCurve*pattern = nullptr,
            double pattern_fraction = 0.0,
            int start = 0,
            int stop = -1,
            bool active = true
    ) : DecayModifier(pattern, start, stop, active) {
#if IMPBFF_VERBOSE
        std::clog << "DecayPattern::DecayPattern" << std::endl;
#endif
        set_pattern_fraction(pattern_fraction);
        set_constant_offset(constant_offset);
    }

};

IMPBFF_END_NAMESPACE


#endif //IMPBFF_DECAYPATTERN_H
