/**
 * \file IMP/bff/DecayLinearization.h
 * \brief Decay modifier add non-linearities
 *
 * \authors Thomas-Otavio Peulen
 * Copyright 2007-2022 IMP Inventors. All rights reserved.
 *
 */

#ifndef IMPBFF_DECAYLINEARIZATION_H
#define IMPBFF_DECAYLINEARIZATION_H

#include <IMP/bff/bff_config.h>

#include <iostream> /* std::cerr */
#include <algorithm> /* std::min */
#include <vector>

#include <IMP/bff/DecayCurve.h>
#include <IMP/bff/DecayModifier.h>

IMPBFF_BEGIN_NAMESPACE


class IMPBFFEXPORT DecayLinearization: public DecayModifier{


public:

    void set_linearization_table(DecayCurve* v) {
        set_data(v);
    }

    DecayCurve* get_linearization_table() {
        return get_data();
    }

    void add(DecayCurve *decay) {
        if ((decay != nullptr) && is_active()) {
            DecayCurve *lt = get_linearization_table();
            lt->resize(decay->size(), 1.0);
            size_t start, stop;
            start = std::min(get_start(), decay->size());
            stop = std::min(get_stop(), decay->size());
            for (size_t i = start; i < stop; i++) {
                decay->y[i] *= lt->y[i];
            }
        }
    }

    DecayLinearization(
            DecayCurve* linearization_table,
            int start,
            int stop,
            bool active = false
    ) : DecayModifier(linearization_table, start, stop, active){
#if IMPBFF_VERBOSE
        std::clog << "DecayLinearization::DecayLinearization(double*, int, bool)" << std::endl;
#endif
    }


};

IMPBFF_END_NAMESPACE


#endif //IMPBFF_DECAYLINEARIZATION_H
