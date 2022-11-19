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

private:

    DecayCurve linearization_table_;

public:

    void set_linearization_table(DecayCurve* v) {
        set_data(v);
    }

    DecayCurve* get_linearization_table() {
        return get_data();
    }

    /// Modify DecayCurve object and apply a linearization
    void add(DecayCurve *out) override {
        if ((out != nullptr) && is_active()) {
            DecayCurve *lt = get_linearization_table();
            lt->resize(out->size(), 1.0);
            size_t start, stop;
            start = std::min(get_start(out), out->size());
            stop = std::min(get_stop(out), out->size());
            for (size_t i = start; i < stop; i++) {
                out->y[i] *= lt->y[i];
            }
        }
    }

    /*!
     * DecayLinearization is a DecayModifier that takes a DecayCurve as an input,
     * computes for the input an (optionally) smoothed linearization_table.
     *
     * Applying DecayLinearization to a DecayCurve using DecayLinearization::add
     * applies the linearization_table to the DecayCurve. The computation of a
     * perturbed fluorescence decay model is a typical application of DecayLinearization.
     */
    DecayLinearization(
            DecayCurve* linearization_table,
            int start, int stop,
            bool active = true,
            int n_window = 5
    ) : DecayModifier(nullptr, start, stop, active){
#if IMPBFF_VERBOSE
        std::clog << "DecayLinearization::DecayLinearization(double*, int, bool)" << std::endl;
#endif
        linearization_table_.set_x(linearization_table->get_x());
        linearization_table_.set_y(linearization_table->get_y());
        linearization_table_.apply_simple_moving_average(
                get_start(linearization_table),
                get_stop(linearization_table),
                n_window,
                true
        );
        set_linearization_table(&linearization_table_);
        set_start(start);
        set_stop(stop);
}


};

IMPBFF_END_NAMESPACE


#endif //IMPBFF_DECAYLINEARIZATION_H
