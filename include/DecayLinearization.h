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

    void set_linearization_table(DecayCurve* v);

    DecayCurve* get_linearization_table();

    /// Modify DecayCurve object and apply a linearization
    void add(DecayCurve *out) override;

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
    );


};

IMPBFF_END_NAMESPACE


#endif //IMPBFF_DECAYLINEARIZATION_H
