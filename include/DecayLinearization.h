/**
 * \file IMP/bff/DecayLinearization.h
 * \brief Decay modifier to add non-linearities
 *
 * This file contains the declaration of the DecayLinearization class, which is a
 * subclass of DecayModifier. DecayLinearization applies a linearization to a
 * DecayCurve object, which represents a decay model. It computes a linearization
 * table and applies it to the input DecayCurve.
 *
 * \authors Thomas-Otavio Peulen
 * Copyright 2007-2022 IMP Inventors. All rights reserved.
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

/**
 * \class DecayLinearization
 * \brief A decay modifier to apply linearization to a DecayCurve
 *
 * The DecayLinearization class is a subclass of DecayModifier. It takes a
 * DecayCurve as input and computes a linearization table for the input curve.
 * The linearization table is then applied to the input DecayCurve using the
 * add() method.
 *
 * Typical applications of DecayLinearization include computing a perturbed
 * fluorescence decay model.
 */
class IMPBFFEXPORT DecayLinearization : public DecayModifier {
private:
    DecayCurve linearization_table_; //!< The linearization table

public:
    /**
     * Set the linearization table for the DecayLinearization object.
     * \param v The linearization table
     */
    void set_linearization_table(DecayCurve* v);

    /**
     * Get the linearization table of the DecayLinearization object.
     * \return The linearization table
     */
    DecayCurve* get_linearization_table();

    /**
     * Apply the linearization table to the input DecayCurve.
     * \param out The output DecayCurve
     */
    void add(DecayCurve* out) override;

    /**
     * Construct a DecayLinearization object.
     * \param linearization_table The linearization table
     * \param start The start index of the linearization window
     * \param stop The stop index of the linearization window
     * \param active Flag indicating if the DecayLinearization is active
     * \param n_window The size of the linearization window
     */
    DecayLinearization(DecayCurve* linearization_table,
                       int start, int stop,
                       bool active = true,
                       int n_window = 5);
};

IMPBFF_END_NAMESPACE

#endif // IMPBFF_DECAYLINEARIZATION_H