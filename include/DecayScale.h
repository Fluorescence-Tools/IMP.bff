/**
 * \file IMP/bff/DecayScale.h
 * \brief Simple Accessible Volume decorator.
 *
 * This file contains the declaration of the DecayScale class, which is a subclass
 * of DecayModifier. DecayScale is used to scale a DecayCurve by a constant factor
 * and subtract a constant background value.
 *
 * \authors Thomas-Otavio Peulen
 * Copyright 2007-2023 IMP Inventors. All rights reserved.
 */

#ifndef IMPBFF_DECAYSCALE_H
#define IMPBFF_DECAYSCALE_H

#include <IMP/bff/bff_config.h>
#include <iostream>
#include <IMP/bff/DecayCurve.h>
#include <IMP/bff/DecayModifier.h>
#include <IMP/bff/DecayRoutines.h>

IMPBFF_BEGIN_NAMESPACE

/**
 * \class DecayScale
 * \brief A class for scaling a DecayCurve by a constant factor and subtracting a constant background value.
 *
 * The DecayScale class is a subclass of DecayModifier and is used to scale a DecayCurve
 * by a constant factor and subtract a constant background value. It provides methods to
 * set and get the constant background value, as well as to set whether the curve should be
 * blanked outside a specified range. The scaling and background subtraction are applied
 * to the DecayCurve when the modify() method is called.
 */
class IMPBFFEXPORT DecayScale : public DecayModifier {
private:
    double _constant_background = 0.0;  //!< A constant that is subtracted from the data
    bool _blank_outside = true;  //!< Flag indicating whether the curve should be blanked outside a specified range

public:
    /**
     * Get the number of photons in the data between the start and stop indices.
     * If the model is scaled to the data, the number of photons is returned.
     * Otherwise, the user-specified number of photons is returned.
     * \return The number of photons in the data.
     */
    double get_number_of_photons();

    /**
     * Get the constant background value.
     * \return The constant background value.
     */
    double get_constant_background() const;

    /**
     * Set the constant background value.
     * \param v The constant background value to be set.
     */
    void set_constant_background(double v);

    /**
     * Set whether the curve should be blanked outside a specified range.
     * \param v Flag indicating whether the curve should be blanked outside a specified range.
     */
    void set_blank_outside(double v);

    /**
     * Get whether the curve should be blanked outside a specified range.
     * \return Flag indicating whether the curve should be blanked outside a specified range.
     */
    bool get_blank_outside();

    /**
     * Set the parameters of the DecayScale object.
     * \param data The DecayCurve to be scaled.
     * \param constant_background The constant background value to be subtracted.
     * \param start The start index of the range to be considered.
     * \param stop The stop index of the range to be considered.
     * \param active Flag indicating whether the DecayScale object is active.
     * \param blank_outside Flag indicating whether the curve should be blanked outside the specified range.
     */
    void set(DecayCurve* data = nullptr, double constant_background = 0.0,
             int start = 0, int stop = -1, bool active = true,
             bool blank_outside = true);

    /**
     * Construct a DecayScale object with the specified parameters.
     * \param data The DecayCurve to be scaled.
     * \param constant_background The constant background value to be subtracted.
     * \param start The start index of the range to be considered.
     * \param stop The stop index of the range to be considered.
     * \param active Flag indicating whether the DecayScale object is active.
     * \param blank_outside Flag indicating whether the curve should be blanked outside the specified range.
     */
    DecayScale(DecayCurve* data = nullptr, double constant_background = 0.0,
               int start = 0, int stop = -1, bool active = true,
               bool blank_outside = true);

    /**
     * Add a DecayCurve to be scaled by the DecayScale object.
     * \param decay The DecayCurve to be added.
     */
    void add(DecayCurve* decay);
};

IMPBFF_END_NAMESPACE

#endif // IMPBFF_DECAYSCALE_H