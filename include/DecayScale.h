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
    double get_number_of_photons();

    double get_constant_background() const;

    void set_constant_background(double v);

    void set_blank_outside(double v);

    bool get_blank_outside();

    void set(
            DecayCurve* data = nullptr,
            double constant_background = 0.0,
            int start = 0, int stop = -1, bool active = true,
            bool blank_outside=true
    );

    DecayScale(
            DecayCurve* data = nullptr,
            double constant_background = 0.0,
            int start=0, int stop=-1, bool active=true,
            bool blank_outside=true
    );

    void add(DecayCurve* decay);

};

IMPBFF_END_NAMESPACE

#endif //IMPBFF_DECAYSCALE_H
