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
    double get_pattern_fraction() const;

    /// Fraction (area) of background pattern
    void set_pattern_fraction(double v);

    void set_pattern(DecayCurve* v);

    DecayCurve* get_pattern();

    void set_constant_offset(double v);

    double get_constant_offset() const;

    void add(DecayCurve* out) override;

    explicit DecayPattern(
            double constant_offset = 0.0,
            DecayCurve*pattern = nullptr,
            double pattern_fraction = 0.0,
            int start = 0,
            int stop = -1,
            bool active = true
    );

};

IMPBFF_END_NAMESPACE


#endif //IMPBFF_DECAYPATTERN_H
