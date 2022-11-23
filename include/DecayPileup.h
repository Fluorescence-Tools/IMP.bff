/**
 * \file IMP/bff/DecayPileup.h
 * \brief Simple Accessible Volume decorator.
 *
 * \authors Thomas-Otavio Peulen
 * Copyright 2007-2022 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_DECAYPILEUP_H
#define IMPBFF_DECAYPILEUP_H

#include <IMP/bff/bff_config.h>

#include <memory> /* shared_ptr */
#include <string>
#include <iostream> /* std::cerr */

#include <IMP/bff/DecayCurve.h>
#include <IMP/bff/DecayModifier.h>
#include <IMP/bff/DecayRoutines.h> /* add_pile_up_to_model */


IMPBFF_BEGIN_NAMESPACE

class IMPBFFEXPORT DecayPileup : public DecayModifier {

private:

    /// Dead time of the instrument in units of the lifetime (usually nanoseconds)
    double instrument_dead_time = std::numeric_limits<double>::epsilon();

    /// Repetition rate of the light source in MHz
    double repetition_rate = 100.0;

    /// Identifier of pile up model
    std::string pile_up_model = "coates";

public:

    void set_pile_up_model(std::string v);

    std::string get_pile_up_model();

    void set_repetition_rate(double v);

    double get_repetition_rate();

    void set_instrument_dead_time(double v);

    double get_instrument_dead_time();

    void add(DecayCurve* out) override;

    DecayPileup(
        DecayCurve* data,
        const char* pile_up_model = "coates",
        double repetition_rate = 100,
        double instrument_dead_time = 120,
        int start = 0, int stop = -1,
        bool active = true
    );

};

IMPBFF_END_NAMESPACE


#endif //IMPBFF_DECAYPILEUP_H
