/**
 * \file IMP/bff/DecayPileup.h
 * \brief Simple Accessible Volume decorator.
 *
 * This file contains the declaration of the DecayPileup class, which is a subclass of DecayModifier.
 * DecayPileup is a decorator that adds pile-up effects to a DecayCurve object.
 *
 * \authors Thomas-Otavio Peulen
 * Copyright 2007-2022 IMP Inventors. All rights reserved.
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

/**
 * \class DecayPileup
 * \brief A decorator that adds pile-up effects to a DecayCurve object.
 *
 * DecayPileup is a subclass of DecayModifier and provides methods to set and get the pile-up model,
 * repetition rate, and instrument dead time. It also overrides the add() method to add pile-up effects
 * to the DecayCurve object.
 */
class IMPBFFEXPORT DecayPileup : public DecayModifier {

private:
    /// Dead time of the instrument in units of the lifetime (usually nanoseconds)
    double instrument_dead_time = std::numeric_limits<double>::epsilon();

    /// Repetition rate of the light source in MHz
    double repetition_rate = 100.0;

    /// Identifier of pile up model
    std::string pile_up_model = "coates";

public:
    /**
     * \brief Set the pile-up model.
     * \param v The pile-up model identifier.
     */
    void set_pile_up_model(std::string v);

    /**
     * \brief Get the pile-up model.
     * \return The pile-up model identifier.
     */
    std::string get_pile_up_model();

    /**
     * \brief Set the repetition rate.
     * \param v The repetition rate in MHz.
     */
    void set_repetition_rate(double v);

    /**
     * \brief Get the repetition rate.
     * \return The repetition rate in MHz.
     */
    double get_repetition_rate();

    /**
     * \brief Set the instrument dead time.
     * \param v The instrument dead time in units of the lifetime.
     */
    void set_instrument_dead_time(double v);

    /**
     * \brief Get the instrument dead time.
     * \return The instrument dead time in units of the lifetime.
     */
    double get_instrument_dead_time();

    /**
     * \brief Add pile-up effects to the DecayCurve object.
     * \param out The DecayCurve object to add pile-up effects to.
     */
    void add(DecayCurve* out) override;

    /**
     * \brief Constructor for DecayPileup.
     * \param data The DecayCurve object to decorate.
     * \param pile_up_model The pile-up model identifier.
     * \param repetition_rate The repetition rate in MHz.
     * \param instrument_dead_time The instrument dead time in units of the lifetime.
     * \param start The start index of the DecayCurve.
     * \param stop The stop index of the DecayCurve.
     * \param active The activity status of the DecayPileup object.
     */
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