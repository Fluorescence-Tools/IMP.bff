/**
 * \file IMP/bff/DecayModifier.h
 * \brief Simple Accessible Volume decorator.
 *
 * \authors Thomas-Otavio Peulen
 * Copyright 2007-2022 IMP Inventors. All rights reserved.
 */

#ifndef IMPBFF_DECAYMODIFIER_H
#define IMPBFF_DECAYMODIFIER_H

#include <IMP/bff/bff_config.h>
#include <memory>
#include <limits>
#include <IMP/bff/DecayRange.h>
#include <IMP/bff/DecayCurve.h>

IMPBFF_BEGIN_NAMESPACE

/**
 * \class DecayModifier
 * \brief A decorator that modifies a DecayCurve within a specified range.
 *
 * The DecayModifier class is a decorator that modifies a DecayCurve object
 * within a specified range. It is used to apply modifications to the decay
 * behavior of a DecayCurve. The modifications are applied by the `add` method,
 * which modifies the input DecayCurve object.
 */
class IMPBFFEXPORT DecayModifier : public DecayRange {

private:
    bool _is_active = true;

protected:
    DecayCurve* data = nullptr;
    DecayCurve* default_data = nullptr;

public:
    /**
     * Set the DecayCurve object to be modified.
     * \param v The DecayCurve object to be modified.
     */
    virtual void set_data(DecayCurve* v);

    /**
     * Get the DecayCurve object being modified.
     * \return The DecayCurve object being modified.
     */
    virtual DecayCurve* get_data();

    /**
     * Check if the DecayModifier is active.
     * \return True if the DecayModifier is active, false otherwise.
     */
    bool is_active() const;

    /**
     * Set the activity status of the DecayModifier.
     * \param v The activity status of the DecayModifier.
     */
    void set_active(bool v);

    /**
     * Set the values of the DecayModifier.
     * \param data The DecayModifier data.
     * \param start The start of the DecayModifier.
     * \param stop The stop of the DecayModifier.
     * \param active If true, the DecayModifier modifies the input decay.
     */
    void set(DecayCurve* data, int start=0, int stop=-1, bool active = true);

    /**
     * Resize the data of the DecayModifier.
     * \param n The new size of the data.
     * \param v The value of the data (if larger than the original size).
     */
    void resize(size_t n, double v = 0.0);

    /**
     * Modify the DecayCurve object.
     * \param out The DecayCurve object to be modified.
     */
    virtual void add(DecayCurve* out) = 0;

    /**
     * Construct a DecayModifier object.
     * \param data The DecayCurve object to be modified.
     * \param start The start of the DecayModifier.
     * \param stop The stop of the DecayModifier.
     * \param active If true, the DecayModifier is active.
     */
    DecayModifier(DecayCurve *data = nullptr,
                  int start = 0, int stop = -1, bool active = true);

    /**
     * Destructor.
     */
    ~DecayModifier() {
        delete default_data;
    }
};

IMPBFF_END_NAMESPACE

#endif // IMPBFF_DECAYMODIFIER_H