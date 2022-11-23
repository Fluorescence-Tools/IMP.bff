/**
 * \file IMP/bff/DecayModifier.h
 * \brief Simple Accessible Volume decorator.
 *
 * \authors Thomas-Otavio Peulen
 * Copyright 2007-2022 IMP Inventors. All rights reserved.
 *
 */

#ifndef IMPBFF_DECAYMODIFIER_H
#define IMPBFF_DECAYMODIFIER_H

#include <IMP/bff/bff_config.h>

#include <memory>
#include <limits>

#include <IMP/bff/DecayRange.h>
#include <IMP/bff/DecayCurve.h>

IMPBFF_BEGIN_NAMESPACE


class IMPBFFEXPORT DecayModifier : public DecayRange {


private:

    bool _is_active = true;

protected:

    DecayCurve* data = nullptr;
    DecayCurve* default_data = nullptr;

public:

    virtual void set_data(DecayCurve* v);

    virtual DecayCurve* get_data();

    bool is_active() const;

    /// Setter to define if DecayModifier is active. An active
    /// DecayModifier modifies a DecayCurve using
    void set_active(bool v);

    /*!
     * Set values of DecayModifier
     * @param data DecayModifier data
     * @param start start of DecayModifier
     * @param stop stop of DecayModifier
     * @param active if active is true DecayModifier::add modifies input decay
     */
    void set(DecayCurve* data, int start=0, int stop=-1, bool active = true);

    /*!
     * Resize data of DecayModifier
     * @param n new size of data
     * @param v value of data (if larger than original size)
     */
    void resize(size_t n, double v = 0.0);

    /*!
     * Modify DecayCurve object
     * @param out DecayCurve that is modified
     */
    virtual void add(DecayCurve* out) = 0;

    /*!
     * DecayModifier modifies DecayCurve in range [start, stop]
     */
    DecayModifier(DecayCurve *data = nullptr,
                  int start = 0, int stop = -1, bool active = true);

    ~DecayModifier() {
        delete default_data;
    }

};

IMPBFF_END_NAMESPACE


#endif //IMPBFF_DECAYMODIFIER_H
