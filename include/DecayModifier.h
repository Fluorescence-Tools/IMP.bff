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

#include <cereal/access.hpp>
#include <cereal/types/base_class.hpp>

IMPBFF_BEGIN_NAMESPACE

IMPBFF_DEPRECATED_HEADER(
    2.25, "The fluorescence-decay classes have moved to tttrlib and will be removed from IMP.bff in the next release; IMP.bff keeps only the structure-related features (AV, PathMap, AVNetworkRestraint). Use tttrlib instead.")

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

    friend class cereal::access;

    /* One serialize(), not a save/load pair: DecayRange already supplies a
       serialize() that this class inherits, and cereal rejects a type that
       offers two candidate functions.

       The two curve pointers differ in ownership. `default_data` is owned (the
       destructor deletes it) and is archived by value. `data` points at a curve
       owned by the caller, so archiving it would silently take a copy and hand
       this object ownership of it; it is cleared on load and must be
       re-attached with set_data(). */
    template<class Archive> void serialize(Archive &ar) {
        ar(cereal::base_class<DecayRange>(this), _is_active);
        bool has_default = (default_data != nullptr);
        ar(has_default);
        if (std::is_base_of<cereal::detail::InputArchiveBase, Archive>::value) {
            delete default_data;
            default_data = has_default ? new DecayCurve() : nullptr;
            data = nullptr;
        }
        if (has_default) ar(*default_data);
    }

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