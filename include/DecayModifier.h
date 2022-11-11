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

    /// Input data potentially by modifier
    DecayCurve* data = nullptr;
    DecayCurve* default_data = nullptr;

public:

    virtual void set_data(DecayCurve* v){
        if(v != nullptr){
            data = v;
        }
    }

    virtual DecayCurve* get_data(){
        auto re = default_data;
        if(data != nullptr){
            re = data;
        }
        return re;
    }

    bool is_active() const {
        return _is_active;
    }

    void set_active(bool v){
        _is_active = v;
    }

    void set(DecayCurve* data, int start=0, int stop=-1, bool active = true){
        DecayRange::set(start, stop);
        set_data(data);
        set_active(active);
    }

    void resize(size_t n, double v = 0.0) {
        default_data->resize(n, v);
        if(data != nullptr){
            data->resize(n, v);
        }
    }

    // Modify a DecayCurve
    virtual void add(DecayCurve* out) = 0;

    DecayModifier(DecayCurve *data = nullptr, int start = 0, int stop = -1, bool active = true) :
            DecayRange(start, stop){
        default_data = new DecayCurve();
        set_data(data);
        set_active(active);
    }

    ~DecayModifier() {
        delete default_data;
    }

};

IMPBFF_END_NAMESPACE


#endif //IMPBFF_DECAYMODIFIER_H
