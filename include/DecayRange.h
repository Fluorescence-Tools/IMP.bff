/**
 * \file IMP/bff/DecayRange.h
 * \brief Defines inspected range of fluorescence decay
 *
 * \authors Thomas-Otavio Peulen
 * Copyright 2007-2022 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_DECAYRANGE_H
#define IMPBFF_DECAYRANGE_H

#include <IMP/bff/bff_config.h>

#include <vector>


IMPBFF_BEGIN_NAMESPACE


class DecayRange{

private:

    /// Start index of the decay range
    size_t _start = 0;

    /// Stop index of the decay range
    size_t _stop = -1;

public:

    void set_start(int v){
        _start = std::max(0, v);
    }

    size_t get_start(DecayCurve* d = nullptr) const {
        size_t re = _start;
        if(d != nullptr){
            re = std::max(0, (int) re);
        }
        return re;
    }

    void set_stop(int v) {
        if(v < 0){
            _stop = std::numeric_limits<size_t>::max();
        } else{
            _stop = v;
        }
    }

    size_t get_stop(DecayCurve* d = nullptr) const {
        size_t re = _stop;
        if(d != nullptr){
            re = (_stop < 0)? d->size() : std::min(d->size(), re);
        }
        return re;
    }

    void set_range(std::vector<int> v){
        if(v.size() > 1){
            set_start(v[0]);
            set_stop(v[1]);
        }
    }

    std::vector<int> get_range(DecayCurve* d = nullptr){
        return std::vector<int>({(int)get_start(d), (int)get_stop(d)});
    }

    void set(int start=0, int stop=-1){
        set_start(start);
        set_stop(stop);
    }

    DecayRange(int start, int stop){
        set_start(start);
        set_stop(stop);
    }

    virtual ~DecayRange() = default;

};

IMPBFF_END_NAMESPACE


#endif //IMPBFF_DECAYRANGE_H
