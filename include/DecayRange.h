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
#include <IMP/bff/DecayCurve.h>

#include <vector>


IMPBFF_BEGIN_NAMESPACE


class IMPBFFEXPORT DecayRange{

private:

    /// Start index of the decay range
    size_t _start = 0;

    /// Stop index of the decay range
    size_t _stop = -1;

public:

    void set_start(int v);

    size_t get_start(DecayCurve* d = nullptr) const;

    void set_stop(int v);

    size_t get_stop(DecayCurve* d = nullptr) const;

    void set_range(std::vector<int> v);

    std::vector<int> get_range(DecayCurve* d = nullptr);

    void set(int start=0, int stop=-1);

    DecayRange(int start, int stop);

    virtual ~DecayRange() = default;

};

IMPBFF_END_NAMESPACE


#endif //IMPBFF_DECAYRANGE_H
