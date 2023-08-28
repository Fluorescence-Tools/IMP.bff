#ifndef IMPBFF_ACCESSIBLEVOLUMESAMPLER_H
#define IMPBFF_ACCESSIBLEVOLUMESAMPLER_H

#include <IMP/bff/bff_config.h>

#include <algorithm>
#include <vector>
#include <map>
#include <random>

#include <IMP/bff/internal/pcg_random.h>

IMPBFF_BEGIN_NAMESPACE


template <typename T> class InverseSampler
{
    // Draw elements using Inverse transform sampling
private:
    // probability is the key and index in the points array
    // is the value. Probabilities may have arbitrary normalization, i.e
    // last element in the map has maximum possible potability as the key.
    std::map<uint32_t, unsigned> map;
    T vec;
    using val_t = typename T::value_type;
    mutable pcg32_fast rng{pcg_extras::seed_seq_from<std::random_device>{}};

public:
    typename T::value_type get_random() const
    {
        // return element index, assuming that elements between the keys
        // in the map are uniformly sampled

        unsigned rnd = rng();
        // find the relevant index range (idxmin..idxmax]
        auto it = map.upper_bound(rnd);
        if (it == map.end()) {
            return vec[map.rbegin()->second];
        }
        uint32_t rmax = it->first;
        unsigned idxmax = it->second;
        --it;
        uint32_t rmin = it->first;
        unsigned idxmin = it->second;

        float relRnd = float(rnd - rmin) / (rmax - rmin);
        unsigned i = idxmin + unsigned(relRnd * (idxmax - idxmin));
        return vec[i];
    }

    template <typename KeyAccessor>
    InverseSampler(T &_vec, KeyAccessor accessor) : vec(_vec){
        // Since entries in the flat_map are only added for unique
        // weights(densities), the map stays small and fast.
        auto adder = [accessor](double sum, const val_t &p) {
            return sum + accessor(p);
        };
        const float sumWeights =
                std::accumulate(vec.begin(), vec.end(), 0.0, adder);
        const float pRatio = rng.max() / sumWeights;

        auto comparator = [accessor](const val_t &v1,
                                     const val_t &v2) -> bool {
            return accessor(v1) > accessor(v2);
        };
        if (not std::is_sorted(vec.begin(), vec.end(), comparator)) {
            std::sort(vec.begin(), vec.end(), comparator);
        }

        map.emplace(0, 0);
        double cumSum = accessor(vec[0]) * pRatio;
        float prevSlope = cumSum;
        for (unsigned i = 1; i < vec.size(); ++i) {
            float slope = accessor(vec[i]) * pRatio;
            if (slope < prevSlope * 0.999) {
                map.emplace(cumSum, i - 1);
                prevSlope = slope;
            }
            cumSum += slope;
        }
        map.emplace(rng.max(), vec.size() - 1);
    }
};

IMPBFF_END_NAMESPACE

#endif //IMP_BFF_ACCESSIBLEVOLUMESAMPLER_H