#include "IMP/bff/DecayPattern.h"

IMPBFF_BEGIN_NAMESPACE


double DecayPattern::get_pattern_fraction() const{
    return pattern_fraction;
}

void DecayPattern::set_pattern_fraction(double v){
    pattern_fraction = std::max(0.0, std::min(1.0, std::abs(v)));
}

void DecayPattern::set_pattern(DecayCurve* v){
    set_data(v);
}

DecayCurve* DecayPattern::get_pattern(){
    return get_data();
}

void DecayPattern::set_constant_offset(double v){
    constant_offset = v;
}

double DecayPattern::get_constant_offset() const{
    return constant_offset;
}

void DecayPattern::add(DecayCurve* out) {
    if(out != nullptr && is_active()) {
        // resize background pattern to size of input decay
        resize(out->size());

        int start = get_start(out);
        int stop = get_stop(out);

        double f = get_pattern_fraction();
        auto bg = get_pattern();


        // Compute pre-factors of decay and added pattern
        double ds, bs;
        if (f > 0.0) {
            ds = std::accumulate(out->y.begin() + start, out->y.begin() + stop, 0.0);
            bs = std::accumulate(bg->y.begin() + start, bg->y.begin() + stop, 0.0);
        } else {
            ds = bs = 1.0;
        }
        double fd = (1. - f);
        double fb = f * ds / bs;

        // mix constant offset, decay, and pattern
        double offset = get_constant_offset();
        for (int i = start; i < stop; i++) {
            out->y[i] = out->y[i] * fd + bg->y[i] * fb + offset;
        }
    }
}

DecayPattern::DecayPattern(
        double constant_offset,
        DecayCurve*pattern,
        double pattern_fraction,
        int start,
        int stop,
        bool active
) : DecayModifier(pattern, start, stop, active) {
#if IMPBFF_VERBOSE
    std::clog << "DecayPattern::DecayPattern" << std::endl;
#endif
    set_pattern_fraction(pattern_fraction);
    set_constant_offset(constant_offset);
}


IMPBFF_END_NAMESPACE