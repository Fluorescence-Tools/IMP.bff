#include "IMP/bff/DecayModifier.h"

IMPBFF_BEGIN_NAMESPACE


void DecayModifier::set_data(DecayCurve* v){
    if(v != nullptr){
        data = v;
    }
}

DecayCurve* DecayModifier::get_data(){
    auto re = default_data;
    if(data != nullptr){
        re = data;
    }
    return re;
}

bool DecayModifier::is_active() const {
    return _is_active;
}

void DecayModifier::set_active(bool v){
    _is_active = v;
}

void DecayModifier::set(DecayCurve* data, int start, int stop, bool active){
    DecayRange::set(start, stop);
    set_data(data);
    set_active(active);
}

void DecayModifier::resize(size_t n, double v) {
    default_data->resize(n, v);
    if(data != nullptr){
        data->resize(n, v);
    }
}

DecayModifier::DecayModifier(DecayCurve *data, int start, int stop, bool active) : DecayRange(start, stop){
    default_data = new DecayCurve();
    set_data(data);
    set_active(active);
}


IMPBFF_END_NAMESPACE