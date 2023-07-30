#include <IMP/bff/DecayRange.h>

IMPBFF_BEGIN_NAMESPACE


void DecayRange::set_start(int v){
    _start = std::max(0, v);
}

size_t DecayRange::get_start(DecayCurve* d) const {
    size_t re = _start;
    if(d != nullptr){
        re = std::max(0, (int) re);
    }
    return re;
}

void DecayRange::set_stop(int v) {
    if(v < 0){
        _stop = -1;
    } else{
        _stop = v;
    }
}

size_t DecayRange::get_stop(DecayCurve* d) const {
    size_t re = _stop;
    if(d != nullptr){
        re = (_stop < 0)? d->size() : std::min(d->size(), re);
    }
    return re;
}

void DecayRange::set_range(std::vector<int> v){
    if(v.size() > 1){
        set_start(v[0]);
        set_stop(v[1]);
    }
}

std::vector<int> DecayRange::get_range(DecayCurve* d){
    return std::vector<int>({(int)get_start(d), (int)get_stop(d)});
}

void DecayRange::set(int start, int stop){
    set_start(start);
    set_stop(stop);
}

DecayRange::DecayRange(int start, int stop){
    set_start(start);
    set_stop(stop);
}

IMPBFF_END_NAMESPACE