#include <IMP/bff/DecayScale.h>

IMPBFF_BEGIN_NAMESPACE

double DecayScale::get_number_of_photons(){
    double re = 0.0;
    DecayCurve* d = get_data();
    for(size_t i = get_start(d); i < get_stop(d); i++)
        re += d->y[i];
    return re;
}

double DecayScale::get_constant_background() const {
    return _constant_background;
}

void DecayScale::set_constant_background(double v){
    _constant_background = v;
}

void DecayScale::set_blank_outside(double v){
    _blank_outside = v;
}

bool DecayScale::get_blank_outside(){
    return _blank_outside;
}

void DecayScale::set(
        DecayCurve* data,
        double constant_background,
        int start, int stop, bool active,
        bool blank_outside
){
    set_data(data);
    set_constant_background(constant_background);
    set_start(start);
    set_stop(stop);
    set_active(active);
    set_blank_outside(blank_outside);
}

DecayScale::DecayScale(
        DecayCurve* data,
        double constant_background,
        int start, int stop, bool active,
        bool blank_outside
) : DecayModifier(data, start, stop, active){
    _blank_outside = blank_outside;
    set_constant_background(constant_background);
}

void DecayScale::add(DecayCurve* decay){
#if IMPBFF_VERBOSE
    std::clog << "DecayScale::add" << std::endl;
#endif
    IMP_USAGE_CHECK(data != nullptr, "Experimental data not set - cannot scale model.");
    if(is_active()){
        auto d = get_data();
        double scale = 0.0;
        double* model = decay->y.data();
        double* data = d->y.data();
        double* squared_data_weights = d->ey.data();
        int n_model = decay->y.size();
        int start = (int) get_start(decay);
        int stop = (int) get_stop(decay);
        auto constant_background = get_constant_background();
        decay_rescale_w_bg(
                model, data,
                squared_data_weights, constant_background,
                &scale, start, stop);
        if(_blank_outside){
            for(int i = 0; i < start; i++){
                model[i] = 0.0;
            }
            for(int i = stop; i < n_model; i++){
                model[i] = 0.0;
            }
        }
    }
}

IMPBFF_END_NAMESPACE