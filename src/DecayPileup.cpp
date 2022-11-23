#include <IMP/bff/DecayPileup.h>


IMPBFF_BEGIN_NAMESPACE

void DecayPileup::set_pile_up_model(std::string v){
    pile_up_model = v;
}

std::string DecayPileup::get_pile_up_model(){
    return pile_up_model;
}

void DecayPileup::set_repetition_rate(double v){
    repetition_rate = v;
}

double DecayPileup::get_repetition_rate(){
    return repetition_rate;
}

void DecayPileup::set_instrument_dead_time(double v){
    instrument_dead_time = std::abs(v);
}

double DecayPileup::get_instrument_dead_time(){
    return instrument_dead_time;
}

void DecayPileup::add(DecayCurve* out) override{
    IMP_USAGE_CHECK(data != nullptr, "Experimental data not set.");
#if IMPBFF_VERBOSE
    std::clog << "DecayPileup::add" << std::endl;
std::clog << "-- data->size(): " << data->size() << std::endl;
std::clog << "-- out->size(): " << out->size() << std::endl;
std::clog << "-- get_repetition_rate(): " << get_repetition_rate() << std::endl;
std::clog << "-- get_instrument_dead_time(): " << get_instrument_dead_time() << std::endl;
std::clog << "-- data->get_acquisition_time(): " << data->get_acquisition_time() << std::endl;
std::clog << "-- data->get_acquisition_time(): " << data->get_acquisition_time() << std::endl;
#endif
    if(out != nullptr && is_active()) {
        out->resize(data->size());
        decay_add_pile_up_to_model(
                out->y.data(), out->y.size(),
                data->y.data(), data->y.size(),
                get_repetition_rate(),
                get_instrument_dead_time(),
                data->get_acquisition_time(),
                pile_up_model.c_str()
        );
    }
}

DecayPileup::DecayPileup(
            DecayCurve* data,
            const char* pile_up_model,
            double repetition_rate,
            double instrument_dead_time,
            int start,
            int stop,
            bool active
    ) : DecayModifier(data, start, stop, active) {
#if IMPBFF_VERBOSE
        std::clog << "DecayPileup::DecayPileup" << std::endl;
std::clog << "-- pile_up_model: " << pile_up_model << std::endl;
std::clog << "-- repetition_rate: " << repetition_rate << std::endl;
std::clog << "-- instrument_dead_time: " << instrument_dead_time << std::endl;
std::clog << "-- start: " << start << std::endl;
std::clog << "-- stop: " << stop << std::endl;
std::clog << "-- active: " << active << std::endl;
#endif
        set_pile_up_model(pile_up_model);
        set_repetition_rate(repetition_rate);
        set_instrument_dead_time(instrument_dead_time);
    }

}

IMPBFF_END_NAMESPACE
