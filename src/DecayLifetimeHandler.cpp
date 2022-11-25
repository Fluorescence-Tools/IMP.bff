#include <IMP/bff/DecayLifetimeHandler.h>

IMPBFF_BEGIN_NAMESPACE

double DecayLifetimeHandler::get_amplitude_threshold(){
    return std::abs(amplitude_threshold);
}

void DecayLifetimeHandler::set_amplitude_threshold(double v){
    amplitude_threshold = v;
}

bool DecayLifetimeHandler::get_use_amplitude_threshold(){
    return use_amplitude_threshold;
}

void DecayLifetimeHandler::set_use_amplitude_threshold(bool v){
    use_amplitude_threshold = v;
}

bool DecayLifetimeHandler::get_abs_lifetime_spectrum() const{
    return abs_lifetime_spectrum;
}

void DecayLifetimeHandler::set_abs_lifetime_spectrum(bool v){
    abs_lifetime_spectrum = v;
}

void DecayLifetimeHandler::set_lifetime_spectrum(std::vector<double> v) {
    _lifetime_spectrum = v;
}

void DecayLifetimeHandler::add_lifetime(double amplitude, double lifetime) {
    _lifetime_spectrum.emplace_back(amplitude);
    _lifetime_spectrum.emplace_back(lifetime);
}

std::vector<double>& DecayLifetimeHandler::get_lifetime_spectrum(){
    lt_ = _lifetime_spectrum;
    if(use_amplitude_threshold){
        discriminate_small_amplitudes(lt_.data(), lt_.size(), amplitude_threshold);
    }
    if (abs_lifetime_spectrum) {
        for(auto &v: lt_) v = std::abs(v);
    }
    return lt_;
}

void DecayLifetimeHandler::get_lifetime_spectrum(double **output_view, int *n_output) {
    auto lt = get_lifetime_spectrum();
    *output_view = lt.data();
    *n_output = (int) lt.size();
}

DecayLifetimeHandler::DecayLifetimeHandler(
        std::vector<double> lifetime_spectrum,
        bool use_amplitude_threshold,
        bool abs_lifetime_spectrum,
        double amplitude_threshold
){
    set_use_amplitude_threshold(use_amplitude_threshold);
    set_abs_lifetime_spectrum(abs_lifetime_spectrum);
    set_amplitude_threshold(amplitude_threshold);
    set_lifetime_spectrum(lifetime_spectrum);
}

IMPBFF_END_NAMESPACE
