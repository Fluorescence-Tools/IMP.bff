#include "IMP/bff/DecayCurve.h"


IMPBFF_BEGIN_NAMESPACE

void DecayCurve::set_shift(double v){
    current_shift = v;
    double* d = _y.data(); int n = _y.size();
    y = shift_array(d, n, current_shift);
}

double DecayCurve::get_shift(){
    return current_shift;
}

std::vector<double> DecayCurve::shift_array(
        double* input, int n_input,
        double shift
){
    auto out = std::vector<double>(n_input, 0.0);
    if (n_input > 0) {
        int s0 = (int) std::floor(shift);
        int s1 = (int) std::ceil(shift);
        int t0 = IMP::bff::mod_p(s0, n_input);
        int t1 = IMP::bff::mod_p(s1, n_input);
        double ts_f = shift - s0;
#if IMPBFF_VERBOSE
        std::clog << "DecayCurve::shift_array" << std::endl;
        std::clog << "-- n_input: " << n_input << std::endl;
        std::clog << "-- shift: " << shift << std::endl;
        std::clog << "-- shift t0: " << t0 << std::endl;
        std::clog << "-- shift t1: " << t1 << std::endl;
        std::clog << "-- shift float: " << ts_f << std::endl;
#endif
        std::vector<double> tmp1(input, input + n_input);
        std::rotate(tmp1.begin(), tmp1.begin() + t0, tmp1.end());

        std::vector<double> tmp2(input, input + n_input);
        std::rotate(tmp2.begin(), tmp2.begin() + t1, tmp2.end());

        for (int i = 0; i < n_input; i++) {
            out[i] = tmp1[i] * (1.0 - ts_f) + tmp2[i] * ts_f;
        }
    }
    return out;
}

void DecayCurve::compute_noise(int noise_model){
    if(noise_model == NOISE_POISSON){
        for (size_t i = 0; i < size(); i++){
            ey[i] = std::max(1.0, std::sqrt(std::abs(y[i])));
        }
    }
}

DecayCurve& DecayCurve::operator+(const DecayCurve& other) const
{
    size_t n_max = std::min(size(), other.size());
    auto d = new DecayCurve();
    d->resize(n_max);
#if IMPBFF_VERBOSE
    std::clog << "DecayCurve::operator+(const DecayCurve& other)" << std::endl;
    std::clog << "-- this->size():" << this->size() << std::endl;
    std::clog << "-- other->size():" << other.size() << std::endl;
    std::clog << "-- new->size():" << d->size() << std::endl;
#endif
    for(size_t i=0; i < n_max; i++){
        d->x[i] = x[i];
        d->_y[i] = y[i] + other.y[i];
        d->y[i] = d->_y[i];

        d->ey[i] = std::sqrt(ey[i]*ey[i] + other.ey[i]*other.ey[i]);
    }
    d->acquisition_time = acquisition_time + other.acquisition_time;
    return *d;
}

DecayCurve& DecayCurve::operator-(const DecayCurve& other) const
{
    size_t n_max = std::min(size(), other.size());
    auto d = new DecayCurve();
    d->resize(n_max);
#if IMPBFF_VERBOSE
    std::clog << "DecayCurve::operator+(const DecayCurve& other)" << std::endl;
std::clog << "-- this->size():" << this->size() << std::endl;
std::clog << "-- other->size():" << other.size() << std::endl;
std::clog << "-- new->size():" << d->size() << std::endl;
#endif
    for(size_t i=0; i < n_max; i++){
        d->x[i] = x[i];
        d->_y[i] = y[i] - other.y[i];
        d->y[i] = d->_y[i];

        d->ey[i] = std::sqrt(ey[i]*ey[i] + other.ey[i]*other.ey[i]);
    }
    d->acquisition_time = acquisition_time + other.acquisition_time;
    return *d;
}

DecayCurve& DecayCurve::operator*(const DecayCurve& other) const
{
    size_t n_max = std::min(size(), other.size());
    auto d = new DecayCurve();
    d->resize(n_max);
#if IMPBFF_VERBOSE
    std::clog << "DecayCurve::operator*(const DecayCurve& other)" << std::endl;
    std::clog << "-- this->size():" << this->size() << std::endl;
    std::clog << "-- other->size():" << other.size() << std::endl;
    std::clog << "-- new->size():" << d->size() << std::endl;
#endif
    for(size_t i=0; i < n_max; i++){
        d->x[i] = x[i];
        d->_y[i] = y[i] * other.y[i];
        d->y[i] = d->_y[i];

        double f1 = other.y[i]*ey[i];
        double f2 = y[i]*other.ey[i];
        d->ey[i] = std::sqrt(f1*f1 + f2*f2);
    }
    d->acquisition_time = acquisition_time + other.acquisition_time;
    return *d;
}

DecayCurve& DecayCurve::operator/(const DecayCurve& other) const
{
    size_t n_max = std::min(size(), other.size());
    auto d = new DecayCurve();
    d->resize(n_max);
#if IMPBFF_VERBOSE
    std::clog << "DecayCurve::operator*(const DecayCurve& other)" << std::endl;
std::clog << "-- this->size():" << this->size() << std::endl;
std::clog << "-- other->size():" << other.size() << std::endl;
std::clog << "-- new->size():" << d->size() << std::endl;
#endif
    for(size_t i=0; i < n_max; i++){
        d->x[i] = x[i];
        d->_y[i] = y[i] / other.y[i];
        d->y[i] = d->_y[i];

        double a_b = y[i] / other.y[i];
        double a_1 = abs(1. + ey[i] / y[i]);
        double b_1 = abs(1. + other.ey[i] / other.y[i]);
        d->ey[i] = a_b * a_1 / b_1;
    }
    d->acquisition_time = acquisition_time + other.acquisition_time;
    return *d;
}

DecayCurve& DecayCurve::operator+(const double v) const {
    auto d = new DecayCurve();
    d->resize(size());
    for (size_t i = 0; i < size(); i++) {
        d->x[i] = x[i];
        d->_y[i] = y[i] + v;
        d->y[i] = d->_y[i];

        d->ey[i] = ey[i];
    }
    return *d;
}

DecayCurve& DecayCurve::operator-(const double v) const {
    auto d = new DecayCurve();
    d->resize(size());
    for (size_t i = 0; i < size(); i++) {
        d->x[i] = x[i];
        d->_y[i] = y[i] - v;
        d->y[i] = d->_y[i];

        d->ey[i] = ey[i];
    }
    return *d;
}

DecayCurve& DecayCurve::operator*(const double v) const {
    auto d = new DecayCurve();
    d->resize(size());
    for (size_t i = 0; i < size(); i++) {
        d->x[i] = x[i];
        d->_y[i] = y[i] * v;
        d->y[i] = d->_y[i];

        d->ey[i] = ey[i] * v;
    }
    return *d;
}

DecayCurve& DecayCurve::operator/(const double v) const {
    auto d = new DecayCurve();
    d->resize(size());
    for (size_t i = 0; i < size(); i++) {
        d->x[i] = x[i];
        d->_y[i] = y[i] / v;
        d->y[i] = d->_y[i];

        d->ey[i] = ey[i] / v;
    }
    return *d;
}

DecayCurve& DecayCurve::operator+=(const double v) {
    for (size_t i = 0; i < size(); i++) {
        y[i] += v;
    }
    return *this;
}

DecayCurve& DecayCurve::operator-=(const double v) {
    for (size_t i = 0; i < size(); i++) {
        y[i] -= v;
    }
    return *this;
}

DecayCurve& DecayCurve::operator*=(const double v) {
#if IMPBFF_VERBOSE
    std::clog << "DecayCurve& DecayCurve::operator*=(const double v)" << std::endl;
    std::clog << "-- v:" << v << std::endl;
#endif
    for (size_t i = 0; i < size(); i++) {
        y[i] *= v;
        ey[i] *= v;
    }
    return *this;
}

DecayCurve& DecayCurve::operator/=(const double v) {
#if IMPBFF_VERBOSE
    std::clog << "DecayCurve& DecayCurve::operator*=(const double v)" << std::endl;
std::clog << "-- v:" << v << std::endl;
#endif
    for (size_t i = 0; i < size(); i++) {
        y[i] /= v;
        ey[i] /= v;
    }
    return *this;
}

DecayCurve& DecayCurve::operator=(const DecayCurve& other){
    // Guard self assignment
    if(this == &other){
        return *this;
    }

    // Assign other to this
    x = other.x;
    y = other.y; _y = other._y;
    ey = other.ey;

    return *this;
}

DecayCurve& DecayCurve::operator<<(double v){
    set_shift(v);
    return *this;
}


void DecayCurve::apply_simple_moving_average(int start, int stop, int n_window, bool normalize){
    double m = std::accumulate(_y.begin() + start, _y.begin() + stop, 0.0);
    m /= std::abs(stop - start);
    for(int i = 0; i < start; i++)
        _y[i] = m;
    for(int i = stop; i < (int) size(); i++)
        _y[i] = m;

    MovingAverage<double, double> ma(n_window);
    for(int i = 0; i < (int) size(); i++){
        _y[i] = ma(_y[i]);
    }

    if(normalize){
        for(size_t i = 0; i < size(); i++)
            _y[i] /= m;
    }
    y = _y;
}


IMPBFF_END_NAMESPACE
