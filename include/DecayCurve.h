/**
 * \file IMP/bff/DecayCurve.h
 * \brief Class for fluorescence decay curves
 *
 * \authors Thomas-Otavio Peulen
 * Copyright 2007-2022 IMP Inventors. All rights reserved.
 *
 */

#ifndef IMPBFF_DECAYCURVE_H
#define IMPBFF_DECAYCURVE_H

#include <IMP/bff/bff_config.h>

#include <cmath>
#include <vector>
#include <numeric>
#include <algorithm>
#include <iostream> /* std::cerr */
#include <cstring> /* strcmp */
#include <limits> /* std::numeric_limits */


IMPBFF_BEGIN_NAMESPACE


enum NoiseModelTypes{
    NOISE_NA, NOISE_POISSON
};



class DecayCurve {

    friend class DecayConvolution;
    friend class DecayPileup;
    friend class DecayLinearization;
    friend class DecayScale;
    friend class DecayScore;
    friend class DecayPattern;

private:

    void compute_noise(int noise_model = NOISE_POISSON);

    int noise_model = NOISE_NA;
    double current_shift = 0.0;
    double acquisition_time = std::numeric_limits<double>::max();
    std::vector<double> dx = {1.0};

    /// x-values
    std::vector<double> x = std::vector<double>();
    /// initial y values
    std::vector<double> _y;
    /// transformed (e.g., shifted) y values
    std::vector<double> y = std::vector<double>();
    /// error of y values
    std::vector<double> ey = std::vector<double>();


public:

    static std::vector<double> shift_array(
            double *input, int n_input,
            double shift,
            bool set_outside = true,
            double outside_value = 0.0
    );

    size_t size() const{
        return x.size();
    }

    bool empty(){
        return x.empty();
    }

    std::vector<double> get_dx(){
        std::vector<double> dx(size(), 0.0);
        if(!empty()){
            for(size_t i = 1; i < x.size(); i++){
                dx[i] = x[i] - x[i - 1];
            }
        }
        return dx;
    }

    void resize(size_t n, double v=0.0, double dx=1.0){
        size_t old_size = size();

        // get last dx to extend the axis linearly
        auto dx_vec = get_dx();
        if(!dx_vec.empty()) dx = dx_vec[dx_vec.size() - 1];

        x.resize(n);
        y.resize(n, v); _y.resize(n, v);
        ey.resize(n, std::numeric_limits<double>::max());

        for(size_t i = old_size; i < n; i++){
            if(i > 0) x[i] = x[i - 1] + dx;
        }
    }

    double get_average_dx(){
        auto dx = get_dx();
        return std::accumulate(dx.begin(), dx.end(), 0.0) / (double) dx.size();
    }

    /*-------*/
    /* x     */
    /*-------*/
    std::vector<double>& get_x(){
        return x;
    }

    void set_x(const std::vector<double>& v){
        resize(v.size());
        x = v;
    }

    void set_x(double *input, int n_input){
        auto vec = std::vector<double>(input, input+n_input);
        set_x(vec);
    }

    /*-------*/
    /* y     */
    /*-------*/
    std::vector<double>& get_y(){
        return y;
    }

    void set_y(std::vector<double>& v){
        resize(v.size());
        y = v; _y = v;
        compute_noise(noise_model);
    }

    void set_y(double *input, int n_input){
        auto vec = std::vector<double>(input, input+n_input);
        set_y(vec);
    }

    /*-------*/
    /* ey    */
    /*-------*/
    std::vector<double>& get_ey(){
        return ey;
    }

    void set_ey(std::vector<double>& v){
        resize(v.size());
        ey = v;
    }

    void set_ey(double *input, int n_input){
        auto vec = std::vector<double>(input, input+n_input);
        set_ey(vec);
    }

    /*--------------------*/
    /* aquisition time    */
    /*--------------------*/
    void set_acquisition_time(double v) {
        if(acquisition_time < 0)
            acquisition_time = std::numeric_limits<double>::max();
        acquisition_time = v;
    }

    double get_acquisition_time() const {
        return acquisition_time;
    }

    /*--------------------*/
    /* time shift         */
    /*--------------------*/
    void set_shift(double v, double outside=0.0, bool shift_set_outside=true){
        current_shift = v;
        double* d = _y.data(); int n = _y.size();
        y = shift_array(d, n, current_shift, shift_set_outside, outside);
    }

    double get_shift(){
        return current_shift;
    }
    
    DecayCurve(
            std::vector<double> x = std::vector<double>(),
            std::vector<double> y  = std::vector<double>(),
            std::vector<double> ey  = std::vector<double>(),
            double acquisition_time = std::numeric_limits<double>::max(),
            int noise_model = NOISE_POISSON,
            int size = -1
    ){
        this->acquisition_time = acquisition_time;
        this->noise_model = noise_model;
        set_ey(ey);
        set_y(y);
        set_x(x);
        if(size > 0) resize(size);
    }

    DecayCurve& operator+(const DecayCurve& other) const;
    DecayCurve& operator*(const DecayCurve& other) const;
    DecayCurve& operator+(double v) const;
    DecayCurve& operator*(double v) const;
    DecayCurve& operator+=(double v);
    DecayCurve& operator*=(double v);
    DecayCurve& operator=(const DecayCurve& other);

};

IMPBFF_END_NAMESPACE


#endif //IMPBFF_DECAYCURVE_H
