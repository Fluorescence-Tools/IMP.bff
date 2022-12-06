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
#include <string>
#include <cstring> /* strcmp */
#include <limits> /* std::numeric_limits */

#include <IMP/bff/internal/MovingAverage.h>
#include <IMP/bff/internal/json.h>
#include <IMP/bff/DecayRoutines.h>



IMPBFF_BEGIN_NAMESPACE


enum NoiseModelTypes{
    NOISE_NA,
    NOISE_POISSON
};



class IMPBFFEXPORT DecayCurve {

    friend class DecayConvolution;
    friend class DecayPileup;
    friend class DecayLinearization;
    friend class DecayScale;
    friend class DecayScore;
    friend class DecayPattern;

private:

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

    void compute_noise(int noise_model = NOISE_POISSON);

public:

    static std::vector<double> shift_array(
            double *input, int n_input,
            double shift
    );

    size_t size() const;

    bool empty();

    std::vector<double> get_dx();

    void resize(size_t n, double v=0.0, double dx=1.0);

    double get_average_dx();

    /*-------*/
    /* x     */
    /*-------*/
    std::vector<double>& get_x();

    void set_x(const std::vector<double>& v);

    void set_x(double *input, int n_input);

    /*-------*/
    /* y     */
    /*-------*/
    std::vector<double>& get_y();

    void set_y(std::vector<double>& v);

    void set_y(double *input, int n_input);

    /*-------*/
    /* ey    */
    /*-------*/
    std::vector<double>& get_ey();

    void set_ey(std::vector<double>& v);

    void set_ey(double *input, int n_input);

    /*--------------------*/
    /* acquisition time   */
    /*--------------------*/
    void set_acquisition_time(double v);

    double get_acquisition_time() const;

    /*--------------------*/
    /* time shift         */
    /*--------------------*/
    void set_shift(double v);
    double get_shift();

    std::string get_json() const;

    /// Read from JSON string
    bool read_json(std::string json_string);


    DecayCurve(
            std::vector<double> x = std::vector<double>(),
            std::vector<double> y  = std::vector<double>(),
            std::vector<double> ey  = std::vector<double>(),
            double acquisition_time = std::numeric_limits<double>::max(),
            int noise_model = NOISE_POISSON,
            int size = -1
    );

    double sum(int start=0, int stop=-1);

    /// Apply a simple moving average (SMA) filter to the data
    void apply_simple_moving_average(int start, int stop, int n_window=5, bool normalize=false);

    DecayCurve& operator+(double v) const;
    DecayCurve& operator-(double v) const;
    DecayCurve& operator*(double v) const;
    DecayCurve& operator/(double v) const;

    DecayCurve& operator+=(double v);
    DecayCurve& operator-=(double v);
    DecayCurve& operator*=(double v);
    DecayCurve& operator/=(double v);

    DecayCurve& operator+(const DecayCurve& other) const;
    DecayCurve& operator-(const DecayCurve& other) const;
    DecayCurve& operator*(const DecayCurve& other) const;
    DecayCurve& operator/(const DecayCurve& other) const;

    DecayCurve& operator=(const DecayCurve& other);

    /// Shift the curve by a float
    DecayCurve& operator<<(double v);

};

IMPBFF_END_NAMESPACE


#endif //IMPBFF_DECAYCURVE_H
