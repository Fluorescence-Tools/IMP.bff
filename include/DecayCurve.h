/**
 * \file IMP/bff/DecayCurve.h
 * \brief Class for fluorescence decay curves
 *
 * \authors Thomas-Otavio Peulen
 * Copyright 2007-2022 IMP Inventors. All rights reserved.
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

/**
 * \class DecayCurve
 * \brief Class for fluorescence decay curves
 *
 * This class represents fluorescence decay curves. It stores the x-values,
 * y-values, and error values of the curve. It also provides various operations
 * and transformations on the curve.
 */
class IMPBFFEXPORT DecayCurve {
    friend class DecayConvolution;
    friend class DecayPileup;
    friend class DecayLinearization;
    friend class DecayScale;
    friend class DecayScore;
    friend class DecayPattern;

private:

    int noise_model = NOISE_NA; /**< The noise model type */
    double current_shift = 0.0; /**< The current shift applied to the curve */
    double acquisition_time = std::numeric_limits<double>::max(); /**< The acquisition time of the curve */
    std::vector<double> dx = {1.0}; /**< The x-value differences */
    std::vector<double> x = std::vector<double>(); /**< The x-values of the curve */
    std::vector<double> _y; /**< The initial y-values of the curve */
    std::vector<double> y = std::vector<double>(); /**< The transformed y-values of the curve */
    std::vector<double> ey = std::vector<double>(); /**< The error values of the curve */

    /**
     * \brief Compute the noise of the curve
     *
     * This function computes the noise of the curve based on the noise model type.
     *
     * \param noise_model The noise model type
     */
    void compute_noise(int noise_model = NOISE_POISSON);

public:

    /**
     * \brief Shift an array by a given value
     *
     * This static function shifts the elements of an array by a given value.
     *
     * \param input The input array
     * \param n_input The size of the input array
     * \param shift The shift value
     * \return The shifted array
     */
    static std::vector<double> shift_array(
            double *input, int n_input,
            double shift
    );
    
    /**
     * \brief Get the size of the curve
     *
     * This function returns the size of the curve, i.e., the number of data points.
     *
     * \return The size of the curve
     */
    size_t size() const;
    
    /**
     * \brief Check if the curve is empty
     *
     * This function checks if the curve is empty, i.e., it has no data points.
     *
     * \return True if the curve is empty, false otherwise
     */
    bool empty();
    
    /**
     * \brief Get the x-value differences
     *
     * This function returns the differences between consecutive x-values.
     *
     * \return The x-value differences
     */
    std::vector<double> get_dx();
    
    /**
     * \brief Resize the curve
     *
     * This function resizes the curve to a given size. If the new size is larger
     * than the current size, the new elements are initialized with a given value.
     *
     * \param n The new size of the curve
     * \param v The value to initialize the new elements with (default: 0.0)
     * \param dx The x-value difference (default: 1.0)
     */
    void resize(size_t n, double v=0.0, double dx=1.0);
    
    /**
     * \brief Get the average x-value difference
     *
     * This function returns the average difference between consecutive x-values.
     *
     * \return The average x-value difference
     */
    double get_average_dx();
    
    /**
     * \brief Get the x-values of the curve
     *
     * This function returns the x-values of the curve.
     *
     * \return The x-values of the curve
     */
    std::vector<double>& get_x();

    void set_x(const std::vector<double>& v);

    /**
     * @brief Sets the x values of the decay curve.
     * @param input The array of x values.
     * @param n_input The number of x values.
     */
    void set_x(double* input, int n_input);

    /**
     * @brief Returns the y values of the decay curve.
     * @return The y values.
     */
    std::vector<double>& get_y();

    /**
     * @brief Sets the y values of the decay curve.
     * @param v The y values.
     */
    void set_y(std::vector<double>& v);

    /**
     * @brief Sets the y values of the decay curve.
     * @param input The array of y values.
     * @param n_input The number of y values.
     */
    void set_y(double* input, int n_input);

    /**
     * @brief Returns the error values of the decay curve.
     * @return The error values.
     */
    std::vector<double>& get_ey();

    /**
     * @brief Sets the error values of the decay curve.
     * @param v The error values.
     */
    void set_ey(std::vector<double>& v);

    /**
     * @brief Sets the error values of the decay curve.
     * @param input The array of error values.
     * @param n_input The number of error values.
     */
    void set_ey(double* input, int n_input);

    /**
     * @brief Sets the acquisition time of the decay curve.
     * @param v The acquisition time.
     */
    void set_acquisition_time(double v);

    /**
     * @brief Returns the acquisition time of the decay curve.
     * @return The acquisition time.
     */
    double get_acquisition_time() const;

    /**
     * @brief Sets the time shift of the decay curve.
     * @param v The time shift.
     */
    void set_shift(double v);

    /**
     * @brief Returns the time shift of the decay curve.
     * @return The time shift.
     */
    double get_shift();

    /**
     * @brief Returns the JSON representation of the decay curve.
     * @return The JSON string.
     */
    std::string get_json() const;

    /**
     * @brief Reads the decay curve from a JSON string.
     * @param json_string The JSON string.
     * @return 0 if successful, -1 otherwise.
     */
    int read_json(std::string json_string);

    DecayCurve(
            std::vector<double> x = std::vector<double>(),
            std::vector<double> y  = std::vector<double>(),
            std::vector<double> ey  = std::vector<double>(),
            double acquisition_time = std::numeric_limits<double>::max(),
            int noise_model = NOISE_POISSON,
            int size = -1
    );

    /**
     * @brief Calculates the sum of y values within a given range.
     * @param start The start index.
     * @param stop The stop index.
     * @return The sum of y values.
     */
    double sum(int start=0, int stop=-1);

     /**
     * @brief Applies a simple moving average (SMA) filter to the data.
     * @param start The start index.
     * @param stop The stop index.
     * @param n_window The window size for the SMA filter.
     * @param normalize Flag indicating whether to normalize the filtered data.
     */
    void apply_simple_moving_average(int start, int stop, int n_window=5, bool normalize=false);

    /**
     * @brief Overloads the addition operator to add a constant value to the decay curve.
     * @param v The constant value to be added.
     * @return A new DecayCurve object representing the result of the addition.
     */
    DecayCurve& operator+(double v) const;

    /**
     * @brief Overloads the subtraction operator to subtract a constant value from the decay curve.
     * @param v The constant value to be subtracted.
     * @return A new DecayCurve object representing the result of the subtraction.
     */
    DecayCurve& operator-(double v) const;

    /**
     * @brief Overloads the multiplication operator to multiply the decay curve by a constant value.
     * @param v The constant value to be multiplied by.
     * @return A new DecayCurve object representing the result of the multiplication.
     */
    DecayCurve& operator*(double v) const;

    /**
     * @brief Overloads the division operator to divide the decay curve by a constant value.
     * @param v The constant value to be divided by.
     * @return A new DecayCurve object representing the result of the division.
     */
    DecayCurve& operator/(double v) const;

    /**
     * @brief Adds a value to the decay curve.
     * @param v The value to be added.
     * @return A reference to the modified DecayCurve object.
     */
    DecayCurve& operator+=(double v);

    /**
     * @brief Subtracts a value from the decay curve.
     * @param v The value to be subtracted.
     * @return A reference to the modified DecayCurve object.
     */
    DecayCurve& operator-=(double v);

    /**
     * @brief Multiplies the decay curve by a value.
     * @param v The value to be multiplied by.
     * @return A reference to the modified DecayCurve object.
     */
    DecayCurve& operator*=(double v);

    /**
     * @brief Divides the decay curve by a value.
     * @param v The value to be divided by.
     * @return A reference to the modified DecayCurve object.
     */
    DecayCurve& operator/=(double v);

    /**
     * @brief Overloads the '+' operator to add two DecayCurve objects.
     * @param other The DecayCurve object to be added.
     * @return A new DecayCurve object that is the result of the addition.
     */
    DecayCurve& operator+(const DecayCurve& other) const;

    /**
     * @brief Overloads the '-' operator to subtract two DecayCurve objects.
     * @param other The DecayCurve object to be subtracted.
     * @return A new DecayCurve object that is the result of the subtraction.
     */
    DecayCurve& operator-(const DecayCurve& other) const;

    /**
     * @brief Overloads the '*' operator to multiply two DecayCurve objects.
     * @param other The DecayCurve object to be multiplied.
     * @return A new DecayCurve object that is the result of the multiplication.
     */
    DecayCurve& operator*(const DecayCurve& other) const;

    /**
     * @brief Overloads the '/' operator to divide two DecayCurve objects.
     * @param other The DecayCurve object to be divided.
     * @return A new DecayCurve object that is the result of the division.
     */
    DecayCurve& operator/(const DecayCurve& other) const;

    /**
     * @brief Assignment operator.
     * @param other The DecayCurve object to be assigned.
     * @return A reference to the assigned DecayCurve object.
     */
    DecayCurve& operator=(const DecayCurve& other);

    /**
     * @brief Shifts the curve by a float value.
     * @param v The value to shift the curve by.
     * @return A reference to the modified DecayCurve object.
     */
    DecayCurve& operator<<(double v);

};

IMPBFF_END_NAMESPACE


#endif //IMPBFF_DECAYCURVE_H
