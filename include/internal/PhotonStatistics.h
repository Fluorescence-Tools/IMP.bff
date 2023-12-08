/**
 * \file IMP/bff/DecayStatistics.h
 * \brief Collection of different statictis for photon data
 *
 * \authors Thomas-Otavio Peulen
 * Copyright 2007-2022 IMP Inventors. All rights reserved.
 *
 */

#ifndef IMPBFF_PHOTONSTATISTICS_H
#define IMPBFF_PHOTONSTATISTICS_H

#include <IMP/bff/bff_config.h>

#include <cmath>
#include <vector>
#include <numeric>
#include <algorithm>
#include <iostream>
#include <cstring> /* strcmp */
#include "omp.h"


IMPBFF_BEGIN_NAMESPACE


/*!
 * Initialize an array containing pre-computed logarithm
 */
void init_fact();

/*!
 * Approximation of log(gamma function). See wikipedia
 *
 * https://en.wikipedia.org/wiki/Gamma_function#The_log-gamma_function
 *
 * @param t input of the gamma function
 * @return approximation of the logarithm of the gamma function
 */
double loggammaf(double t);

/*!
 * log-likelihood w(C|m) for Cp + 2Cs
 *
 * @param C number of counts in channel
 * @param m model function
 * @return log-likelihood w(C|m) for Cp + 2Cs
 */
double wcm(int C, double m);

/*!
 * Compute the -log-likelihood for Cp + 2Cs of a single micro time channel.
 *
 * Compute score of model counts in a parallel and perpendicular detection channel
 * and the experimental counts for a micro time channel.
 *
 * This function computes a score for the experimental counts (C) in a channel
 * where the experimental counts were computed by the sum of the counts in the
 * parallel (P) and the perpendicular (S) channel by the equation C = P + 2 S.
 *
 * This function considers that the number of counts C = P + 2S is not Poissonian.
 * The score relates to a maximum likelihood function.
 *
 * @param C number of experimental counts (P + 2 S) in a micro time channel
 * @param mp number of counts of the model in parallel detection channel
 * @param ms number of counts of the model in the perpendicular detection channel
 * @return
 */
double wcm_p2s(int C, double mp, double ms);

/*!
 * Compute the overall -log-likelihood for Cp + 2Cs for all micro time channels
 *
 * @param C array of experimental counts in Jordi format
 * @param M array model function in Jordi format
 * @param Nchannels number of micro time channels in parallel and perpendicular
 * (half the number of elements in C and M).
 * @return -log-likelihood for Cp + 2Cs for all micro time channels
 */
double Wcm_p2s(int* C, double* M, int Nchannels);

/*!
 * Compute overall 2I* for Cp + 2Cs
 *
 * This function computes the overall 2I* for the model function Cp + 2Cs that
 * is computed by parallel signal (Cp) and the perpendicular signal (Cs). For
 * the definition of function 2I* see "An Experimental Comparison of the
 * Maximum Likelihood Estimation and Nonlinear Least-Squares Fluorescence Lifetime
 * Analysis of Single Molecules, Michael Maus, Mircea Cotlet, Johan Hofkens,
 * Thomas Gensch, Frans C. De Schryver, J. Schaffer, and C. A. M. Seidel, Anal.
 * Chem. 2001, 73, 9, 2078–2086".
 *
 * @param C array of experimental counts in Jordi format
 * @param M array model function in Jordi format
 * @param Nchannels number of micro time channels in parallel and perpendicular
 * (half the number of elements in C and M).
 * @return 2I* for Cp + 2Cs
 */
double twoIstar_p2s(int* C, double* M, int Nchannels);


// overall log-likelihood w(C,M) for 1 ch data
double twoIstar_1ch(int* C, double* M, int Ndata);


/*!
 * Compute overall 2I* for Cp & Cs
 *
 * This function computes 2I* for Cp and Cs. Cp and Cs are the model signals in
 * the parallel and perpendicular channel. Contrary to twoIstar_p2s the overall
 * 2I* is the sum of 2I* for Cp and Cs.
 *
 * @param C array of experimental counts in Jordi format
 * @param M array model function in Jordi format
 * @param Nchannels number of micro time channels in parallel and perpendicular
 * (half the number of elements in C and M).
 * @return 2I* for Cp & Cs
 */
double twoIstar(int* C, double* M, int Nchannels);

/*!
 * Compute overall -log-likelihood for Cp & Cs
 *
 * @param C array of experimental counts in Jordi format
 * @param M array model function in Jordi format
 * @param Nchannels number of micro time channels in parallel and perpendicular
 * (half the number of elements in C and M).
 * @return -log-likelihood for Cp & Cs
 */
double Wcm(int* C, double* M, int Nchannels);


namespace statistics{


    /**
     * @brief Calculates the Neyman chi-squared statistic.
     * 
     * The Neyman chi-squared statistic is used to compare a model to observed data. It measures the difference between the expected model values and the observed data values, taking into account the uncertainties in the observed data.
     * 
     * @param data An array of observed data values.
     * @param model An array of expected model values.
     * @param start The starting index of the data and model arrays to consider.
     * @param stop The stopping index of the data and model arrays to consider.
     * @return The Neyman chi-squared statistic.
     */
    inline double neyman(double* data, double *model, int start, int stop){
        double chi2 = 0.0;
        for(int i = start; i < stop; i++){
            double mu = model[i];
            double m = std::max(1., data[i]);
            chi2 += (mu - m) * (mu - m) / m;
        }
        return chi2;
    }

    /**
     * Calculates the chi-squared value for comparing Poisson distributed data to a model.
     *
     * @param data   Pointer to an array of observed data.
     * @param model  Pointer to an array of model predictions.
     * @param start  The starting index of the data and model arrays.
     * @param stop   The stopping index of the data and model arrays.
     *
     * @return The chi-squared value.
     */
    inline double poisson(double* data, double *model, int start, int stop){
        double chi2 = 0.0;
        for(int i = start; i < stop; i++){
            double mu = model[i];
            double m = data[i];
            chi2 += 2 * std::abs(mu);
            chi2 -= 2 * m * (1 + log(std::max(1.0, mu) / std::max(1.0, m)));
        }
        return chi2;
    }

    /**
     * @brief Calculates the Pearson chi-square statistic for a given set of data and model.
     * 
     * The Pearson chi-square statistic is a measure of the difference between observed data and expected model values.
     * This function calculates the chi-square statistic by comparing the data and model values from the specified start index to the stop index.
     * 
     * @param data An array of double values representing the observed data.
     * @param model An array of double values representing the expected model values.
     * @param start The start index of the data and model arrays to consider.
     * @param stop The stop index (exclusive) of the data and model arrays to consider.
     * @return The calculated Pearson chi-square statistic.
     */
    inline double pearson(double* data, double *model, int start, int stop){
        double chi2 = 0.0;
        for(int i = start; i < stop; i++){
            double m = model[i];
            double d = data[i];
            if (m > 0) {
                chi2 += (m-d) / m;
            }
        }
        return chi2;
    }

    /**
     * @brief Calculates the chi-squared value for a given data set and model.
     * 
     * This function calculates the chi-squared value for a given data set and model.
     * It iterates over the data array from the start index to the stop index (exclusive),
     * and calculates the chi-squared value using the formula:
     * 
     * chi2 += (mu - m) * (mu - m) / mu + log(mu/mu_p) - (mu_p - m) * (mu_p - m) / mu_p
     * 
     * where:
     * - `mu` is the model value at the current index
     * - `m` is the data value at the current index
     * - `mu_p` is the square root of (0.25 + m * m) - 0.5
     * 
     * If `mu_p` is less than or equal to 1.e-12, the iteration continues to the next index.
     * 
     * Reference:
     * Xiangpan Ji, Wenqiang Gu, Xin Qian, Hanyu Wei, Chao Zhang.
     * Combined Neyman–Pearson chi-square: An improved approximation to the Poisson-likelihood chi-square.
     * Nuclear Instruments and Methods in Physics Research Section A: Accelerators, Spectrometers, Detectors and Associated Equipment, Volume 967, 2020, 163677.
     * https://doi.org/10.1016/j.nima.2020.163677
     * 
     * @param data The data array.
     * @param model The model array.
     * @param start The start index of the iteration.
     * @param stop The stop index of the iteration (exclusive).
     * @return The calculated chi-squared value.
     */
    inline double gauss(double* data, double *model, int start, int stop){
        double chi2 = 0.0;
        for(int i = start; i < stop; i++){
            double mu = model[i];
            double m = data[i];
            double mu_p = std::sqrt(.25 + m * m) - 0.5;
            if(mu_p <= 1.e-12) continue;
            chi2 += (mu - m) * (mu - m) / mu + std::log(mu/mu_p) - (mu_p - m) * (mu_p - m) / mu_p;
        }
        return chi2;
    }

    /**
     * @brief Calculates the chi-squared value for a given set of data and model.
     * 
     * The cnp function calculates the chi-squared value for a given set of data and model.
     * It iterates over the data array from the start index to the stop index (exclusive),
     * and for each element, it calculates the chi-squared value using the formula:
     * 
     *     chi2 += (mu - m) * (mu - m) / (3. / (1./m + 2./mu))
     * 
     * where m is the data value at the current index, mu is the model value at the current index,
     * and chi2 is the accumulated chi-squared value.
     * 
     * If the data value at the current index is less than or equal to 1e-12, it is skipped.
     * 
     * This function is based on the method described in the paper "Combined Neyman–Pearson chi-square:
     * An improved approximation to the Poisson-likelihood chi-square" by Xiangpan Ji, Wenqiang Gu, Xin Qian,
     * Hanyu Wei, and Chao Zhang.
     * 
     * Reference:
     * Xiangpan Ji, Wenqiang Gu, Xin Qian, Hanyu Wei, Chao Zhang.
     * Combined Neyman–Pearson chi-square: An improved approximation to the Poisson-likelihood chi-square.
     * Nuclear Instruments and Methods in Physics Research Section A: Accelerators, Spectrometers, Detectors and Associated Equipment, Volume 967, 2020, 163677.
     * https://doi.org/10.1016/j.nima.2020.163677
     * 
     * @param data The array of data values.
     * @param model The array of model values.
     * @param start The start index of the iteration (inclusive).
     * @param stop The stop index of the iteration (exclusive).
     * @return The calculated chi-squared value.
     */
    inline double cnp(double* data, double *model, int start, int stop){
        double chi2 = 0.0;
        for(int i = start; i < stop; i++){
            double m = data[i];
            double mu = model[i];
            if(m <= 1e-12) continue;
            chi2 += (mu - m) * (mu - m) / (3. / (1./m + 2./mu));
        }
        return chi2;
    }

    /// Sum of squared weighted residuals
    inline double sswr(double* data, double *model, double *data_noise, int start, int stop){
        double chi2 = 0.0;
        for(int i=start;i<stop;i++){
            double d = (data[i] - model[i]) / data_noise[i];
            chi2 += d * d;
        }
        return chi2;
    }

    /**
     * @brief Calculates the chi2 measure for counting data.
     *
     * This function calculates the chi2 measure for counting data using the specified method.
     * The chi2 measure is a statistical measure used to compare observed data with a model.
     *
     * @param data The observed data.
     * @param model The model data.
     * @param weights The weights for each data point.
     * @param x_min The minimum value of the data range to consider. Default is -1.
     * @param x_max The maximum value of the data range to consider. Default is -1.
     * @param type The type of chi2 measure to calculate. Default is "neyman".
     * @return The calculated chi2 measure.
     *
     * @see Xiangpan Ji, Wenqiang Gu, Xin Qian, Hanyu Wei, Chao Zhang. "Combined Neyman–Pearson Chi-square: 
     * An Improved Approximation to the Poisson-likelihood Chi-square." arXiv preprint arXiv:1903.07185 (2019).
     * @see https://arxiv.org/pdf/1903.07185.pdf
     */
    double chi2_counting(
            std::vector<double> &data,
            std::vector<double> &model,
            std::vector<double> &weights,
            int x_min = -1,
            int x_max = -1,
            const char* type="neyman"
    );
}

IMPBFF_END_NAMESPACE

#endif //IMPBFF_PHOTONSTATISTICS_H
