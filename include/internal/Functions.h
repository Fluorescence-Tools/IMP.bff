#ifndef IMPBFF_FUNCTIONS_H
#define IMPBFF_FUNCTIONS_H

#include <IMP/bff/bff_config.h>

#include <vector>
#include <algorithm>  /* std::min, std::max */
#include <math.h>
#include <iostream>
#include <chrono>
#include <bson.h>

IMPBFF_BEGIN_NAMESPACE

namespace Functions {

    /*!
     * Shifts a vector by a floating number.
     *
     * This function shifts the y-axis and consider non-integer values by
     * determining the integer part of the shift and the
     * floating part of the shift value, e.g., for the shift 4.66
     * the integer part is 4 and the floating part is 0.66.
     * Next the array is shifted by the integer part and a copy
     * of the array is shifted by the integer part + 1. Finally,
     * the weighted sum of the both arrays is calculated.
     * @param value
     * @param x
     */
    void shift(double value, std::vector<double> &x);


    /*!
     * Rolls a vector by an integer
     *
     * @param value
     * @param y
     */
    void roll(int value, std::vector<double> &y);


    /*!
     * Allocates memory on a pointer to an array and copies the content
     * of a vector to the newly allocated memory.
     *
     * This function is mainly used for numpy array SWIGs
     *
     * @param v
     * @param out
     * @param nout
     */
    void copy_vector_to_array(std::vector<double> &v, double **out, int *nout);

    /*!
     * Copy an array and copies the content
     * of a vector to already allocated memory.
     *
     * This function is mainly used for numpy array SWIGs
     *
     * @param v
     * @param out
     * @param nout
     */
    void copy_vector_to_array(std::vector<double> &v, double *out, int nout);

    /*!
     *
     * @param in
     * @param nin
     * @param v
     */
    void copy_array_to_vector(double *in, int nin, std::vector<double> &v);

    /*!
     * This function copies two vectors of equal size to an interleaved array
     *
     * @param v1
     * @param v2
     * @param out
     * @param nout
     */
    void copy_two_vectors_to_interleaved_array(
            std::vector<double> &v1,
            std::vector<double> &v2,
            double **out, int *nout
    );

    /*!
     * Maps the array of values @param in place to values to the interval
     * (min, max) determined by the parameters @param bound_1 and @param bound_2.
     * The values are mapped to the interval by:
     *
     * \f$
     * max(bound_1, bound_2) - abs(bound_1-bound_2)/(exp(value / abs(bound_1-bound_2))+1)
     * \f$
     *
     * @tparam T The type of the values
     * @param values (InOut)
     * @param n_values (In)
     * @param lower_bound (In)
     * @param upper_bound (In)
     */
    template <typename T>
    void value2internal(
            T *values,
            int n_values,
            double lower_bound,
            double upper_bound
            )
    {
        T delta = upper_bound - lower_bound;
        for(int i = 0; i<n_values; i++)
        {
            values[i] = upper_bound - delta / (exp(values[i]/delta) + 1.0);
        }
    }

    template <typename T>
    void internal2value(T *values, int n_values, double lower_bound, double upper_bound){
        T delta = upper_bound - lower_bound;
        for(int i = 0; i<n_values; i++)
        {
            values[i] = delta * log(delta / (upper_bound - values[i]) - 1.0);
        }
    }

    template <typename T>
    void bound_values(T *values, int n_values, double lower_bound, double upper_bound){
        for(int i = 0; i<n_values; i++)
        {
          values[i] = (T) std::min(std::max((double)values[i], lower_bound), upper_bound);
        }
    }

    /*!
     * Calculates the discrete difference for an vector
     *
     * @param v
     * @return
     */
    std::vector<double> diff(std::vector<double> v);

    /*!
     * Returns the current time in milliseconds as a uint64
     * @return
     */
    uint64_t get_time();

    /*!
     * Adds the content in the bson_t document src to the document dst omitting the keys
     * provided by the vector skip.
     */
    void add_documents(bson_t *src, bson_t *dst, std::vector<std::string> skip);

    /*!
     * Returns true if the key associated to @param iter is in the list of vectors @param skip
     *
     * @param iter pointer to a bson_iter_t
     * @param skip vector of strings containing keys that are skipped by iter
     */
    bool bson_iter_skip(bson_iter_t *iter, std::vector<std::string> *skip);

    /*!
     * Returns a vector with a size that is is min(a.size(), b.size())
     * @tparam T The type of the vectors
     * @param a
     * @param b
     * @return a vector of size min(a.size(), b.size())
     */
    template <typename T>
    std::vector<T> get_vector_of_min_size(std::vector<T> a, std::vector<T>b){
        size_t size_min = std::min(a.size(), b.size());
        size_t size_max = std::max(a.size(), b.size());
        if(size_max != size_min){
            std::cerr << "The vectors have differ in length: " << size_min << "<" << size_max << std::endl;
        }
        std::vector<T> result_vector;
        result_vector.resize(size_min);
        return  result_vector;
    }
}

IMPBFF_END_NAMESPACE

#endif //IMPBFF_FUNCTIONS_H
