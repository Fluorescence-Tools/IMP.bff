/**
 * \file IMP/bff/Histogram.h
 * \brief Simple histogram class
 *
 * \authors Thomas-Otavio Peulen
 * Copyright 2007-2022 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_HISTOGRAM_H
#define IMPBFF_HISTOGRAM_H

#include <IMP/bff/bff_config.h>

#include <algorithm>
#include <vector>
#include <cstdio>
#include <string>
#include <map>
#include <cmath>


IMPBFF_BEGIN_NAMESPACE

/// Histogram axis types
typedef enum{
    AXIS_LIN,         /// Linear spaced axis
    AXIS_LOG10,         /// Log spaced axis
    AXIS_ARB          /// Arbitrary spaced axis
} HistogramAxisTypes;


/// Linear spaced bin edges
template <typename T>
inline void linspace(double start, double stop, T *bin_edges, int n_bins){
    double bin_width = (stop - start) / n_bins;
    for(int i=0; i<n_bins; i++){
        bin_edges[i] = start + ((T) i) * bin_width;
    }
}


/// Log spaced bin edges
template <typename T>
inline void logspace(double start, double stop, T *bin_edges, int n_bins){
    linspace(std::log(start), std::log(stop), bin_edges, n_bins);
    for(int i=0; i<n_bins; i++){
        bin_edges[i] = std::pow(10.0, bin_edges[i]);
    }
}


/*! Searches for the bin index of a value within a list of bin edges
 *
 * If a value is inside the bounds find the bin.
 * The search partitions the bin_edges in upper and lower ranges and
 * adapts the edge for the upper and lower range depending if the target
 * value is bigger or smaller than the bin in the middle.

 * @tparam T type of the histogram values
 * @param value value that is searched
 * @param bin_edges array of bin edges
 * @param n_bins number of bins
 * @return negative value if the search value is out of the bounds. Otherwise the bin number
 * is returned.
 */
template <typename T>
inline int search_bin_idx(T value, T *bin_edges, int n_bins){
    int b, e, m;

    // ignore values outside of the bounds
    if ((value < bin_edges[0]) || (value > bin_edges[n_bins - 2])) {
        return -1;
    }

    b = 0;
    e = n_bins;
    do {
        m = (e - b) / 2 + b;
        if (value > bin_edges[m]) {
            b = m;
        } else {
            e = m;
        }
    } while ((value < bin_edges[m]) || (value >= bin_edges[m + 1]));
    return m;
}


/*!
 * Calculates for a linear axis the bin index for a particular value.
 *
 * @tparam T Type of the values in the histogram
 * @param begin value of first bin
 * @param bin_width width of histogram bins
 * @param value value the bin is searched for
 * @return
 */
template <typename T>
inline int calc_bin_idx(T begin, T bin_width, T value){
    return  ((value - begin) / bin_width);
}


/// Template class for Histogram axis
template<class T>
class HistogramAxis{

private:

    /// Name of the axis
    std::string name;

    /// First value of the axis
    double begin;

    /// Last value of the axis
    double end;

    /// Number of bins of the axis
    int n_bins;

    /// Axis bin width
    double bin_width; // for logarithmic spacing the bin_width in logarithms

    /// Vector of bin edges
    std::vector<T> bin_edges;

    /// Type of axis
    int axis_type;

public:

    /*!
     * Recalculates the bin edges of the axis
    */
    void update(){
        bin_edges.resize(n_bins);
        switch (HistogramAxis::axis_type){
            case 0:
                // linear axis
                bin_width = (end - begin) / n_bins;
                linspace(begin, end, bin_edges.data(), bin_edges.size());
                break;
            case 1:
                // logarithmic axis
                bin_width = (log(end) - log(begin)) / n_bins;
                logspace(begin, end, bin_edges.data(), bin_edges.size());
                break;
        }
    }

    void setAxisType(const int axis_type) {
        HistogramAxis::axis_type = axis_type;
    }

    int getNumberOfBins(){
        return n_bins;
    }

    int getBinIdx(T value){
        switch (axis_type){
            case AXIS_LIN:
                // linear
                return calc_bin_idx(begin, bin_width, value);
            case AXIS_LOG10:
                // logarithm
                return calc_bin_idx(begin, bin_width, std::log10(value));
            default:
                return search_bin_idx(value, bin_edges.data(), bin_edges.size());
        }
    }

    T* getBins(){
        return bin_edges.data();
    }

    void getBins(T* bin_edges, int n_bins){
        for(int i = 0; i < n_bins; i++){
            bin_edges[i] = this->bin_edges[i];
        }
    }

    const std::string &getName() const {
        return name;
    }

    void setName(const std::string &name) {
        HistogramAxis::name = name;
    }

    HistogramAxis() = default;
    HistogramAxis(std::string name, T begin, T end, int n_bins, int axis_type = AXIS_LIN)
    {

        // make sure that begin < end
        if(begin > end){
            T temp = begin;
            begin = end;
            end = temp;
        }
        HistogramAxis::begin = begin;
        HistogramAxis::end = end;
        HistogramAxis::n_bins = n_bins;
        HistogramAxis::setAxisType(axis_type);
        HistogramAxis::name = name;
        HistogramAxis::update();
    }
};


template<class T>
class Histogram {

private:

    std::map<size_t , HistogramAxis<T>> axes;

    std::vector<T> histogram; // A 1D array of that contains the histogram
    int number_of_axis;
    int n_total_bins;
    size_t getAxisDimensions(){
        return axes.size();
    }

public:

    void update(T *data, int n_rows_data, int n_cols_data){
        int axis_index;
        int global_bin_idx;
        int global_bin_offset;
        int n_axis;
        int current_bin_idx, current_n_bins;
        HistogramAxis<T> *current_axis;
        T data_value;

        // update the axes
        for(const auto& p : axes){
            axes[p.first].update();
        }

        // initialize a new empty histogram
        // clear the memory of the old histogram
        histogram.clear();
        // 1. count the total number of bins
        n_total_bins = 1;
        n_axis = 0;
        for(const auto& p : axes){
            axis_index = p.first;
            n_axis += 1;
            n_total_bins *= axes[axis_index].getNumberOfBins();
        }
        // 2. fill the histogram with zeros
        histogram.resize(n_total_bins);
        for(global_bin_idx=0; global_bin_idx < n_total_bins; global_bin_idx++){
            histogram[global_bin_idx] = 0.0;
        }

        // Fill the histogram
        // Very instructive for multi-dimensional array indexing
        // https://eli.thegreenplace.net/2015/memory-layout-of-multi-dimensional-arrays/
        bool is_inside;
        for(int i_row=0; i_row<n_rows_data; i_row++){
            // in this loop the position within the 1D array is calculated
            global_bin_offset = 0;
            is_inside = true;
            for(const auto& p : axes){
                axis_index = p.first;
                current_axis = &axes[axis_index];

                data_value = data[i_row*n_axis + axis_index];
                current_bin_idx = current_axis->getBinIdx(data_value);
                current_n_bins = current_axis->getNumberOfBins();

                if( (current_bin_idx < 0) || (current_bin_idx >= current_n_bins) ){
                    is_inside = false;
                    break;
                }
                global_bin_offset = current_bin_idx + current_n_bins * global_bin_offset;
            }
            if(is_inside){
                histogram[global_bin_offset] += 1.0;
            }
        }
    }

    void get_histogram(T** hist, int* dim){
        *hist = histogram;
        *dim = n_total_bins;
    }

    void set_axis(size_t data_column, HistogramAxis<T> &new_axis){
        axes[data_column] = new_axis;
    }

    void set_axis(
            size_t data_column,
            std::string name,
            T begin, T end, int n_bins,
            std::string axis_type
            ){
        HistogramAxis<T> new_axis(name, begin, end, n_bins, axis_type);
        set_axis(data_column, new_axis);
    }

    HistogramAxis<T> get_axis(size_t axis_index){
        return axes[axis_index];
    }

    Histogram() = default;

    ~Histogram() = default;

};


/*!
 * Compute a 1D Histogram for a set of data points
 *
 * @tparam T Type of the values in the histogram
 * @param[in] data array of data points
 * @param[in] n_data number of data points
 * @param[in] weights weight of each data point
 * @param[in] n_weights number of weights
 * @param[in] bin_edges contains the edges of the histogram in ascending order (from small to large)
 * @param[in] n_bins the number of bins in the histogram
 * @param[out] hist array that will contain the histogram as an output
 * @param[out] n_hist number of elements in the output array
 * @param[in] axis_type either a linear, logarithmic or arbitrary axis
 * @param[in] use_weights if true uses weights
 */
template<typename T>
void histogram1D(
        T *data, int n_data,
        double *weights, int n_weights,
        T *bin_edges, int n_bins,
        double *hist, int n_hist,
        const char axis_type,
        bool use_weights
) {
    T v; // stores the data value in iterations
    int i, bin_idx;
    T lower, upper, bin_width;
    bool is_log10 = (axis_type == AXIS_LOG10);
    bool is_lin = (axis_type == AXIS_LIN);

    if (is_lin || is_log10) {
        if (is_log10) {
            lower = std::log10(bin_edges[0]);
            upper = std::log10(bin_edges[n_bins - 1]);
        } else {
            lower = bin_edges[0];
            upper = bin_edges[n_bins - 1];
        }
        bin_width = (upper - lower) / (n_bins - 1);

        for (i = 0; i < n_data; i++) {
            v = data[i];
            if(is_log10){
                if(v == 0){
                    continue;
                } else {
                    v = std::log10(v);
                }
            }
            bin_idx = calc_bin_idx(lower, bin_width, v);
            // ignore values outside bounds
            if ((bin_idx <= n_bins) && (bin_idx >= 0)){
                hist[bin_idx] += (use_weights) ? weights[i] : 1;
            }
        }
    } else {
        for (i = 0; i < n_data; i++) {
            v = data[i];
            bin_idx = search_bin_idx(v, bin_edges, n_bins);
            if(bin_idx > 0)
                hist[bin_idx] += (use_weights) ? weights[i] : 1;
        }
    }

}

void bincount1D(int* data, int n_data, int* bins, int n_bins){
    for(int j=0; j < n_data; j++)
    {
        int v = data[j];
        if( (v >= 0) && (v < n_bins) )
            bins[v]++;
    }
}

IMPBFF_END_NAMESPACE


#endif //IMPBFF_HISTOGRAM_H
