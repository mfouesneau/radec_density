/**
 * This module implements an online 2d histogram.
 *
 * The data can be given from a streaming input of (x, y, weight) tuples.
 * and it generates on the fly the un-normalized density and the corresponding
 * normalization factor.
 *
 */
#ifndef HISTOGRAM_HH
#define HISTOGRAM_HH
#include <stdlib.h>
#include <random>
#include <iostream>
#include <iomanip>
#include <vector>       // std::vector
#include <valarray>     // std::valarray, std::slice
#include <sstream>      // std::stringstream

/** 2D version of the valarray
 *  convenient for instanciation and code clarity
 *
 *  @todo: make it a template
 */
typedef std::valarray<std::valarray<double> > valarray2d;

/** 
 * Make a string representation of an array
 *
 * @param arr 	array to print
 * @return sstr stream instance of the array
 */
template <typename T> std::string toString(std::valarray<T> arr);

/** 
 * Make a string representation of an array
 *
 * @param arr 	array to print
 * @return sstr stream instance of the array
 */
std::string toString(valarray2d arr);

/**
 * Generate `num` evenly sampled values from `start` to `stop`
 * @param start  first value
 * @param stop last value
 * @param num number of points to generate
 * @return x  sequence of values
 */
std::valarray<double> linspace(double start, double stop, size_t num);

/**
 * Generate evenly sampled values from `start` to `stop` with step `dx`
 * @param start  first value
 * @param stop last value
 * @param dx step size 
 * @return x  sequence of values
 */
std::valarray<double> arange(double start, double stop, double dx);


/** Main histogram class */
class Histogram2d
{
public:
    /* constructors */
    Histogram2d (const Histogram2d* h);
    Histogram2d (double xmin, double xmax, double ymin, double ymax, double dx, double dy);

    /* adding data */
    void add_value(double x, double y, double w);
    void add_value(double x, double y);

    /* convenient representation */
    std::string toString();

    /* data */
    double dx;                    //< step along x, pixel width
    double dy;                    //< step along y, pixel width
    double norm;                  //< normalization factor: sum of all pixels
    valarray2d density;           //< array of the density map
    std::valarray<double> x_bins; //< array of bin edges along the x-axis
    std::valarray<double> y_bins; //< array of bin edges along the y-axis
    std::valarray<double> extent; //< [xmin, xmax, ymin, ymax]
    std::valarray<size_t> shape;  //< density shape (Ny, Nx) pixels

private:

    /* sets the internal structure */
    void Setup(double x_min, double x_max, double y_min, double y_max, 
            double dx, double dy);

};
#endif
