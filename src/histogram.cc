#include "histogram.hh"

/** Make a string representation of an array
 *
 */
template <typename T> std::string toString(std::valarray<T> arr){
    std::stringstream t;
    size_t nx = arr.size();
    for (size_t j=0; j < nx; j++){
        t << std::setw(12) << std::to_string(arr[j]) << "   ";
    }
    t << std::endl;
    return t.str();
}

/** Make a string representation of an array
 *
 */
std::string toString(valarray2d arr){
    std::stringstream t;
    size_t nx = arr.size();
    size_t ny = arr[0].size();
    for (size_t j=0; j < nx; j++){
        for (size_t i=0; i < ny; i++){
            t << std::setw(12) << std::to_string(arr[j][i]) << "   ";
        }
        t << std::endl;
    }
    return t.str();
}


/**
 * Generate `num` evenly sampled values from `start` to `stop`
 * @param start  first value
 * @param stop last value
 * @param num number of points to generate
 * @return x  sequence of values
 */
std::valarray<double> linspace(double start, double stop, size_t num){
    double delta = stop - start;
    double dx = std::ceil(delta / (double)(num - 1));
    std::valarray<double> x(num);
    
    for (size_t i=0; i < num; i++){
        x[i] = start + dx * i;
    }
    return x;
}


/**
 * Generate evenly sampled values from `start` to `stop` with step `dx`
 * @param start  first value
 * @param stop last value
 * @param dx step size 
 * @return x  sequence of values
 */
std::valarray<double> arange(double start, double stop, double dx){
    size_t n = std::ceil((stop - start) / dx);
    std::valarray<double> x(n);
    
    for (size_t i=0; i < n; i++){
        x[i] = start + dx * i;
    }
    return x;
}
	
/**
 * Stream constructor. Initialize an empty histogram ready to be filled
 * @param x_min  minimum x value
 * @param x_max  maximum x value
 * @param y_min  minimum y value
 * @param y_max  maximum y value
 * @param dx     step along x
 * @param dy     step along y
 */
Histogram2d::Histogram2d(double xmin, double xmax, double ymin, double ymax,
        double dx, double dy){
    this->Setup(xmin, xmax, ymin, ymax, dx, dy);
}

/**
 * Stream constructor. Initialize an empty histogram ready to be filled
 * @param h      histogram of reference
 */
Histogram2d::Histogram2d(const Histogram2d* h){
    std::valarray<double> e = h->extent;
		this->Setup(e[0], e[1], e[2], e[3], h->dx, h->dy);
}


/** Prepare the object
 * @param x_min  minimum x value
 * @param x_max  maximum x value
 * @param y_min  minimum y value
 * @param y_max  maximum y value
 * @param dx     step along x
 * @param dy     step along y
 */
void Histogram2d::Setup(double x_min, double x_max, double y_min, double y_max,
        double dx, double dy){

    this->dx = dx;
    this->dy = dy;
    this->x_bins = arange(x_min, x_max + 0.1 * dx, dx);
    this->y_bins = arange(y_min, y_max + 0.1 * dy, dy);
    this->density = valarray2d(std::valarray<double>(this->x_bins.size() - 1), this->y_bins.size() -1);
    this->shape = std::valarray<size_t>({x_bins.size() - 1, y_bins.size() - 1});
    this->extent = std::valarray<double>({x_min, x_max, y_min, y_max});
    this->norm = 0.;

    // put histogram at 0
    for (size_t i=0; i < this->shape[0]; i++){
        for (size_t j=0; j < this->shape[1]; j++){
            this->density[j][i] = 0.;
        }
    }

}

/** Update the histogram with a new value
 * @param x  value along x
 * @param y  value along y
 * @param w  weight
 */
void Histogram2d::add_value(double x, double y, double w){
    if ((x <= this->extent[1]) && (x >= this->extent[0]) &&
        (y <= this->extent[3]) && (y >= this->extent[2])){

        size_t xbin = std::ceil((x - this->extent[0]) / this->dx) - 1; // starts at 0
        size_t ybin = std::ceil((y - this->extent[2]) / this->dy) - 1;
//        std::cout << " x, y, xb, yb "
//            << std::setw(12) << x
//            << std::setw(12) << y
//            << std::setw(12) << xbin
//            << std::setw(12) << ybin << std::endl;
        this->density[ybin][xbin] += w;
        this->norm += w;
    }
}

/** Update the histogram with a new value
 * @param x  value along x
 * @param y  value along y
 */
void Histogram2d::add_value(double x, double y){
    this->add_value(x,  y,  1.);
}


/** Make a string representation of the histogram
 *
 */
std::string Histogram2d::toString(){
    std::stringstream t;
    for (size_t j=0; j < this->shape[1]; j++){
        t << std::setw(12) << 0.5 * (this->y_bins[j+1] + this->y_bins[j]) << " | ";
        for (size_t i=0; i < this->shape[0]; i++){
            t << std::setw(12) << std::to_string(this->density[j][i]) << "   ";
        }
        t << std::endl;
    }
    t << "            " << "   ";
    for (size_t i=0; i < this->shape[0]; i++){
        t << "------------" << "   ";
    }
    t << std::endl;
    t << "            " << "   ";
    for (size_t i=0; i < this->shape[0]; i++){
        t << std::setw(12) << 0.5 * (this->x_bins[i+1] + this->x_bins[i]) << "   ";
    }
    t << std::endl;
    return t.str();
}


///** testing histogram code with a Gaussian */
//int main(int argc, char *argv[])
//{
//	
//    // xmin, xmax, ymin, ymax, dx, dy
//    Histogram2d H(0., 2., 0., 2, 0.01, 0.01);
//
//	// define random generator
//    std::random_device rd;
//    std::mt19937 generator(rd());
//    std::normal_distribution<double> randn(0, 1);
//
//
//    size_t N = 400000; // number of data points
//    double x = 0.;
//    double y = 0.;
//    for (int i=0 ; i < N; i++){
//        x = 0.1 * randn(generator) + 0.5;
//        y = 0.1 * randn(generator) + 1;
//        H.add_value(x, y);
//    }
//    
//    // print only the image
//    std::cout << toString(H.density) << std::endl;
//    // better repr for stdout readings
//    // std::cout << H.toString() << std::endl;
//    return 0;
//}

// vim: expandtab:ts=4:softtabstop=4:shiftwidth=4 
