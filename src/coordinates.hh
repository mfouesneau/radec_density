/**
 * Tools for coordinate transformation. 
 * 
 * Using this class you can transform one set of coordinates into another
 * system. Available transformations are
 * 
 *  +----------+----------+----------+
 *  |  name    |   from   |    to    |
 *  +----------+----------+----------+
 *  | GAL2ICRS | galactic | ICRS     |
 *  | ICRS2GAL | ICRS     | galactic | 
 *  | ECL2ICRS | ecliptic | ICRS     |
 *  | ICRS2ECL | ICRS     | ecliptic |
 *  | GAL2ECL  | galactic | ecliptic |
 *  | ECL2GAL  | ecliptic | galactic |
 *  +----------+----------+----------+
 *  
 *  Transformations are based on
 *  	Hipparcos Explanatory Vol 1 section 1.5
 *      Murray, 1983, section 10.2
 *  	van Altena et al. 2012, Chapter 4.5 of "Astrometry for Astrophysics"
 *  
 *  Galactic pole coordinates are based on J2000 values.
 *
 *  @author fouesneau
 * 
 */
#pragma once

#include <stdlib.h>
#include <random>
#include <iostream>
#include <iomanip>
#include <map>          // std::map
#include <vector>       // std::vector
#include <valarray>     // std::valarray, std::slice
#include <sstream>      // std::stringstream
#include <algorithm>
#include <string>
#include <cmath>

/** 2D version of the valarray
 *  convenient for instanciation and code clarity
 *
 *  @todo: make it a template
 */
typedef std::valarray<std::valarray<double> > valarray2d;

/**
 * Split a string according to a given delimiter.
 *
 * @param  s          string to split
 * @param  delim      delimiter
 * @param  skip_empty set to discard empty elements
 * @return elements   vector of strings
 */
std::vector<std::string> split_string(std::string &s, char delim, bool skip_empty);


/** 
 * Return a valarray from a vector
 *
 * @param v  vector to convert
 * @return a new array
 */
template <typename T>
std::valarray<T> vector2valarray( std::vector<T> v );


namespace coordinates {
    const double PI = std::atan(1)*4; 

	/**
	 * Convert degrees to radians.
	 *
	 * @param deg  value in degrees
	 * @return rad value in radians
	 */
	double radians(double deg);

	/**
	 * Convert radians to degrees.
     *
	 * @param rad value in radians
	 * @return deg  value in degrees
	 */
	double degrees(double rad);

	/**
	 * Convert arcsec to degrees.
     *
	 * @param angle value in arcseconds
	 * @return deg  value in degrees
	 */
	double arcsec2degrees(double angle);

	/**
	 * Convert arcsec to degrees.
     *
	 * @param angle value in arcminutes
	 * @return deg  value in degrees
	 */
	double arcmin2degrees(double angle);

	/**
	 * Convert arcsec to radians.
     *
	 * @param angle value in arcseconds
	 * @return rad  value in radians
	 */
    double arcsec2radians(double angle);

	/**
	 * Convert arcsec to radians.
     *
	 * @param angle value in arcminutes
	 * @return rad  value in degrees
	 */
    double arcmin2radians(double angle);

	/**
	 * Transform dms angle into degrees.
	 * 
	 * @param dms   string angle
	 * @param delimiter  delimiter between values
	 * @return	deg  converted angle
	 */
	double parseDmsToDegrees(std::string dms, char delimiter);

	/**
	 * Transform hms angle into degrees.
	 * 
	 * @param hms   string angle
	 * @param delimiter  delimiter between values
	 * @return	deg  converted angle
	 */
	double parseHmsToDegrees(std::string hms, char delimiter);


	/** 
     * Obliquity of the Ecliptic (arcsec).
	 */
	const double obliquityOfEcliptic = radians(84381.41100 / 3600.0);

	/** Galactic pole in ICRS coordinates.
	 * (see Hipparcos Explanatory Vol 1 section 1.5, and Murray, 1983, section
	 * 10.2) J2000
	 */
	const double alpha_galactic_pole = radians(192.85948);
	const double delta_galactic_pole = radians(27.12825);

	/** The galactic longitude of the ascending node of the galactic plane on the
	 * equator of ICRS. 
     * (see Hipparcos Explanatory Vol 1 section 1.5, and Murray, 1983, section
     * 10.2)
	 */
	const double omega = radians(32.93192);
	
    
    /** 
     * Vector dot product  z = transpose(x) * y.
     *
     * @param x vector 1
     * @param y vector 2
     * @return z 
     */
    double dot(std::valarray<double> x, std::valarray<double> y);

    /** 
     * Transpose matrix C = transpose(A).
     *
     * @param A  matrix to transpose
     * @return C
     */
    valarray2d transpose(valarray2d A);
	
    /** 
     * Matrix-Matrix multiplication.
     *
     * @param A  matrix 1
     * @param B  matrix 2
     * @return C
     */
    valarray2d dot(valarray2d A, valarray2d B);

    /**
     * Matrix-Vector multiplication y = A * x.
     *
     * @param A matrix
     * @param x vector
     * @return y 
     */
    std::valarray<double> dot(valarray2d A, std::valarray<double> x);

    /**
     * Vector-matrix multiplication y = transpose(x) * A.
     *
     * @param x vector
     * @param A matrix
     * @return y 
     */
    std::valarray<double> dot(std::valarray<double> x, valarray2d A);
 
	/**
	 * Construct the rotation matrix associated with the rotation of given angle
	 * along given x, y, z, axis-vector.
	 * 
	 * @param axis Axis around which to rotate ("x", "y", or "z")
	 * @param angle the rotation angle
	 * @return mat  rotation matrix (3x3)
	 * @throws Exception Wrong axis 
	 */
    valarray2d elementaryRotationMatrix( std::string axis, double angle);

	/** Rotation matrix for the transformation from ICRS to Galactic
	 *  coordinates. See equation (4.25) in chapter 4.5 of "Astrometry for
	 *  Astrophysics", 2012, van Altena et al.
	 */
	const valarray2d _matA = elementaryRotationMatrix("z", PI / 2.0 + alpha_galactic_pole);
	const valarray2d _matB = elementaryRotationMatrix("x", PI / 2.0 - delta_galactic_pole);
	const valarray2d _matC = elementaryRotationMatrix("z", -omega);
	const valarray2d _rotationMatrixIcrsToGalactic = dot(_matC, dot(_matB, _matA));
	
	/** Rotation matrix for the transformation from Galactic to ICRS coordinates. */
	const valarray2d _rotationMatrixGalacticToIcrs = transpose(_rotationMatrixIcrsToGalactic);

	/** Rotation matrix for the transformation from Ecliptic to ICRS coordinates. */
    const valarray2d _rotationMatrixEclipticToIcrs = elementaryRotationMatrix("x", -1 * (obliquityOfEcliptic));

	/** Rotation matrix for the transformation from ICRS to Ecliptic coordinates. */
    const valarray2d _rotationMatrixIcrsToEcliptic = transpose(_rotationMatrixEclipticToIcrs);

	/** Rotation matrix for the transformation from Galactic to Ecliptic coordinates. */
    const valarray2d _rotationMatrixGalacticToEcliptic = dot(_rotationMatrixIcrsToEcliptic, 
                                                             _rotationMatrixGalacticToIcrs); 

	/** Rotation matrix for the transformation from Ecliptic to Galactic coordinates. */
	const valarray2d _rotationMatrixEclipticToGalactic = transpose(_rotationMatrixGalacticToEcliptic);

	/** 
     * Mappings between transformations and matrices.
     *
     * @param name name of the transformation
     * @return mat transformation matrix appropriately
     */
    valarray2d _get_rotationMatrix(std::string name);

	/** 
     * Convert spherical coordinates to Cartesian ones.
	 * Note that the angle coordinates follow the astronomical convention of
	 * using elevation (declination, latitude) rather than its complement
	 * (pi/2-elevation), where the latter is commonly used in the mathematical
	 * treatment of spherical coordinates.
	 *  
	 * @param r  length of input Cartesian vector.
	 * @param phi longitude-like angle (e.g., right ascension, ecliptic longitude) in radians
	 * @param theta latitude-like angle (e.g., declination, ecliptic latitude) in radians
	 * @return (x,y,z) coordinates
	 */
    std::valarray<double> sphericalToCartesian(double r, double phi, double theta);

	/** 
     * Convert Cartesian to spherical coordinates. 
	 * Note that the angle coordinates follow the astronomical convention of
	 * using elevation (declination, latitude) rather than its complement
	 * (pi/2-elevation), which is commonly used in the mathematical treatment
	 * of spherical coordinates.
	 *
	 *  @param x   Cartesian vector component along the X-axis
	 *  @param y   Cartesian vector component along the Y-axis
	 *  @param z   Cartesian vector component along the Z-axis
	 *  @return [r, theta, phi] the vector in spherical coordinates.
	 */
    std::valarray<double> cartesianToSpherical(double x, double y, double z);

	/**
	 * apply transformation coordinates
	 * 
	 * @param name transformation name
	 * @param phi first coordinate in radians (or degrees if use_degrees)
	 * @param theta coordinate in radians (or degrees if use_degrees)
	 * @param use_degrees  if set, takes input angles in degrees and produces outputs also in degrees
	 * @return (a,b)  transformed coordinates in radians (or degrees if use_degrees)
	 */
    std::valarray<double> apply_transformation(std::string name, double phi, 
			double theta, bool use_degrees);

	/** Compute the angular distance between 2 points on the sphere
	 * @param ra1    first longitude  in radians
	 * @param dec1   first latitude   in radians
	 * @param ra2    second latitude  in radians
	 * @param dec2   second longitude in radians
	 * @return distance in radians
	 */
	double spherical_distance_radians(double ra1, double dec1, double ra2, double dec2);

	/** Compute the angular distance between 2 points on the sphere in degrees
	 * @param ra1    first longitude  in degrees
	 * @param dec1   first latitude   in degrees
	 * @param ra2    second latitude  in degrees
	 * @param dec2   second longitude in degrees
	 * @return distance in degrees
	 */
    double spherical_distance_degrees(double ra1, double dec1, double ra2, double dec2);
} // end of namespace coordinates


// vim: expandtab:ts=4:softtabstop=4:shiftwidth=4
