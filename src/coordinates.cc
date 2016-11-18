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
#include "coordinates.hh"

/**
 * Split a string according to a given delimiter.
 *
 * @param  s          string to split
 * @param  delim      delimiter
 * @param  skip_empty set to discard empty elements
 * @return elements   vector of strings
 */
std::vector<std::string> split_string(std::string &s, char delim, 
        bool skip_empty=true) {
    std::stringstream ss(s);
    std::string item;
    std::vector<std::string> elements;

    while (std::getline(ss, item, delim)) {
        if (item.size() > 0) {
            elements.push_back(item);
        } else if (not skip_empty) {
            elements.push_back(item);
        }
    }
    return elements;
}

/** 
 * Return a valarray from a vector
 *
 * @param v  vector to convert
 * @return a new array
 */
template <typename T>
std::valarray<T> vector2valarray( std::vector<T> v ){
    std::valarray<T> a(&(v[0]), v.size());
    return a;
}


namespace coordinates {

	/**
	 * Convert degrees to radians.
     *
	 * @param deg  value in degrees
	 * @return rad value in radians
	 */
	double radians(double deg){
		return PI / 180. * deg;
	}

	/**
	 * Convert radians to degrees.
     *
	 * @param rad value in radians
	 * @return deg  value in degrees
	 */
	double degrees(double rad){
		return 180. / PI * rad;
	}

	/**
	 * Convert arcsec to degrees.
     *
	 * @param angle value in arcseconds
	 * @return deg  value in degrees
	 */
	double arcsec2degrees(double angle){
		return angle / 3600.;
	}

	/**
	 * Convert arcsec to degrees.
     *
	 * @param angle value in arcminutes
	 * @return deg  value in degrees
	 */
	double arcmin2degrees(double angle){
		return angle / 60.;
	}

	/**
	 * Convert arcsec to radians.
     *
	 * @param angle value in arcseconds
	 * @return rad  value in radians
	 */
    double arcsec2radians(double angle){
		return radians(arcsec2degrees(angle));
	}

	/**
	 * Convert arcsec to radians.
     *
	 * @param angle value in arcminutes
	 * @return rad  value in degrees
	 */
    double arcmin2radians(double angle){
		return radians(arcmin2degrees(angle));
	}

	/**
	 * Transform dms angle into degrees.
	 * 
	 * @param dms   string angle
	 * @param delimiter  delimiter between values
	 * @return	deg  converted angle
	 * @throws Exception  if the format is not recognized
	 */
	double parseDmsToDegrees(std::string dms, char delimiter=':') {

        std::vector<std::string> elements = split_string(dms, delimiter); 
		if (elements.size() > 3) throw std::length_error("Wrong input format");

		double deg = std::stod(elements[0])
                     + (std::stod(elements[1])
                        + std::stod(elements[2]) / 60. ) / 60. ;
		return deg;
	}

	/**
	 * Transform hms angle into degrees.
	 * 
	 * @param hms   string angle
	 * @param delimiter  delimiter between values
	 * @return	deg  converted angle
	 * @throws Exception  if the format is not recognized
	 */
	double parseHmsToDegrees(std::string hms, char delimiter = ':') {
		return parseDmsToDegrees(hms, delimiter) * 15.;   // 15 = 360 / 24.
	}


    /** 
     * Vector dot product  z = transpose(x) * y.
     *
     * @param x vector 1
     * @param y vector 2
     * @return z 
     */
    double dot(std::valarray<double> x, std::valarray<double> y) {
        if (x.size() != y.size()) 
            throw std::length_error("Illegal vector dimensions.");

        double sum = 0.0;
        for (size_t i = 0; i < x.size(); ++i)
            sum += x[i] * y[i];
        return sum;
    }

    /** 
     * Transpose matrix C = transpose(A).
     *
     * @param A  matrix to transpose
     * @return C
     */
    valarray2d transpose(valarray2d A) {
        size_t m = A.size();
        size_t n = A[0].size();
        valarray2d C = valarray2d(std::valarray<double>(m), n);
        for (size_t i = 0; i < m; ++i)
            for (size_t j = 0; j < n; ++j)
                C[j][i] = A[i][j];
        return C;
    }
	
    /** 
     * Matrix-Matrix multiplication.
     *
     * @param A  matrix 1
     * @param B  matrix 2
     * @return C
     */
    valarray2d dot(valarray2d A, valarray2d B) {
        size_t mA = A.size();
        size_t nA = A[0].size();
        size_t mB = B.size();
        size_t nB = B[0].size();
        if (nA != mB) throw std::length_error("Illegal matrix dimensions.");
        valarray2d C = valarray2d(std::valarray<double>(mA), nB);
        for (size_t i = 0; i < mA; i++)
            for (size_t j = 0; j < nB; j++)
                for (size_t k = 0; k < nA; k++)
                    C[i][j] += A[i][k] * B[k][j];
        return C;
    }

    /**
     * Matrix-Vector multiplication y = A * x.
     *
     * @param A matrix
     * @param x vector
     * @return y 
     */
    std::valarray<double> dot(valarray2d A, std::valarray<double> x) {
        size_t m = A.size();
        size_t n = A[0].size();
        if (x.size() != n) {
        	throw std::length_error("Illegal matrix dimensions.");
        }
        std::valarray<double> y(m);
        for (size_t i = 0; i < m; i++)
            for (size_t j = 0; j < n; j++)
                y[i] += A[i][j] * x[j];
        return y;
    }


    /**
     * Vector-matrix multiplication y = transpose(x) * A.
     *
     * @param x vector
     * @param A matrix
     * @return y 
     */
    std::valarray<double> dot(std::valarray<double> x, valarray2d A) {
        size_t m = A.size();
        size_t n = A[0].size();
        if (x.size() != m) {
        	throw std::length_error("Illegal matrix dimensions.");
        }
        std::valarray<double> y(n);
        for (size_t j = 0; j < n; j++)
            for (size_t i = 0; i < m; i++)
                y[j] += A[i][j] * x[i];
        return y;
    }
 
	/**
	 * Construct the rotation matrix associated with the rotation of given angle
	 * along given x, y, z, axis-vector.
	 * 
	 * @param axis Axis around which to rotate ("x", "y", or "z")
	 * @param angle the rotation angle
	 * @return mat  rotation matrix (3x3)
	 * @throws Exception Wrong axis 
	 */
    valarray2d elementaryRotationMatrix( std::string axis, double angle)
    {
        std::vector<double> r1, r2, r3;

        valarray2d mat = valarray2d(std::valarray<double>(3), 3);

		if ((axis == "x") or (axis == "X")) {
            mat[0][0] = 1.;
            mat[0][1] = 0;
            mat[0][2] = 0;
            mat[1][0] = 0;
            mat[1][1] = cos(angle);
            mat[1][2] = sin(angle);
            mat[2][0] = 0;
            mat[2][1] = -sin(angle);
            mat[2][2] = cos(angle);
		} else if ((axis == "y") or (axis == "Y")) {
            mat[0][0] = cos(angle);
            mat[0][1] = 1;
            mat[0][2] = -sin(angle);
            mat[1][0] = 0.;
            mat[1][1] = 1;
            mat[1][2] = 0;
            mat[2][0] = sin(angle);
            mat[2][1] = 0.;
            mat[2][2] = cos(angle);
		} else if ((axis == "z") or (axis == "Z")) {
            mat[0][0] = cos(angle);
            mat[0][1] = sin(angle);
            mat[0][2] = 0.;
            mat[1][0] = -sin(angle);
            mat[1][1] = cos(angle);
            mat[1][2] = 0;
            mat[2][0] = 0.;
            mat[2][1] = 0.;
            mat[2][2] = 1.;
		}
    
		return mat;
	}


	/** 
     * Mappings between transformations and matrices.
     *
     * @param name name of the transformation
     * @return mat transformation matrix appropriately
     */
    valarray2d _get_rotationMatrix(std::string name){
        std::string NAME = name;
        std::transform( NAME.begin(), NAME.end(), NAME.begin(), toupper);
        if (NAME == "GAL2ICRS")
            return _rotationMatrixGalacticToIcrs;
        else if (NAME == "ICRS2GAL")
            return _rotationMatrixIcrsToGalactic;
        else if (NAME == "ECL2ICRS")
            return _rotationMatrixEclipticToIcrs;
        else if (NAME == "ICRS2ECL")
            return _rotationMatrixIcrsToEcliptic;
        else if (NAME == "GAL2ECL")
            return _rotationMatrixGalacticToEcliptic;
        else if (NAME == "ECL2GAL")
            return _rotationMatrixEclipticToGalactic;
        else
            throw std::runtime_error("Connot find this transformation name");
    }

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
    std::valarray<double> sphericalToCartesian(double r, double phi, double theta){
	    double ctheta = cos(theta);
	    double x = r * cos(phi) * ctheta;
	    double y = r * sin(phi) * ctheta;
	    double z = r * sin(theta);
        
        std::valarray<double> v(3);
        v[0] = x;
        v[1] = y;
        v[2] = z;
        return v;
	}


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
    std::valarray<double> cartesianToSpherical(double x, double y, double z) {
	    double rCylSq = x * x + y * y;
		double r = sqrt(rCylSq + z * z);
		if (r == 0.){ throw std::runtime_error("Error: point is at distance zero."); }
        std::valarray<double> v(3);
        v[0] = r;
        v[1] = atan2(y, x); 
        v[2] = atan2(z, sqrt(rCylSq));
        return v;
	}


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
			double theta, bool use_degrees=true) {

		double r = 1.;
        std::valarray<double> xyz;
		if (use_degrees){
			xyz = sphericalToCartesian(r, radians(phi), radians(theta));
		} else {
			xyz = sphericalToCartesian(r, phi, theta);
		}
        valarray2d mat = _get_rotationMatrix(name);
        
        std::valarray<double> xyzrot = dot(mat, xyz);
        std::valarray<double> rab = cartesianToSpherical(xyzrot[0], xyzrot[1], xyzrot[2]);
        std::valarray<double> m(2);
	    if (use_degrees){
            m[0] = degrees(rab[1]);
            m[1] = degrees(rab[2]);
	    } else {
            m[0] = rab[1];
            m[1] = rab[2];
	    }
        return m;
	}


	/** Compute the angular distance between 2 points on the sphere
	 * @param ra1    first longitude  in radians
	 * @param dec1   first latitude   in radians
	 * @param ra2    second latitude  in radians
	 * @param dec2   second longitude in radians
	 * @return distance in radians
	 */
	double spherical_distance_radians(double ra1, double dec1, double ra2, double dec2){
		return 2 * asin(
				sqrt( 
						sin( (dec1 - dec2) / 2)  * sin( (dec1 - dec2) / 2) +
						cos(dec1) * cos(dec2) * ( 
								sin( (ra1 - ra2) / 2) * sin( (ra1 - ra2) / 2) )
						));
	}


	/** Compute the angular distance between 2 points on the sphere in degrees
	 * @param ra1    first longitude  in degrees
	 * @param dec1   first latitude   in degrees
	 * @param ra2    second latitude  in degrees
	 * @param dec2   second longitude in degrees
	 * @return distance in degrees
	 */
    double spherical_distance_degrees(double ra1, double dec1, double ra2, double dec2){
		return degrees(spherical_distance_radians(
				radians(ra1), radians(dec1), radians(ra2), radians(dec2))
				);
	
	}

	/**
	 * Return the matrix of precession between two epochs (IAU 1976, FK5).
	 *
	 * Though the matrix method itself is rigorous, the precession
	 * angles are expressed through canonical polynomials which are
	 * valid only for a limited time span.  There are also known
	 * errors in the IAU precession rate.  
	 * The absolute accuracy of the present formulation is 
	 * better than 0.1 arcsec from 1960AD to 2040AD, 
	 * better than 1 arcsec from 1640AD to 2360AD,
	 * and remains below 3 arcsec for the whole of the period 500BC to 3000AD. 
	 * 
	 * The errors exceed 10 arcsec outside the range 1200BC to 3900AD, exceed
	 * 100 arcsec outside 4200BC to 5600AD and exceed 1000 arcsec outside 6800BC
	 * to 8200AD. 
 	 *
	 * References:
	 *  Lieske,J.H., 1979. Astron.Astrophys.,73,282. equations (6) & (7), p283.
	 *  Kaplan,G.H., 1981. USNO circular no. 163, pA2.
	 *  
	 * @param begEpoch beginning Julian epoch (e.g. 2000 for J2000)
	 * @param endEpoch ending Julian epoch
	 * @return 
	 * @return the precession matrix as a 3x3 matrix, 
	 *              where pos(endEpoch) = rotMat * pos(begEpoch)
	 */
    valarray2d get_FK5PrecessMatrix(double begEpoch, double endEpoch){

	    // Interval between basic epoch J2000.0 and beginning epoch (JC)
	    double t0 = (begEpoch - 2000.0) / 100.0;                  // origin J2000

	    // Interval over which precession required (JC)
	    double dt = (endEpoch - begEpoch) / 100.0;               // centuries

	    // Euler angles
	    double RadPerDeg  = PI / 180.0;                          // radians per degree
	    double ArcSecPerDeg = 60.0 * 60.0;                       // arcseconds per degree
	    double RadPerArcSec = RadPerDeg / ArcSecPerDeg;          // radians per arcsec
	    double tas2r = dt * RadPerArcSec;
	    double w = 2306.2181 + (1.39656 - 0.000139 * t0) * t0;   // arcsec

	    double zeta = (w + ((0.30188 - 0.000344 * t0) + 0.017998 * dt) * dt) * tas2r;
	    double z = (w + ((1.09468 + 0.000066 * t0) + 0.018203 * dt) * dt) * tas2r;
	    double theta = ((2004.3109 + ( - 0.85330 - 0.000217 * t0) * t0)
	                    + (( - 0.42665 - 0.000217 * t0) - 0.041833 * dt) * dt) * tas2r;
	    
	    // Rotation matrix
	    valarray2d mat = dot(dot(elementaryRotationMatrix("z", z), 
	                             elementaryRotationMatrix("y", theta)), 
	                         elementaryRotationMatrix("z", zeta));
	    return mat;
	}

    /**
	 * Apply precession and proper motion of stars.
	 * 
	 *  References:
	 *  "The Astronomical Almanac" for 1987, page B39
	 *  P.T. Wallace's prec routine
	 *  
	 * @param pos initial mean FK5 cartesian position (au)
	 * @param pm  initial mean FK5 cartesian velocity (au per Jul. year) 
	 * @param fromDate date of initial coordinates (Julian epoch)
	 * @param toDate   date to which to precess (Julian epoch)
	 * @return ( final mean FK5 cartesian position (au)
	 *           final mean FK5 cartesian velocity (au per Julian year)
	 *           )
	 */
	valarray2d fk5xyz_applyPrecession(std::valarray<double>& pos, std::valarray<double>& pm, 
            double fromDate, double toDate){
        
	    // compute new precession constants
	    valarray2d rotMat = get_FK5PrecessMatrix(fromDate, toDate);
	    // correct for velocity (proper motion and radial velocity)
        std::valarray<double> tempP(3);

	    double dt = toDate - fromDate;
	    for (size_t k=0; k < 3; ++k){
	    	tempP[k] = pos[k] + pm[k] * dt;
	    }

	    // precess positions and velocity
        std::valarray<double> newpos = dot(rotMat, tempP);
        std::valarray<double> newpm = dot(rotMat, pm);
        valarray2d res = {newpos, newpm};
	    return res;
	}

    /**
	 * Apply precession and proper motion of stars
	 * 
	 *  References:
	 *  "The Astronomical Almanac" for 1987, page B39
	 *  P.T. Wallace's prec routine
	 * 
	 * @param phi first coordinate in radians (or degrees if use_degrees)
	 * @param theta coordinate in radians (or degrees if use_degrees)
	 * @param use_degrees  if set, takes input angles in degrees and produces outputs also in degrees
	 * @return (a,b)  transformed coordinates in radians (or degrees if use_degrees)
	 * @throws Exception Something went wrong in the spherical transformation
	 */
    valarray3d apply_precession(std::valarray<double>& phi, std::valarray<double>& theta, 
            std::valarray<double>& muphi, std::valarray<double>& mutheta, 
			double fromDate, double toDate,
			bool use_degrees){

        valarray2d newpos = valarray2d(std::valarray<double>(2), phi.size());
        valarray2d newvel = valarray2d(std::valarray<double>(2), phi.size());
		double r = 1.;
        std::valarray<double> xyz;
		std::valarray<double> vxyz;
		std::valarray<double> xyzrot;
		std::valarray<double> vxyzrot;
		std::valarray<double> rab;
		std::valarray<double> vab;

		// compute new precession constants
		valarray2d rotMat = get_FK5PrecessMatrix(fromDate, toDate);
		double dt = toDate - fromDate;
        std::valarray<double> tempP(3);

		for (size_t k=0; k<phi.size(); ++k){
			if (use_degrees){
				xyz = sphericalToCartesian(r, radians(phi[k]), radians(theta[k]));
				vxyz = sphericalToCartesian(r, radians(muphi[k]), radians(mutheta[k]));
			} else {
				xyz = sphericalToCartesian(r, phi[k], theta[k]);
				vxyz = sphericalToCartesian(r, muphi[k], mutheta[k]);
			}
			// correct for velocity (proper motion and radial velocity)
			for (size_t j=0; j < 3; ++j){
				tempP[j] = xyz[j] + vxyz[j] * dt;
			}
			// precess position and velocity
			xyzrot = dot(rotMat, tempP);
			vxyzrot = dot(rotMat, vxyz);

			rab = cartesianToSpherical(xyzrot[0], xyzrot[1], xyzrot[2]);
			vab = cartesianToSpherical(vxyzrot[0], vxyzrot[1], vxyzrot[2]);
			if (use_degrees){
				newpos[k][0] = degrees(rab[1]);
				newpos[k][1] = degrees(rab[2]);
				newvel[k][0] = degrees(vab[1]);
				newvel[k][1] = degrees(vab[2]);
			} else {
				newpos[k][0] = rab[1];
				newpos[k][1] = rab[2];
				newvel[k][1] = vab[1];
				newvel[k][1] = vab[2];
			}
		}
		valarray3d res = {newpos, newvel};
		return res;
	}
} // end of namespace coordinates


/*-----------------------------------------------------------------------------
 *  TESTING
 *-----------------------------------------------------------------------------*/
//int main(int argc, char *argv[])
//{
//    
//    for (size_t i=0; i<360; ++i){
//
//        double k = i - 180.;
//        std::valarray<double> r1 = coordinates::apply_transformation("ICRS2GAL", k, 0., true);
//        std::valarray<double> r2 = coordinates::apply_transformation("GAL2ICRS", r1[0], r1[1], true);
//
//        std::cout << "GAL(" << k << ", 0)"
//                  << "-> IRCS(" << r1[0] << ", " << r1[1] << ") "
//                  << "->  GAL(" << r2[0] << ", " << r2[1] << ") "
//                  << std::endl;
//    }
//    return 0;
//}

// vim: expandtab:ts=4:softtabstop=4:shiftwidth=4
