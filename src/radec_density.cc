#include <stdlib.h>
#include <iostream>
#include <vector>       // std::vector
#include <valarray>     // std::valarray, std::slice
#include <sstream>      // std::stringstream
#include <algorithm>
#include <string>

#include "coordinates.hh"  // astro coordinates
#include "mfio.hh"         // ascii input
#include "histogram.hh"    // histogram making
#include "pbar.hh"

#include "mfopts.hh"      // option parser

/**
 * Split a string according to a given delimiter.
 *
 * @param  s          string to split
 * @param  delim      delimiter
 * @param  skip_empty set to discard empty elements
 * @return elements   vector of strings
 */
std::vector<std::string> split(std::string &s, char delim, 
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
 * Counts the number of lines in a file.
 *
 * @param fname        input filename
 * @param buffer_size  the size of the buffer (too large/small is not good.
 *                     16K = 16384 is what GNU wc command uses
 *                     64K is the pipe buffer in Linux systems
 * @return nlines      the number of lines found in stream
 */
size_t count_lines(std::string fname, const int buffer_size = 8 * 1024 * 1024)
{
    std::FILE* stream = std::fopen(fname.c_str(), "r");  //fopen works with char*

    std::vector<char> buffer(buffer_size);
    size_t size;
    size_t line_count = 0;

    while ((size = std::fread(buffer.data(), sizeof(char), buffer_size, stream)) > 0) {
       line_count += std::count_if(buffer.begin(), 
                                   buffer.begin() + size, 
                                   [](char ch) { return ch == '\n'; });
    }
    std::fclose(stream);
    return line_count;
}


/**
 * Make a density map from the input file.
 *
 * @param fname          file to read from
 * @param cols           sequence of columns of RA, DEC, Weights
 * @param dra            ra bin size in degrees
 * @param ddec           dec bin size in degrees
 * @param transform      coordinate transformation (see coordinates)
 * @param comment        string indicating comment line
 * @param delimiter      delimiter between columns
 */
Histogram2d make_density_map(std::string fname,
		std::vector<size_t> cols,
		double dra=2, double ddec=2,
		std::string transform="",
		std::string comment="#",
		std::string delimiter=" ",
        bool radians=false)
{

//    std::cout << "# Checking file size: ";
//    size_t nlines = count_lines(fname);
//    std::cout << fname << " contains " << nlines << " lines" << std::endl;  

    AsciiData data(fname);
    data.usecols = cols;
    data.comment = comment;
    data.delim   = delimiter;

    Histogram2d H(-180, 180, -90, 90, dra, ddec);
    
    // precompute condition
    bool coord_transformation = (transform.size() > 0);

    PBar pb;
    pb.set_description("Histogram");
    pb.start();

    // check if we have ra, dec or ra, dec, weight
    // and work accordingly
    if (cols.size() == 2){
        double ra;
        double dec;
        for(size_t lineno=1; ++data; ++lineno, ++pb){
            ra = data.next_as<double>();
            dec = data.next_as<double>();
            if (coord_transformation){
                std::valarray<double> r = coordinates::apply_transformation(
                        transform, ra, dec, true); // true = use_degrees
                ra = r[0];
                dec = r[1];
            }
            if (radians){
                ra  = coordinates::degrees(ra);
                dec = coordinates::degrees(dec);
            }
            if (ra > 180) ra -= 360;
            H.add_value(ra, dec);
        }
    } else if (cols.size() == 3) {
        double ra;
        double dec;
        double w;
        for(size_t lineno=1; ++data; ++lineno, ++pb){
            ra = data.next_as<double>();
            dec = data.next_as<double>();
            w = data.next_as<double>();
            if (coord_transformation){
                std::valarray<double> r = coordinates::apply_transformation(
                        transform, ra, dec, true); // true = use_degrees
                ra = r[0];
                dec = r[1];
            }
            if (radians){
                ra  = coordinates::degrees(ra);
                dec = coordinates::degrees(dec);
            }
            if (ra > 180) ra -= 360;
            H.add_value(ra, dec, w);
        }
    }
    pb.finish();
    return H;
}

int main(int argc, char *argv[])
{
    // default values
    std::string fname = "-";
    std::string oname = "tmp.hist";
    double dra = 2;
    double ddec = 2;
    std::string transform = "";
    std::string comment = "#";
    std::string delimiter = " ";
    std::vector<size_t> usecols = {0, 1};
    bool radians = false;

    // parse options
    // =============
    
    mfopts::Options  opt(argv[0], "");

    opt.add_option("-h,--help", "Display help message");
    opt.add_option("-i,--input", "input data file");
    opt.add_option("-o,--output", "output histogram to file");
    opt.add_option("-c,--columns" , "column indexes to use (RA, DEC [,Weight]) (starts at 0)");
    opt.add_option("--dra" , "RA bin width in degrees");
    opt.add_option("--ddec" , "RA bin width in degrees");
    opt.add_option("-t,--transform", "coordinate transformation");
    opt.add_option("--comment", "comment line character in input file");
    opt.add_option("-d,--delim", "delimiter character in input file");
    opt.add_option("-r,--radians", "set if input ra,dec are in radians");

    opt.parse_options(argc, argv);

    if (opt.has_option("--help"))
    {
        std::cout << opt.help() << std::endl;
        exit(0);
    }

    if (opt.has_option("--input"))
        fname = opt.get_option<std::string>("--input");
    if (opt.has_option("--output"))
        oname = opt.get_option<std::string>("--output");
    if (opt.has_option("--transform"))
        transform = opt.get_option<std::string>("--transform");
    if (opt.has_option("--comment"))
        comment = opt.get_option<std::string>("--comment");
    if (opt.has_option("--delim"))
        delimiter = opt.get_option<std::string>("--delim");
    if (opt.has_option("--dra"))
        dra = opt.get_option<double>("--dra");
    if (opt.has_option("--ddec"))
        ddec = opt.get_option<double>("--ddec");
    if (opt.has_option("--radians"))
        radians = true;

    // parse columns separated by commas
    if (opt.has_option("--columns")){
        std::string s = opt.get_option<std::string>("--columns");
        std::vector<std::string> cstr = split(s, ',');
        usecols.clear();
        for (size_t i=0; i < cstr.size(); ++i){
            size_t k = std::atoi(cstr[i].c_str());
            usecols.push_back(k);
        }
    }

    Histogram2d H = make_density_map(fname, usecols, dra, ddec, 
                                     transform, comment, delimiter,
                                     radians);

    // output result
    if (oname.size() > 0){
        std::ofstream outputFile( oname );
        outputFile << toString(H.density) << std::endl;
        outputFile.close();
    } else {
        std::cout << toString(H.density) << std::endl;
    }



    std::cout << "# done" << std::endl;
    return 0;
}

// vim: expandtab:ts=4:softtabstop=4:shiftwidth=4

