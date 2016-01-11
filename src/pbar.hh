/**
 * One-line refreshing progress bar inspired by wget that shows ETA (time remaining).
 * 
 *   description [#---------] k/n  10% [time: 01:23:45s, eta: 01d 00:12:34, 2.7 iters/sec]
 *  
 *
 *  .. code:: c
 *
 *
 *		size_t n = 10000;
 *		PBar pb(n);
 *      pb.set_description("short loop");
 *
 *      std::cout << "test short loop" << std::endl;
 *      // neat to also increment the pb during the for declaration!!
 *		for(size_t i=0; i <= n; ++i, ++pb) {
 *			usleep(200);
 *		}
 */

#ifndef PBAR_HH
#define PBAR_HH

#include <iostream>
#include <string>
#include <sstream>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <cassert>
#include <algorithm>
#include <cmath>

/**
 * Main class 
 */
class PBar {

    private:
        typedef std::chrono::system_clock system_clock;
        typedef std::chrono::duration<size_t, system_clock::period> duration;

        size_t n;                    //< total number of iterations (<=0 for unknown)
        size_t cur;                  //< current number of iterations done.
        unsigned char width;         //< How many chars the entire line can be.
        std::string desc="";         //< description
        std::string units = "iters"; //< units of iteration
        size_t current_length = 0;   //< width of the last print out

        std::chrono::system_clock::time_point startTime; //< time of start
        std::chrono::system_clock::time_point lastCheck; //< time of last update

        /*  private methods */
        std::string format_interval(duration t);
        std::string build_str_meter(std::string desc, size_t n, size_t total, duration dt);
        void print_status(std::string s);

    public:
        PBar(size_t _n=0) : n(_n), cur(0), width(10), desc(""), units("iters") {}
        void reset( size_t _n );
        void start();
        void finish();
        void set_description(std::string desc);
        void set_units(std::string units);
        void operator++();
        void operator+=( const size_t d );
};

#endif // PBAR_HH
// vim: expandtab:ts=4:softtabstop=4:shiftwidth=4

