/**
 * Trying to make IO simpler.
 *
 * AsciiData implements a proxy to reading line by line the content of an ASCII
 * file or stdin and allows one to parse value by value into the relevant
 * variables.
 *
 * @author M. Fouesneau
 *
 * example
 * -------
 *
 *  // open data file 
 *  std::string fname = "sample.cat";
 *  AsciiData data(fname);
 *
 *  // select columns to use
 *  data.usecols = {0, 1, 3};
 *
 *  // define data field values and types
 *  double ra;
 *  double dec;
 *  double gmag;
 *
 *  while(data++){
 *      data.next_value() >> ra;
 *      data.next_value() >> dec;
 *      data.next_value() >> gmag;
 *  }
 */
#pragma once

#include <iostream>
#include <iomanip>
#include <istream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <vector>
#include <memory>


/*-----------------------------------------------------------------------------
 *  Ascii Data class definition
 *-----------------------------------------------------------------------------*/
class AsciiData
{
public:
	/* methods */
	AsciiData (std::string fname);
    std::string next_as_string();           /**< get next value as string */
    template <typename T> T next_as();      /**< get next value on the line */
	bool operator++();                      /**< alias to buffer_line */
	bool isEOF();                           /**< check if EOF reached */

	/* data */
	std::string comment = "#";              /**< comment character */
	std::string delim = " ";                /**< comment character */
	std::vector<size_t> usecols = {};       /**< column to use only */

private:
	/* data */
	// std::ifstream is;                       /**< input stream */
    std::shared_ptr<std::istream> is;       /**< pointer to stream */
	std::string fname = "";                 /**< filename of the source */
	bool eof = false;                       /**< eof reached */
	std::string current_line;               /**< current line from is */
	size_t current_column = 0;              /**< current column position */
	/* methods */
	std::string::const_iterator current_i;  /**< current character in line */
	bool buffer_line();                     /**< load a new line in buffer */
};

// vim: expandtab:ts=4:softtabstop=4:shiftwidth=4
