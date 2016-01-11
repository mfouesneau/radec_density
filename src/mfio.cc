/**
 * Trying to make IO simpler.
 *
 * AsciiData implements a proxy to reading line by line the content of an ASCII
 * file or stdin and allows one to parse value by value into the relevant variables
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
 *  while(data++){
 *		double ra   = data.next_as<double>();
 *		double dec  = data.next_as<double>();
 *		double gmag = data.next_as<double>();
 *  }
 */
#include "mfio.hh"

/** 
 * Constructor
 *
 * @param fname filename
 */
AsciiData::AsciiData(std::string fname){
    if ((fname.size() == 0) or (fname.compare("-") == 0)){
        // use stdin
        this->is.reset(&std::cin, [](...){}); // do nothing for delete
        this->fname = "stdin";
    } else {
        this->fname = fname;
        this->is.reset(new std::ifstream(fname));
        if (!*this->is){
             throw std::runtime_error("\nERROR: Cannot open file " + fname);
        }
    }
}

/** 
 * Check if EOF reached 
 *
 * @return eof true if EOF has been reached
 */
bool AsciiData::isEOF(){
	return this->eof;
}

/**
 * Load a line into the buffer
 *
 * @return eof true if EOF has been reached
 */
bool AsciiData::buffer_line(){
	this->current_line.clear();
	this->current_column = 0;
	if (this->eof) {return false;}
	// EOF empty entry if nothing to read
	if ((not getline(*this->is, this->current_line))){
		this->eof = true;
	}
	// empty line or comment line, but not EOF -> go to next line
	if ((this->current_line.size() == 0) 
	     or (this->current_line.substr(0, this->comment.size()) == this->comment)) {
		this->buffer_line();
	}
	this->current_i = this->current_line.begin();
	return (not this->eof);
}

/**
 * alias to buffer_line 
 */
bool AsciiData::operator++(){
	return this->buffer_line();
}

/** 
 * Read the next value into a string.
 *
 * @return value
 */
std::string AsciiData::next_as_string(){
	if (this->eof) {
		return std::string("");
	} else if (this->current_line.size() == 0) {
		this->buffer_line();
	}

	std::string line = this->current_line;

	// now we have data
	std::string strnum;
	strnum.clear();
	// scan all characters
	for (std::string::const_iterator i = this->current_i;
			i != line.end(); ++i){

		// if i is not a delimiter, then append it to strnum
		if (this->delim.find(*i) == std::string::npos){ 
			strnum += *i; 
			// if this is not the last character of the line
			// continue gathering
			if (i + 1 != line.end()) { continue; }
		}

		// if strnum is empty, it means the previous character is also
		// a delimiter (multiple delimiters together). Ignore
		if (strnum.empty()) { continue; }

        // now we only keep it if the column is considered in usecols
        bool found = false;
        if (this->usecols.size() > 0){
            for (size_t key=0; key < this->usecols.size(); ++key){
                if (this->usecols[key] == this->current_column)
                    found = true;
            }
        }
		if (not found){ 
			this->current_column ++;
			strnum.clear();
			continue; 
		}

		// if we reach this point, we have some data to take care of
		this->current_column ++;
		this->current_i = i;
		break;
	}
    return strnum;
}

/**
 * <template>
 * return the next value as given type.
 *
 * @param T type to use
 * @return val the value
 */
template <typename T> T AsciiData::next_as(){
    T val;
	std::istringstream(this->next_as_string()) >> val;
    return val;
}

/* specialization definitions */
template <> double AsciiData::next_as<double>(){
    double val;
	std::istringstream(this->next_as_string()) >> val;
    return val;
}
template <> float AsciiData::next_as<float>(){
    float val;
	std::istringstream(this->next_as_string()) >> val;
    return val;
}
template <> size_t AsciiData::next_as<size_t>(){
    size_t val;
	std::istringstream(this->next_as_string()) >> val;
    return val;
}
template <> int AsciiData::next_as<int>(){
    int val;
	std::istringstream(this->next_as_string()) >> val;
    return val;
}
template <> long AsciiData::next_as<long>(){
    long val;
	std::istringstream(this->next_as_string()) >> val;
    return val;
}
template <> std::string AsciiData::next_as<std::string>(){
    return this->next_as_string();
}
template <> bool AsciiData::next_as<bool>(){
    bool val;
	std::istringstream(this->next_as_string()) >> val;
    return val;
}

/* explicit instanciation of the templates for double, float, ... */
template double AsciiData::next_as<double>();
template float AsciiData::next_as<float>();
template size_t AsciiData::next_as<size_t>();
template int AsciiData::next_as<int>();
template long AsciiData::next_as<long>();
template std::string AsciiData::next_as<std::string>();
template bool AsciiData::next_as<bool>();



/*-----------------------------------------------------------------------------
 *  Testing 
 *-----------------------------------------------------------------------------*/
//int main(int argc, char *argv[])
//{
//	std::string fname = "sample.cat";
//
//	AsciiData data(fname);
//	data.usecols = {0, 1, 3};
//
//	// define data field values
//	double ra;
//	double dec;
//	// std::string id;
//	double gmag;
//
//	// print some header
//	std::cout << std::string("").append( 6 + 3 * (12 + 1), '=') << std::endl;
//	std::cout << " Filename: " << fname << std::endl;
//	std::cout << std::string("").append( 6 + 3 * (12 + 1), '=') << std::endl;
//	std::cout << std::setw(6) << "lineno" << " "
//		  << std::setw(12) << "RA" << " "
//	          << std::setw(12) << "DEC" << " "
//	          << std::setw(12) << "Gmag"
//		  << std::endl;
//	std::cout << std::string("").append( 6 + 3 * (12 + 1), '-') << std::endl;
//
//	// print content with lineno
//	// a while loop would work too
//	// the for below makes use that data++ or ++data
//	// returns if EOF is reach as well as read the next line in the datafile
//	for(size_t lineno=1; ++data; ++lineno){
//		ra  = data.next_as<double>();
//		dec = data.next_as<double>();
//		// id = data.next_as<std::string>();
//		gmag = data.next_as<double>();
//
//		std::cout << std::setw(6)  << lineno << " "
//			  << std::setw(12) << ra << " " 
//		          << std::setw(12) << dec << " "
//		          << std::setw(12) << gmag << " "
//		          // << id << " "
//			  << std::endl;
//	}
//	std::cout << std::string("").append( 6 + 3 * (12 + 1), '=') << std::endl;
//
//	return 0;
//}
// vim: expandtab:ts=4:softtabstop=4:shiftwidth=4
