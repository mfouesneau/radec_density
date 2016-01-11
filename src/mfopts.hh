/**
 * A dead simple option parser.
 *
 *
 * This is a very simplified parser for command line options
 * It allows to declare options and access them. In addition, it generates the
 * usual help message.
 *
 * What this does not do:
 *
 *    * it does not check that an option used was not declare.
 *    * it does not filter the command line to leave only residuals.
 *
 * @author M. Fouesneau
 *
 * example::
 *
 *   int main(int argc, char * argv[])
 *   {
 *       Options opt(argv[0]);
 *       opt.add_option("-h,--help", "display help message");
 *       opt.add_option("-i,--input", "input file");
 *
 *       opt.parse_options(argc, argv);
 *
 *       if (opt.has_option("-h")){
 *           std::cout<< opt.help();
 *           exit(0);
 *       }
 *
 *       if (opt.has_option("--input")){
 *           std::cout << opt.get_option<std::string>("--input") << std::endl;
 *           std::cout << "input = " << opt.get_option<double>("--input") << std::endl;
 *           std::cout << "input = " << opt.get_option<float>("--input") << std::endl;
 *       }
 *
 *       return 0;
 *  }
 *
 */
#pragma once

#include <stdlib.h>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <sstream>

namespace mfopts{

/**
 * Split a string according to a given delimiter.
 *
 * @param  s          string to split
 * @param  delim      delimiter
 * @param  skip_empty set to discard empty elements
 * @return elements   vector of strings
 */
std::vector<std::string> split(const std::string &s, char delim, 
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

/*-----------------------------------------------------------------------------
 * Class Option
 *
 * It collects all aliases and the corresponding value of one option.
 *-----------------------------------------------------------------------------*/
class Option
{
public:
    Option (const std::string& code, const std::string& helpmsg);
    Option (const std::string& code, const std::string& helpmsg, const std::string& def);
    std::string has_option(char** begin, char** end);
    std::string has_option(int argc, char** begin);
    bool has_short();
    void set_value(const std::string& value);
    std::string get_value();
    bool name_is(const std::string& name);
    std::vector<std::string> codes;   /**< option name and aliases             */
    std::string helpmsg;              /**< description / help message          */

private:
    std::string value = "";           /**< value of the option                 */
    bool _short = false;               /**< true if has short option (e.g., -h) */
};


/**
 * Constructor
 *
 * @param code      option codes separated by commas. (do no use spaces)
 * @param helpmsg   help message description
 */
Option::Option(const std::string& code, const std::string& helpmsg){
    this->helpmsg = helpmsg;
    this->codes   = split(code, ',');
    for( size_t i=0; i < this->codes.size(); ++i){
        if ((codes[i][0] == '-') && (codes[i][1] != '-')){
            this->_short = true;
        }
    }
}

/**
 * Constructor with default value
 *
 * @param code      option codes separated by commas. (do no use spaces)
 * @param helpmsg   help message description
 * @param def       default value
 */
Option::Option(const std::string& code, const std::string& helpmsg, 
        const std::string& def){
    this->helpmsg = helpmsg;
    this->codes   = split(code, ',');
    this->set_value(def);
    for( size_t i=0; i < this->codes.size(); ++i){
        if ((codes[i][0] == '-') && (codes[i][1] != '-')){
            this->_short = true;
        }
    }
}

/**
 * return true is there is a short option in the list
 */
bool Option::has_short(){
    return this->_short;
}

/**
 * check if the option exists in the command line
 *
 * @param begin  argv
 * @param end    argv + argc
 * @return name  the name under-which the option was used
 */
std::string Option::has_option(char** begin, char** end){
    for (size_t i = 0; i < this->codes.size(); ++i) {
        if (std::find(begin, end, this->codes[i]) != end)
            return codes[i];
    }
    return "";
}

/**
 * check if the option exists in the command line.
 *
 * @param begin  argv
 * @param argc   the size of argv
 * @return name  the name under-which the option was used
 */
std::string Option::has_option(int argc, char** begin){
    char** end = begin + argc; 
    return this->has_option(begin, end);
}


/**
 * Set the value of the option
 *
 * @param value     value to set
 */
void Option::set_value(const std::string& value){
    this->value = value;
}

/**
 * Get the value of the option
 *
 * @return value     value to set
 */
std::string Option::get_value(){
    return this->value;
}

/** 
 * return the current option has a given name in its aliases.
 *
 * @param name    name to check
 */
bool Option::name_is(const std::string& name){
    for (size_t i = 0; i < this->codes.size(); ++i) {
        if (codes[i] == name)
            return true;
    }
    return false;
}


/*-----------------------------------------------------------------------------
 * Class Options is the manager of Option objects
 *-----------------------------------------------------------------------------*/
class Options
{
public:
    Options (const std::string& cmdname, const std::string& args);
    void parse_options(int _argc, char** _argv);
    bool has_option(const std::string& option);
    std::string find_option(const std::string& option);
    template <typename T> T get_option(const std::string& option);
    void add_option(std::string code, std::string helpmsg);
    std::string help();

private:
    std::vector<Option> _opts;         /**< declared options */
    std::string cmdname = "<Command>"; /**< program/code name */
    std::string args = "<args>";       /**< string for arguments */
    char** _argv;                      /**< argv content when parsing */
    int   _argc;                       /**< argv length when parsing */
};

/** 
 * Constructor.
 *
 * <cmdname> [OPTION...] <args>
 *
 * @param cmdname   command name
 * @param args      argument string
 */
Options::Options (const std::string& cmdname, const std::string& args="<args>"){
    this->cmdname = cmdname;
    this->args = args;

}


/**
 *  Check is one option is found.
 *
 *  @param option    option name
 *  @return b        true if found
 */
bool Options::has_option(const std::string& option){
    char** begin = this->_argv;
    char** end = begin + _argc; 
    for (size_t i=0; i< this->_opts.size(); ++i){
        if (_opts[i].name_is(option)) {
            std::string optname = _opts[i].has_option(begin, end);
            if (optname.size() > 0) { return true; }
        }
    }
    return false;
}


/**
 * Add one option to the declarations.
 * 
 * @param code      option codes (e.g. -h,--help)
 * @param helpmsg   help message
 */
void Options::add_option(std::string code, std::string helpmsg){
    this->_opts.push_back(Option(code, helpmsg));
}


/**
 * Parse the command line.
 *
 * @param _argc  size of _argv
 * @param _argv  _argv command line
 */
void Options::parse_options(int _argc, char** _argv){
    this->_argc = _argc;
    this->_argv = _argv;

    char** begin = this->_argv;
    char** end = begin + _argc; 
    for (size_t i=0; i<this->_opts.size(); ++i){
        std::string optname = _opts[i].has_option(begin, end);
        if (optname.size() > 0){
            char ** itr = std::find(begin, end, optname);
            if (itr != end && ++itr != end)
            {
                _opts[i].set_value(*itr);
            }
        }
    }
}

/**
 * <template>
 * Get a given option and parse to a given type.
 *
 * @param option  option name to get
 * @param T       type to infer
 * @return val    value
 */
template <typename T> 
T Options::get_option(const std::string& option){
    for (size_t i=0; i<this->_opts.size(); ++i){
        if (this->_opts[i].name_is(option)){
            T val;
            std::istringstream(this->_opts[i].get_value()) >> val;
            return val;
        }
    }
    throw std::runtime_error("options not found");
}

/* explicit specializations of the templates for double, float, ... */

template <> double Options::get_option<double>(const std::string& option){
    std::string val = this->get_option<std::string>(option);
    double v = 0;
    if (val.size() > 0)
        std::istringstream(val) >> v;
    return v;
}

template <> float Options::get_option<float>(const std::string& option){
    std::string val = this->get_option<std::string>(option);
    float v = 0;
    if (val.size() > 0)
        std::istringstream(val) >> v;
    return v;
}

template <> size_t Options::get_option<size_t>(const std::string& option){
    std::string val = this->get_option<std::string>(option);
    size_t v = 0;
    if (val.size() > 0)
        std::istringstream(val) >> v;
    return v;
}

template <> int Options::get_option<int>(const std::string& option){
    std::string val = this->get_option<std::string>(option);
    int v = 0;
    if (val.size() > 0)
        std::istringstream(val) >> v;
    return v;
}

/* explicit instanciation of the templates for double, float, ... */

template double Options::get_option<double>(const std::string& option);
template float Options::get_option<float>(const std::string& option);
template int Options::get_option<int>(const std::string& option);
template size_t Options::get_option<size_t>(const std::string& option);


/** 
 * Construct and retur the help message
 *
 * @return str   help message
 */
std::string Options::help(){

    std::stringstream msg;

    msg << "Usage:" << std::endl;
    msg << "  " << this->cmdname << " [OPTION...] "<< this->args 
        << std::endl << std::endl;

    for (size_t i=0; i<this->_opts.size() ; ++i){
        Option o = _opts[i];
        size_t length = o.codes[0].size();
        if (o.has_short()){
            msg << "  " << o.codes[0];
        } else {
            msg << "      " << o.codes[0];  // add -x, equiv. space
            length += 4;
        }
        

        for (size_t j=1; j < o.codes.size(); ++j){
            msg << ", " << o.codes[j];
            length += o.codes[j].size() + 2;
        }
        std::string padding = "";
        padding.append(20 - length, ' ');
        msg << padding << o.helpmsg << std::endl;
    }
    return msg.str();
}

} // end namespace
// vim: expandtab:ts=4:softtabstop=4:shiftwidth=4
