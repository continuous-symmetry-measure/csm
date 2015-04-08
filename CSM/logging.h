/*
 * This file contains a very simple logging interface.
 *
 * The logging is eventually done by boost::log, but it is somewhat hidden away
 *
 * Developed by Itay Zandbank
 */

#ifndef LOGGING_H
#define LOGGING_H

#ifdef __APPLE__
/*
 * boost::log requires the Boost libraries, which require everything to be compiled with C++ 11 on
 * OS X. So instead we're using a stupid log replacement, that will eventually be replaced with the Pythong
 * logging library.
 */

#include <iostream>
#include <fstream>
#include <string>

enum LoggingSeverities
{
    fatal, error, info, debug
};

extern std::ofstream null_output;

//#define LOG(sev) ((sev==fatal || sev==error) ? std::cerr : (sev==info ? std::cout : (std::ostream)null_output))

#define LOG(sev) ((sev==fatal || sev==error) ? std::cerr : (sev==info ? std::cout : null_output))

#else

#include <boost/log/trivial.hpp>

#define LOG(sev) BOOST_LOG_TRIVIAL(sev)

#endif

void init_logging();  // Initialize the logging subsystem
void set_file_logging(std::string path);

#endif