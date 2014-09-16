/*
 * This file contains a very simple logging interface.
 *
 * The logging is eventually done by boost::log, but it is somewhat hidden away
 *
 * Developed by Itay Zandbank
 */

#ifndef LOGGING_H
#define LOGGING_H

#include <boost/log/trivial.hpp>

#define LOG(sev) BOOST_LOG_TRIVIAL(sev)

void init_logging();  // Initialize the logging subsystem
void set_debug_logging(bool enabled);  // Enable or disable debug logging

#endif