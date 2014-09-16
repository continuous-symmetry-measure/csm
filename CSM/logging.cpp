/*
 * Implementation of the thin wrapper around boost::log
 *
 * Developed by Itay Zandbank
 */

#include "logging.h"

#include <boost/log/utility/setup/console.hpp>
#include <boost/log/utility/setup/file.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>
#include <boost/log/sources/logger.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/trivial.hpp>

using namespace boost::log;

void init_logging()
{
	add_common_attributes();
	register_simple_formatter_factory< trivial::severity_level, char >("Severity");
	add_console_log(std::cerr,
		keywords::format = "[%TimeStamp%]  [%Severity%]: %Message%",
		keywords::filter = trivial::severity >= trivial::severity_level::error);
	add_console_log(std::cout, keywords::format = "%Message%");
	add_file_log("csm.log",
		keywords::format = "[%TimeStamp%]  [%Severity%]: %Message%");

	set_debug_logging(false);
}

void set_debug_logging(bool enable)
{
	if (enable)
		core::get()->set_filter(trivial::severity >= trivial::severity_level::trace);
	else
		core::get()->set_filter(trivial::severity >= trivial::severity_level::info);
}
