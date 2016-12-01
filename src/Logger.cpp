/** 
* @file Logger.cpp
* @synopsis implemention of log init
* @author Yongsheng Chen
* @date 2012-03-28
*/

#include "Logger.h"

log4cplus::Logger g_logger = log4cplus::Logger::getInstance("Global Logger");

bool init_logger(const std::string &appenders, const std::string &log_out_pathname )
{
	SharedAppenderPtr append;
	if(appenders=="file")
	{
		if(!log_out_pathname.empty())
		{
			append = SharedAppenderPtr(new FileAppender(log_out_pathname));
		}
		else
		{
			cerr<<"Please input the log pathname!\n";
		}
	}
	else
	{
		append = SharedAppenderPtr(new ConsoleAppender());
	}	
	//string pattern = "%D{%Y-%m-%d %H:%M:%S} %-5p - %m [%l]%n";
	string pattern = "%D{%Y-%m-%d %H:%M:%S} %-5p - %m%n";
	auto_ptr<Layout> layout(new PatternLayout(pattern));
	//std::auto_ptr layout(new TTCCLayout());
	append->setLayout(layout);
	g_logger.addAppender(append);
	return true;
};

bool init_logger()
{
	SharedAppenderPtr append = SharedAppenderPtr(new ConsoleAppender());
	string pattern = "%D{%Y-%m-%d %H:%M:%S} %-5p - %m\t#L%n";
	auto_ptr<Layout> layout(new PatternLayout(pattern));
	append->setLayout(layout);
	g_logger.addAppender(append);
	return true;
};

