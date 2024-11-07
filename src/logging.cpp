/*******************************************************************************

Copyright (C) The University of Auckland

OpenCOR is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

OpenCOR is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <https://gnu.org/licenses>.

*******************************************************************************/


#include "logging.h"
#include <iostream>
#include <sstream>

#ifdef ENABLE_LOGGING
#include <spdlog/cfg/env.h> // support for loading levels from the environment variable
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>
#endif

namespace BG {
#ifdef ENABLE_LOGGING


namespace{
    class BGLogger{
        public:
         BGLogger(){
             spdlog::cfg::load_env_levels();
         };

         ~BGLogger(){
             //Cleanup at close

         };
    };
    BGLogger logger;
}

void _logCritical(std::string msg){
    spdlog::critical(msg);
}

void _logDebug(std::string msg){
    spdlog::debug(msg);
}

void _logInfo(std::string msg){
    spdlog::info(msg);
}

void _logWarn(std::string msg){
    spdlog::warn(msg);
}

void setLogLevelInfo()
{
    spdlog::set_level(spdlog::level::info); // Set global log level to info
}

void setLogLevelWarn()
{
    spdlog::set_level(spdlog::level::warn); // Set global log level to warn
}

void setLogLevelOff()
{
    spdlog::set_level(spdlog::level::off); // Set global log level to off
}

void setLogLevelDebug()
{
    spdlog::set_level(spdlog::level::debug); // Set global log level to debug
}

void setLogLevelCritical()
{
    spdlog::set_level(spdlog::level::critical); // Set global log level to critical
}

void logToFile(std::string path)
{
    try{
        auto log = spdlog::basic_logger_mt("file", path);
        spdlog::set_default_logger(log);
    }catch(std::exception& ex){
        std::cerr << "Failed to create file logger " << ex.what() <<" file name "<<path<< std::endl;
    }
}


void logToConsole()
{
    try {
        auto log = spdlog::stdout_color_mt("console");
        spdlog::set_default_logger(log);
    } catch (std::exception &ex) {
        std::cerr << "Failed to create console logger " << ex.what() << std::endl;
    }
}

#else
void setLogLevelInfo(){}

void setLogLevelWarn(){}

void setLogLevelOff(){}

void setLogLevelDebug(){}

void setLogLevelCritical(){}

#endif

}