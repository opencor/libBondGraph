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

#pragma once

#include "export.h"
#include <memory>
#include <sstream>

namespace BG {
//Logging related setups

EXPORTED void setLogLevelInfo();

EXPORTED void setLogLevelWarn();

EXPORTED void setLogLevelOff();

EXPORTED void setLogLevelDebug();

EXPORTED void setLogLevelCritical();

//Templated function should not be defined in cpp as the compiler is not able to find them when the library is included downstream
#ifndef ENABLE_LOGGING
template<typename Arg, typename... Args>
void logCritical(Arg&& arg, Args&&... args){}

template<typename Arg, typename... Args>
void logDebug(Arg &&arg, Args &&...args){}

template<typename Arg, typename... Args>
void logInfo(Arg &&arg, Args &&...args){}

template<typename Arg, typename... Args>
void logWarn(Arg &&arg, Args &&...args){}

#else

void _logCritical(std::string msg);
void _logDebug(std::string msg);
void _logInfo(std::string msg);
void _logWarn(std::string msg);

template<typename Arg, typename... Args>
void logCritical(Arg &&arg, Args &&...args)
{
    std::ostringstream out;
    out << std::forward<Arg>(arg);
    using expander = int[];
    (void)expander {0, (void(out << " " << std::forward<Args>(args)), 0)...};
    std::string msg = out.str();
    _logCritical(msg);
}

template<typename Arg, typename... Args>
void logDebug(Arg &&arg, Args &&...args)
{
    std::ostringstream out;
    out << std::forward<Arg>(arg);
    using expander = int[];
    (void)expander {0, (void(out << " " << std::forward<Args>(args)), 0)...};
    std::string msg = out.str();
    _logDebug(msg);
}

template<typename Arg, typename... Args>
void logInfo(Arg &&arg, Args &&...args)
{
    std::ostringstream out;
    out << std::forward<Arg>(arg);
    using expander = int[];
    (void)expander {0, (void(out << " " << std::forward<Args>(args)), 0)...};
    std::string msg = out.str();
    _logInfo(msg);
}

template<typename Arg, typename... Args>
void logWarn(Arg &&arg, Args &&...args)
{
    std::ostringstream out;
    out << std::forward<Arg>(arg);
    using expander = int[];
    (void)expander {0, (void(out << " " << std::forward<Args>(args)), 0)...};
    std::string msg = out.str();
    _logWarn(msg);
}

#endif

}

