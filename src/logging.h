#ifndef LOGGING__H__
#define LOGGING__H__
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

#endif