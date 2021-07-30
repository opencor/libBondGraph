#ifndef EXCEPTIONS_H_
#define EXCEPTIONS_H_
#include <stdexcept>
#include <string>
#include "export.h"
#include "logging.h"

namespace BG {

	class EXPORTED BGException : public std::runtime_error {
	public:
		BGException(const std::string & inMessage) : std::runtime_error(inMessage) {
#ifdef ENABLE_LOGGING
				logCritical(inMessage);
#endif
        }
	};
}

#endif