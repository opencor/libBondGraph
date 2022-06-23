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

#include <stdexcept>
#include <string>
#include "export.h"
#include "logging.h"

namespace BG {

	class  BGException : public std::runtime_error {
	public:
		BGException(const std::string & inMessage) : std::runtime_error(inMessage) {
#ifdef ENABLE_LOGGING
				logCritical(inMessage);
#endif
        }
	};
}
