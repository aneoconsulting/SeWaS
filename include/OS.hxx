/*
  SeWaS
  Copyright (C) 2018  ANEO

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU Affero General Public License as published
  by the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Affero General Public License for more details.

  You should have received a copy of the GNU Affero General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
==============================================================================*/

#pragma once

#include <unistd.h>
#include <limits.h>
#include <string>
#include <iomanip>
#include <sstream>

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX // prevent windows from redefining min and max macros
#endif
#endif

namespace SWS
{
  inline std::string getDateTime()
  {
    std::string date;

    auto t = std::time(nullptr);
    auto lt = std::localtime(&t);

    std::ostringstream oss;
    oss << std::put_time(lt, "%d-%m-%Y_%H-%M-%S");
    date=oss.str();

    return date;
  }

  inline std::string getHostName()
  {
    static char hostname[HOST_NAME_MAX];

#ifdef _WIN32
    GetComputerName(hostname, HOST_NAME_MAX);
#else
    gethostname(hostname, HOST_NAME_MAX);
#endif

    return hostname;
  }

  inline auto getPID()
  {
#ifdef _WIN32
    return static_cast<int>(::GetCurrentProcessId());
#else
    return static_cast<int>(::getpid());
#endif
  }
};
