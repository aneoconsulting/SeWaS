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

#include <string>
#include <vector>

#include <spdlog/spdlog.h>
#include <spdlog/async.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>

#include "OS.hxx"
#include "Constants.hxx"

class LogManager{
public:
  static LogManager * getInstance();
  static void releaseInstance();


#ifdef VERBOSE
#define LOG(LEVEL, ...)                                 \
  LogManager::getInstance()->log<LEVEL>(__VA_ARGS__);
#else
#define LOG(LEVEL, ...)                                 \
  if constexpr (LEVEL <= SWS::LOG_ERROR){                   \
    LogManager::getInstance()->log<LEVEL>(__VA_ARGS__); \
  }
#endif


  inline auto & getLogger()
  {
    return logger_;
  }

  template<SWS::LogLevels level, typename... LogMsg>
  inline void log(LogMsg &&... msg)
  {
    switch(level){
    case SWS::LOG_TRACE:
      logger_->trace(msg...);
      break;
    case SWS::LOG_DEBUG:
      logger_->debug(msg...);
      break;
    case SWS::LOG_INFO:
      logger_->info(msg...);
      break;
    case SWS::LOG_WARN:
      logger_->warn(msg...);
      break;
    case SWS::LOG_ERROR:
      logger_->error(msg...);
      break;
    case SWS::LOG_CRITICAL:
      logger_->critical(msg...);
      break;
    default:
      logger_->warn("The selected log level ({}) is not supported. Fall back to INFO level.", level);
      logger_->info(msg...);
      break;
    }
  }

private:
  LogManager();
  ~LogManager();

  void setVerbosityLevel();

  static inline LogManager * pInstance_ = nullptr;

  std::string logFileName_;

  /* Thread pool to be used by spdlog engine for asynchronous logging */
  std::shared_ptr<spdlog::details::thread_pool> threadPool_;

  std::shared_ptr<spdlog::logger> logger_;
};
