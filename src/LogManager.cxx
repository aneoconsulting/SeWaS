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

#include "LogManager.hxx"
#include "Version.hxx"

LogManager * LogManager::pInstance_ = nullptr;

LogManager * LogManager::getInstance()
{
  if (nullptr == pInstance_){
    pInstance_ = new LogManager();
    return pInstance_;
  }
  else{
    return pInstance_;
  }
}

void LogManager::releaseInstance()
{
  if (pInstance_){
    delete pInstance_;
    pInstance_ = nullptr;
  }
}

LogManager::LogManager()
{
  logFileName_ = "sewas_" + SWS::getDateTime() + "_" + SWS::getHostName() + "_" + std::to_string(SWS::getPID()) + ".log";

  threadPool_ = std::make_shared<spdlog::details::thread_pool>(8192, 1);

  auto consoleSink=std::make_shared<spdlog::sinks::ansicolor_stdout_sink_mt>();
  auto fileSink=std::make_shared<spdlog::sinks::basic_file_sink_mt>(logFileName_);

  /* Build the unique logger for dumping messages in both console and logfile */
  // logger_ = std::make_shared<spdlog::logger>("logger", std::initializer_list<spdlog::sink_ptr>{consoleSink, fileSink});
  logger_ = std::make_shared<spdlog::async_logger>("logger",
                                                   std::initializer_list<spdlog::sink_ptr>{consoleSink, fileSink},
                                                   threadPool_,
                                                   spdlog::async_overflow_policy::block);

  spdlog::register_logger(logger_);

  setVerbosityLevel();

  spdlog::set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%^%l%$] [process %P] [thread %t] %v");

  logger_->info("SeWaS {}.{}.{}", SEWAS_VERSION_MAJOR, SEWAS_VERSION_MINOR, SEWAS_VERSION_PATCH);
}

void LogManager::setVerbosityLevel()
{
#ifdef VERBOSE
  const char * vl=getenv("VERBOSITY_LEVEL");
  if (vl){
    if (strcmp("OFF", vl) == 0){
      logger_->set_level(spdlog::level::off);
    }
    else if (strcmp("CRITICAL", vl) == 0){
      logger_->set_level(spdlog::level::critical);
    }
    else if (strcmp("ERROR", vl) == 0){
      logger_->set_level(spdlog::level::err);
    }
    else if (strcmp("WARN", vl) == 0){
      logger_->set_level(spdlog::level::warn);
    }
    else if (strcmp("INFO", vl) == 0){
      logger_->set_level(spdlog::level::info);
    }
    else if (strcmp("DEBUG", vl) == 0){
      logger_->set_level(spdlog::level::debug);
    }
    else if (strcmp("TRACE", vl) == 0){
      logger_->set_level(spdlog::level::trace);
    }
    else{
      spdlog::warn("Unrecognized token {} when selecting verbosity level. Accepted values are: {}. Falling back to INFO mode", vl, "OFF, CRITICAL, ERROR, WARN, INFO, DEBUG, TRACE");
      logger_->set_level(spdlog::level::info);
    }
  }
  else{
    logger_->set_level(spdlog::level::info);
  }
#else
  logger_->set_level(spdlog::level::err);
#endif
}

LogManager::~LogManager()
{
  spdlog::drop_all();
  spdlog::shutdown();
}
