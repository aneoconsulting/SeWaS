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

#include "MetricsManager.hxx"

#include <cassert>

MetricsManager * MetricsManager::pInstance_ = nullptr;

MetricsManager * MetricsManager::getInstance(const std::string applicationName)
{
  if (nullptr == pInstance_){
    pInstance_ = new MetricsManager(applicationName);
    return pInstance_;
  }
  else{
    return pInstance_;
  }
}

MetricsManager * MetricsManager::getInstance()
{
  assert(nullptr != pInstance_);
  return pInstance_;
}

void MetricsManager::releaseInstance()
{
  if (pInstance_){
    delete pInstance_;
    pInstance_ = nullptr;
  }
}

MetricsManager::MetricsManager(const std::string applicationName) : applicationName_(applicationName)
{
}

MetricsManager::~MetricsManager()
{
}
