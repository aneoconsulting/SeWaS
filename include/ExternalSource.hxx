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

#include "Config.hxx"

class ExternalSource
{
public:
  static ExternalSource* getInstance(const std::string eSourceForceFile = "");
  static void releaseInstance();

  const SWS::SpatialField& operator()(const SWS::Directions& d) const;

  int intializeSourceForce();

private:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW // To force an object of this class to be allocated as aligned

  ExternalSource(const std::string eSourceForceFile);
  ~ExternalSource();

  static ExternalSource* pInstance_;

  std::string eSourceForceFile_;

  SWS::SourceForce F_;
};
