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

#include <fstream>
#include <iostream>
#include <string>

#include "CartesianMesh3D.hxx"
#include "Config.hxx"
#include "Mesh3DPartitioning.hxx"
#include "SEWASParameterManager.hxx"

class DataSet
{
public:
  static DataSet* getInstance(const SEWASParameterManager& pm);
  static DataSet* getInstance();

  static void releaseInstance();

  inline const auto& lambda(const SWS::StressFieldComponents& sc) const { return lambda_(sc); }

  inline const auto& mu(const SWS::StressFieldComponents& sc) const { return mu_(sc); }

  inline const auto& rho() const { return rho_; }

  inline const auto& b(const SWS::Directions& d) const { return b_(d); }

  void initialize(const int ii, const int jj, const int kk);

private:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW // To force an object of this class to be allocated as aligned

  DataSet(const SEWASParameterManager& pm);
  ~DataSet();

  const auto getLayer(const int k, const int ii, const int jj, const int kk) const;

  static DataSet* pInstance_;

  const SEWASParameterManager& pm_;

  SWS::Elasticity lambda_;
  SWS::Elasticity mu_;
  SWS::MaterialDensity rho_;
  SWS::Buoyancy b_;
};
