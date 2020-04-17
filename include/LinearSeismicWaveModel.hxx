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

#include <iostream>

#include <Eigen/Core>

#include "CentralFDOperator.hxx"
#include "Config.hxx"
#include "HaloManager.hxx"
#include "SEWASParameterManager.hxx"

class LinearSeismicWaveModel
{
public:
  static LinearSeismicWaveModel* getInstance(const CentralFDOperator& fdo,
                                             const SEWASParameterManager& pm);
  static LinearSeismicWaveModel* getInstance();

  static void releaseInstance();

  int propagate() noexcept;

  inline const auto& v(const SWS::Directions& d) const { return v_(d); }

  inline const auto& sigma(const SWS::StressFieldComponents& sc) const { return sigma_(sc); }

  inline const auto& nt() const { return nt_; }

  /* Entry points called by the PaRSEC runtime when executing tasks' BODY */
  static int initializeFieldsWrapper(const int ii, const int jj, const int kk);

  static int computeVelocityWrapper(const int d, const int ts, const int ii, const int jj, const int kk);

  static int computeStressWrapper(const int sc, const int ts, const int ii, const int jj, const int kk);

private:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW // To force an object of this class to be allocated as aligned

  LinearSeismicWaveModel(const CentralFDOperator& fdo,
                         const SEWASParameterManager& pm);
  ~LinearSeismicWaveModel();

  void initialize(const int ii, const int jj, const int kk);

  void setVelocitySourceLocations();

  void addVelocitySource(const SWS::Directions& d, const int ts, const int ii, const int jj, const int kk);

  int computeVelocity(const SWS::Directions& d, const int& ts, const int& ii, const int& jj, const int& kk);

  int computeStress(const SWS::StressFieldComponents& sc,
                    const int& ts,
                    const int& ii,
                    const int& jj,
                    const int& kk);

  static LinearSeismicWaveModel* pInstance_;

  const CentralFDOperator& fdo_;
  const SEWASParameterManager& pm_;

  int nt_;
  float tmax_;

  SWS::RealType dt_;

  SWS::Velocity v_;

  SWS::StressField sigma_;
};
