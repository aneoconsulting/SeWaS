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

#include "Config.hxx"

class CartesianMesh3D
{
public:
  static CartesianMesh3D* getInstance(const int nx = 1,
                                      const int ny = 1,
                                      const int nz = 1,
                                      const SWS::RealType ds = 1.);
  static void releaseInstance();

  inline const int& nx() const { return nx_; }
  inline const int& ny() const { return ny_; }
  inline const int& nz() const { return nz_; }

  inline const SWS::RealType& dx() const { return dx_; }
  inline const SWS::RealType& dy() const { return dy_; }
  inline const SWS::RealType& dz() const { return dz_; }

  inline const SWS::RealType& ds() const { return ds_; }

private:
  CartesianMesh3D(const int nx,
                  const int ny,
                  const int nz,
                  const SWS::RealType ds);
  ~CartesianMesh3D();

  static CartesianMesh3D* pInstance_;

  int nx_;
  int ny_;
  int nz_;

  SWS::RealType lx_;
  SWS::RealType ly_;
  SWS::RealType lz_;

  SWS::RealType ds_;

  SWS::RealType dx_;
  SWS::RealType dy_;
  SWS::RealType dz_;
};
