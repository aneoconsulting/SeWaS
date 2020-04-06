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

#include "CartesianMesh3D.hxx"

#include "LogManager.hxx"

CartesianMesh3D* CartesianMesh3D::pInstance_ = nullptr;

CartesianMesh3D*
CartesianMesh3D::getInstance(const int nx, const int ny, const int nz, const SWS::RealType ds)
{
  if (nullptr == pInstance_) {
    pInstance_ = new CartesianMesh3D(nx, ny, nz, ds);
    return pInstance_;
  } else {
    return pInstance_;
  }
}

void
CartesianMesh3D::releaseInstance()
{
  if (pInstance_) {
    delete pInstance_;
    pInstance_ = nullptr;
  }
}

CartesianMesh3D::CartesianMesh3D(const int nx, const int ny, const int nz, const SWS::RealType ds)
{
  ds_ = ds;

  dx_ = ds_;
  dy_ = ds_;
  dz_ = ds_;

  nx_ = nx;
  ny_ = ny;
  nz_ = nz;

  lx_ = nx * dx_;
  ly_ = ny * dy_;
  lz_ = nz * dz_;

  LOG(SWS::LOG_INFO, "(dx, dy, dz)=({}, {}, {}), (lx, ly, lz)=({}, {}, {})", dx_, dy_, dz_, lx_, ly_, lz_);
}

CartesianMesh3D::~CartesianMesh3D() {}
