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
CartesianMesh3D::getInstance(const SEWASParameterManager& pm)
{
  if (nullptr == pInstance_) {
    pInstance_ = new CartesianMesh3D(pm);
    return pInstance_;
  } else {
    return pInstance_;
  }
}

CartesianMesh3D*
CartesianMesh3D::getInstance()
{
  assert(nullptr != pInstance_);
  return pInstance_;
}

void
CartesianMesh3D::releaseInstance()
{
  if (pInstance_) {
    delete pInstance_;
    pInstance_ = nullptr;
  }
}

CartesianMesh3D::CartesianMesh3D(const SEWASParameterManager& pm)
{
  ds_ = pm.ds();

  dx_ = ds_;
  dy_ = ds_;
  dz_ = ds_;

  nx_ = pm.nx();
  ny_ = pm.ny();
  nz_ = pm.nz();

  lx_ = pm.lx();
  ly_ = pm.ly();
  lz_ = pm.lz();
}

CartesianMesh3D::~CartesianMesh3D() {}
