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

#include <iostream>
#include <string>

#include "CartesianMesh3D.hxx"
#include "ExternalSource.hxx"
#include "Mesh3DPartitioning.hxx"

ExternalSource* ExternalSource::pInstance_ = nullptr;

ExternalSource*
ExternalSource::getInstance(const std::string eSourceForceFile)
{
  if (nullptr == pInstance_) {
    pInstance_ = new ExternalSource(eSourceForceFile);
    return pInstance_;
  } else {
    return pInstance_;
  }
}

void
ExternalSource::releaseInstance()
{
  if (pInstance_) {
    delete pInstance_;
    pInstance_ = nullptr;
  }
}

const SWS::SpatialField&
ExternalSource::operator()(const SWS::Directions& d) const
{
  return F_(d);
}

ExternalSource::ExternalSource(const std::string eSourceForceFile)
  : eSourceForceFile_(eSourceForceFile)
{
  for (auto d : { SWS::X, SWS::Y, SWS::Z }) {
    Mesh3DPartitioning::getInstance()->buildSpatialField(F_(d));
  }

  // Fill the Force field with input data
  intializeSourceForce();
}

int
ExternalSource::intializeSourceForce()
{

  const auto& lnxx = Mesh3DPartitioning::getInstance()->lnxx();
  const auto& lnyy = Mesh3DPartitioning::getInstance()->lnyy();
  const auto& lnzz = Mesh3DPartitioning::getInstance()->lnzz();

  for (int ii = 0; ii < lnxx; ii++) {
    for (int jj = 0; jj < lnyy; jj++) {
      for (int kk = 0; kk < lnzz; kk++) {
        F_(SWS::X)(ii, jj, kk) = 1.0;
        F_(SWS::Y)(ii, jj, kk) = 2.0;
        F_(SWS::Z)(ii, jj, kk) = 3.0;
      }
    }
  }

  return 0;
}

ExternalSource::~ExternalSource() {}
