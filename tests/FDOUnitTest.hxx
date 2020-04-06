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

#include "CartesianMesh3D.hxx"
#include "CentralFDOperator.hxx"

class FDOTest : public ::testing::Test
{
protected:
  void SetUp() override
  {
    ASSERT_TRUE((std::is_same<SWS::RealType, float>::value || std::is_same<SWS::RealType, double>::value));

    nx_ = 4;
    ny_ = 4;
    nz_ = 4;

    int hnx = CentralFDOperator::hnx();
    int hny = CentralFDOperator::hny();
    int hnz = CentralFDOperator::hnz();

    // 4x4x4 tile with halos
    cx_ = hnx + 4 + hnx;
    cy_ = hny + 4 + hny;
    cz_ = hnz + 4 + hnz;

    ds_ = 0.1;

    // Create the computational domain
    ASSERT_NE(nullptr, CartesianMesh3D::getInstance(nx_, ny_, nz_, ds_));

    auto pCartesianMesh = CartesianMesh3D::getInstance();

    // Create a finite difference operator
    fdo_ =
      std::make_shared<CentralFDOperator>(pCartesianMesh->dx(), pCartesianMesh->dy(), pCartesianMesh->dz());

    fijk_.resize(cx_, cy_, cz_);
    fijk_.setHaloSize(hnx, hny, hnz);
    fijk_ = 0.0;

    /* k == 2
       ------

      4 | 8 | 12 | 16
      ---------------
      3 | 7 | 11 | 15
      ---------------
      2 | 6 | 10 | 14
      ---------------
      1 | 5 |  9 | 13

      ---> i

      /////////////////

      (i == 2, j == 2)
      ----------------
      59
      --
      43
      --
      27
      --
      11
    */

    for (int k = fijk_.kStart(); k < fijk_.kEnd(); k++) {
      SWS::RealType val = 16 * (k - fijk_.kStart()) + 1.0;
      for (int i = fijk_.iStart(); i < fijk_.iEnd(); i++) {
        for (int j = fijk_.jStart(); j < fijk_.jEnd(); j++) {
          fijk_(i, j, k) = val;
          val += 1.0;
        }
      }
    }

    // fijk_(i0,j0,j0) == 43
    i0_ = fijk_.iStart() + 2;
    j0_ = fijk_.jStart() + 2;
    k0_ = fijk_.kStart() + 2;

    rijk_.resize(cx_, cy_, cz_);
    fijk_.setHaloSize(hnx, hny, hnz);
    rijk_ = 0.0;
  }

  void TearDown() override { CartesianMesh3D::releaseInstance(); }

  const SWS::RealType c1_ = -1. / 24;
  const SWS::RealType c2_ = +1. / 24;
  const SWS::RealType c3_ = +9. / 8;
  const SWS::RealType c4_ = -9. / 8;

  int nx_;
  int ny_;
  int nz_;

  int cx_;
  int cy_;
  int cz_;

  SWS::RealType ds_;

  std::shared_ptr<CentralFDOperator> fdo_;

  SWS::SpatialBlockField<SWS::RealType> fijk_;
  SWS::SpatialBlockField<SWS::RealType> rijk_;

  int i0_;
  int j0_;
  int k0_;
};
