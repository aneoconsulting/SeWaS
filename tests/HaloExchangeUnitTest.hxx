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

#include "CartesianMesh3D.hxx"
#include "CentralFDOperator.hxx"
#include "HaloManager.hxx"
#include "Mesh3DPartitioning.hxx"
#include "SEWASParameterManager.hxx"

class HaloExchangeTest : public ::testing::Test
{
protected:
  void SetUp() override
  {
    ASSERT_TRUE((std::is_same<SWS::RealType, float>::value || std::is_same<SWS::RealType, double>::value));

    nx_ = 8;
    ny_ = 8;
    nz_ = 8;

    char** argv = (char**)calloc(argc_, sizeof(char*));
    for (int i = 0; i < argc_; i++) {
      argv[i] = (char*)argv_[i];
    }

    auto pm = *std::make_unique<SEWASParameterManager>(&argc_, &argv);

    int hnx = CentralFDOperator::hnx();
    int hny = CentralFDOperator::hny();
    int hnz = CentralFDOperator::hnz();

    ASSERT_NE(nullptr, CartesianMesh3D::getInstance(pm));
    ASSERT_NE(nullptr, Mesh3DPartitioning::getInstance(pm, hnx, hny, hnz));

    auto init = [&](SWS::SpatialBlockField<SWS::RealType>& fijk) {
      for (int k = fijk.kStart(); k < fijk.kEnd(); k++) {
        SWS::RealType val = 16 * (k - fijk.kStart()) + 1.0;
        for (int i = fijk.iStart(); i < fijk.iEnd(); i++) {
          for (int j = fijk.jStart(); j < fijk.jEnd(); j++) {
            fijk(i, j, k) = val;
            val += 1.0;
          }
        }
      }
    };

    for (auto d : { SWS::X, SWS::Y, SWS::Z }) {
      Mesh3DPartitioning::getInstance()->buildSpatialField(v_(d));

      ASSERT_EQ(2, v_(d).dimension(SWS::X));
      ASSERT_EQ(2, v_(d).dimension(SWS::Y));
      ASSERT_EQ(2, v_(d).dimension(SWS::Z));

      const auto nxx = v_(d).dimension(SWS::X);
      const auto nyy = v_(d).dimension(SWS::Y);
      const auto nzz = v_(d).dimension(SWS::Z);

      for (int kk = 0; kk < nzz; kk++) {
        for (int ii = 0; ii < nxx; ii++) {
          for (int jj = 0; jj < nyy; jj++) {
            init(v_(d)(ii, jj, kk));
          } // kk
        }   // jj
      }     // ii
    }       // d

    for (auto sc : { SWS::XX, SWS::YY, SWS::ZZ, SWS::XY, SWS::XZ, SWS::YZ }) {
      Mesh3DPartitioning::getInstance()->buildSpatialField(sigma_(sc));

      ASSERT_EQ(2, sigma_(sc).dimension(SWS::X));
      ASSERT_EQ(2, sigma_(sc).dimension(SWS::Y));
      ASSERT_EQ(2, sigma_(sc).dimension(SWS::Z));

      const auto nxx = sigma_(sc).dimension(SWS::X);
      const auto nyy = sigma_(sc).dimension(SWS::Y);
      const auto nzz = sigma_(sc).dimension(SWS::Z);

      for (int kk = 0; kk < nzz; kk++) {
        for (int ii = 0; ii < nxx; ii++) {
          for (int jj = 0; jj < nyy; jj++) {
            init(sigma_(sc)(ii, jj, kk));
          } // kk
        }   // jj
      }     // ii
    }       // sc

    ASSERT_NE(nullptr, HaloManager::getInstance(sigma_, v_, hnx, hny, hnz));

    /* k == 2
       ------

      4 | 8 | 12 | 16 |  | 4 | 8 | 12 | 16
      --------------- |  | ---------------
      3 | 7 | 11 | 15 |  | 3 | 7 | 11 | 15
      --------------- |  | ---------------
      2 | 6 | 10 | 14 |  | 2 | 6 | 10 | 14
      --------------- |  | ---------------
      1 | 5 |  9 | 13 |  | 1 | 5 |  9 | 13
      ----------------|  |----------------

      ----------------|  |----------------
      4 | 8 | 12 | 16 |  | 4 | 8 | 12 | 16
      --------------- |  | ---------------
      3 | 7 | 11 | 15 |  | 3 | 7 | 11 | 15
      --------------- |  | ---------------
      2 | 6 | 10 | 14 |  | 2 | 6 | 10 | 14
      --------------- |  | ---------------
      1 | 5 |  9 | 13 |  | 1 | 5 |  9 | 13


          ii == 0              ii == 1

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

    const auto lnxx = Mesh3DPartitioning::getInstance()->lnxx();
    const auto lnyy = Mesh3DPartitioning::getInstance()->lnyy();
    const auto lnzz = Mesh3DPartitioning::getInstance()->lnzz();

    for (auto d : { SWS::X, SWS::Y, SWS::Z }) {
      vH_(d).resize(lnxx, lnyy, lnzz);

      for (int lkk = 0; lkk < lnxx; lkk++) {
        for (int lii = 0; lii < lnyy; lii++) {
          for (int ljj = 0; ljj < lnzz; ljj++) {
            buildArenas(vH_(d)(lii, ljj, lkk), lii, ljj, lkk);
          } // lii
        }   // ljj
      }     // lkk
    }

    for (auto sc : { SWS::XX, SWS::YY, SWS::ZZ, SWS::XY, SWS::XZ, SWS::YZ }) {
      sigmaH_(sc).resize(lnxx, lnyy, lnzz);

      for (int lkk = 0; lkk < lnxx; lkk++) {
        for (int lii = 0; lii < lnyy; lii++) {
          for (int ljj = 0; ljj < lnzz; ljj++) {
            buildArenas(sigmaH_(sc)(lii, ljj, lkk), lii, ljj, lkk);
          } // lii
        }   // ljj
      }     // lkk
    }
  }

  void TearDown() override
  {
    const auto lnxx = Mesh3DPartitioning::getInstance()->lnxx();
    const auto lnyy = Mesh3DPartitioning::getInstance()->lnyy();
    const auto lnzz = Mesh3DPartitioning::getInstance()->lnzz();

    for (auto l : { SWS::Locations::LEFT,
                    SWS::Locations::RIGHT,
                    SWS::Locations::BACKWARD,
                    SWS::Locations::FORWARD,
                    SWS::Locations::BOTTOM,
                    SWS::Locations::TOP }) {

      for (int lkk = 0; lkk < lnxx; lkk++) {
        for (int lii = 0; lii < lnyy; lii++) {
          for (int ljj = 0; ljj < lnzz; ljj++) {

            for (auto d : { SWS::X, SWS::Y, SWS::Z }) {
              auto& halo = vH_(d)(lii, lkk, lkk)(l);
              delete[] halo;
              halo = nullptr;
            }

            for (auto sc : { SWS::XX, SWS::YY, SWS::ZZ, SWS::XY, SWS::XZ, SWS::YZ }) {
              auto& halo = sigmaH_(sc)(lii, lkk, lkk)(l);
              delete[] halo;
              halo = nullptr;
            }

          } // ljj
        }   // lii
      }     // lkk

    } // l

    HaloManager::releaseInstance();
    Mesh3DPartitioning::releaseInstance();
    CartesianMesh3D::releaseInstance();
  }

private:
  void buildArenas(SWS::Halo& halo, const int lii, const int ljj, const int lkk)
  {
    for (auto l : { SWS::Locations::LEFT,
                    SWS::Locations::RIGHT,
                    SWS::Locations::BACKWARD,
                    SWS::Locations::FORWARD,
                    SWS::Locations::BOTTOM,
                    SWS::Locations::TOP }) {

      size_t hsize = 0;

      const auto cx = Mesh3DPartitioning::getInstance()->ccx()[lii];
      const auto cy = Mesh3DPartitioning::getInstance()->ccy()[ljj];
      const auto cz = Mesh3DPartitioning::getInstance()->ccz()[lkk];

      switch (l) {
        case SWS::Locations::LEFT:
        case SWS::Locations::RIGHT:
          hsize = 2 * cy * cz;
          break;
        case SWS::Locations::BACKWARD:
        case SWS::Locations::FORWARD:
          hsize = cx * 2 * cz;
          break;
        case SWS::Locations::BOTTOM:
        case SWS::Locations::TOP:
          hsize = cx * cy * 2;
          break;
        default:
          break;
      }

      ASSERT_EQ(hsize, HaloManager::getInstance()->getHaloSize(l, lii, ljj, lkk));

      halo(l) = new SWS::RealType[hsize];
    }
  }

protected:
  static inline int argc_ = 17;
  static inline const char* argv_[] = { "sewas", 
                                        "--cx", "4", "--cy", "4", "--cz", "4", "--P", "1", "--Q", "1", "--R", "1", "--nthreads", "1", "--dfile", "TestX.json" };

  int nx_;
  int ny_;
  int nz_;

  /* The halo exchange test is performed for a given time-step */
  const int ts_ = 0;

  SWS::Velocity v_;
  SWS::StressField sigma_;

  /* Halos */
  SWS::VelocityHalo vH_;
  SWS::StressFieldHalo sigmaH_;
};
