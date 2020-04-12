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

#include <type_traits>

#include <gtest/gtest.h>

#include "HaloUpdateUnitTest.hxx"
#include "Indexer.hxx"

TEST_F(HaloUpdateTest, UpdateLeftHalo)
{
  const auto d = SWS::X;

  const auto cx = Mesh3DPartitioning::getInstance()->ccx()[0];
  const auto cy = Mesh3DPartitioning::getInstance()->ccy()[0];
  const auto cz = Mesh3DPartitioning::getInstance()->ccz()[0];

  const auto hnx = CentralFDOperator::hnx();
  const auto hny = CentralFDOperator::hny();
  const auto hnz = CentralFDOperator::hnz();

  const auto iEnd = v_(0)(0, 0, 0).iEnd();

  HaloManager::extractVelocityHaloWrapper(SWS::RIGHT, d, ts_, 0, 0, 0, vH_(d)(0, 0, 0)(SWS::RIGHT));

  memcpy(vH_(d)(1, 0, 0)(SWS::LEFT),
         vH_(d)(0, 0, 0)(SWS::RIGHT),
         HaloManager::getInstance()->getHaloSize(SWS::RIGHT, 0, 0, 0) * sizeof(SWS::RealType));

  int iStartH = 0, iEndH = hnx;
  int jStartH = 0, jEndH = cy;
  int kStartH = 0, kEndH = cz;

  int iShift = 0;
  int jShift = 0;
  int kShift = 0;

  Indexer<SWS::Ordering> indexer(iEndH, jEndH, kEndH);

  HaloManager::updateVelocityWrapper(SWS::LEFT, d, ts_, 1, 0, 0, vH_(d)(1, 0, 0)(SWS::LEFT));

  for (int i = iStartH; i < iEndH; i++) {
    for (int j = jStartH; j < jEndH; j++) {
      for (int k = kStartH; k < kEndH; k++) {

        auto index = indexer(i, j, k);

        const auto& v = v_(d)(1, 0, 0)(i + iShift, j + jShift, k + kShift);

        if constexpr (std::is_same<SWS::RealType, float>::value) {
          EXPECT_FLOAT_EQ(v, vH_(d)(1, 0, 0)(SWS::LEFT)[index]);
        } else if constexpr (std::is_same<SWS::RealType, double>::value) {
          EXPECT_DOUBLE_EQ(v, vH_(d)(1, 0, 0)(SWS::LEFT)[index]);
        } else {
          std::cerr << "Unsupported element type for SWS::RealType : " << typeid(SWS::RealType).name()
                    << "\n";
          exit(SWS::INVALID_REAL_TYPE);
        }
      }
    }
  }
}

TEST_F(HaloUpdateTest, UpdateRightHalo)
{
  const auto d = SWS::X;

  const auto cx = Mesh3DPartitioning::getInstance()->ccx()[0];
  const auto cy = Mesh3DPartitioning::getInstance()->ccy()[0];
  const auto cz = Mesh3DPartitioning::getInstance()->ccz()[0];

  const auto hnx = CentralFDOperator::hnx();
  const auto hny = CentralFDOperator::hny();
  const auto hnz = CentralFDOperator::hnz();

  const auto iEnd = v_(0)(0, 0, 0).iEnd();

  HaloManager::extractVelocityHaloWrapper(SWS::LEFT, d, ts_, 1, 0, 0, vH_(d)(1, 0, 0)(SWS::LEFT));

  memcpy(vH_(d)(0, 0, 0)(SWS::RIGHT),
         vH_(d)(1, 0, 0)(SWS::LEFT),
         HaloManager::getInstance()->getHaloSize(SWS::RIGHT, 0, 0, 0) * sizeof(SWS::RealType));

  int iStartH = 0, iEndH = hnx;
  int jStartH = 0, jEndH = cy;
  int kStartH = 0, kEndH = cz;

  int iShift = iEnd;
  int jShift = 0;
  int kShift = 0;

  Indexer<SWS::Ordering> indexer(iEndH, jEndH, kEndH);

  HaloManager::updateVelocityWrapper(SWS::RIGHT, d, ts_, 0, 0, 0, vH_(d)(0, 0, 0)(SWS::RIGHT));

  for (int i = iStartH; i < iEndH; i++) {
    for (int j = jStartH; j < jEndH; j++) {
      for (int k = kStartH; k < kEndH; k++) {

        auto index = indexer(i, j, k);

        const auto& v = v_(d)(0, 0, 0)(i + iShift, j + jShift, k + kShift);

        if constexpr (std::is_same<SWS::RealType, float>::value) {
          EXPECT_FLOAT_EQ(v, vH_(d)(0, 0, 0)(SWS::RIGHT)[index]);
        } else if constexpr (std::is_same<SWS::RealType, double>::value) {
          EXPECT_DOUBLE_EQ(v, vH_(d)(0, 0, 0)(SWS::RIGHT)[index]);
        } else {
          std::cerr << "Unsupported element type for SWS::RealType : " << typeid(SWS::RealType).name()
                    << "\n";
          exit(SWS::INVALID_REAL_TYPE);
        }
      }
    }
  }
}

TEST_F(HaloUpdateTest, UpdateBackwardHalo)
{
  const auto d = SWS::X;

  const auto cx = Mesh3DPartitioning::getInstance()->ccx()[0];
  const auto cy = Mesh3DPartitioning::getInstance()->ccy()[0];
  const auto cz = Mesh3DPartitioning::getInstance()->ccz()[0];

  const auto hnx = CentralFDOperator::hnx();
  const auto hny = CentralFDOperator::hny();
  const auto hnz = CentralFDOperator::hnz();

  const auto jEnd = v_(0)(0, 0, 0).jEnd();

  HaloManager::extractVelocityHaloWrapper(SWS::FORWARD, d, ts_, 0, 0, 0, vH_(d)(0, 0, 0)(SWS::FORWARD));

  memcpy(vH_(d)(0, 1, 0)(SWS::BACKWARD),
         vH_(d)(0, 0, 0)(SWS::FORWARD),
         HaloManager::getInstance()->getHaloSize(SWS::FORWARD, 0, 0, 0) * sizeof(SWS::RealType));

  int iStartH = 0, iEndH = cx;
  int jStartH = 0, jEndH = hny;
  int kStartH = 0, kEndH = cz;

  int iShift = 0;
  int jShift = 0;
  int kShift = 0;

  Indexer<SWS::Ordering> indexer(iEndH, jEndH, kEndH);

  HaloManager::updateVelocityWrapper(SWS::BACKWARD, d, ts_, 0, 1, 0, vH_(d)(0, 1, 0)(SWS::BACKWARD));

  for (int i = iStartH; i < iEndH; i++) {
    for (int j = jStartH; j < jEndH; j++) {
      for (int k = kStartH; k < kEndH; k++) {

        auto index = indexer(i, j, k);

        const auto& v = v_(d)(0, 1, 0)(i + iShift, j + jShift, k + kShift);

        if constexpr (std::is_same<SWS::RealType, float>::value) {
          EXPECT_FLOAT_EQ(v, vH_(d)(0, 1, 0)(SWS::BACKWARD)[index]);
        } else if constexpr (std::is_same<SWS::RealType, double>::value) {
          EXPECT_DOUBLE_EQ(v, vH_(d)(0, 1, 0)(SWS::BACKWARD)[index]);
        } else {
          std::cerr << "Unsupported element type for SWS::RealType : " << typeid(SWS::RealType).name()
                    << "\n";
          exit(SWS::INVALID_REAL_TYPE);
        }
      }
    }
  }
}

TEST_F(HaloUpdateTest, UpdateForwardHalo)
{
  const auto d = SWS::X;

  const auto cx = Mesh3DPartitioning::getInstance()->ccx()[0];
  const auto cy = Mesh3DPartitioning::getInstance()->ccy()[0];
  const auto cz = Mesh3DPartitioning::getInstance()->ccz()[0];

  const auto hnx = CentralFDOperator::hnx();
  const auto hny = CentralFDOperator::hny();
  const auto hnz = CentralFDOperator::hnz();

  const auto jEnd = v_(0)(0, 0, 0).jEnd();

  HaloManager::extractVelocityHaloWrapper(SWS::BACKWARD, d, ts_, 0, 1, 0, vH_(d)(0, 1, 0)(SWS::BACKWARD));

  memcpy(vH_(d)(0, 0, 0)(SWS::FORWARD),
         vH_(d)(1, 0, 0)(SWS::BACKWARD),
         HaloManager::getInstance()->getHaloSize(SWS::FORWARD, 0, 0, 0) * sizeof(SWS::RealType));

  int iStartH = 0, iEndH = cx;
  int jStartH = 0, jEndH = hny;
  int kStartH = 0, kEndH = cz;

  int iShift = 0;
  int jShift = jEnd;
  int kShift = 0;

  Indexer<SWS::Ordering> indexer(iEndH, jEndH, kEndH);

  HaloManager::updateVelocityWrapper(SWS::FORWARD, d, ts_, 0, 0, 0, vH_(d)(0, 0, 0)(SWS::FORWARD));

  for (int i = iStartH; i < iEndH; i++) {
    for (int j = jStartH; j < jEndH; j++) {
      for (int k = kStartH; k < kEndH; k++) {

        auto index = indexer(i, j, k);

        const auto& v = v_(d)(0, 0, 0)(i + iShift, j + jShift, k + kShift);

        if constexpr (std::is_same<SWS::RealType, float>::value) {
          EXPECT_FLOAT_EQ(v, vH_(d)(0, 0, 0)(SWS::FORWARD)[index]);
        } else if constexpr (std::is_same<SWS::RealType, double>::value) {
          EXPECT_DOUBLE_EQ(v, vH_(d)(0, 0, 0)(SWS::FORWARD)[index]);
        } else {
          std::cerr << "Unsupported element type for SWS::RealType : " << typeid(SWS::RealType).name()
                    << "\n";
          exit(SWS::INVALID_REAL_TYPE);
        }
      }
    }
  }
}

TEST_F(HaloUpdateTest, UpdateBottomHalo)
{
  const auto d = SWS::X;

  const auto cx = Mesh3DPartitioning::getInstance()->ccx()[0];
  const auto cy = Mesh3DPartitioning::getInstance()->ccy()[0];
  const auto cz = Mesh3DPartitioning::getInstance()->ccz()[0];

  const auto hnx = CentralFDOperator::hnx();
  const auto hny = CentralFDOperator::hny();
  const auto hnz = CentralFDOperator::hnz();

  const auto kEnd = v_(0)(0, 0, 0).kEnd();

  HaloManager::extractVelocityHaloWrapper(SWS::TOP, d, ts_, 0, 0, 0, vH_(d)(0, 0, 0)(SWS::TOP));

  memcpy(vH_(d)(0, 0, 1)(SWS::BOTTOM),
         vH_(d)(0, 0, 0)(SWS::TOP),
         HaloManager::getInstance()->getHaloSize(SWS::TOP, 0, 0, 0) * sizeof(SWS::RealType));

  int iStartH = 0, iEndH = cx;
  int jStartH = 0, jEndH = cy;
  int kStartH = 0, kEndH = hnx;

  int iShift = 0;
  int jShift = 0;
  int kShift = 0;

  Indexer<SWS::Ordering> indexer(iEndH, jEndH, kEndH);

  HaloManager::updateVelocityWrapper(SWS::BOTTOM, d, ts_, 0, 0, 1, vH_(d)(0, 0, 1)(SWS::BOTTOM));

  for (int i = iStartH; i < iEndH; i++) {
    for (int j = jStartH; j < jEndH; j++) {
      for (int k = kStartH; k < kEndH; k++) {

        auto index = indexer(i, j, k);

        const auto& v = v_(d)(0, 0, 1)(i + iShift, j + jShift, k + kShift);

        if constexpr (std::is_same<SWS::RealType, float>::value) {
          EXPECT_FLOAT_EQ(v, vH_(d)(0, 0, 1)(SWS::BOTTOM)[index]);
        } else if constexpr (std::is_same<SWS::RealType, double>::value) {
          EXPECT_DOUBLE_EQ(v, vH_(d)(0, 0, 1)(SWS::BOTTOM)[index]);
        } else {
          std::cerr << "Unsupported element type for SWS::RealType : " << typeid(SWS::RealType).name()
                    << "\n";
          exit(SWS::INVALID_REAL_TYPE);
        }
      }
    }
  }
}

TEST_F(HaloUpdateTest, UpdateTopHalo)
{
  const auto d = SWS::X;

  const auto cx = Mesh3DPartitioning::getInstance()->ccx()[0];
  const auto cy = Mesh3DPartitioning::getInstance()->ccy()[0];
  const auto cz = Mesh3DPartitioning::getInstance()->ccz()[0];

  const auto hnx = CentralFDOperator::hnx();
  const auto hny = CentralFDOperator::hny();
  const auto hnz = CentralFDOperator::hnz();

  const auto kEnd = v_(0)(0, 0, 0).kEnd();

  HaloManager::extractVelocityHaloWrapper(SWS::BOTTOM, d, ts_, 0, 0, 1, vH_(d)(0, 0, 1)(SWS::BOTTOM));

  memcpy(vH_(d)(0, 0, 0)(SWS::TOP),
         vH_(d)(0, 0, 1)(SWS::BOTTOM),
         HaloManager::getInstance()->getHaloSize(SWS::TOP, 0, 0, 0) * sizeof(SWS::RealType));

  int iStartH = 0, iEndH = cx;
  int jStartH = 0, jEndH = cy;
  int kStartH = 0, kEndH = hnz;

  int iShift = 0;
  int jShift = 0;
  int kShift = kEnd;

  Indexer<SWS::Ordering> indexer(iEndH, jEndH, kEndH);

  HaloManager::updateVelocityWrapper(SWS::TOP, d, ts_, 0, 0, 0, vH_(d)(0, 0, 0)(SWS::TOP));

  for (int i = iStartH; i < iEndH; i++) {
    for (int j = jStartH; j < jEndH; j++) {
      for (int k = kStartH; k < kEndH; k++) {

        auto index = indexer(i, j, k);

        const auto& v = v_(d)(0, 0, 0)(i + iShift, j + jShift, k + kShift);

        if constexpr (std::is_same<SWS::RealType, float>::value) {
          EXPECT_FLOAT_EQ(v, vH_(d)(0, 0, 0)(SWS::TOP)[index]);
        } else if constexpr (std::is_same<SWS::RealType, double>::value) {
          EXPECT_DOUBLE_EQ(v, vH_(d)(0, 0, 0)(SWS::TOP)[index]);
        } else {
          std::cerr << "Unsupported element type for SWS::RealType : " << typeid(SWS::RealType).name()
                    << "\n";
          exit(SWS::INVALID_REAL_TYPE);
        }
      }
    }
  }
}
