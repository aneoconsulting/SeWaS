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

#include "HaloExtractUnitTest.hxx"
#include "Indexer.hxx"

TEST_F(HaloExtractTest, ExtractLeftHalo)
{
  const auto d = SWS::X;
  const auto l = SWS::Locations::LEFT;

  const auto cx = Mesh3DPartitioning::getInstance()->ccx()[0];
  const auto cy = Mesh3DPartitioning::getInstance()->ccy()[0];
  const auto cz = Mesh3DPartitioning::getInstance()->ccz()[0];

  const auto hnx = CentralFDOperator::hnx();
  const auto hny = CentralFDOperator::hnx();
  const auto hnz = CentralFDOperator::hnx();

  int iStartH = 0, iEndH = hnx;
  int jStartH = 0, jEndH = cy;
  int kStartH = 0, kEndH = cz;

  int iShift = hnx;
  int jShift = 0;
  int kShift = 0;

  Indexer<SWS::Ordering> indexer(iEndH, jEndH, kEndH);

  auto& vH = vH_(d)(0, 0, 0)(l);

  HaloManager::extractVelocityHaloWrapper(l, d, ts_, 0, 0, 0, vH);

  for (int i = iStartH; i < iEndH; i++) {
    for (int j = jStartH; j < jEndH; j++) {
      for (int k = kStartH; k < kEndH; k++) {

        auto index = indexer(i, j, k);

        const auto& v = v_(d)(0, 0, 0)(i + iShift, j + jShift, k + kShift);

        if constexpr (std::is_same<SWS::RealType, float>::value) {
          EXPECT_FLOAT_EQ(v, vH[index]);
        } else if constexpr (std::is_same<SWS::RealType, double>::value) {
          EXPECT_DOUBLE_EQ(v, vH[index]);
        } else {
          std::cerr << "Unsupported element type for SWS::RealType : " << typeid(SWS::RealType).name()
                    << "\n";
          exit(SWS::INVALID_REAL_TYPE);
        }
      }
    }
  }
}

TEST_F(HaloExtractTest, ExtractRightHalo)
{
  const auto d = SWS::X;
  const auto l = SWS::Locations::RIGHT;

  const auto cx = Mesh3DPartitioning::getInstance()->ccx()[0];
  const auto cy = Mesh3DPartitioning::getInstance()->ccy()[0];
  const auto cz = Mesh3DPartitioning::getInstance()->ccz()[0];

  const auto hnx = CentralFDOperator::hnx();
  const auto hny = CentralFDOperator::hnx();
  const auto hnz = CentralFDOperator::hnx();

  const auto iEnd = v_(0)(0, 0, 0).iEnd();

  int iStartH = 0, iEndH = hnx;
  int jStartH = 0, jEndH = cy;
  int kStartH = 0, kEndH = cz;

  int iShift = iEnd - hnx;
  int jShift = 0;
  int kShift = 0;

  Indexer<SWS::Ordering> indexer(iEndH, jEndH, kEndH);

  auto& vH = vH_(d)(0, 0, 0)(l);

  HaloManager::extractVelocityHaloWrapper(l, d, ts_, 0, 0, 0, vH);

  for (int i = iStartH; i < iEndH; i++) {
    for (int j = jStartH; j < jEndH; j++) {
      for (int k = kStartH; k < kEndH; k++) {

        auto index = indexer(i, j, k);

        const auto& v = v_(d)(0, 0, 0)(i + iShift, j + jShift, k + kShift);

        if constexpr (std::is_same<SWS::RealType, float>::value) {
          EXPECT_FLOAT_EQ(v, vH[index]);
        } else if constexpr (std::is_same<SWS::RealType, double>::value) {
          EXPECT_DOUBLE_EQ(v, vH[index]);
        } else {
          std::cerr << "Unsupported element type for SWS::RealType : " << typeid(SWS::RealType).name()
                    << "\n";
          exit(SWS::INVALID_REAL_TYPE);
        }
      }
    }
  }
}

TEST_F(HaloExtractTest, ExtractBackwardHalo)
{
  const auto d = SWS::X;
  const auto l = SWS::Locations::BACKWARD;

  const auto cx = Mesh3DPartitioning::getInstance()->ccx()[0];
  const auto cy = Mesh3DPartitioning::getInstance()->ccy()[0];
  const auto cz = Mesh3DPartitioning::getInstance()->ccz()[0];

  const auto hnx = CentralFDOperator::hnx();
  const auto hny = CentralFDOperator::hnx();
  const auto hnz = CentralFDOperator::hnx();

  int iStartH = 0, iEndH = cx;
  int jStartH = 0, jEndH = hny;
  int kStartH = 0, kEndH = cz;

  int iShift = 0;
  int jShift = hny;
  int kShift = 0;

  Indexer<SWS::Ordering> indexer(iEndH, jEndH, kEndH);

  auto& vH = vH_(d)(0, 0, 0)(l);

  HaloManager::extractVelocityHaloWrapper(l, d, ts_, 0, 0, 0, vH);

  for (int i = iStartH; i < iEndH; i++) {
    for (int j = jStartH; j < jEndH; j++) {
      for (int k = kStartH; k < kEndH; k++) {

        auto index = indexer(i, j, k);

        const auto& v = v_(d)(0, 0, 0)(i + iShift, j + jShift, k + kShift);

        if constexpr (std::is_same<SWS::RealType, float>::value) {
          EXPECT_FLOAT_EQ(v, vH[index]);
        } else if constexpr (std::is_same<SWS::RealType, double>::value) {
          EXPECT_DOUBLE_EQ(v, vH[index]);
        } else {
          std::cerr << "Unsupported element type for SWS::RealType : " << typeid(SWS::RealType).name()
                    << "\n";
          exit(SWS::INVALID_REAL_TYPE);
        }
      }
    }
  }
}

TEST_F(HaloExtractTest, ExtractForwardHalo)
{
  const auto d = SWS::X;
  const auto l = SWS::Locations::FORWARD;

  const auto cx = Mesh3DPartitioning::getInstance()->ccx()[0];
  const auto cy = Mesh3DPartitioning::getInstance()->ccy()[0];
  const auto cz = Mesh3DPartitioning::getInstance()->ccz()[0];

  const auto hnx = CentralFDOperator::hnx();
  const auto hny = CentralFDOperator::hnx();
  const auto hnz = CentralFDOperator::hnx();

  const auto jEnd = v_(0)(0, 0, 0).jEnd();

  int iStartH = 0, iEndH = cx;
  int jStartH = 0, jEndH = hny;
  int kStartH = 0, kEndH = cz;

  int iShift = 0;
  int jShift = jEnd - hny;
  int kShift = 0;

  Indexer<SWS::Ordering> indexer(iEndH, jEndH, kEndH);

  auto& vH = vH_(d)(0, 0, 0)(l);

  HaloManager::extractVelocityHaloWrapper(l, d, ts_, 0, 0, 0, vH);

  for (int i = iStartH; i < iEndH; i++) {
    for (int j = jStartH; j < jEndH; j++) {
      for (int k = kStartH; k < kEndH; k++) {

        auto index = indexer(i, j, k);

        const auto& v = v_(d)(0, 0, 0)(i + iShift, j + jShift, k + kShift);

        if constexpr (std::is_same<SWS::RealType, float>::value) {
          EXPECT_FLOAT_EQ(v, vH[index]);
        } else if constexpr (std::is_same<SWS::RealType, double>::value) {
          EXPECT_DOUBLE_EQ(v, vH[index]);
        } else {
          std::cerr << "Unsupported element type for SWS::RealType : " << typeid(SWS::RealType).name()
                    << "\n";
          exit(SWS::INVALID_REAL_TYPE);
        }
      }
    }
  }
}

TEST_F(HaloExtractTest, ExtractBottomHalo)
{
  const auto d = SWS::X;
  const auto l = SWS::Locations::BOTTOM;

  const auto cx = Mesh3DPartitioning::getInstance()->ccx()[0];
  const auto cy = Mesh3DPartitioning::getInstance()->ccy()[0];
  const auto cz = Mesh3DPartitioning::getInstance()->ccz()[0];

  const auto hnx = CentralFDOperator::hnx();
  const auto hny = CentralFDOperator::hnx();
  const auto hnz = CentralFDOperator::hnx();

  int iStartH = 0, iEndH = cx;
  int jStartH = 0, jEndH = cy;
  int kStartH = 0, kEndH = hnz;

  int iShift = 0;
  int jShift = 0;
  int kShift = hnz;

  Indexer<SWS::Ordering> indexer(iEndH, jEndH, kEndH);

  auto& vH = vH_(d)(0, 0, 0)(l);

  HaloManager::extractVelocityHaloWrapper(l, d, ts_, 0, 0, 0, vH);

  for (int i = iStartH; i < iEndH; i++) {
    for (int j = jStartH; j < jEndH; j++) {
      for (int k = kStartH; k < kEndH; k++) {

        auto index = indexer(i, j, k);

        const auto& v = v_(d)(0, 0, 0)(i + iShift, j + jShift, k + kShift);

        if constexpr (std::is_same<SWS::RealType, float>::value) {
          EXPECT_FLOAT_EQ(v, vH[index]);
        } else if constexpr (std::is_same<SWS::RealType, double>::value) {
          EXPECT_DOUBLE_EQ(v, vH[index]);
        } else {
          std::cerr << "Unsupported element type for SWS::RealType : " << typeid(SWS::RealType).name()
                    << "\n";
          exit(SWS::INVALID_REAL_TYPE);
        }
      }
    }
  }
}

TEST_F(HaloExtractTest, ExtractTopHalo)
{
  const auto d = SWS::X;
  const auto l = SWS::Locations::TOP;

  const auto cx = Mesh3DPartitioning::getInstance()->ccx()[0];
  const auto cy = Mesh3DPartitioning::getInstance()->ccy()[0];
  const auto cz = Mesh3DPartitioning::getInstance()->ccz()[0];

  const auto hnx = CentralFDOperator::hnx();
  const auto hny = CentralFDOperator::hnx();
  const auto hnz = CentralFDOperator::hnx();

  const auto kEnd = v_(0)(0, 0, 0).kEnd();

  int iStartH = 0, iEndH = cx;
  int jStartH = 0, jEndH = cy;
  int kStartH = 0, kEndH = hnz;

  int iShift = 0;
  int jShift = 0;
  int kShift = kEnd - hnz;

  Indexer<SWS::Ordering> indexer(iEndH, jEndH, kEndH);

  auto& vH = vH_(d)(0, 0, 0)(l);

  HaloManager::extractVelocityHaloWrapper(l, d, ts_, 0, 0, 0, vH);

  for (int i = iStartH; i < iEndH; i++) {
    for (int j = jStartH; j < jEndH; j++) {
      for (int k = kStartH; k < kEndH; k++) {

        auto index = indexer(i, j, k);

        const auto& v = v_(d)(0, 0, 0)(i + iShift, j + jShift, k + kShift);

        if constexpr (std::is_same<SWS::RealType, float>::value) {
          EXPECT_FLOAT_EQ(v, vH[index]);
        } else if constexpr (std::is_same<SWS::RealType, double>::value) {
          EXPECT_DOUBLE_EQ(v, vH[index]);
        } else {
          std::cerr << "Unsupported element type for SWS::RealType : " << typeid(SWS::RealType).name()
                    << "\n";
          exit(SWS::INVALID_REAL_TYPE);
        }
      }
    }
  }
}
