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

#include "FDOUnitTest.hxx"

TEST_F(FDOTest, PartialDerivativeOverX)
{
  rijk_(i0_, j0_) += fdo_->apply<SWS::X>(fijk_, i0_, j0_, SWS::DK);

  for (int k = fijk_.kStart(); k < fijk_.kEnd(); k++) {

    const auto r =
      (c1_ * fijk_(i0_ + 1, j0_, k) + c2_ * fijk_(i0_ - 2, j0_, k) +
       c3_ * fijk_(i0_, j0_, k) + c4_ * fijk_(i0_ - 1, j0_, k)) /
      ds_;

    if constexpr (std::is_same<SWS::RealType, float>::value) {
      EXPECT_FLOAT_EQ(rijk_(i0_, j0_, k), r) << " k = " << k;
    } else if constexpr (std::is_same<SWS::RealType, double>::value) {
      EXPECT_DOUBLE_EQ(rijk_(i0_, j0_, k), r) << " k = " << k;
    } else {
      std::cerr << "Unsupported element type for SWS::RealType : "
                << typeid(SWS::RealType).name() << "\n";
      exit(SWS::INVALID_REAL_TYPE);
    }
  }
}

TEST_F(FDOTest, PartialDerivativeOverY)
{
  rijk_(i0_, j0_) += fdo_->apply<SWS::Y>(fijk_, i0_, j0_, SWS::DK);

  for (int k = fijk_.kStart(); k < fijk_.kEnd(); k++) {

    const auto r =
      (c1_ * fijk_(i0_, j0_ + 1, k) + c2_ * fijk_(i0_, j0_ - 2, k) +
       c3_ * fijk_(i0_, j0_, k) + c4_ * fijk_(i0_, j0_ - 1, k)) /
      ds_;

    if constexpr (std::is_same<SWS::RealType, float>::value) {
      EXPECT_FLOAT_EQ(rijk_(i0_, j0_, k), r) << " k = " << k;
    } else if constexpr (std::is_same<SWS::RealType, double>::value) {
      EXPECT_DOUBLE_EQ(rijk_(i0_, j0_, k), r) << " k = " << k;
    } else {
      std::cerr << "Unsupported element type for SWS::RealType : "
                << typeid(SWS::RealType).name() << "\n";
      exit(SWS::INVALID_REAL_TYPE);
    }
  }
}

TEST_F(FDOTest, PartialDerivativeOverZ)
{
  rijk_(i0_, j0_) += fdo_->apply<SWS::Z>(fijk_, i0_, j0_, SWS::DK);

  for (int k = fijk_.kStart(); k < fijk_.kEnd(); k++) {

    const auto r =
      (c1_ * fijk_(i0_, j0_, k + 1) + c2_ * fijk_(i0_, j0_, k - 2) +
       c3_ * fijk_(i0_, j0_, k) + c4_ * fijk_(i0_, j0_, k - 1)) /
      ds_;

    if constexpr (std::is_same<SWS::RealType, float>::value) {
      EXPECT_FLOAT_EQ(rijk_(i0_, j0_, k), r) << " k = " << k;
    } else if constexpr (std::is_same<SWS::RealType, double>::value) {
      EXPECT_DOUBLE_EQ(rijk_(i0_, j0_, k), r) << " k = " << k;
    } else {
      std::cerr << "Unsupported element type for SWS::RealType : "
                << typeid(SWS::RealType).name() << "\n";
      exit(SWS::INVALID_REAL_TYPE);
    }
  }
}
