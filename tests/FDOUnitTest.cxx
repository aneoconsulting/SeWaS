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

#include <cmath>
#include <type_traits>

#include <gtest/gtest.h>

#include "FDOUnitTest.hxx"

TEST_F(FDOTest, PartialDerivativeOverX)
{
  rijk_(i0_, j0_) += fdo_->apply<SWS::X>(fijk_, i0_, j0_, SWS::DK);

  for (int k = fijk_.kStart(); k < fijk_.kEnd(); k++) {

    const auto r = (c1_ * fijk_(i0_ + 1, j0_, k) + c2_ * fijk_(i0_ - 2, j0_, k) + c3_ * fijk_(i0_, j0_, k) +
                    c4_ * fijk_(i0_ - 1, j0_, k)) /
                   ds_;

    if constexpr (std::is_same<SWS::RealType, float>::value) {
      EXPECT_FLOAT_EQ(rijk_(i0_, j0_, k), r) << " k = " << k;
    } else if constexpr (std::is_same<SWS::RealType, double>::value) {
      EXPECT_DOUBLE_EQ(rijk_(i0_, j0_, k), r) << " k = " << k;
    } else {
      std::cerr << "Unsupported element type for SWS::RealType : " << typeid(SWS::RealType).name() << "\n";
      exit(SWS::INVALID_REAL_TYPE);
    }
  }
}

TEST_F(FDOTest, PartialDerivativeOverY)
{
  rijk_(i0_, j0_) += fdo_->apply<SWS::Y>(fijk_, i0_, j0_, SWS::DK);

  for (int k = fijk_.kStart(); k < fijk_.kEnd(); k++) {

    const auto r = (c1_ * fijk_(i0_, j0_ + 1, k) + c2_ * fijk_(i0_, j0_ - 2, k) + c3_ * fijk_(i0_, j0_, k) +
                    c4_ * fijk_(i0_, j0_ - 1, k)) /
                   ds_;

    if constexpr (std::is_same<SWS::RealType, float>::value) {
      EXPECT_FLOAT_EQ(rijk_(i0_, j0_, k), r) << " k = " << k;
    } else if constexpr (std::is_same<SWS::RealType, double>::value) {
      EXPECT_DOUBLE_EQ(rijk_(i0_, j0_, k), r) << " k = " << k;
    } else {
      std::cerr << "Unsupported element type for SWS::RealType : " << typeid(SWS::RealType).name() << "\n";
      exit(SWS::INVALID_REAL_TYPE);
    }
  }
}

TEST_F(FDOTest, PartialDerivativeOverZ)
{
  rijk_(i0_, j0_) += fdo_->apply<SWS::Z>(fijk_, i0_, j0_, SWS::DK);

  for (int k = fijk_.kStart(); k < fijk_.kEnd(); k++) {

    const auto r = (c1_ * fijk_(i0_, j0_, k + 1) + c2_ * fijk_(i0_, j0_, k - 2) + c3_ * fijk_(i0_, j0_, k) +
                    c4_ * fijk_(i0_, j0_, k - 1)) /
                   ds_;

    if constexpr (std::is_same<SWS::RealType, float>::value) {
      EXPECT_FLOAT_EQ(rijk_(i0_, j0_, k), r) << " k = " << k;
    } else if constexpr (std::is_same<SWS::RealType, double>::value) {
      EXPECT_DOUBLE_EQ(rijk_(i0_, j0_, k), r) << " k = " << k;
    } else {
      std::cerr << "Unsupported element type for SWS::RealType : " << typeid(SWS::RealType).name() << "\n";
      exit(SWS::INVALID_REAL_TYPE);
    }
  }
}

TEST_F(FDOTest, PartialDerivativeOverZ_LegacyKernelVariant)
{
  rijk_(i0_, j0_) += fdo_->applyZ<SWS::KernelVariant::Legacy>(fijk_, i0_, j0_, SWS::DK);

  for (int k = fijk_.kStart(); k < fijk_.kEnd(); k++) {

    const auto r = (c1_ * fijk_(i0_, j0_, k + 1) + c2_ * fijk_(i0_, j0_, k - 2) + c3_ * fijk_(i0_, j0_, k) +
                    c4_ * fijk_(i0_, j0_, k - 1)) /
                   ds_;

    if constexpr (std::is_same<SWS::RealType, float>::value) {
      EXPECT_FLOAT_EQ(rijk_(i0_, j0_, k), r) << " k = " << k;
    } else if constexpr (std::is_same<SWS::RealType, double>::value) {
      EXPECT_DOUBLE_EQ(rijk_(i0_, j0_, k), r) << " k = " << k;
    } else {
      std::cerr << "Unsupported element type for SWS::RealType : " << typeid(SWS::RealType).name() << "\n";
      exit(SWS::INVALID_REAL_TYPE);
    }
  }
}

TEST_F(FDOTest, PartialDerivativeOverZ_UnalignedKernelVariant)
{
  rijk_(i0_, j0_) += fdo_->applyZ<SWS::KernelVariant::Unaligned>(fijk_, i0_, j0_, SWS::DK);

  for (int k = fijk_.kStart(); k < fijk_.kEnd(); k++) {

    const auto r = (c1_ * fijk_(i0_, j0_, k + 1) + c2_ * fijk_(i0_, j0_, k - 2) + c3_ * fijk_(i0_, j0_, k) +
                    c4_ * fijk_(i0_, j0_, k - 1)) /
                   ds_;

    if constexpr (std::is_same<SWS::RealType, float>::value) {
      EXPECT_FLOAT_EQ(rijk_(i0_, j0_, k), r) << " k = " << k;
    } else if constexpr (std::is_same<SWS::RealType, double>::value) {
      EXPECT_DOUBLE_EQ(rijk_(i0_, j0_, k), r) << " k = " << k;
    } else {
      std::cerr << "Unsupported element type for SWS::RealType : " << typeid(SWS::RealType).name() << "\n";
      exit(SWS::INVALID_REAL_TYPE);
    }
  }
}

/* CentralFDOperator::tap<d>() (raw kernel variant) must reproduce apply<d>/applyZ for every
 * direction and both stencil shifts, on a non-polynomial field that also fills the halos. */
TEST_F(FDOTest, RawTapMatchesApplyAllDirectionsAndShifts)
{
  for (int i = 0; i < cx_; i++)
    for (int j = 0; j < cy_; j++)
      for (int k = 0; k < cz_; k++)
        fijk_(i, j, k) =
          static_cast<SWS::RealType>(std::sin(0.7 * i + 0.3 * j) * std::cos(0.5 * k) + 0.1 * i * j);

  const long o = fijk_.index(0, 0, 0);
  const long sx = fijk_.index(1, 0, 0) - o, sy = fijk_.index(0, 1, 0) - o, sz = fijk_.index(0, 0, 1) - o;
  const SWS::RealType* p = fijk_.data().data();

  const double tol = std::is_same<SWS::RealType, float>::value ? 1e-4 : 1e-11;

  for (short shift = 0; shift <= 1; shift++) {
    for (int i = fijk_.iStart(); i < fijk_.iEnd(); i++) {
      for (int j = fijk_.jStart(); j < fijk_.jEnd(); j++) {
        SWS::SpatialBlockField1D rx = fdo_->apply<SWS::X>(fijk_, i + shift, j, SWS::DK);
        SWS::SpatialBlockField1D ry = fdo_->apply<SWS::Y>(fijk_, i, j + shift, SWS::DK);
        SWS::SpatialBlockField1D rz =
          fdo_->applyZ<SWS::KernelVariant::Unaligned>(fijk_, i, j, SWS::DK + shift);

        for (int k = fijk_.kStart(); k < fijk_.kEnd(); k++) {
          const auto c = fijk_.index(i, j, k);
          const auto kk = k - fijk_.kStart();
          EXPECT_NEAR(fdo_->tap<SWS::X>(p + c, sx, shift), rx(kk), tol * std::abs(rx(kk)) + tol)
            << "X shift " << shift << " at " << i << "," << j << "," << k;
          EXPECT_NEAR(fdo_->tap<SWS::Y>(p + c, sy, shift), ry(kk), tol * std::abs(ry(kk)) + tol)
            << "Y shift " << shift << " at " << i << "," << j << "," << k;
          EXPECT_NEAR(fdo_->tap<SWS::Z>(p + c, sz, shift), rz(kk), tol * std::abs(rz(kk)) + tol)
            << "Z shift " << shift << " at " << i << "," << j << "," << k;
        }
      }
    }
  }
}
