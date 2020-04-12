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

#include "ComputeVelocityUnitTest.hxx"

TEST_F(ComputeVelocityTest, ComputeVelocityX)
{
  auto d = SWS::X;

  int ts = 2;

  int ii = 0;
  int jj = 0;
  int kk = 0;

  // LinearSeismicWaveModel::computeVelocityWrapper(d, ts, ii, jj, kk);
}
