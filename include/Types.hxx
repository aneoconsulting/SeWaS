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

#include <type_traits>

#ifdef SEWAS_WITH_PARSEC
#include <parsec/parsec_config.h>
#include <parsec/datatype.h>
#endif

#include <Eigen/Core>

namespace SWS
{
  using RealType=float;

#ifdef SEWAS_WITH_PARSEC
  constexpr auto MPIRealType=parsec_datatype_float_t;
#endif

  // Data alignment
  constexpr int Alignment=64;

  // (k)
  typedef Eigen::Array<RealType, 1, Eigen::Dynamic, Eigen::RowMajor|Eigen::AutoAlign> SpatialBlockField1D;
}
