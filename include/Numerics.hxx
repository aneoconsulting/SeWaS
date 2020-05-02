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

#include <cmath>
#include <execution>
#include <functional>

#include "Config.hxx"
#include "Tiles.hxx"

namespace SWS {
namespace Numerics {
static inline auto
Norm(const SWS::SpatialField& f)
{
  Tiles tiles(f.dimension(SWS::X), f.dimension(SWS::Y), f.dimension(SWS::Z));

  auto norm2 = std::transform_reduce(
    std::execution::par, tiles().begin(), tiles().end(), 0.0, std::plus<>(), [&](const Tile& tile) {
      auto [ii, jj, kk] = tile;
      return f(ii, jj, kk).norm2();
    });

  return std::sqrt(norm2);
}
}
}