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

#include <tuple>
#include <vector>

namespace SWS {
using Tile = std::tuple<int, int, int>;

class Tiles
{
public:
  Tiles(const int nxx, const int nyy, const int nzz)
    : tiles_(std::vector<Tile>(nxx * nyy * nzz))
  {
    int index = 0;
    for (auto ii = 0; ii < nxx; ii++) {
      for (auto jj = 0; jj < nyy; jj++) {
        for (auto kk = 0; kk < nzz; kk++) {
          tiles_[index++] = std::make_tuple(ii, jj, kk);
        }
      }
    }
  }

  const auto& operator()() const { return tiles_; }

private:
  std::vector<Tile> tiles_;
};
}
