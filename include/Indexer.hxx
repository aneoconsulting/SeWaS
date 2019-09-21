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

enum Orderings
{
    X_MAJOR,
    Y_MAJOR,
    Z_MAJOR
};

template<short Ordering>
class Indexer
{
public:
    Indexer(const int &nx, const int &ny, const int &nz) : nx_(nx),
                                                           ny_(ny),
                                                           nz_(nz)
    {
    }

    ~Indexer()
    {
    }
    
    inline auto operator()(const int i, const int j, const int k) const
    {
        static_assert(Orderings::X_MAJOR == Ordering || Orderings::Y_MAJOR == Ordering || Orderings::Z_MAJOR == Ordering,
                      "Unrecognized storage ordering");

        if constexpr (Orderings::X_MAJOR == Ordering)
        {
            return (j * nx_ + i) * nz_ + k; // x-major
        }
        else if constexpr (Orderings::Y_MAJOR == Ordering)
        {
            return (i * ny_ + j) * nz_ + k; // y-major
        }
        else if constexpr (Orderings::Z_MAJOR == Ordering)
        {
            return (k * ny_ + j) * nx_ + i; // z-major
        }
        else
        {
            return -1;
        }
    }

private:
    const int &nx_;
    const int &ny_;
    const int &nz_;
};