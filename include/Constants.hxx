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

namespace SWS
{
  enum Directions{
    X=0,
    Y,
    Z,
    DIM
  };

  enum StressFieldComponents{
    XX=0,
    YY,
    ZZ,
    XY,
    XZ,
    YZ,
    NB_STRESS_FIELD_COMPONENTS
  };

  enum Locations{
    LEFT=0,
    RIGHT,
    BACKWARD,
    FORWARD,
    BOTTOM,
    TOP,
    NB_LOCATIONS
  };


  enum DIJK {DI=0, DJ=0, DK=0}; // Identifiers to be used when applying the finite-difference operator

  enum ClusterNodes{MP_0, MP_1, MP_2, MP_3, MP_4, MP_5, MP_6, MP_7, MP_8, MP_9, MP_10, MP_11, MP_12, MP_13, MP_14, MP_15, MP_16, MP_17, NB_CLUSTER_NODES}; // Cluster nodes
  enum Cores{C0, C1, C2, C3, NB_CORES}; // CPU cores

  enum ColorizationStrategies{CORE, TIME_STEP};
}
