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

#include <deque>
#include <vector>

#include <Eigen/Core>
#include <unsupported/Eigen/CXX11/Tensor>

#ifdef USE_VTK
#include <vtkActor.h>
#endif

#include "Constants.hxx"
#include "SpatialBlockField.hxx"
#include "Types.hxx"

#include "OS.hxx"

namespace SWS {
constexpr int packetSize = Eigen::internal::packet_traits<RealType>::size;

// (i, j, k)
// TODO since Eigen::Tensor is unstable, using our internal SpatialBlockField
// class for the moment typedef Eigen::Tensor<RealType, DIM, Eigen::AutoAlign>
// SpatialBlockField;

// (ii, jj, kk)(i, j, k)
typedef Eigen::Tensor<SpatialBlockField<RealType>, DIM, Eigen::AutoAlign> SpatialField;

// // (ts, ii, jj, kk)(i, j, k)
// typedef Eigen::Tensor<SpatialBlockField<RealType>, 1+DIM, Eigen::AutoAlign>
// TimeSpatialField;

// (x|y|z)(ii, jj, kk)(i, j, k)
typedef Eigen::Array<SpatialField, DIM, 1, Eigen::ColMajor> Velocity;

// (xx|yy|zz|xy|xz|yz)(ii, jj, kk)(i, j, k)
typedef Eigen::Array<SpatialField, NB_STRESS_FIELD_COMPONENTS, 1, Eigen::ColMajor> StressField;

// (left|right|bottom|top|backward|forward)(index)
typedef Eigen::Array<RealType*, NB_LOCATIONS, 1, Eigen::ColMajor | Eigen::AutoAlign> Halo;
// typedef Map<Eigen::Array<RealType, NB_LOCATIONS, Eigen::Dynamic,
// Eigen::ColMajor|Eigen::AutoAlign>, Eigen::Aligned64> Halo;

// (ii,jj,kk)(left|right|bottom|top|backward|forward)(index)
typedef Eigen::Tensor<Halo, DIM, Eigen::AutoAlign> SpatialFieldHalo;

// (x|y|z)(ii,jj,kk)(left|right|bottom|top|backward|forward)(index)
typedef Eigen::Array<SpatialFieldHalo, DIM, 1, Eigen::ColMajor> VelocityHalo;

// (xx|yy|zz|xy|xz|yz)(ii,jj,kk)(left|right|bottom|top|backward|forward)(index)
typedef Eigen::Array<SpatialFieldHalo, NB_STRESS_FIELD_COMPONENTS, 1, Eigen::ColMajor> StressFieldHalo;

// (x|y|z)(ii, jj, kk)(i, j, k)
typedef Eigen::Array<SpatialField, DIM, 1, Eigen::ColMajor> SourceForce;

// (x|y|z)(ii, jj, kk)(i, j, k)
typedef Eigen::Array<SpatialField, DIM, 1, Eigen::ColMajor> Buoyancy;

// (xx|yy|zz|xy|xz|yz)(ii, jj, kk)(i, j, k)
typedef Eigen::Array<SpatialField, NB_STRESS_FIELD_COMPONENTS, 1, Eigen::ColMajor> Elasticity;

// (ii, jj, kk)(i, j, k)
typedef SpatialField MaterialDensity;

// (ii, jj, kk)
typedef Eigen::Tensor<int, DIM, Eigen::AutoAlign> PriorityField;

// (taskType)(ii, jj, kk)
typedef Eigen::Array<PriorityField, 1, NB_TASK_TYPES, Eigen::RowMajor | Eigen::AutoAlign> TaskPriorities;

// (taskType)
typedef Eigen::Array<int, 1, NB_TASK_TYPES, Eigen::RowMajor | Eigen::AutoAlign> TaskGroups;

#ifdef USE_VTK
// Data structure for storing VTK actors
// (ii,jj,kk)(actor)
typedef Eigen::Tensor<vtkActor*, DIM, Eigen::AutoAlign> VTKSpatialActors;

// (x|y|z)(ii,jj,kk)(actor)
typedef Eigen::Array<VTKSpatialActors, DIM, 1, Eigen::ColMajor> VTKVelocityActors;

// (xx|yy|zz|xy|xz|yz)(ii,jj,kk)(actor)
typedef Eigen::Array<VTKSpatialActors, NB_STRESS_FIELD_COMPONENTS, 1, Eigen::ColMajor> VTKStressActors;

// (node,core)(actor)
typedef Eigen::Array<vtkActor*, NB_CLUSTER_NODES, NB_CORES, Eigen::ColMajor> VTKClusterActors;
#endif
}
