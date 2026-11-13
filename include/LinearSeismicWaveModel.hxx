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

#include <iostream>

#include <Eigen/Core>

#include "CentralFDOperator.hxx"
#include "Config.hxx"
#include "HaloManager.hxx"
#include "SEWASParameterManager.hxx"

class LinearSeismicWaveModel
{
public:
  static LinearSeismicWaveModel* getInstance(const CentralFDOperator& fdo,
                                             const SEWASParameterManager& pm);
  static LinearSeismicWaveModel* getInstance();

  static void releaseInstance();

  int propagate() noexcept;

  inline const auto& v(const SWS::Directions& d) const { return v_(d); }

  inline const auto& sigma(const SWS::StressFieldComponents& sc) const { return sigma_(sc); }

  inline const auto& nt() const { return nt_; }

  /* Entry points called by the PaRSEC runtime when executing tasks' BODY */
  static int initializeFieldsWrapper(const int ii, const int jj, const int kk);

  static int computeVelocityWrapper(const int d, const int ts, const int ii, const int jj, const int kk);

  static int computeStressWrapper(const int sc, const int ts, const int ii, const int jj, const int kk);

  /* Tile-level entry points used by the sequential and StarPU engines: update every
   * velocity (resp. stress) component of tile (ii, jj, kk). With --fuse-fields, all
   * components are updated in a single (i, j) sweep (and the diagonal stresses share
   * Dx(vX), Dy(vY), Dz(vZ)); otherwise this loops over computeVelocity/computeStress, one
   * full-tile sweep per component. PaRSEC keeps the per-component wrappers above. */
  static int computeVelocityTileWrapper(const int ts, const int ii, const int jj, const int kk);

  static int computeStressTileWrapper(const int ts, const int ii, const int jj, const int kk);

private:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW // To force an object of this class to be allocated as aligned

  LinearSeismicWaveModel(const CentralFDOperator& fdo,
                         const SEWASParameterManager& pm);
  ~LinearSeismicWaveModel();

  void initialize(const int ii, const int jj, const int kk);

  void setVelocitySourceLocations();

  void addVelocitySource(const SWS::Directions& d, const int ts, const int ii, const int jj, const int kk);

  /* Templated on the Dz kernel variant (see CentralFDOperator::applyZ<Variant>()) so the
   * --kernel-variant choice is resolved once per (component, tile, phase) call in
   * computeVelocityWrapper()/computeStressWrapper() below, not per grid cell. */
  template<SWS::KernelVariant Variant>
  int computeVelocity(const SWS::Directions& d, const int& ts, const int& ii, const int& jj, const int& kk);

  template<SWS::KernelVariant Variant>
  int computeStress(const SWS::StressFieldComponents& sc,
                    const int& ts,
                    const int& ii,
                    const int& jj,
                    const int& kk);

  /* Single-sweep variants of the above over all components of one tile (Phase 1, item 0 of
   * sewas_transformation_plan_en.md). Only defined with BLOCKWISE_FDO. */
  template<SWS::KernelVariant Variant>
  int computeVelocityFused(const int ts, const int ii, const int jj, const int kk);

  template<SWS::KernelVariant Variant>
  int computeStressFused(const int ts, const int ii, const int jj, const int kk);

  /* `raw` kernel variant (Phase 1.1b): single-sweep raw-pointer tile kernels over all
   * components, built on CentralFDOperator::tap<d>(). Independent of BLOCKWISE_FDO and
   * --fuse-fields. */
  int computeVelocityRaw(const int ts, const int ii, const int jj, const int kk);

  int computeStressRaw(const int ts, const int ii, const int jj, const int kk);

  template<SWS::KernelVariant Variant>
  int computeVelocityTile(const int ts, const int ii, const int jj, const int kk);

  template<SWS::KernelVariant Variant>
  int computeStressTile(const int ts, const int ii, const int jj, const int kk);

  static LinearSeismicWaveModel* pInstance_;

  const CentralFDOperator& fdo_;
  const SEWASParameterManager& pm_;

  int nt_;
  float tmax_;

  SWS::RealType dt_;

  SWS::Velocity v_;

  SWS::StressField sigma_;
};
