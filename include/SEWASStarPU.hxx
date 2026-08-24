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

#ifdef SEWAS_WITH_STARPU
#include <vector>

#include <starpu.h>
#include <starpu_mpi.h>

#include "Config.hxx"
#include "SEWASParameterManager.hxx"

/* SEWASStarPU drives the same two-phase-per-timestep algorithm as
 * SEWASSequential (see SEWASSequential.hxx), but submits one StarPU task per
 * local tile per phase (parallel across tiles within a node, via
 * --nthreads), and performs cross-rank halo exchange with the explicit
 * starpu_mpi_isend_detached/irecv_detached API, which -- unlike PaRSEC and
 * StarPU's own starpu_mpi_task_insert API -- does not exhibit a cross-rank
 * deadlock under this access pattern. Same-rank neighbor tiles are handled
 * with a plain memcpy, exactly like SEWASSequential::sendreceive(). */
class SEWASStarPU
{
public:
  static SEWASStarPU* getInstance(const int nt = 1, const int nxx = 1, const int nyy = 1, const int nzz = 1);
  static void releaseInstance();

  static void init(SEWASParameterManager& pm);
  static void finalize();

  int run();

private:
  enum DataTransferStages
  {
    SEND,
    RECEIVE
  };

  SEWASStarPU(const int nt, const int nxx, const int nyy, const int nzz);
  ~SEWASStarPU();

  void buildVelocityArenas();
  void buildStressArenas();
  void destroyVelocityArenas();
  void destroyStressArenas();

  void exchangeVelocityHalo(const SWS::Directions d, const int ts);
  void exchangeStressHalo(const SWS::StressFieldComponents sc, const int ts);

  size_t haloHandleIndex(const int lii, const int ljj, const int lkk, const SWS::Locations l) const noexcept;

  starpu_mpi_tag_t makeTag(const bool isStress,
                           const int comp,
                           const int ts,
                           const int ii,
                           const int jj,
                           const int kk,
                           const SWS::Locations l) const noexcept;

  bool neighborCoords(const int ii,
                      const int jj,
                      const int kk,
                      const SWS::Locations l,
                      int& nii,
                      int& njj,
                      int& nkk) const noexcept;

  static SWS::Locations opposite(const SWS::Locations l) noexcept;

  static void initializeFieldsCodelet(void* buffers[], void* cl_arg);
  static void phaseAVelocityCodelet(void* buffers[], void* cl_arg);
  static void phaseBStressCodelet(void* buffers[], void* cl_arg);

  static SEWASStarPU* pInstance_;

  int world_;
  int rank_;

  int nt_;
  int nxx_;
  int nyy_;
  int nzz_;

  int lnxx_;
  int lnyy_;
  int lnzz_;

  /* Persistent send/receive halo buffers, sized to the local tile grid
   * (mirrors SEWASSequential's own arenas, and HaloManager's own sigmaH_/vH_
   * sizing): sigmaH_[u](sc)(lii,ljj,lkk)(l), vH_[u](d)(lii,ljj,lkk)(l). */
  SWS::StressFieldHalo sigmaH_[2];
  SWS::VelocityHalo vH_[2];

  /* StarPU handles wrapping the buffers above, flattened per (component,
   * local tile, location); only used for cross-rank transfers -- same-rank
   * neighbors are handled with a plain memcpy and never touch these. */
  std::vector<starpu_data_handle_t> vSendHandles_[SWS::DIM];
  std::vector<starpu_data_handle_t> vRecvHandles_[SWS::DIM];
  std::vector<starpu_data_handle_t> sSendHandles_[SWS::NB_STRESS_FIELD_COMPONENTS];
  std::vector<starpu_data_handle_t> sRecvHandles_[SWS::NB_STRESS_FIELD_COMPONENTS];
};
#endif
