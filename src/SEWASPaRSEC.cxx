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

#ifdef SEWAS_WITH_PARSEC
#include <mpi.h>

#include <parsec/parsec_config.h>
#include <parsec/datatype.h>
#include <parsec/arena.h>

#include "Config.hxx"
#include "Mesh3DPartitioning.hxx"
#include "SEWASPaRSEC.hxx"
#include "LinearSeismicWaveModel.hxx"
#include "HaloManager.hxx"
#ifdef USE_VTK
#include "VisualizationManager.hxx"
#endif

#include "sewas.h"

extern parsec_context_t * g_parsec;

SEWASPaRSEC * SEWASPaRSEC::pInstance_ = nullptr;

SEWASPaRSEC * SEWASPaRSEC::getInstance(const int nt,
				     const int nxx, const int nyy, const int nzz)
{
  if (nullptr == pInstance_){
    pInstance_ = new SEWASPaRSEC(nt, nxx, nyy, nzz);
    return pInstance_;
  }
  else{
    return pInstance_;
  }
}

void SEWASPaRSEC::releaseInstance()
{
  if (pInstance_){
    delete pInstance_;
    pInstance_ = nullptr;
  }
}

int SEWASPaRSEC::run()
{
  int status=0;

  MPI_Barrier(MPI_COMM_WORLD);

  status=parsec_context_start(g_parsec);
  PARSEC_CHECK_ERROR(status, "parsec_context_start");

  MPI_Barrier(MPI_COMM_WORLD);

  status=parsec_context_wait(g_parsec);
  PARSEC_CHECK_ERROR(status, "parsec_context_wait");

  return status;
}

void SEWASPaRSEC::buildEngine()
{
  // TODO create and initialize the PaRSEC engine here
}

void SEWASPaRSEC::buildDataDescriptor()
{
  pDDesc_ = (parsec_data_collection_t *) malloc(sizeof(parsec_data_collection_t));

  parsec_data_collection_init(pDDesc_, world_, rank_);

  pDDesc_->rank_of  = Mesh3DPartitioning::rank_of;

  pDDesc_->vpid_of  = Mesh3DPartitioning::vpid_of;

  pDDesc_->data_of  = Mesh3DPartitioning::data_of;
  pDDesc_->data_key = Mesh3DPartitioning::data_key;
}

void SEWASPaRSEC::buildDAG()
{
  pDAG_ = (parsec_sewas_taskpool_t *) parsec_sewas_new(nt_,
                                                     nxx_, nyy_, nzz_,
                                                     (void *) &LinearSeismicWaveModel::initializeFieldsWrapper,
                                                     (void *) &LinearSeismicWaveModel::computeVelocityWrapper,
                                                     (void *) &HaloManager::extractVelocityHaloWrapper,
                                                     (void *) &HaloManager::updateVelocityWrapper,
                                                     (void *) &LinearSeismicWaveModel::computeStressWrapper,
                                                     (void *) &HaloManager::extractStressHaloWrapper,
                                                     (void *) &HaloManager::updateStressWrapper,
#ifdef USE_VTK
                                                     (void *) &VisualizationManager::displayVelocityWrapper,
                                                     (void *) &VisualizationManager::displayStressWrapper,
#else
                                                     nullptr,
                                                     nullptr,
#endif
                                                     pDDesc_);
  assert(nullptr != pDAG_);
}

void SEWASPaRSEC::enqueueDAG()
{
  int status = parsec_enqueue(g_parsec, (parsec_taskpool_t *) pDAG_);
  PARSEC_CHECK_ERROR(status, "parsec_enqueue");
}

void SEWASPaRSEC::addArena(const short arena_idx, const SWS::Locations l)
{
  parsec_datatype_t oldtype = SWS::MPIRealType;
  parsec_datatype_t newtype;
  ptrdiff_t lb, extent;

  size_t asize=1;
  if (SWS::NB_LOCATIONS != l){
    asize = HaloManager::getInstance()->getHaloSize(l, 0, 0, 0);
  }

  // TODO we assume that each spatial block of cells have the same halo size
  parsec_type_create_contiguous(asize, oldtype, &newtype);
  parsec_type_extent(newtype, &lb, &extent);
  parsec_arena_construct(pDAG_->arenas[arena_idx],
			 extent,
			 SWS::Alignment,
			 newtype);
}

SEWASPaRSEC::SEWASPaRSEC(const int nt,
		       const int nxx, const int nyy, const int nzz):nt_(nt),
								    nxx_(nxx),
								    nyy_(nyy),
								    nzz_(nzz)
{
  MPI_Comm_size(MPI_COMM_WORLD, &world_);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank_);

  buildEngine();
  buildDataDescriptor();
  buildDAG();

  /* Default halo arena */
  addArena(PARSEC_sewas_DEFAULT_ARENA);

  /* Stress field halo arenas */
  addArena(PARSEC_sewas_SLEFT_HALO_ARENA, SWS::LEFT);
  addArena(PARSEC_sewas_SRIGHT_HALO_ARENA, SWS::RIGHT);
  addArena(PARSEC_sewas_SBACKWARD_HALO_ARENA, SWS::BACKWARD);
  addArena(PARSEC_sewas_SFORWARD_HALO_ARENA, SWS::FORWARD);
  addArena(PARSEC_sewas_SBOTTOM_HALO_ARENA, SWS::BOTTOM);
  addArena(PARSEC_sewas_STOP_HALO_ARENA, SWS::TOP);

  /* Velocity halo arenas */
  addArena(PARSEC_sewas_VLEFT_HALO_ARENA, SWS::LEFT);
  addArena(PARSEC_sewas_VRIGHT_HALO_ARENA, SWS::RIGHT);
  addArena(PARSEC_sewas_VBACKWARD_HALO_ARENA, SWS::BACKWARD);
  addArena(PARSEC_sewas_VFORWARD_HALO_ARENA, SWS::FORWARD);
  addArena(PARSEC_sewas_VBOTTOM_HALO_ARENA, SWS::BOTTOM);
  addArena(PARSEC_sewas_VTOP_HALO_ARENA, SWS::TOP);

  enqueueDAG();
}

SEWASPaRSEC::~SEWASPaRSEC()
{
  parsec_taskpool_free((parsec_taskpool_t *) pDAG_);
  parsec_data_collection_destroy(pDDesc_);
}
#endif
