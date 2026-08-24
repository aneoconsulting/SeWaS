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

#include "SEWASStarPU.hxx"

#ifdef SEWAS_WITH_STARPU

#include <algorithm>
#include <cstring>

#include "ExecutionContext.hxx"
#include "HaloManager.hxx"
#include "LinearSeismicWaveModel.hxx"
#include "LogManager.hxx"
#include "Mesh3DPartitioning.hxx"

namespace {
struct starpu_codelet
makeCodelet(void (*func)(void*[], void*), const char* name)
{
  struct starpu_codelet cl;
  memset(&cl, 0, sizeof(cl));
  cl.cpu_funcs[0] = func;
  cl.cpu_funcs_name[0] = name;
  cl.nbuffers = 0;
  return cl;
}

/* Clamp into StarPU's currently valid priority range: the active scheduling
 * policy may support only a single priority level (or a narrower range than
 * MinimumCommunicationPriorityEvaluator's raw output), and starpu_task_insert
 * gives no other feedback for an out-of-range value. */
int
clampPriority(const int priority) noexcept
{
  const int minPrio = starpu_sched_get_min_priority();
  const int maxPrio = starpu_sched_get_max_priority();
  return (std::max)(minPrio, (std::min)(maxPrio, priority));
}
}

SEWASStarPU* SEWASStarPU::pInstance_ = nullptr;

SEWASStarPU*
SEWASStarPU::getInstance(const int nt, const int nxx, const int nyy, const int nzz)
{
  if (nullptr == pInstance_) {
    pInstance_ = new SEWASStarPU(nt, nxx, nyy, nzz);
  }
  return pInstance_;
}

void
SEWASStarPU::releaseInstance()
{
  if (pInstance_) {
    delete pInstance_;
    pInstance_ = nullptr;
  }
}

void
SEWASStarPU::init(SEWASParameterManager& pm)
{
  struct starpu_conf conf;
  starpu_conf_init(&conf);
  conf.ncpus = pm.nthreads();

  int status = starpu_init(&conf);
  if (status != 0) {
    LOG(SWS::LOG_CRITICAL, "starpu_init failed with status {}", status);
    exit(SWS::OBJECT_CREATION_FAILURE);
  }

  /* MPI is already initialized by ExecutionContext::init(); StarPU-MPI must not
   * re-initialize it -- same ownership model SEWASPaRSEC::init() follows for PaRSEC. */
  status = starpu_mpi_init_comm(&pm.argc(), &pm.argv(), 0, MPI_COMM_WORLD);
  if (status != 0) {
    LOG(SWS::LOG_CRITICAL, "starpu_mpi_init_comm failed with status {}", status);
    exit(SWS::OBJECT_CREATION_FAILURE);
  }
}

void
SEWASStarPU::finalize()
{
  starpu_mpi_shutdown();
  starpu_shutdown();
}

size_t
SEWASStarPU::haloHandleIndex(const int lii, const int ljj, const int lkk, const SWS::Locations l) const noexcept
{
  return (((size_t)lii * lnyy_ + ljj) * lnzz_ + lkk) * SWS::NB_LOCATIONS + l;
}

starpu_mpi_tag_t
SEWASStarPU::makeTag(const bool isStress,
                     const int comp,
                     const int ts,
                     const int ii,
                     const int jj,
                     const int kk,
                     const SWS::Locations l) const noexcept
{
  int64_t t = isStress ? 1 : 0;
  t = t * 8 + comp;
  t = t * SWS::NB_LOCATIONS + l;
  t = t * nxx_ + ii;
  t = t * nyy_ + jj;
  t = t * nzz_ + kk;
  t = t * (nt_ + 1) + ts;
  return (starpu_mpi_tag_t)t;
}

bool
SEWASStarPU::neighborCoords(const int ii,
                            const int jj,
                            const int kk,
                            const SWS::Locations l,
                            int& nii,
                            int& njj,
                            int& nkk) const noexcept
{
  nii = ii;
  njj = jj;
  nkk = kk;

  switch (l) {
    case SWS::LEFT:
      if (ii <= 0) return false;
      nii = ii - 1;
      break;
    case SWS::RIGHT:
      if (ii >= nxx_ - 1) return false;
      nii = ii + 1;
      break;
    case SWS::BACKWARD:
      if (jj <= 0) return false;
      njj = jj - 1;
      break;
    case SWS::FORWARD:
      if (jj >= nyy_ - 1) return false;
      njj = jj + 1;
      break;
    case SWS::BOTTOM:
      if (kk <= 0) return false;
      nkk = kk - 1;
      break;
    case SWS::TOP:
      if (kk >= nzz_ - 1) return false;
      nkk = kk + 1;
      break;
    default:
      LOG(SWS::LOG_ERROR,
          "Unknown location {} requested within SEWASStarPU::neighborCoords()",
          static_cast<int>(l));
      return false;
  }

  return true;
}

SWS::Locations
SEWASStarPU::opposite(const SWS::Locations l) noexcept
{
  switch (l) {
    case SWS::LEFT:
      return SWS::RIGHT;
    case SWS::RIGHT:
      return SWS::LEFT;
    case SWS::BACKWARD:
      return SWS::FORWARD;
    case SWS::FORWARD:
      return SWS::BACKWARD;
    case SWS::BOTTOM:
      return SWS::TOP;
    case SWS::TOP:
      return SWS::BOTTOM;
    default:
      LOG(SWS::LOG_ERROR, "Unknown location {} requested within SEWASStarPU::opposite()", static_cast<int>(l));
      return SWS::NB_LOCATIONS;
  }
}

void
SEWASStarPU::buildVelocityArenas()
{
  auto* pHaloManager = HaloManager::getInstance();

  for (auto d : { SWS::X, SWS::Y, SWS::Z }) {
    for (auto u : { SEND, RECEIVE }) {
      vH_[u](d).resize(lnxx_, lnyy_, lnzz_);
    }

    vSendHandles_[d].assign((size_t)lnxx_ * lnyy_ * lnzz_ * SWS::NB_LOCATIONS, nullptr);
    vRecvHandles_[d].assign((size_t)lnxx_ * lnyy_ * lnzz_ * SWS::NB_LOCATIONS, nullptr);

    for (int lii = 0; lii < lnxx_; lii++) {
      for (int ljj = 0; ljj < lnyy_; ljj++) {
        for (int lkk = 0; lkk < lnzz_; lkk++) {
          for (auto l : { SWS::LEFT, SWS::RIGHT, SWS::BACKWARD, SWS::FORWARD, SWS::BOTTOM, SWS::TOP }) {
            const auto hsize = pHaloManager->getHaloSize(l, lii, ljj, lkk);

            auto* sendBuf = new SWS::RealType[hsize];
            memset(sendBuf, 0, hsize * sizeof(SWS::RealType));
            vH_[SEND](d)(lii, ljj, lkk)(l) = sendBuf;

            auto* recvBuf = new SWS::RealType[hsize];
            memset(recvBuf, 0, hsize * sizeof(SWS::RealType));
            vH_[RECEIVE](d)(lii, ljj, lkk)(l) = recvBuf;

            const auto idx = haloHandleIndex(lii, ljj, lkk, l);
            starpu_vector_data_register(
              &vSendHandles_[d][idx], STARPU_MAIN_RAM, (uintptr_t)sendBuf, hsize, sizeof(SWS::RealType));
            starpu_vector_data_register(
              &vRecvHandles_[d][idx], STARPU_MAIN_RAM, (uintptr_t)recvBuf, hsize, sizeof(SWS::RealType));
          }
        }
      }
    }
  }
}

void
SEWASStarPU::buildStressArenas()
{
  auto* pHaloManager = HaloManager::getInstance();

  for (auto sc : { SWS::XX, SWS::YY, SWS::ZZ, SWS::XY, SWS::XZ, SWS::YZ }) {
    for (auto u : { SEND, RECEIVE }) {
      sigmaH_[u](sc).resize(lnxx_, lnyy_, lnzz_);
    }

    sSendHandles_[sc].assign((size_t)lnxx_ * lnyy_ * lnzz_ * SWS::NB_LOCATIONS, nullptr);
    sRecvHandles_[sc].assign((size_t)lnxx_ * lnyy_ * lnzz_ * SWS::NB_LOCATIONS, nullptr);

    for (int lii = 0; lii < lnxx_; lii++) {
      for (int ljj = 0; ljj < lnyy_; ljj++) {
        for (int lkk = 0; lkk < lnzz_; lkk++) {
          for (auto l : { SWS::LEFT, SWS::RIGHT, SWS::BACKWARD, SWS::FORWARD, SWS::BOTTOM, SWS::TOP }) {
            const auto hsize = pHaloManager->getHaloSize(l, lii, ljj, lkk);

            auto* sendBuf = new SWS::RealType[hsize];
            memset(sendBuf, 0, hsize * sizeof(SWS::RealType));
            sigmaH_[SEND](sc)(lii, ljj, lkk)(l) = sendBuf;

            auto* recvBuf = new SWS::RealType[hsize];
            memset(recvBuf, 0, hsize * sizeof(SWS::RealType));
            sigmaH_[RECEIVE](sc)(lii, ljj, lkk)(l) = recvBuf;

            const auto idx = haloHandleIndex(lii, ljj, lkk, l);
            starpu_vector_data_register(
              &sSendHandles_[sc][idx], STARPU_MAIN_RAM, (uintptr_t)sendBuf, hsize, sizeof(SWS::RealType));
            starpu_vector_data_register(
              &sRecvHandles_[sc][idx], STARPU_MAIN_RAM, (uintptr_t)recvBuf, hsize, sizeof(SWS::RealType));
          }
        }
      }
    }
  }
}

void
SEWASStarPU::destroyVelocityArenas()
{
  for (auto d : { SWS::X, SWS::Y, SWS::Z }) {
    for (int lii = 0; lii < lnxx_; lii++) {
      for (int ljj = 0; ljj < lnyy_; ljj++) {
        for (int lkk = 0; lkk < lnzz_; lkk++) {
          for (auto l : { SWS::LEFT, SWS::RIGHT, SWS::BACKWARD, SWS::FORWARD, SWS::BOTTOM, SWS::TOP }) {
            const auto idx = haloHandleIndex(lii, ljj, lkk, l);

            starpu_data_unregister(vSendHandles_[d][idx]);
            starpu_data_unregister(vRecvHandles_[d][idx]);

            delete[] vH_[SEND](d)(lii, ljj, lkk)(l);
            delete[] vH_[RECEIVE](d)(lii, ljj, lkk)(l);
          }
        }
      }
    }
  }
}

void
SEWASStarPU::destroyStressArenas()
{
  for (auto sc : { SWS::XX, SWS::YY, SWS::ZZ, SWS::XY, SWS::XZ, SWS::YZ }) {
    for (int lii = 0; lii < lnxx_; lii++) {
      for (int ljj = 0; ljj < lnyy_; ljj++) {
        for (int lkk = 0; lkk < lnzz_; lkk++) {
          for (auto l : { SWS::LEFT, SWS::RIGHT, SWS::BACKWARD, SWS::FORWARD, SWS::BOTTOM, SWS::TOP }) {
            const auto idx = haloHandleIndex(lii, ljj, lkk, l);

            starpu_data_unregister(sSendHandles_[sc][idx]);
            starpu_data_unregister(sRecvHandles_[sc][idx]);

            delete[] sigmaH_[SEND](sc)(lii, ljj, lkk)(l);
            delete[] sigmaH_[RECEIVE](sc)(lii, ljj, lkk)(l);
          }
        }
      }
    }
  }
}

void
SEWASStarPU::exchangeVelocityHalo(const SWS::Directions d, const int ts)
{
  auto* pHaloManager = HaloManager::getInstance();
  auto& vH_S = vH_[SEND];
  auto& vH_R = vH_[RECEIVE];

  for (int ii = 0; ii < nxx_; ii++) {
    for (int jj = 0; jj < nyy_; jj++) {
      for (int kk = 0; kk < nzz_; kk++) {
        if (Mesh3DPartitioning::rank_of(ii, jj, kk) != rank_) continue;

        const int lii = ii % lnxx_;
        const int ljj = jj % lnyy_;
        const int lkk = kk % lnzz_;

        for (auto l : { SWS::LEFT, SWS::RIGHT, SWS::BACKWARD, SWS::FORWARD, SWS::BOTTOM, SWS::TOP }) {
          int nii, njj, nkk;
          if (!neighborCoords(ii, jj, kk, l, nii, njj, nkk)) {
            continue;
          }

          const int nrank = Mesh3DPartitioning::rank_of(nii, njj, nkk);
          const auto hsize = pHaloManager->getHaloSize(l, lii, ljj, lkk);
          const auto lOpp = opposite(l);

          if (nrank == rank_) {
            const int nlii = nii % lnxx_, nljj = njj % lnyy_, nlkk = nkk % lnzz_;
            memcpy(vH_R(d)(nlii, nljj, nlkk)(lOpp), vH_S(d)(lii, ljj, lkk)(l), hsize * sizeof(SWS::RealType));
          } else {
            const auto idx = haloHandleIndex(lii, ljj, lkk, l);

            const auto outTag = makeTag(false, d, ts, ii, jj, kk, l);
            starpu_mpi_isend_detached(vSendHandles_[d][idx], nrank, outTag, MPI_COMM_WORLD, nullptr, nullptr);

            const auto inTag = makeTag(false, d, ts, nii, njj, nkk, lOpp);
            starpu_mpi_irecv_detached(vRecvHandles_[d][idx], nrank, inTag, MPI_COMM_WORLD, nullptr, nullptr);
          }
        }
      }
    }
  }
}

void
SEWASStarPU::exchangeStressHalo(const SWS::StressFieldComponents sc, const int ts)
{
  auto* pHaloManager = HaloManager::getInstance();
  auto& sigmaH_S = sigmaH_[SEND];
  auto& sigmaH_R = sigmaH_[RECEIVE];

  for (int ii = 0; ii < nxx_; ii++) {
    for (int jj = 0; jj < nyy_; jj++) {
      for (int kk = 0; kk < nzz_; kk++) {
        if (Mesh3DPartitioning::rank_of(ii, jj, kk) != rank_) continue;

        const int lii = ii % lnxx_;
        const int ljj = jj % lnyy_;
        const int lkk = kk % lnzz_;

        for (auto l : { SWS::LEFT, SWS::RIGHT, SWS::BACKWARD, SWS::FORWARD, SWS::BOTTOM, SWS::TOP }) {
          int nii, njj, nkk;
          if (!neighborCoords(ii, jj, kk, l, nii, njj, nkk)) {
            continue;
          }

          const int nrank = Mesh3DPartitioning::rank_of(nii, njj, nkk);
          const auto hsize = pHaloManager->getHaloSize(l, lii, ljj, lkk);
          const auto lOpp = opposite(l);

          if (nrank == rank_) {
            const int nlii = nii % lnxx_, nljj = njj % lnyy_, nlkk = nkk % lnzz_;
            memcpy(
              sigmaH_R(sc)(nlii, nljj, nlkk)(lOpp), sigmaH_S(sc)(lii, ljj, lkk)(l), hsize * sizeof(SWS::RealType));
          } else {
            const auto idx = haloHandleIndex(lii, ljj, lkk, l);

            const auto outTag = makeTag(true, sc, ts, ii, jj, kk, l);
            starpu_mpi_isend_detached(sSendHandles_[sc][idx], nrank, outTag, MPI_COMM_WORLD, nullptr, nullptr);

            const auto inTag = makeTag(true, sc, ts, nii, njj, nkk, lOpp);
            starpu_mpi_irecv_detached(sRecvHandles_[sc][idx], nrank, inTag, MPI_COMM_WORLD, nullptr, nullptr);
          }
        }
      }
    }
  }
}

void
SEWASStarPU::initializeFieldsCodelet(void* /*buffers*/[], void* cl_arg)
{
  int ii, jj, kk;
  starpu_codelet_unpack_args(cl_arg, &ii, &jj, &kk);
  LinearSeismicWaveModel::initializeFieldsWrapper(ii, jj, kk);
}

void
SEWASStarPU::phaseAVelocityCodelet(void* /*buffers*/[], void* cl_arg)
{
  int ts, ii, jj, kk, lii, ljj, lkk;
  starpu_codelet_unpack_args(cl_arg, &ts, &ii, &jj, &kk, &lii, &ljj, &lkk);

  auto& sigmaH_R = pInstance_->sigmaH_[RECEIVE];
  auto& vH_S = pInstance_->vH_[SEND];

  for (auto l : { SWS::LEFT, SWS::RIGHT, SWS::BACKWARD, SWS::FORWARD, SWS::BOTTOM, SWS::TOP }) {
    for (auto sc : { SWS::XX, SWS::YY, SWS::ZZ, SWS::XY, SWS::XZ, SWS::YZ }) {
      HaloManager::updateStressWrapper(l, sc, ts - 1, ii, jj, kk, sigmaH_R(sc)(lii, ljj, lkk)(l));
    }
  }

  for (auto d : { SWS::X, SWS::Y, SWS::Z }) {
    LinearSeismicWaveModel::computeVelocityWrapper(d, ts, ii, jj, kk);

    for (auto l : { SWS::LEFT, SWS::RIGHT, SWS::BACKWARD, SWS::FORWARD, SWS::BOTTOM, SWS::TOP }) {
      HaloManager::extractVelocityHaloWrapper(l, d, ts, ii, jj, kk, vH_S(d)(lii, ljj, lkk)(l));
    }
  }
}

void
SEWASStarPU::phaseBStressCodelet(void* /*buffers*/[], void* cl_arg)
{
  int ts, ii, jj, kk, lii, ljj, lkk;
  starpu_codelet_unpack_args(cl_arg, &ts, &ii, &jj, &kk, &lii, &ljj, &lkk);

  auto& vH_R = pInstance_->vH_[RECEIVE];
  auto& sigmaH_S = pInstance_->sigmaH_[SEND];

  for (auto l : { SWS::LEFT, SWS::RIGHT, SWS::BACKWARD, SWS::FORWARD, SWS::BOTTOM, SWS::TOP }) {
    for (auto d : { SWS::X, SWS::Y, SWS::Z }) {
      HaloManager::updateVelocityWrapper(l, d, ts, ii, jj, kk, vH_R(d)(lii, ljj, lkk)(l));
    }
  }

  for (auto sc : { SWS::XX, SWS::YY, SWS::ZZ, SWS::XY, SWS::XZ, SWS::YZ }) {
    LinearSeismicWaveModel::computeStressWrapper(sc, ts + 1, ii, jj, kk);

    for (auto l : { SWS::LEFT, SWS::RIGHT, SWS::BACKWARD, SWS::FORWARD, SWS::BOTTOM, SWS::TOP }) {
      HaloManager::extractStressHaloWrapper(l, sc, ts + 1, ii, jj, kk, sigmaH_S(sc)(lii, ljj, lkk)(l));
    }
  }
}

int
SEWASStarPU::run()
{
  int status = 0;

  static struct starpu_codelet initCl = makeCodelet(initializeFieldsCodelet, "initializeFieldsCodelet");
  static struct starpu_codelet phaseACl = makeCodelet(phaseAVelocityCodelet, "phaseAVelocityCodelet");
  static struct starpu_codelet phaseBCl = makeCodelet(phaseBStressCodelet, "phaseBStressCodelet");

  auto* pPriorityManager = Mesh3DPartitioning::getInstance()->getTaskPriorityManager();

  auto forEachLocalTile = [this](auto&& fn) {
    for (int ii = 0; ii < nxx_; ii++) {
      for (int jj = 0; jj < nyy_; jj++) {
        for (int kk = 0; kk < nzz_; kk++) {
          if (Mesh3DPartitioning::rank_of(ii, jj, kk) != rank_) continue;
          fn(ii, jj, kk, ii % lnxx_, jj % lnyy_, kk % lnzz_);
        }
      }
    }
  };

  forEachLocalTile([&](int ii, int jj, int kk, int /*lii*/, int /*ljj*/, int /*lkk*/) {
    const auto prio = clampPriority(pPriorityManager->getPriority(INITIALIZE_FIELDS, 0, ii, jj, kk));
    starpu_task_insert(&initCl,
                       STARPU_VALUE, &ii, sizeof(ii),
                       STARPU_VALUE, &jj, sizeof(jj),
                       STARPU_VALUE, &kk, sizeof(kk),
                       STARPU_PRIORITY, (unsigned long long)prio,
                       0);
  });
  starpu_task_wait_for_all();

  for (int ts = 2; ts <= nt_ - 2; ts += 2) {
    LOG(SWS::LOG_TRACE, "[start] Processing time-step {}", ts);

    forEachLocalTile([&](int ii, int jj, int kk, int lii, int ljj, int lkk) {
      /* Bundles UPDATE_STRESS + COMPUTE_VELOCITY + EXTRACT_VELOCITY_HALO into one
       * task; take the most urgent of the three sub-steps' priorities. */
      const auto prio = clampPriority((std::max)({
        pPriorityManager->getPriority(UPDATE_STRESS, ts, ii, jj, kk),
        pPriorityManager->getPriority(COMPUTE_VELOCITY, ts, ii, jj, kk),
        pPriorityManager->getPriority(EXTRACT_VELOCITY_HALO, ts, ii, jj, kk),
      }));
      starpu_task_insert(&phaseACl,
                         STARPU_VALUE, &ts, sizeof(ts),
                         STARPU_VALUE, &ii, sizeof(ii),
                         STARPU_VALUE, &jj, sizeof(jj),
                         STARPU_VALUE, &kk, sizeof(kk),
                         STARPU_VALUE, &lii, sizeof(lii),
                         STARPU_VALUE, &ljj, sizeof(ljj),
                         STARPU_VALUE, &lkk, sizeof(lkk),
                         STARPU_PRIORITY, (unsigned long long)prio,
                         0);
    });
    starpu_task_wait_for_all();

    for (auto d : { SWS::X, SWS::Y, SWS::Z }) {
      exchangeVelocityHalo(d, ts);
    }
    starpu_mpi_wait_for_all(MPI_COMM_WORLD);

    forEachLocalTile([&](int ii, int jj, int kk, int lii, int ljj, int lkk) {
      /* Bundles UPDATE_VELOCITY + COMPUTE_STRESS + EXTRACT_STRESS_HALO into one
       * task; take the most urgent of the three sub-steps' priorities. */
      const auto prio = clampPriority((std::max)({
        pPriorityManager->getPriority(UPDATE_VELOCITY, ts, ii, jj, kk),
        pPriorityManager->getPriority(COMPUTE_STRESS, ts, ii, jj, kk),
        pPriorityManager->getPriority(EXTRACT_STRESS_HALO, ts, ii, jj, kk),
      }));
      starpu_task_insert(&phaseBCl,
                         STARPU_VALUE, &ts, sizeof(ts),
                         STARPU_VALUE, &ii, sizeof(ii),
                         STARPU_VALUE, &jj, sizeof(jj),
                         STARPU_VALUE, &kk, sizeof(kk),
                         STARPU_VALUE, &lii, sizeof(lii),
                         STARPU_VALUE, &ljj, sizeof(ljj),
                         STARPU_VALUE, &lkk, sizeof(lkk),
                         STARPU_PRIORITY, (unsigned long long)prio,
                         0);
    });
    starpu_task_wait_for_all();

    for (auto sc : { SWS::XX, SWS::YY, SWS::ZZ, SWS::XY, SWS::XZ, SWS::YZ }) {
      exchangeStressHalo(sc, ts);
    }
    starpu_mpi_wait_for_all(MPI_COMM_WORLD);

    LOG(SWS::LOG_TRACE, "[stop] Processing time-step {}", ts);
  }

  return status;
}

SEWASStarPU::SEWASStarPU(const int nt, const int nxx, const int nyy, const int nzz)
  : nt_(nt)
  , nxx_(nxx)
  , nyy_(nyy)
  , nzz_(nzz)
{
  world_ = ExecutionContext::world();
  rank_ = ExecutionContext::rank();

  auto* pMesh = Mesh3DPartitioning::getInstance();
  lnxx_ = pMesh->lnxx();
  lnyy_ = pMesh->lnyy();
  lnzz_ = pMesh->lnzz();

  /* Evaluate task priorities: currently unused by this phase-barrier design (there is
   * no cross-phase overlap left to schedule around), kept for parity with SEWASPaRSEC
   * and as a hook for future, finer-grained pipelining. */
  pMesh->getTaskPriorityManager()->evaluate();

  buildVelocityArenas();
  buildStressArenas();
}

SEWASStarPU::~SEWASStarPU()
{
  destroyVelocityArenas();
  destroyStressArenas();
}

#endif // SEWAS_WITH_STARPU
