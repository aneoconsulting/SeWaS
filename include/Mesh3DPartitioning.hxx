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

#include <cassert>
#include <cstdarg>
#include <iostream>
#include <vector>

#ifdef SEWAS_WITH_PARSEC
#include <parsec/data_distribution.h>
#include <parsec/parsec_config.h>
#endif

#include "Config.hxx"
#include "DataSet.hxx"
#include "LinearSeismicWaveModel.hxx"
#include "MinimumCommunicationPriorityEvaluator.hxx"

class Mesh3DPartitioning
{
public:
  static Mesh3DPartitioning* getInstance(const int cx,
                                         const int cy,
                                         const int cz,
                                         const int hnx,
                                         const int hny,
                                         const int hnz,
                                         const int P,
                                         const int Q,
                                         const int R);
  static Mesh3DPartitioning* getInstance();

  static void releaseInstance();

  void buildSpatialField(SWS::SpatialField& f);

  inline const auto& ccx() const { return ccx_; }
  inline const auto& ccy() const { return ccy_; }
  inline const auto& ccz() const { return ccz_; }

  inline const auto& nxx() const { return nxx_; }
  inline const auto& nyy() const { return nyy_; }
  inline const auto& nzz() const { return nzz_; }

  inline const auto& P() const { return P_; }
  inline const auto& Q() const { return Q_; }
  inline const auto& R() const { return R_; }

  inline const auto& lnxx() const { return lnxx_; }
  inline const auto& lnyy() const { return lnyy_; }
  inline const auto& lnzz() const { return lnzz_; }

  inline const auto& tileSize() const { return tileSize_; }

  inline const auto& lncells() const { return lncells_; }

  inline const auto lii(const int ii) const { return ii % lnxx_; }
  inline const auto ljj(const int jj) const { return jj % lnyy_; }
  inline const auto lkk(const int kk) const { return kk % lnzz_; }

  inline const auto lp(const int ii) const { return ii / lnxx_; }
  inline const auto lq(const int jj) const { return jj / lnyy_; }
  inline const auto lr(const int kk) const { return kk / lnzz_; }

  static inline int rank_of(const int ii, const int jj, const int kk)
  {
    int rank = 0;

    int P = pInstance_->P();
    int Q = pInstance_->Q();

    int lp = pInstance_->lp(ii);
    int lq = pInstance_->lq(jj);
    int lr = pInstance_->lr(kk);

    rank = P * Q * lr + P * lq + lp;

    return rank;
  }

#ifdef SEWAS_WITH_PARSEC
  static unsigned int rank_of(parsec_data_collection_t* desc, ...);

  static int vpid_of(parsec_data_collection_t* desc, ...);

  static parsec_data_t* data_of(parsec_data_collection_t* desc, ...);

  static inline parsec_data_key_t data_key(parsec_data_collection_t* desc, ...)
  {
    int k;
    va_list ap;
    (void)desc;
    va_start(ap, desc);
    k = va_arg(ap, int);
    va_end(ap);
    return (uint64_t)k;
  }
#endif

  class TaskPriorityManager
  {
  public:
    TaskPriorityManager() {}

    ~TaskPriorityManager() {}

    inline void evaluate()
    {
      auto pMeshPartitioning = Mesh3DPartitioning::getInstance();

      auto nxx = pMeshPartitioning->nxx();
      auto nyy = pMeshPartitioning->nyy();
      auto nzz = pMeshPartitioning->nzz();

      auto lnxx = pMeshPartitioning->lnxx();
      auto lnyy = pMeshPartitioning->lnyy();
      auto lnzz = pMeshPartitioning->lnzz();

      auto priorityEvaluator =
        std::make_unique<MinimumCommunicationPriorityEvaluator>(nxx, nyy, nzz, lnxx, lnyy, lnzz);

      priorityEvaluator->evaluate<Mesh3DPartitioning>(taskPriorities_);
    }

    static inline int getPriorityWrapper(const int taskType,
                                         const int ts,
                                         const int ii,
                                         const int jj,
                                         const int kk)
    {
      return Mesh3DPartitioning::getInstance()->getTaskPriorityManager()->getPriority(
        taskType, ts, ii, jj, kk);
    }

    inline int getPriority(const int taskType, const int ts, const int ii, const int jj, const int kk)
    {
      return taskPriorities_(taskType)(ii, jj, kk);
    }

  private:
    SWS::TaskPriorities taskPriorities_;
  };

  inline TaskPriorityManager* getTaskPriorityManager() { return pTaskPriorityManager_; }

private:
  Mesh3DPartitioning(const int cx,
                     const int cy,
                     const int cz,
                     const int hnx,
                     const int hny,
                     const int hnz,
                     const int P,
                     const int Q,
                     const int R);
  ~Mesh3DPartitioning();

  static Mesh3DPartitioning* pInstance_;

  /* Priority manager */
  TaskPriorityManager* pTaskPriorityManager_;

  /* Size of a single sub-block */
  const int cx_;
  const int cy_;
  const int cz_;

  /* Halo size */
  const int hnx_;
  const int hny_;
  const int hnz_;

  /* Sizes of all sub-blocks according to their coordinates */
  std::vector<int> ccx_;
  std::vector<int> ccy_;
  std::vector<int> ccz_;

  /* Number of processes on each of the three axes */
  const int P_;
  const int Q_;
  const int R_;

  /* Number of sub-blocks on each of the three axes */
  int nxx_;
  int nyy_;
  int nzz_;

  /* Local number of cells belonging to the current process */
  int lncells_;

  /* Number of cells within a single tile */
  int tileSize_;

  /* Local number of sub-blocks on each of the three axes */
  int lnxx_;
  int lnyy_;
  int lnzz_;
};
