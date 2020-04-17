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

#include "Mesh3DPartitioning.hxx"
#include "CartesianMesh3D.hxx"
#include "HaloManager.hxx"

#ifdef SEWAS_WITH_PARSEC
#include <parsec/parsec_config.h>
//#include <parsec/vpmap.h>
#endif

Mesh3DPartitioning* Mesh3DPartitioning::pInstance_ = nullptr;

Mesh3DPartitioning*
Mesh3DPartitioning::getInstance(const SEWASParameterManager& pm, const int hnx, const int hny, const int hnz)
{
  if (nullptr == pInstance_) {
    pInstance_ = new Mesh3DPartitioning(pm, hnx, hny, hnz);
    return pInstance_;
  } else {
    return pInstance_;
  }
}

Mesh3DPartitioning*
Mesh3DPartitioning::getInstance()
{
  assert(nullptr != pInstance_);
  return pInstance_;
}

void
Mesh3DPartitioning::releaseInstance()
{
  if (pInstance_) {
    delete pInstance_;
    pInstance_ = nullptr;
  }
}

void
Mesh3DPartitioning::buildSpatialField(SWS::SpatialField& f)
{
  f.resize(lnxx_, lnyy_, lnzz_);
  for (int ii = 0; ii < lnxx_; ii++) {
    for (int jj = 0; jj < lnyy_; jj++) {
      for (int kk = 0; kk < lnzz_; kk++) {
        f(ii, jj, kk).resize(ccx_[ii], ccy_[jj], ccz_[kk]);
        f(ii, jj, kk).setHaloSize(hnx_, hny_, hnz_);
      }
    }
  }
}

#ifdef SEWAS_WITH_PARSEC
unsigned int
Mesh3DPartitioning::rank_of(parsec_data_collection_t* desc, ...)
{
  unsigned int rank;

  int ii, jj, kk;
  va_list ap;
  (void)desc;
  va_start(ap, desc);
  ii = va_arg(ap, int);
  jj = va_arg(ap, int);
  kk = va_arg(ap, int);
  va_end(ap);

  rank = rank_of(ii, jj, kk);

  return rank;
}

int
Mesh3DPartitioning::vpid_of(parsec_data_collection_t* desc, ...)
{
  int vpid = 0;

  /*
    int ii, jj, kk;
    va_list ap;
    (void)desc;
    va_start(ap, desc);
    ii = va_arg(ap, int);
    jj = va_arg(ap, int);
    kk = va_arg(ap, int);
    va_end(ap);

    int lii=pInstance_->lii(ii);

    int nb_vp=vpmap_get_nb_vp();

    vpid=lii/(pInstance_->lnxx() / nb_vp);

    assert(vpid < nb_vp);
  */

  return vpid;
}

parsec_data_t*
Mesh3DPartitioning::data_of(parsec_data_collection_t* desc, ...)
{
  int ii, jj, kk;
  va_list ap;
  (void)desc;
  va_start(ap, desc);
  ii = va_arg(ap, int);
  jj = va_arg(ap, int);
  kk = va_arg(ap, int);
  va_end(ap);

  int lii = pInstance_->lii(ii);
  int ljj = pInstance_->ljj(jj);
  int lkk = pInstance_->lkk(kk);

  assert(pInstance_->rank_of(desc, ii, jj, kk) == desc->myrank);

  const auto& rho = DataSet::getInstance()->rho();

  return (parsec_data_t*)&rho(lii, ljj, lkk);
}
#endif

Mesh3DPartitioning::Mesh3DPartitioning(const SEWASParameterManager& pm,
                                       const int hnx,
                                       const int hny,
                                       const int hnz)
  : pm_(pm)
  , cx_(pm.cx())
  , cy_(pm.cy())
  , cz_(pm.cz())
  , hnx_(hnx)
  , hny_(hny)
  , hnz_(hnz)
  , P_(pm.P())
  , Q_(pm.Q())
  , R_(pm.R())
{
  lnxx_ = pm.lnxx();
  lnyy_ = pm.lnyy();
  lnzz_ = pm.lnzz();

  // For each dimension, we are adding halo (hnx, hny, hnz) to the orginal size of the Sub-block of cells

  ccx_.resize(lnxx_);
  for (int ii = 0; ii < lnxx_; ii++)
    ccx_[ii] = hnx + pm.cx() + hnx;

  ccy_.resize(lnyy_);
  for (int jj = 0; jj < lnyy_; jj++)
    ccy_[jj] = hny + pm.cy() + hny;

  ccz_.resize(lnzz_);
  for (int kk = 0; kk < lnzz_; kk++)
    ccz_[kk] = hnz + pm.cz() + hnz;

  tileSize_ = ccx_[0] * ccy_[0] * ccz_[0];

  lncells_ = (lnxx_ * lnyy_ * lnzz_) * tileSize_;

  /* Create an instance of the task priority manager */
  pTaskPriorityManager_ = new TaskPriorityManager();
}

Mesh3DPartitioning::~Mesh3DPartitioning()
{
  delete pTaskPriorityManager_;
  pTaskPriorityManager_ = nullptr;
}
