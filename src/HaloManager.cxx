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

#include <iostream>
#include <cassert>

#include "HaloManager.hxx"
#include "CartesianMesh3D.hxx"
#include "Mesh3DPartitioning.hxx"
#include "MetricsManager.hxx"
#include "LogManager.hxx"

HaloManager * HaloManager::pInstance_ = nullptr;

HaloManager * HaloManager::getInstance(SWS::StressField & sigma, SWS::Velocity & v,
				       const int & hnx, const int & hny, const int & hnz)
{
  if (nullptr == pInstance_){
    pInstance_ = new HaloManager(sigma, v, hnx, hny, hnz);
  }
  return pInstance_;
}

HaloManager * HaloManager::getInstance()
{
  assert(nullptr != pInstance_);
  return pInstance_;
}

void HaloManager::releaseInstance()
{
  if (nullptr != pInstance_){
    delete pInstance_;
    pInstance_ = nullptr;
  }
}

int HaloManager::extractVelocityHaloWrapper(const int l,
					    const int d,
					    const int ts,
					    const int ii, const int jj, const int kk,
					    void * vH)
{
  auto * const pMesh=Mesh3DPartitioning::getInstance();

  const int lii=pMesh->lii(ii);
  const int ljj=pMesh->ljj(jj);
  const int lkk=pMesh->lkk(kk);

  pInstance_->setVelocityHalo((const SWS::Locations) l, (const SWS::Directions) d, lii, ljj, lkk, vH);
  return pInstance_->extractVelocityHalo((const SWS::Locations) l, (const SWS::Directions) d, ts, lii, ljj, lkk);
}

int HaloManager::updateVelocityWrapper(const int l,
				       const int d,
				       const int ts,
				       const int ii, const int jj, const int kk,
				       void * vH)
{
  auto * const pMesh=Mesh3DPartitioning::getInstance();

  const int lii=pMesh->lii(ii);
  const int ljj=pMesh->ljj(jj);
  const int lkk=pMesh->lkk(kk);

  pInstance_->setVelocityHalo((const SWS::Locations) l, (const SWS::Directions) d, lii, ljj, lkk, vH);
  return pInstance_->updateVelocity((const SWS::Locations) l, (const SWS::Directions) d, ts, lii, ljj, lkk);
}

int HaloManager::extractStressHaloWrapper(const int l,
					  const int sc,
					  const int ts,
					  const int ii, const int jj, const int kk,
					  void * sigmaH)
{
  auto * const pMesh=Mesh3DPartitioning::getInstance();

  const int lii=pMesh->lii(ii);
  const int ljj=pMesh->ljj(jj);
  const int lkk=pMesh->lkk(kk);

  pInstance_->setStressHalo((const SWS::Locations) l, (const SWS::StressFieldComponents) sc, lii, ljj, lkk, sigmaH);
  return pInstance_->extractStressHalo((const SWS::Locations) l, (const SWS::StressFieldComponents) sc, ts, lii, ljj, lkk);
}

int HaloManager::updateStressWrapper(const int l,
				     const int sc,
				     const int ts,
				     const int ii, const int jj, const int kk,
				     void * sigmaH)
{
  auto * const pMesh=Mesh3DPartitioning::getInstance();

  const int lii=pMesh->lii(ii);
  const int ljj=pMesh->ljj(jj);
  const int lkk=pMesh->lkk(kk);

  pInstance_->setStressHalo((const SWS::Locations) l, (const SWS::StressFieldComponents) sc, lii, ljj, lkk, sigmaH);
  return pInstance_->updateStress((const SWS::Locations) l, (const SWS::StressFieldComponents) sc, ts, lii, ljj, lkk);
}

size_t HaloManager::getHaloSize(const int ii, const int jj, const int kk) const
{
  size_t hs=0;

  for (auto l : {SWS::LEFT, SWS::RIGHT, SWS::BACKWARD, SWS::FORWARD, SWS::BOTTOM, SWS::TOP}){
    hs+=getHaloSize(l, ii, jj, kk);
  }

  return hs;
}

size_t HaloManager::getHaloSize(const SWS::Locations l,
                                const int ii, const int jj, const int kk) const
{
  size_t hs=0;

  const auto cx=Mesh3DPartitioning::getInstance()->ccx()[ii];
  const auto cy=Mesh3DPartitioning::getInstance()->ccy()[jj];
  const auto cz=Mesh3DPartitioning::getInstance()->ccz()[kk];

  switch(l){
  case SWS::LEFT:
  case SWS::RIGHT:
    hs=hnx_*cy*cz;
    break;
  case SWS::BACKWARD:
  case SWS::FORWARD:
    hs=cx*hny_*cz;
    break;
  case SWS::BOTTOM:
  case SWS::TOP:
    hs=cx*cy*hnz_;
    break;
  default:
    LOG(SWS::ERROR, "Unknown location {} requested within HaloManager::getHaloSize()", l);
    break;
  }

  return hs;
}

int HaloManager::updateStress(const SWS::Locations l,
			      const SWS::StressFieldComponents sc,
			      const int ts,
			      const int ii, const int jj, const int kk) noexcept
{
  // MetricsManager::getInstance()->start("UpdateStress");

  int iStartH, iEndH;
  int jStartH, jEndH;
  int kStartH, kEndH;

  int iShift;
  int jShift;
  int kShift;

  setUpdateOffsets(l,
                   ii, jj, kk,
                   iStartH, iEndH,
                   jStartH, jEndH,
                   kStartH, kEndH,
                   iShift, jShift, kShift);

  size_t index=0;
  for (int i=iStartH; i<iEndH; i++){
    for (int j=jStartH; j<jEndH; j++){
      for (int k=kStartH; k<kEndH; k++){
	index=k*iEndH*jEndH+j*iEndH+i;
        assert(index < getHaloSize(l, ii, jj, kk));
	sigma_(sc)(ii,jj,kk)(i+iShift,j+jShift,k+kShift)=sigmaH_(sc)(ii,jj,kk)(l)[index];
      }
    }
  }

  // MetricsManager::getInstance()->stop("UpdateStress");

  return 0;
}

int HaloManager::updateVelocity(const SWS::Locations l,
				const SWS::Directions d,
				const int ts,
				const int ii, const int jj, const int kk) noexcept
{
  // MetricsManager::getInstance()->start("UpdateVelocity");

  int iStartH, iEndH;
  int jStartH, jEndH;
  int kStartH, kEndH;

  int iShift;
  int jShift;
  int kShift;

  setUpdateOffsets(l,
                   ii, jj, kk,
                   iStartH, iEndH,
                   jStartH, jEndH,
                   kStartH, kEndH,
                   iShift, jShift, kShift);

  size_t index=0;
  for (int i=iStartH; i<iEndH; i++){
    for (int j=jStartH; j<jEndH; j++){
      for (int k=kStartH; k<kEndH; k++){
	index=k*iEndH*jEndH+j*iEndH+i;
        assert(index < getHaloSize(l, ii, jj, kk));

        v_(d)(ii,jj,kk)(i+iShift,j+jShift,k+kShift)=vH_(d)(ii,jj,kk)(l)[index];
      }
    }
  }

  // MetricsManager::getInstance()->stop("UpdateVelocity");

  return 0;
}

int HaloManager::extractStressHalo(const SWS::Locations l,
				   const SWS::StressFieldComponents sc,
				   const int ts,
				   const int ii, const int jj, const int kk) noexcept
{
  // MetricsManager::getInstance()->start("ExtractStressHalo");

  int iStartH, iEndH;
  int jStartH, jEndH;
  int kStartH, kEndH;

  int iShift;
  int jShift;
  int kShift;

  setExtractOffsets(l,
                    ii, jj, kk,
                    iStartH, iEndH,
                    jStartH, jEndH,
                    kStartH, kEndH,
                    iShift, jShift, kShift);

  // TODO the following copy can be optimized by viewing the Halo data pointer as an Eigen's Array
  size_t index=0;
  for (int i=iStartH; i<iEndH; i++){
    for (int j=jStartH; j<jEndH; j++){
      for (int k=kStartH; k<kEndH; k++){
	index=k*iEndH*jEndH+j*iEndH+i;
        assert(index < getHaloSize(l, ii, jj, kk));
	sigmaH_(sc)(ii,jj,kk)(l)[index]=sigma_(sc)(ii,jj,kk)(i+iShift,j+jShift,k+kShift);
      }
    }
  }

  // MetricsManager::getInstance()->stop("ExtractStressHalo");

  return 0;
}

int HaloManager::extractVelocityHalo(const SWS::Locations l,
				     const SWS::Directions d,
				     const int ts,
				     const int ii, const int jj, const int kk) noexcept
{
  // MetricsManager::getInstance()->start("ExtractVelocityHalo");

  int iStartH, iEndH;
  int jStartH, jEndH;
  int kStartH, kEndH;

  int iShift;
  int jShift;
  int kShift;

  setExtractOffsets(l,
                    ii, jj, kk,
                    iStartH, iEndH,
                    jStartH, jEndH,
                    kStartH, kEndH,
                    iShift, jShift, kShift);

  size_t index=0;
  for (int i=iStartH; i<iEndH; i++){
    for (int j=jStartH; j<jEndH; j++){
      for (int k=kStartH; k<kEndH; k++){
	index=k*iEndH*jEndH+j*iEndH+i;

        assert(index < getHaloSize(l, ii, jj, kk));

        vH_(d)(ii,jj,kk)(l)[index]=v_(d)(ii,jj,kk)(i+iShift,j+jShift,k+kShift);
      }
    }
  }

  // MetricsManager::getInstance()->stop("ExtractVelocityHalo");

  return 0;
}

void HaloManager::setExtractOffsets(const SWS::Locations l,
                                    const int ii, const int jj, const int kk,
                                    int & iStartH, int & iEndH,
                                    int & jStartH, int & jEndH,
                                    int & kStartH, int & kEndH,
                                    int & iShift, int & jShift, int & kShift) noexcept
{
  const auto cx=Mesh3DPartitioning::getInstance()->ccx()[ii];
  const auto cy=Mesh3DPartitioning::getInstance()->ccy()[jj];
  const auto cz=Mesh3DPartitioning::getInstance()->ccz()[kk];

  const auto iEnd=v_(0)(ii,jj,kk).iEnd();
  const auto jEnd=v_(0)(ii,jj,kk).jEnd();
  const auto kEnd=v_(0)(ii,jj,kk).kEnd();

  iStartH=0, iEndH=cx;
  jStartH=0, jEndH=cy;
  kStartH=0, kEndH=cz;

  iShift=0;
  jShift=0;
  kShift=0;

  switch(l){
  case SWS::LEFT:
    iStartH=0;
    iEndH=hnx_;
    iShift=hnx_;
    break;
  case SWS::RIGHT:
    iStartH=0;
    iEndH=hnx_;
    iShift=iEnd-hnx_;
    break;
  case SWS::BACKWARD:
    jStartH=0;
    jEndH=hny_;
    jShift=hny_;
    break;
  case SWS::FORWARD:
    jStartH=0;
    jEndH=hny_;
    jShift=jEnd-hny_;
    break;
  case SWS::BOTTOM:
    kStartH=0;
    kEndH=hnz_;
    kShift=hnz_;
    break;
  case SWS::TOP:
    kStartH=0;
    kEndH=hnz_;
    kShift=kEnd-hnz_;
    break;
  default:
    LOG(SWS::ERROR, "Unknown halo location {} requested within HaloManager::setExtractOffsets()", l);
    break;
  }
}

void HaloManager::setUpdateOffsets(const SWS::Locations l,
                                   const int ii, const int jj, const int kk,
                                   int & iStartH, int & iEndH,
                                   int & jStartH, int & jEndH,
                                   int & kStartH, int & kEndH,
                                   int & iShift, int & jShift, int & kShift) noexcept
{
  const auto cx=Mesh3DPartitioning::getInstance()->ccx()[ii];
  const auto cy=Mesh3DPartitioning::getInstance()->ccy()[jj];
  const auto cz=Mesh3DPartitioning::getInstance()->ccz()[kk];

  const auto iEnd=v_(0)(ii,jj,kk).iEnd();
  const auto jEnd=v_(0)(ii,jj,kk).jEnd();
  const auto kEnd=v_(0)(ii,jj,kk).kEnd();

  iStartH=0, iEndH=cx;
  jStartH=0, jEndH=cy;
  kStartH=0, kEndH=cz;

  iShift=0;
  jShift=0;
  kShift=0;

  switch(l){
  case SWS::LEFT:
    iStartH=0;
    iEndH=hnx_;
    iShift=0;
    break;
  case SWS::RIGHT:
    iStartH=0;
    iEndH=hnx_;
    iShift=iEnd;
    break;
  case SWS::BACKWARD:
    jStartH=0;
    jEndH=hny_;
    jShift=0;
    break;
  case SWS::FORWARD:
    jStartH=0;
    jEndH=hny_;
    jShift=jEnd;
    break;
  case SWS::BOTTOM:
    kStartH=0;
    kEndH=hnz_;
    kShift=0;
    break;
  case SWS::TOP:
    kStartH=0;
    kEndH=hnz_;
    kShift=kEnd;
    break;
  default:
    LOG(SWS::ERROR, "Unknown halo location {} requested within HaloManager::setUpdateOffsets()", l);
    break;
  }
}

HaloManager::HaloManager(SWS::StressField & sigma, SWS::Velocity & v,
			 const int & hnx, const int & hny, const int & hnz):sigma_(sigma),
									    v_(v),
									    hnx_(hnx),
									    hny_(hny),
									    hnz_(hnz)
{

  const auto & lnxx=Mesh3DPartitioning::getInstance()->lnxx();
  const auto & lnyy=Mesh3DPartitioning::getInstance()->lnyy();
  const auto & lnzz=Mesh3DPartitioning::getInstance()->lnzz();

  for (auto sc : {SWS::XX, SWS::YY, SWS::ZZ, SWS::XY, SWS::XZ, SWS::YZ}){
    sigmaH_(sc).resize(lnxx,lnyy,lnzz);
  }

  for (auto d : {SWS::X, SWS::Y, SWS::Z}){
    vH_(d).resize(lnxx,lnyy,lnzz);
  }

}

HaloManager::HaloManager():sigma_(pInstance_->sigma()),
			   v_(pInstance_->v()),
			   hnx_(pInstance_->hnx()),
			   hny_(pInstance_->hny()),
			   hnz_(pInstance_->hnz())
{
  sigmaH_=pInstance_->sigmaH();
  vH_=pInstance_->vH();
}

HaloManager::~HaloManager()
{
}
