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
#define _USE_MATH_DEFINES // for C++
#include <cmath>

#include "ExecutionContext.hxx"

#ifdef SEWAS_WITH_PARSEC
#include <parsec/parsec_config.h>
#include <parsec/data_distribution.h>

#include "sewas.h"
#endif

#include "Config.hxx"
#include "LogManager.hxx"
#include "CentralFDOperator.hxx"
#include "CartesianMesh3D.hxx"
#include "DataSet.hxx"
#include "ExternalSource.hxx"
#include "LinearSeismicWaveModel.hxx"
#include "HaloManager.hxx"
#include "IOManager.hxx"
#include "Mesh3DPartitioning.hxx"
#include "SEWASSequential.hxx"
#include "SEWASPaRSEC.hxx"
#include "MetricsManager.hxx"


LinearSeismicWaveModel * LinearSeismicWaveModel::pInstance_ = nullptr;

LinearSeismicWaveModel * LinearSeismicWaveModel::getInstance(const CentralFDOperator & fdo,
                                                             const SEWASParameterManager & pm,
							     const int nt, const float tmax){
  if (nullptr == pInstance_){
    pInstance_ = new LinearSeismicWaveModel(fdo, pm, nt, tmax);
    return pInstance_;
  }
  else{
    return pInstance_;
  }
}

LinearSeismicWaveModel * LinearSeismicWaveModel::getInstance(){
  assert(nullptr != pInstance_);
  return pInstance_;
}

void LinearSeismicWaveModel::releaseInstance(){
  if (pInstance_){
    delete pInstance_;
    pInstance_ = nullptr;
  }
}

int LinearSeismicWaveModel::propagate() noexcept {

  const int & nxx=Mesh3DPartitioning::getInstance()->nxx();
  const int & nyy=Mesh3DPartitioning::getInstance()->nyy();
  const int & nzz=Mesh3DPartitioning::getInstance()->nzz();

#ifdef SEWAS_WITH_PARSEC
  LOG(SWS::LOG_INFO, "Starting the PaRSEC-based runner");

  SEWASPaRSEC * sewasPaRSEC=SEWASPaRSEC::getInstance(nt_,
                                                     nxx, nyy, nzz);
  LOG(SWS::LOG_INFO, "PaRSEC-based runner started");


  LOG(SWS::LOG_INFO, "Starting execution of the core seismic simulation");
  sewasPaRSEC->run();
  LOG(SWS::LOG_INFO, "Seismic simulation completed");

  SEWASPaRSEC::releaseInstance();
#else
  LOG(SWS::LOG_INFO, "Starting the sequential runner");
  SEWASSequential * sewasSequential=SEWASSequential::getInstance(nt_,
                                                                 nxx, nyy, nzz);
  LOG(SWS::LOG_INFO, "Sequential runner started");


  LOG(SWS::LOG_INFO, "Starting execution of the core seismic simulation");
  sewasSequential->run();
  LOG(SWS::LOG_INFO, "Seismic simulation completed");

  SEWASSequential::releaseInstance();
#endif

  return 0;
}

int LinearSeismicWaveModel::initializeFieldsWrapper(const int ii, const int jj, const int kk)
{
  LogManager::getInstance()->log<SWS::LOG_TRACE>("Initializing data on tile ({}, {}, {})", ii, jj, kk);

  /* Initialize rho, mu and lambda */
  DataSet::getInstance()->initialize(ii, jj, kk);

  /* Initialize velocity and stress fields */
  pInstance_->initialize(ii, jj, kk);

  LogManager::getInstance()->log<SWS::LOG_TRACE>("Completed initializing data on tile ({}, {}, {})", ii, jj, kk);

  return 0;
}

int LinearSeismicWaveModel::computeVelocityWrapper(const int d,
						   const int ts,
						   const int ii, const int jj, const int kk)
{
  auto * pMesh=Mesh3DPartitioning::getInstance();

  const int lii=pMesh->lii(ii);
  const int ljj=pMesh->ljj(jj);
  const int lkk=pMesh->lkk(kk);

  return pInstance_->computeVelocity((SWS::Directions) d, ts, lii, ljj, lkk);
}

int LinearSeismicWaveModel::computeVelocity(const SWS::Directions & d,
					    const int & ts,
					    const int & ii, const int & jj, const int & kk){

  MetricsManager::getInstance()->start("ComputeVelocity");

  IOManager::getInstance().start<COMPUTE_VELOCITY>(ts, ii, jj, kk);

  const auto & bx=DataSet::getInstance()->b(SWS::X)(ii,jj,kk);
  const auto & by=DataSet::getInstance()->b(SWS::Y)(ii,jj,kk);
  const auto & bz=DataSet::getInstance()->b(SWS::Z)(ii,jj,kk);

  auto & vX=v_(SWS::X)(ii,jj,kk);
  auto & vY=v_(SWS::Y)(ii,jj,kk);
  auto & vZ=v_(SWS::Z)(ii,jj,kk);

  const auto & sigmaXX=sigma_(SWS::XX)(ii,jj,kk);
  const auto & sigmaYY=sigma_(SWS::YY)(ii,jj,kk);
  const auto & sigmaZZ=sigma_(SWS::ZZ)(ii,jj,kk);
  const auto & sigmaXY=sigma_(SWS::XY)(ii,jj,kk);
  const auto & sigmaXZ=sigma_(SWS::XZ)(ii,jj,kk);
  const auto & sigmaYZ=sigma_(SWS::YZ)(ii,jj,kk);

  const auto iStart=vX.iStart();
  const auto iEnd=vX.iEnd();

  const auto jStart=vX.jStart();
  const auto jEnd=vX.jEnd();

  const auto kStart=vX.kStart();
  const auto kEnd=vX.kEnd();

  const auto ds=CartesianMesh3D::getInstance()->ds();
  const auto dt=2*dt_;

  constexpr SWS::RealType c2=1./24;
  constexpr SWS::RealType c3=9./8;

  switch(d){
  case SWS::X:{

    LOG(SWS::LOG_TRACE, "[start] Computing Vx at time-step {} on tile ({}, {}, {})", ts, ii, jj, kk);

#ifdef BLOCKWISE_FDO
    for (int i=iStart; i<iEnd; i++){
      for (int j=jStart; j<jEnd; j++){
        vX(i,j)+=(fdo_.apply<SWS::X>(sigmaXX, i, j, SWS::DK)
                + fdo_.apply<SWS::Y>(sigmaXY, i, j, SWS::DK)
                + fdo_.apply<SWS::Z>(sigmaXZ, i, j, SWS::DK))*dt*bx(i,j);
      }
    }
#else
    for (int i=iStart; i<iEnd; i++){
      for (int j=jStart; j<jEnd; j++){
        for (int k=kStart; k<kEnd; k++){

          SWS::RealType vx=0.0;

          vx+=(c3*(sigmaXX(i,j,k)-sigmaXX(i-1,j,k)) - c2*(sigmaXX(i+1,j,k)-sigmaXX(i-2,j,k)))/ds; // Dx(sigmaXX)
          vx+=(c3*(sigmaXY(i,j,k)-sigmaXY(i,j-1,k)) - c2*(sigmaXY(i,j+1,k)-sigmaXY(i,j-2,k)))/ds; // Dy(sigmaXY)
          vx+=(c3*(sigmaXZ(i,j,k)-sigmaXZ(i,j,k-1)) - c2*(sigmaXZ(i,j,k+1)-sigmaXZ(i,j,k-2)))/ds; // Dz(sigmaXZ)
          vx*=dt*bx(i,j,k);

          vX(i,j,k)+=vx;

        }
      }
    }
#endif

    LOG(SWS::LOG_TRACE, "[stop] Computing Vx at time-step {} on tile ({}, {}, {})", ts, ii, jj, kk);

    addVelocitySource(SWS::X, ts, ii, jj, kk);

    LOG(SWS::LOG_TRACE, "||Vx({},{},{},{})||^2 = {}", ts, ii, jj, kk, vX.norm2());

    break;
  }
  case SWS::Y:{

    LOG(SWS::LOG_TRACE, "[start] Computing Vy at time-step {} on tile ({}, {}, {})", ts, ii, jj, kk);

#ifdef BLOCKWISE_FDO
    for (int i=iStart; i<iEnd; i++){
      for (int j=jStart; j<jEnd; j++){
        vY(i,j)+=(fdo_.apply<SWS::X>(sigmaXY, i+1, j,   SWS::DK)
                + fdo_.apply<SWS::Y>(sigmaYY, i,   j+1, SWS::DK)
                + fdo_.apply<SWS::Z>(sigmaYZ, i,   j,   SWS::DK))*dt*by(i,j);
      }
    }
#else
    for (int i=iStart; i<iEnd; i++){
      for (int j=jStart; j<jEnd; j++){
        for (int k=kStart; k<kEnd; k++){

          SWS::RealType vy=0.0;

          vy+=(c3*(sigmaXY(i+1,j,k)-sigmaXY(i,j,k)) - c2*(sigmaXY(i+2,j,k)-sigmaXY(i-1,j,k)))/ds; // Dx(sigmaXY)
          vy+=(c3*(sigmaYY(i,j+1,k)-sigmaYY(i,j,k)) - c2*(sigmaYY(i,j+2,k)-sigmaYY(i,j-1,k)))/ds; // Dy(sigmaYY)
          vy+=(c3*(sigmaYZ(i,j,k)-sigmaYZ(i,j,k-1)) - c2*(sigmaYZ(i,j,k+1)-sigmaYZ(i,j,k-2)))/ds; // Dz(sigmaYZ)
          vy*=dt*by(i,j,k);

          vY(i,j,k)+=vy;

        }
      }
    }
#endif

    LOG(SWS::LOG_TRACE, "[stop] Computing Vy at time-step {} on tile ({}, {}, {})", ts, ii, jj, kk);

    addVelocitySource(SWS::Y, ts, ii, jj, kk);

    LOG(SWS::LOG_TRACE, "||Vy({},{},{},{})||^2 = {}", ts, ii, jj, kk, vY.norm2());

    break;
  }
  case SWS::Z:{

    LOG(SWS::LOG_TRACE, "[start] Computing Vz at time-step {} on tile ({}, {}, {})", ts, ii, jj, kk);

#ifdef BLOCKWISE_FDO
    for (int i=iStart; i<iEnd; i++){
      for (int j=jStart; j<jEnd; j++){
        vZ(i,j)+=(fdo_.apply<SWS::X>(sigmaXZ, i+1, j, SWS::DK)
                + fdo_.apply<SWS::Y>(sigmaYZ, i,   j, SWS::DK)
                + fdo_.apply<SWS::Z>(sigmaZZ, i,   j, SWS::DK+1))*dt*bz(i,j);
      }
    }
#else
    for (int i=iStart; i<iEnd; i++){
      for (int j=jStart; j<jEnd; j++){
        for (int k=kStart; k<kEnd; k++){

          SWS::RealType vz=0.0;

          vz+=(c3*(sigmaXZ(i+1,j,k)-sigmaXZ(i,j,k)) - c2*(sigmaXZ(i+2,j,k)-sigmaXZ(i-1,j,k)))/ds; // Dx(sigmaXZ)
          vz+=(c3*(sigmaYZ(i,j,k)-sigmaYZ(i,j-1,k)) - c2*(sigmaYZ(i,j+1,k)-sigmaYZ(i,j-2,k)))/ds; // Dy(sigmaYZ)
          vz+=(c3*(sigmaZZ(i,j,k+1)-sigmaZZ(i,j,k)) - c2*(sigmaZZ(i,j,k+2)-sigmaZZ(i,j,k-1)))/ds; // Dz(sigmaZZ)
          vz*=dt*bz(i,j,k);

          vZ(i,j,k)+=vz;

        }
      }
    }
#endif

    LOG(SWS::LOG_TRACE, "[stop] Computing Vz at time-step {} on tile ({}, {}, {})", ts, ii, jj, kk);

    LOG(SWS::LOG_TRACE, "||Vz({},{},{},{})||^2 = {}", ts, ii, jj, kk, vZ.norm2());

    break;
  }
  default:
    LOG(SWS::LOG_ERROR, "Unknown spatial direction {} requested within LinearSeismicWaveModel::computeVelocity()", d);
    break;
  }

  IOManager::getInstance().stop<COMPUTE_VELOCITY>(ts, ii, jj, kk);

  MetricsManager::getInstance()->stop("ComputeVelocity");

  return 0;

}

int LinearSeismicWaveModel::computeStressWrapper(const int sc,
						 const int ts,
						 const int ii, const int jj, const int kk)
{
  auto * pMesh=Mesh3DPartitioning::getInstance();

  const int lii=pMesh->lii(ii);
  const int ljj=pMesh->ljj(jj);
  const int lkk=pMesh->lkk(kk);

  return pInstance_->computeStress((SWS::StressFieldComponents) sc, ts, lii, ljj, lkk);
}

int LinearSeismicWaveModel::computeStress(const SWS::StressFieldComponents & sc,
					  const int & ts,
					  const int & ii, const int & jj, const int & kk){

  MetricsManager::getInstance()->start("ComputeStress");

  IOManager::getInstance().start<COMPUTE_STRESS>(ts, ii, jj, kk);

  const auto & lambda=DataSet::getInstance()->lambda(sc)(ii,jj,kk);
  const auto & mu=DataSet::getInstance()->mu(sc)(ii,jj,kk);

  const auto & vX=v_(SWS::X)(ii,jj,kk);
  const auto & vY=v_(SWS::Y)(ii,jj,kk);
  const auto & vZ=v_(SWS::Z)(ii,jj,kk);

  auto & sigmaXX=sigma_(SWS::XX)(ii,jj,kk);
  auto & sigmaYY=sigma_(SWS::YY)(ii,jj,kk);
  auto & sigmaZZ=sigma_(SWS::ZZ)(ii,jj,kk);
  auto & sigmaXY=sigma_(SWS::XY)(ii,jj,kk);
  auto & sigmaXZ=sigma_(SWS::XZ)(ii,jj,kk);
  auto & sigmaYZ=sigma_(SWS::YZ)(ii,jj,kk);

  const auto iStart=vX.iStart();
  const auto iEnd=vX.iEnd();

  const auto jStart=vX.jStart();
  const auto jEnd=vX.jEnd();

  const auto kStart=vX.kStart();
  const auto kEnd=vX.kEnd();

  const auto ds=CartesianMesh3D::getInstance()->ds();
  const auto dt=2*dt_;

  constexpr SWS::RealType c3=9./8;
  constexpr SWS::RealType c2=1./24;

  switch(sc){
  case SWS::XX:

    LOG(SWS::LOG_TRACE, "[start] Computing Sxx at time-step {} on tile ({}, {}, {})", ts, ii, jj, kk);

#ifdef BLOCKWISE_FDO
    for (int i=iStart; i<iEnd; i++){
      for (int j=jStart; j<jEnd; j++){
        sigmaXX(i,j)+=((fdo_.apply<SWS::Y>(vY, i, j, SWS::DK) + fdo_.apply<SWS::Z>(vZ, i, j, SWS::DK))*lambda(i,j)
                     + (lambda(i,j)+2.0*mu(i,j))*fdo_.apply<SWS::X>(vX, i+1, j, SWS::DK))*dt;
      }
    }
#else
    for (int i=iStart; i<iEnd; i++){
      for (int j=jStart; j<jEnd; j++){
        for (int k=kStart; k<kEnd; k++){

          SWS::RealType sxx=0.0;

          sxx+=(c3*(vY(i,j,k)-vY(i,j-1,k)) - c2*(vY(i,j+1,k)-vY(i,j-2,k)))/ds; // Dy(vY)
          sxx+=(c3*(vZ(i,j,k)-vZ(i,j,k-1)) - c2*(vZ(i,j,k+1)-vZ(i,j,k-2)))/ds; // Dz(vZ)
          sxx*=lambda(i,j,k);
          sxx+=(lambda(i,j,k)+2.0*mu(i,j,k))*(c3*(vX(i+1,j,k)-vX(i,j,k)) - c2*(vX(i+2,j,k)-vX(i-1,j,k)))/ds;
          sxx*=dt;

          sigmaXX(i,j,k)+=sxx;

        }
      }
    }
#endif

    LOG(SWS::LOG_TRACE, "[stop] Computing Sxx at time-step {} on tile ({}, {}, {})", ts, ii, jj, kk);

    LOG(SWS::LOG_TRACE, "||sigmaXX({},{},{},{})||^2 = {}", ts, ii, jj, kk, sigmaXX.norm2());

    break;
  case SWS::YY:

    LOG(SWS::LOG_TRACE, "[start] Computing Syy at time-step {} on tile ({}, {}, {})", ts, ii, jj, kk);

#ifdef BLOCKWISE_FDO
    for (int i=iStart; i<iEnd; i++){
      for (int j=jStart; j<jEnd; j++){
        sigmaYY(i,j)+=((fdo_.apply<SWS::X>(vX, i+1, j, SWS::DK) + fdo_.apply<SWS::Z>(vZ, i, j, SWS::DK))*lambda(i,j)
                     + (lambda(i,j)+2.0*mu(i,j))*fdo_.apply<SWS::Y>(vY, i, j, SWS::DK))*dt;
      }
    }
#else
    for (int i=iStart; i<iEnd; i++){
      for (int j=jStart; j<jEnd; j++){
        for (int k=kStart; k<kEnd; k++){

          SWS::RealType syy=0.0;

          syy+=(c3*(vX(i+1,j,k)-vX(i,j,k)) - c2*(vX(i+2,j,k)-vX(i-1,j,k)))/ds; // Dx(vX)
          syy+=(c3*(vZ(i,j,k)-vZ(i,j,k-1)) - c2*(vZ(i,j,k+1)-vZ(i,j,k-2)))/ds; // Dz(vZ)
          syy*=lambda(i,j,k);
          syy+=(lambda(i,j,k)+2.0*mu(i,j,k))*(c3*(vY(i,j,k)-vY(i,j-1,k)) - c2*(vY(i,j+1,k)-vY(i,j-2,k)))/ds;
          syy*=dt;

          sigmaYY(i,j,k)+=syy;

        }
      }
    }
#endif

    LOG(SWS::LOG_TRACE, "[stop] Computing Syy at time-step {} on tile ({}, {}, {})", ts, ii, jj, kk);

    LOG(SWS::LOG_TRACE, "||sigmaYY({},{},{},{})||^2 = {}", ts, ii, jj, kk, sigmaYY.norm2());

    break;
  case SWS::ZZ:

    LOG(SWS::LOG_TRACE, "[start] Computing Szz at time-step {} on tile ({}, {}, {})", ts, ii, jj, kk);

#ifdef BLOCKWISE_FDO
    for (int i=iStart; i<iEnd; i++){
      for (int j=jStart; j<jEnd; j++){
        sigmaZZ(i,j)+=((fdo_.apply<SWS::X>(vX, i+1, j, SWS::DK) + fdo_.apply<SWS::Y>(vY, i, j, SWS::DK))*lambda(i,j)
                     + (lambda(i,j)+2.0*mu(i,j))*fdo_.apply<SWS::Z>(vZ, i, j, SWS::DK))*dt;
      }
    }
#else
    for (int i=iStart; i<iEnd; i++){
      for (int j=jStart; j<jEnd; j++){
        for (int k=kStart; k<kEnd; k++){

          SWS::RealType szz=0.0;

          szz+=(c3*(vX(i+1,j,k)-vX(i,j,k)) - c2*(vX(i+2,j,k)-vX(i-1,j,k)))/ds; // Dx(vX)
          szz+=(c3*(vY(i,j,k)-vY(i,j-1,k)) - c2*(vY(i,j+1,k)-vY(i,j-2,k)))/ds; // Dy(vY)
          szz*=lambda(i,j,k);
          szz+=(lambda(i,j,k)+2.0*mu(i,j,k))*(c3*(vZ(i,j,k)-vZ(i,j,k-1)) - c2*(vZ(i,j,k+1)-vZ(i,j,k-2)))/ds;
          szz*=dt;

          sigmaZZ(i,j,k)+=szz;

        }
      }
    }
#endif

    LOG(SWS::LOG_TRACE, "[stop] Computing Szz at time-step {} on tile ({}, {}, {})", ts, ii, jj, kk);

    LOG(SWS::LOG_TRACE, "||sigmaZZ({},{},{},{})||^2 = {}", ts, ii, jj, kk, sigmaZZ.norm2());

    break;
  case SWS::XY:

    LOG(SWS::LOG_TRACE, "[start] Computing Sxy at time-step {} on tile ({}, {}, {})", ts, ii, jj, kk);

#ifdef BLOCKWISE_FDO
    for (int i=iStart; i<iEnd; i++){
      for (int j=jStart; j<jEnd; j++){
        sigmaXY(i,j)+= (fdo_.apply<SWS::Y>(vX, i, j+1, SWS::DK) + fdo_.apply<SWS::X>(vY, i, j, SWS::DK))*mu(i,j)*dt;
      }
    }
#else
    for (int i=iStart; i<iEnd; i++){
      for (int j=jStart; j<jEnd; j++){
        for (int k=kStart; k<kEnd; k++){

          SWS::RealType sxy=0.0;

          sxy+=(c3*(vX(i,j+1,k)-vX(i,j,k)) - c2*(vX(i,j+2,k)-vX(i,j-1,k)))/ds; // Dy(vX)
          sxy+=(c3*(vY(i,j,k)-vY(i-1,j,k)) - c2*(vY(i+1,j,k)-vY(i-2,j,k)))/ds; // Dx(vY)
          sxy*=mu(i,j,k)*dt;

          sigmaXY(i,j,k)+=sxy;

        }
      }
    }
#endif

    LOG(SWS::LOG_TRACE, "[stop] Computing Sxy at time-step {} on tile ({}, {}, {})", ts, ii, jj, kk);
    
    LOG(SWS::LOG_TRACE, "||sigmaXY({},{},{},{})||^2 = {}", ts, ii, jj, kk, sigmaXY.norm2());

    break;
  case SWS::XZ:

    LOG(SWS::LOG_TRACE, "[start] Computing Sxz at time-step {} on tile ({}, {}, {})", ts, ii, jj, kk);

#ifdef BLOCKWISE_FDO
    for (int i=iStart; i<iEnd; i++){
      for (int j=jStart; j<jEnd; j++){
        sigmaXZ(i,j)+= (fdo_.apply<SWS::Z>(vX, i, j, SWS::DK+1) + fdo_.apply<SWS::X>(vZ, i, j, SWS::DK))*mu(i,j)*dt;
      }
    }
#else
    for (int i=iStart; i<iEnd; i++){
      for (int j=jStart; j<jEnd; j++){
        for (int k=kStart; k<kEnd; k++){

          SWS::RealType sxz=0.0;

          sxz+=(c3*(vX(i,j,k+1)-vX(i,j,k)) - c2*(vX(i,j,k+2)-vX(i,j,k-1)))/ds; // Dz(vX)
          sxz+=(c3*(vZ(i,j,k)-vZ(i-1,j,k)) - c2*(vZ(i+1,j,k)-vZ(i-2,j,k)))/ds; // Dx(vZ)
          sxz*=mu(i,j,k)*dt;

          sigmaXZ(i,j,k)+=sxz;

        }
      }
    }
#endif

    LOG(SWS::LOG_TRACE, "[stop] Computing Sxz at time-step {} on tile ({}, {}, {})", ts, ii, jj, kk);
    
    LOG(SWS::LOG_TRACE, "||sigmaXZ({},{},{},{})||^2 = {}", ts, ii, jj, kk, sigmaXZ.norm2());

    break;
  case SWS::YZ:

    LOG(SWS::LOG_TRACE, "[start] Computing Syz at time-step {} on tile ({}, {}, {})", ts, ii, jj, kk);

#ifdef BLOCKWISE_FDO
    for (int i=iStart; i<iEnd; i++){
      for (int j=jStart; j<jEnd; j++){
        sigmaYZ(i,j)+= (fdo_.apply<SWS::Z>(vY, i, j, SWS::DK+1) + fdo_.apply<SWS::Y>(vZ, i, j+1, SWS::DK))*mu(i,j)*dt;
      }
    }
#else
    for (int i=iStart; i<iEnd; i++){
      for (int j=jStart; j<jEnd; j++){
        for (int k=kStart; k<kEnd; k++){

          SWS::RealType syz=0.0;

          syz+=(c3*(vY(i,j,k+1)-vY(i,j,k)) - c2*(vY(i,j,k+2)-vY(i,j,k-1)))/ds; // Dz(vY)
          syz+=(c3*(vZ(i,j+1,k)-vZ(i,j,k)) - c2*(vZ(i,j+2,k)-vZ(i,j-1,k)))/ds; // Dy(vZ)
          syz*=mu(i,j,k)*dt;

          sigmaYZ(i,j,k)+=syz;

        }
      }
    }
#endif

    LOG(SWS::LOG_TRACE, "[stop] Computing Syz at time-step {} on tile ({}, {}, {})", ts, ii, jj, kk);

    LOG(SWS::LOG_TRACE, "||sigmaYZ({},{},{},{})||^2 = {}", ts, ii, jj, kk, sigmaYZ.norm2());

    break;
  default:
    LOG(SWS::LOG_ERROR, "Unknown stress component {} requested within LinearSeismicWaveModel::computeStress()", sc);
    break;
  }

  IOManager::getInstance().stop<COMPUTE_STRESS>(ts, ii, jj, kk);

  MetricsManager::getInstance()->stop("ComputeStress");

  return 0;
}

void LinearSeismicWaveModel::addVelocitySource(const SWS::Directions & d,
                                               const int ts,
                                               const int ii, const int jj, const int kk)
{
  if (!v_(d)(ii,jj,kk).hasSource()) return;

  /* Get the source coordinates */
  const auto & is=v_(d)(ii,jj,kk).is();
  const auto & js=v_(d)(ii,jj,kk).js();
  const auto & ks=v_(d)(ii,jj,kk).ks();

  const auto & rho=DataSet::getInstance()->rho()(ii,jj,kk);

  auto & vX=v_(SWS::X)(ii,jj,kk);
  auto & vY=v_(SWS::Y)(ii,jj,kk);

  const auto nsrc=pm_.sources().size();

  const auto dt=2*dt_;
  const auto l=ts/2;

  switch(d){
  case SWS::X:{

    LOG(SWS::LOG_TRACE, "[start] Updating velocity source for Vx at time-step {} on tile ({}, {}, {})", ts, ii, jj, kk);

    static const SWS::RealType a=sin(135.*M_PI/180.)*10.e6*2.0*M_PI*M_PI*7.*7.;
    for (auto s=0; s<nsrc; ++s){
      auto i=is[s];
      auto j=js[s];
      auto k=ks[s];
      vX(i,j,k)-=a*((l-1)*dt-1.2/7.)*exp(-M_PI*M_PI*7.*7.*((l-1)*dt-1.2/7.)*((l-1)*dt-1.2/7.))*dt/rho(i,j,k);
    }

    LOG(SWS::LOG_TRACE, "[stop] Updating velocity source for Vx at time-step {} on tile ({}, {}, {})", ts, ii, jj, kk);

    break;
  }
  case SWS::Y:{

    LOG(SWS::LOG_TRACE, "[start] Updating velocity source for Vy at time-step {} on tile ({}, {}, {})", ts, ii, jj, kk);

    static const SWS::RealType b=cos(135.*M_PI/180.)*10.e6*2.0*M_PI*M_PI*7.*7.;
    for (auto s=0; s<nsrc; ++s){
      auto i=is[s];
      auto j=js[s];
      auto k=ks[s];
      vY(i,j,k)-=b*((l-1)*dt-1.2/7.)*exp(-M_PI*M_PI*7.*7.*((l-1)*dt-1.2/7.)*((l-1)*dt-1.2/7.))*dt/rho(i,j,k);
    }

    LOG(SWS::LOG_TRACE, "[stop] Updating velocity source for Vy at time-step {} on tile ({}, {}, {})", ts, ii, jj, kk);

    break;
  }
  case SWS::Z:
    break;
  default:
    LOG(SWS::LOG_ERROR, "Unknown spatial direction {} requested within LinearSeismicWaveModel::addVelocitySource()", d);
    break;
  }
}

void LinearSeismicWaveModel::initialize(const int ii, const int jj, const int kk)
{
  auto * pMesh=Mesh3DPartitioning::getInstance();

  const int lii=pMesh->lii(ii);
  const int ljj=pMesh->ljj(jj);
  const int lkk=pMesh->lkk(kk);

  LOG(SWS::LOG_DEBUG, "[start] Initializing velocity and stress fields on tile ({}, {}, {})", ii, jj, kk);

  // Velocity
  for (auto d : {SWS::X, SWS::Y, SWS::Z}){
    v_(d)(lii,ljj,lkk)=0.0;
  }

  // Stress field
  for (auto sc : {SWS::XX, SWS::YY, SWS::ZZ, SWS::XY, SWS::XZ, SWS::YZ}){
    sigma_(sc)(lii,ljj,lkk)=0.0;
  }

  LOG(SWS::LOG_DEBUG, "[stop] Initializing velocity and stress fields on tile ({}, {}, {})", ii, jj, kk);
}

void LinearSeismicWaveModel::setVelocitySourceLocations()
{
  auto pMesh=Mesh3DPartitioning::getInstance();

  const auto ds=CartesianMesh3D::getInstance()->ds();

  LOG(SWS::LOG_DEBUG, "[start] Evaluating velocity source locations");

  // FIXME for the moment we assume that all SpatialBlocks have the same size
  const auto cx=pMesh->ccx()[0];
  const auto cy=pMesh->ccy()[0];
  const auto cz=pMesh->ccz()[0];

  for (auto & source : pm_.sources()){
    const auto xs=pm_.xs().at(source);
    const auto ys=pm_.ys().at(source);
    const auto zs=pm_.zs().at(source);

    // Find the coordinates of the spatial cell to which belongs the source s
    const int i=floor(xs/ds)+1;
    const int j=floor(ys/ds)+1;
    const int k=floor(zs/ds)+1;

    // Find the coordinates of the spatial block to which belongs the source s
    const auto ii=i/cx;
    const auto jj=j/cy;
    const auto kk=k/cz;

    // Relocating the source location to take into account the halo
    const auto is=i%cx + v_(0)(ii,jj,kk).iStart();
    const auto js=j%cy + v_(0)(ii,jj,kk).jStart();
    const auto ks=k%cz + v_(0)(ii,jj,kk).kStart();

    auto rank=ExecutionContext::rank();

    if (rank == Mesh3DPartitioning::getInstance()->rank_of(ii, jj, kk)){
      // The location (is,js,ks) contains a source
      LOG(SWS::LOG_INFO, "Found velocity source at cell ({}, {}, {}) within the tile ({}, {}, {})", is, js, ks, ii, jj, kk);
      for (auto d : {SWS::X, SWS::Y, SWS::Z}){
        v_(d)(ii,jj,kk).addSource(is, js, ks);
      }
    }
  } // sources

  LOG(SWS::LOG_DEBUG, "[stop] Evaluating velocity source locations");
}

LinearSeismicWaveModel::LinearSeismicWaveModel(const CentralFDOperator & fdo,
                                               const SEWASParameterManager & pm,
					       const int nt, const float tmax):fdo_(fdo),
                                                                               pm_(pm)
{
  nt_=2*nt; // Using of a staggered-grid
  tmax_=tmax;

  dt_=tmax_/nt_;

  LOG(SWS::LOG_INFO, "Adjusted time-step: Dt={}", dt_);
  LOG(SWS::LOG_INFO, "Adjusted time-step count: Nt={}", nt_);

  auto pMesh = Mesh3DPartitioning::getInstance();

  // Velocity
  for (auto d : {SWS::X, SWS::Y, SWS::Z}){
    pMesh->buildSpatialField(v_(d));
  }
  setVelocitySourceLocations();

  // Stress field
  for (auto sc : {SWS::XX, SWS::YY, SWS::ZZ, SWS::XY, SWS::XZ, SWS::YZ}){
    pMesh->buildSpatialField(sigma_(sc));
  }

  if (nullptr == HaloManager::getInstance(sigma_, v_, fdo_.hnx(), fdo_.hny(), fdo_.hnz())){
    LOG(SWS::LOG_CRITICAL, "Unable to create an instance of HaloManager. Exiting...");
    exit(SWS::OBJECT_CREATION_FAILURE);
  }

    // Build the IO manager
  IOManager::getInstance(nt_, pMesh->lnxx(), pMesh->lnyy(), pMesh->lnzz()).init();
}

LinearSeismicWaveModel::~LinearSeismicWaveModel(){
  IOManager::getInstance().finalize();

  HaloManager::releaseInstance();
}
