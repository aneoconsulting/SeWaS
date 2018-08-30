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

#include <algorithm>

#include "Config.hxx"
#include "HaloManager.hxx"

class SEWASSequential{
public:
  static SEWASSequential * getInstance(const int nt=1,
                                      const int nxx=1, const int nyy=1, const int nzz=1);
  static void releaseInstance();

  int run();

private:
  SEWASSequential(const int nt,
                 const int nxx, const int nyy, const int nzz);
  ~SEWASSequential();

  enum DataTransferStages { SEND, RECEIVE };

  template<typename H>
  inline void buildArenas(H & halo, const int u) noexcept
  {
    halo(u).resize(nxx_,nyy_,nzz_);
    for (int ii=0; ii<nxx_; ii++){
      for (int jj=0; jj<nyy_; jj++){
        for (int kk=0; kk<nzz_; kk++){

          for (auto l : {SWS::LEFT, SWS::RIGHT, SWS::BACKWARD, SWS::FORWARD, SWS::BOTTOM, SWS::TOP}){
            const auto hsize=HaloManager::getInstance()->getHaloSize(l, ii, jj, kk);
            halo(u)(ii,jj,kk)(l) = new SWS::RealType[hsize];
            memset(halo(u)(ii,jj,kk)(l), 0, hsize);
          }

        } // kk
      } // jj
    } // ii
  }

  template<typename H>
  inline void destroyArenas(H & halo, const int u) noexcept
  {
    halo(u).resize(nxx_,nyy_,nzz_);
    for (int ii=0; ii<nxx_; ii++){
      for (int jj=0; jj<nyy_; jj++){
        for (int kk=0; kk<nzz_; kk++){

          for (auto l : {SWS::LEFT, SWS::RIGHT, SWS::BACKWARD, SWS::FORWARD, SWS::BOTTOM, SWS::TOP}){
            if (halo(u)(ii,jj,kk)(l)){
              delete[] halo(u)(ii,jj,kk)(l);
              halo(u)(ii,jj,kk)(l)=nullptr;
            }
          }

        } // kk
      } // jj
    } // ii
  }


  // FIXME use metaprogramming to simplify the following sendreceive overloaded methods

  inline int sendreceive(const SWS::VelocityHalo & vH_S,
                         const int d,
                         const int ii, const int jj, const int kk,
                         const SWS::Locations l) noexcept
  {
    int status=0;

    auto & vH_R=vH_[RECEIVE];

    const auto hsize=HaloManager::getInstance()->getHaloSize(l, ii, jj, kk)*sizeof(SWS::RealType);

    switch(l){
    case SWS::LEFT:
      if (ii > 0){
        memcpy(vH_R(d)(ii-1,jj,kk)(SWS::RIGHT), vH_S(d)(ii,jj,kk)(l), hsize);
      }
      break;
    case SWS::RIGHT:
      if (ii < nxx_-1){
        memcpy(vH_R(d)(ii+1,jj,kk)(SWS::LEFT), vH_S(d)(ii,jj,kk)(l), hsize);
      }
      break;
    case SWS::BACKWARD:
      if (jj > 0){
        memcpy(vH_R(d)(ii,jj-1,kk)(SWS::FORWARD), vH_S(d)(ii,jj,kk)(l), hsize);
      }
      break;
    case SWS::FORWARD:
      if (jj < nyy_-1){
        memcpy(vH_R(d)(ii,jj+1,kk)(SWS::BACKWARD), vH_S(d)(ii,jj,kk)(l), hsize);
      }
      break;
    case SWS::BOTTOM:
      if (kk > 0){
        memcpy(vH_R(d)(ii,jj,kk-1)(SWS::TOP), vH_S(d)(ii,jj,kk)(l), hsize);
      }
      break;
    case SWS::TOP:
      if (kk < nzz_-1){
        memcpy(vH_R(d)(ii,jj,kk+1)(SWS::BOTTOM), vH_S(d)(ii,jj,kk)(l), hsize);
      }
      break;
    default:
      fprintf(stdout, "send(%d, %d, %d, %d) - Unknown location (%d)\n", ii, jj, kk, l, l);
    }

    return status;
  }

  inline int sendreceive(const SWS::StressFieldHalo & sigmaH_S,
                         const int sc,
                         const int ii, const int jj, const int kk,
                         const SWS::Locations l) noexcept
  {
    int status=0;

    auto & sigmaH_R=sigmaH_[RECEIVE];

    const auto hsize=HaloManager::getInstance()->getHaloSize(l, ii, jj, kk)*sizeof(SWS::RealType);

    switch(l){
    case SWS::LEFT:
      if (ii > 0){
        memcpy(sigmaH_R(sc)(ii-1,jj,kk)(SWS::RIGHT), sigmaH_S(sc)(ii,jj,kk)(l), hsize);
      }
      break;
    case SWS::RIGHT:
      if (ii < nxx_-1){
        memcpy(sigmaH_R(sc)(ii+1,jj,kk)(SWS::LEFT), sigmaH_S(sc)(ii,jj,kk)(l), hsize);
      }
      break;
    case SWS::BACKWARD:
      if (jj > 0){
        memcpy(sigmaH_R(sc)(ii,jj-1,kk)(SWS::FORWARD), sigmaH_S(sc)(ii,jj,kk)(l), hsize);
      }
      break;
    case SWS::FORWARD:
      if (jj < nyy_-1){
        memcpy(sigmaH_R(sc)(ii,jj+1,kk)(SWS::BACKWARD), sigmaH_S(sc)(ii,jj,kk)(l), hsize);
      }
      break;
    case SWS::BOTTOM:
      if (kk > 0){
        memcpy(sigmaH_R(sc)(ii,jj,kk-1)(SWS::TOP), sigmaH_S(sc)(ii,jj,kk)(l), hsize);
      }
      break;
    case SWS::TOP:
      if (kk < nzz_-1){
        memcpy(sigmaH_R(sc)(ii,jj,kk+1)(SWS::BOTTOM), sigmaH_S(sc)(ii,jj,kk)(l), hsize);
      }
      break;
    default:
      fprintf(stdout, "send(%d, %d, %d, %d) - Unknown location (%d)\n", ii, jj, kk, l, l);
    }

    return status;
  }


  static SEWASSequential * pInstance_;

  int nt_;

  int nxx_;
  int nyy_;
  int nzz_;


  // Allocate send/receive temporary buffers
  SWS::StressFieldHalo sigmaH_[2];
  SWS::VelocityHalo vH_[2];
};
