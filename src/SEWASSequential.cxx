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

#include "HaloManager.hxx"
#include "LinearSeismicWaveModel.hxx"
#include "Mesh3DPartitioning.hxx"
#include "SEWASSequential.hxx"

SEWASSequential* SEWASSequential::pInstance_ = nullptr;

SEWASSequential*
SEWASSequential::getInstance(const int nt, const int nxx, const int nyy, const int nzz)
{
  if (nullptr == pInstance_) {
    pInstance_ = new SEWASSequential(nt, nxx, nyy, nzz);
    return pInstance_;
  } else {
    return pInstance_;
  }
}

void
SEWASSequential::releaseInstance()
{
  if (pInstance_) {
    delete pInstance_;
    pInstance_ = nullptr;
  }
}

int
SEWASSequential::run()
{
  int status = 0;

  for (int ii = 0; ii < nxx_; ii++) {
    for (int jj = 0; jj < nyy_; jj++) {
      for (int kk = 0; kk < nzz_; kk++) {
        LinearSeismicWaveModel::initializeFieldsWrapper(ii, jj, kk);
      }
    }
  }

  auto& vH_S = vH_[SEND];
  auto& vH_R = vH_[RECEIVE];

  auto& sigmaH_S = sigmaH_[SEND];
  auto& sigmaH_R = sigmaH_[RECEIVE];

  for (int ts = 2; ts <= nt_ - 2; ts += 2) {

    LOG(SWS::LOG_TRACE, "[start] Processing time-step {}", ts);

    for (int ii = 0; ii < nxx_; ii++) {
      for (int jj = 0; jj < nyy_; jj++) {
        for (int kk = 0; kk < nzz_; kk++) {

          for (auto l : { SWS::LEFT, SWS::RIGHT, SWS::BACKWARD, SWS::FORWARD, SWS::BOTTOM, SWS::TOP }) {
            for (auto sc : { SWS::XX, SWS::YY, SWS::ZZ, SWS::XY, SWS::XZ, SWS::YZ }) {
              HaloManager::updateStressWrapper(l, sc, ts - 1, ii, jj, kk, sigmaH_R(sc)(ii, jj, kk)(l));
            }
          }

          for (auto d : { SWS::X, SWS::Y, SWS::Z }) {
            LinearSeismicWaveModel::computeVelocityWrapper(d, ts, ii, jj, kk);

            for (auto l : { SWS::LEFT, SWS::RIGHT, SWS::BACKWARD, SWS::FORWARD, SWS::BOTTOM, SWS::TOP }) {
              HaloManager::extractVelocityHaloWrapper(l, d, ts, ii, jj, kk, vH_S(d)(ii, jj, kk)(l));

              sendreceive(vH_S, d, ii, jj, kk, l);
            }
          }

        } // kk
      }   // jj
    }     // ii

    for (int ii = 0; ii < nxx_; ii++) {
      for (int jj = 0; jj < nyy_; jj++) {
        for (int kk = 0; kk < nzz_; kk++) {

          for (auto l : { SWS::LEFT, SWS::RIGHT, SWS::BACKWARD, SWS::FORWARD, SWS::BOTTOM, SWS::TOP }) {
            for (auto d : { SWS::X, SWS::Y, SWS::Z }) {
              HaloManager::updateVelocityWrapper(l, d, ts, ii, jj, kk, vH_R(d)(ii, jj, kk)(l));
            }
          }

          for (auto sc : { SWS::XX, SWS::YY, SWS::ZZ, SWS::XY, SWS::XZ, SWS::YZ }) {
            LinearSeismicWaveModel::computeStressWrapper(sc, ts + 1, ii, jj, kk);

            for (auto l : { SWS::LEFT, SWS::RIGHT, SWS::BACKWARD, SWS::FORWARD, SWS::BOTTOM, SWS::TOP }) {
              HaloManager::extractStressHaloWrapper(l, sc, ts + 1, ii, jj, kk, sigmaH_S(sc)(ii, jj, kk)(l));

              sendreceive(sigmaH_S, sc, ii, jj, kk, l);
            }
          }

        } // kk
      }   // jj
    }     // ii

    LOG(SWS::LOG_TRACE, "[stop] Processing time-step {}", ts);

  } // ts

  return status;
}

SEWASSequential::SEWASSequential(const int nt, const int nxx, const int nyy, const int nzz)
  : nt_(nt)
  , nxx_(nxx)
  , nyy_(nyy)
  , nzz_(nzz)
{

  for (auto sc : { SWS::XX, SWS::YY, SWS::ZZ, SWS::XY, SWS::XZ, SWS::YZ }) {
    for (auto u : { SEND, RECEIVE }) {
      buildArenas<decltype(sigmaH_[u])>(sigmaH_[u], sc);
    }
  }

  for (auto d : { SWS::X, SWS::Y, SWS::Z }) {
    for (auto u : { SEND, RECEIVE }) {
      buildArenas<decltype(vH_[u])>(vH_[u], d);
    }
  }
}

SEWASSequential::~SEWASSequential()
{

  for (auto sc : { SWS::XX, SWS::YY, SWS::ZZ, SWS::XY, SWS::XZ, SWS::YZ }) {
    for (auto u : { SEND, RECEIVE }) {
      destroyArenas<decltype(sigmaH_[u])>(sigmaH_[u], sc);
    }
  }

  for (auto d : { SWS::X, SWS::Y, SWS::Z }) {
    for (auto u : { SEND, RECEIVE }) {
      destroyArenas<decltype(vH_[u])>(vH_[u], d);
    }
  }
}
