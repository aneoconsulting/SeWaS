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

#include "DataSet.hxx"
#include "LogManager.hxx"
#include "Mesh3DPartitioning.hxx"

DataSet* DataSet::pInstance_ = nullptr;

DataSet*
DataSet::getInstance(const SEWASParameterManager& pm)
{
  if (nullptr == pInstance_) {
    pInstance_ = new DataSet(pm);
  }
  return pInstance_;
}

DataSet*
DataSet::getInstance()
{
  assert(nullptr != pInstance_);
  return pInstance_;
}

void
DataSet::releaseInstance()
{
  if (pInstance_) {
    delete pInstance_;
    pInstance_ = nullptr;
  }
}

DataSet::DataSet(const SEWASParameterManager& pm)
  : pm_(pm)
{
  // Lamé parameters
  for (auto sc = 0; sc < SWS::NB_STRESS_FIELD_COMPONENTS; sc++) {
    Mesh3DPartitioning::getInstance()->buildSpatialField(lambda_(sc));
    Mesh3DPartitioning::getInstance()->buildSpatialField(mu_(sc));
  }

  // Density
  Mesh3DPartitioning::getInstance()->buildSpatialField(rho_);

  // Buoyancy
  for (auto d : { SWS::X, SWS::Y, SWS::Z }) {
    Mesh3DPartitioning::getInstance()->buildSpatialField(b_(d));
  }
}

DataSet::~DataSet() {}

const auto
DataSet::getLayer(const int k, const int ii, const int jj, const int kk) const
{
  /* Get the layer to which belongs the plane of order k within the tile (ii,jj,kk) */

  // FIXME we assume a uniform tile size
  const auto cz = pm_.cz(); /* Actual z-size of the tile without halo and padding */

  const auto lz = (kk * cz + k) * pm_.ds();

  for (auto& layer : pm_.layers()) {
    const auto start = pm_.start().at(layer)[SWS::Z];
    const auto end = pm_.end().at(layer)[SWS::Z];

    if (start <= lz && lz < end)
      return layer;
  }

  LOG(SWS::LOG_CRITICAL,
      "The plane {} from the tile {} does not belong to any layer. Check the topology file configuration. "
      "Exiting...");

  exit(SWS::BAD_TOPOLOGY_FILE_CONFIGURATION);
}

void
DataSet::initialize(const int ii, const int jj, const int kk)
{
  LOG(SWS::LOG_INFO, "[start] Initializing Lamé parameters on tile ({}, {}, {})", ii, jj, kk);

  auto pMesh = Mesh3DPartitioning::getInstance();

  const int lii = pMesh->lii(ii);
  const int ljj = pMesh->ljj(jj);
  const int lkk = pMesh->lkk(kk);

  auto& _rho = rho_(lii, ljj, lkk);

  for (int k = _rho.kStart(); k < _rho.kEnd(); k++) {

    const auto layer = getLayer(k - _rho.kStart(), ii, jj, kk);

    const auto& vp = pm_.Vp().at(layer);
    const auto& vs = pm_.Vs().at(layer);
    const auto& rho = pm_.rho().at(layer);

    /*
      vs=sqrt(mu/rho)
      vp=sqrt((lambda+2*mu)/rho)
    */
    auto vp2 = vp * vp;
    auto vs2 = vs * vs;

    for (int i = _rho.iStart(); i < _rho.iEnd(); i++) {
      for (int j = _rho.jStart(); j < _rho.jEnd(); j++) {

        _rho(i, j, k) = rho;

        // Initialize Lamé parameters for the current layer
        for (auto sc = 0; sc < SWS::NB_STRESS_FIELD_COMPONENTS; sc++) {
          auto& _mu = mu_(sc)(lii, ljj, lkk);
          auto& _lambda = lambda_(sc)(lii, ljj, lkk);

          _mu(i, j, k) = _rho(i, j, k) * vs2;
          _lambda(i, j, k) = _rho(i, j, k) * vp2 - 2 * _mu(i, j, k);
        }

        for (auto d = 0; d < SWS::DIM; d++) {
          b_(d)(lii, ljj, lkk)(i, j, k) = 1. / rho_(lii, ljj, lkk)(i, j, k);
        }

      } // j
    }   // i

  } // k

  LOG(SWS::LOG_INFO, "[stop] Initializing Lamé parameters on tile ({}, {}, {})", ii, jj, kk);
}
