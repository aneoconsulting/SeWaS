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
#include <memory>
#include <string>
#include <thread>

#include "CartesianMesh3D.hxx"
#include "Config.hxx"
#include "DataSet.hxx"
#include "ExecutionContext.hxx"
#include "ExternalSource.hxx"
#include "LinearSeismicWaveModel.hxx"
#include "LogManager.hxx"
#include "Mesh3DPartitioning.hxx"
#include "SEWASParameterManager.hxx"
#ifdef VISUALIZE_EXECUTION
#include "VisualizationManager.hxx"
#endif
#include "MetricsManager.hxx"

#ifdef VISUALIZE_EXECUTION
bool bContinueRendering = true;
std::thread tRenderer;

void
render()
{
  if (nullptr == VisualizationManager::getInstance()) {
    LOG("Visualization Manager is not yet instantiated in render(). Exiting...");
    exit(SWS::INSTANCE_ACCESS_VIOLATION);
  }

  while (bContinueRendering) {
    VisualizationManager::getInstance()->render();
  }
}
#endif

int
main(int argc, char* argv[])
{
  int status = 0;

  // Create the logger
  if (nullptr == LogManager::getInstance()) {
    std::cerr << "Failed to create the logger. Exiting...\n";
    return SWS::OBJECT_CREATION_FAILURE;
  }
  auto pLogger = LogManager::getInstance();

  /* Create a generic metrics manager */
  if (nullptr ==
      MetricsManager::getInstance(
        "SEWAS", { "Global", "Core simulation", "Initialization", "ComputeVelocity", "ComputeStress" })) {
    LOG(SWS::LOG_CRITICAL, "Unable to create the metrics manager. Exiting...");
    return SWS::OBJECT_CREATION_FAILURE;
  }
  auto mm = MetricsManager::getInstance();

  mm->start("Global");

  mm->start("Initialization");

  /* Intialize application parameters */
  auto pm = *std::make_unique<SEWASParameterManager>(&argc, &argv);

  const auto tmax = pm.tmax();
  const auto nt = pm.nt();
  const auto ds = pm.ds();
  const auto nx = pm.nx(), ny = pm.ny(), nz = pm.nz();
  const auto cx = pm.cx(), cy = pm.cy(), cz = pm.cz();
  const auto P = pm.P(), Q = pm.Q(), R = pm.R();
  const auto nthreads = pm.nthreads();

  ExecutionContext::init(pm);

  // Create the computational domain
  if (nullptr == CartesianMesh3D::getInstance(nx, ny, nz, ds)) {
    LOG(SWS::LOG_CRITICAL, "Unable to create the mesh. Exiting...");
    return SWS::OBJECT_CREATION_FAILURE;
  }
  auto pCartesianMesh = CartesianMesh3D::getInstance();

  // Create a finite difference operator
  CentralFDOperator fdo(pCartesianMesh->dx(), pCartesianMesh->dy(), pCartesianMesh->dz());

  // Create the domain partitioning
  if (nullptr == Mesh3DPartitioning::getInstance(cx, cy, cz, fdo.hnx(), fdo.hny(), fdo.hnz(), P, Q, R)) {
    LOG(SWS::LOG_CRITICAL, "Unable to partition the mesh. Exiting...");
    return SWS::OBJECT_CREATION_FAILURE;
  }
  auto pMeshPartitioning = Mesh3DPartitioning::getInstance();

  // Create a dataset:
  /* - Lamé parameters (lambda, mu)
     - Material density (rho)
     - velocities of primary and seconday waves
  */
  if (nullptr == DataSet::getInstance(pm)) {
    LOG(SWS::LOG_CRITICAL, "Unable to create the dataset. Exiting...");
    return SWS::OBJECT_CREATION_FAILURE;
  }
  auto pDataSet = DataSet::getInstance();

  // Create the seismic wave model
  if (nullptr == LinearSeismicWaveModel::getInstance(fdo, pm, nt, tmax)) {
    LOG(SWS::LOG_CRITICAL, "Unable to create the linear seismic wave model. Exiting...");
    return SWS::OBJECT_CREATION_FAILURE;
  }
  auto pSEWAS = LinearSeismicWaveModel::getInstance();

  auto rank = ExecutionContext::rank();

#ifdef VISUALIZE_EXECUTION
  if (0 == rank) {
    // Create and start the visualization engine : only the master process
    if (nullptr == VisualizationManager::getInstance(pSEWAS->nt(),
                                                     pMeshPartitioning->nxx(),
                                                     pMeshPartitioning->nyy(),
                                                     pMeshPartitioning->nzz(),
                                                     nthreads)) {
      LOG(SWS::CRITICAL, "Unable to create the visualization engine. Exiting...");
      return SWS::OBJECT_CREATION_FAILURE;
    }
    tRenderer = std::thread(render);
  }
#endif

  ExecutionContext::barrier();

  mm->stop("Initialization");

  mm->start("Core simulation");

  // Let's rock!
  pSEWAS->propagate();

  mm->stop("Core simulation");

  LOG(SWS::LOG_INFO, "MPI process {} completed the simulation", rank);

#ifdef VISUALIZE_EXECUTION
  if (0 == rank) {
    bContinueRendering = false;
    tRenderer.join();
  }
#endif

  ExecutionContext::barrier();

  if (0 == rank) {
    LOG(SWS::LOG_INFO, "||Vx(0,0,0)||^2 = {}", pSEWAS->v(SWS::X)(0, 0, 0).norm2());
    LOG(SWS::LOG_INFO, "||Vy(0,0,0)||^2 = {}", pSEWAS->v(SWS::Y)(0, 0, 0).norm2());
    LOG(SWS::LOG_INFO, "||Vz(0,0,0)||^2 = {}", pSEWAS->v(SWS::Z)(0, 0, 0).norm2());
  }

  // Release memory
  CartesianMesh3D::releaseInstance();
  DataSet::releaseInstance();
  ExternalSource::releaseInstance();
  LinearSeismicWaveModel::releaseInstance();
#ifdef VISUALIZE_EXECUTION
  if (0 == rank) {
    VisualizationManager::releaseInstance();
  }
#endif

  mm->stop("Global");

  if (0 == rank) {
    mm->show();
  }

  MetricsManager::releaseInstance();

  ExecutionContext::finalize();

  LogManager::releaseInstance();

  return status;
}
