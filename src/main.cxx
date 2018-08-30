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
#include <string>
#include <memory>
#include <thread>

#include <mpi.h>

#ifdef SEWAS_WITH_PARSEC
#include <parsec/parsec_config.h>
#include <parsec.h>
#include <parsec/data_distribution.h>
#include <parsec/datatype.h>
#include <parsec/utils/mca_param.h>
#include <parsec/scheduling.h>

#include "sewas.h"
#endif

#include "Config.hxx"
#include "SEWASParameterManager.hxx"
#include "CartesianMesh3D.hxx"
#include "DataSet.hxx"
#include "ExternalSource.hxx"
#include "LinearSeismicWaveModel.hxx"
#include "Mesh3DPartitioning.hxx"
#ifdef VISUALIZE_EXECUTION
#include "VisualizationManager.hxx"
#endif
#include "MetricsManager.hxx"


#ifdef VISUALIZE_EXECUTION
bool bContinueRendering = true;
std::thread tRenderer;

void render()
{
  if (nullptr == VisualizationManager::getInstance()){
    std::cerr << "Visualization Manager is not initialized. Exiting...\n";
    exit(-1);
  }

  while(bContinueRendering){
    VisualizationManager::getInstance()->render();
  }
}
#endif

#ifdef SEWAS_WITH_PARSEC
parsec_context_t * g_parsec;
#endif

int main (int argc, char* argv[])
{
  int status=0;

  /* Create a generic metrics manager */
  if (nullptr == MetricsManager::getInstance("SEWAS")){
    std::cerr << "Unable to create the metrics manager. Exiting...\n";
    return -1;
  }
  auto mm=MetricsManager::getInstance();

  mm->setSampledEvents({"Global", "Core simulation", "Initialization", "ComputeVelocity", "ComputeStress"});

  mm->start("Global");

  mm->start("Initialization");

  /* Intialize application parameters */
  auto pm=*std::make_unique<SEWASParameterManager>();
  pm.parse(argc, argv);

  const auto tmax=pm.tmax();
  const auto nt=pm.nt();
  const auto ds=pm.ds();
  const auto nx=pm.nx(), ny=pm.ny(), nz=pm.nz();
  const auto cx=pm.cx(), cy=pm.cy(), cz=pm.cz();
  const auto P=pm.P(), Q=pm.Q(), R=pm.R();
  const auto nthreads=pm.nthreads();

  /* Start the MPI runtime */
  {
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_SERIALIZED, &provided);
    if (MPI_THREAD_SERIALIZED != provided)
      std::cerr << "WARNING: The thread level supported by the MPI implementation is " << provided << "\n";
  }

#ifdef SEWAS_WITH_PARSEC
  /* PaRSEC initialization */

  int parsec_argc=0;
  char ** parsec_argv=(char **) calloc(argc, sizeof(char *));

  parsec_argv[parsec_argc++]=argv[0]; /* App name */

  for (int i=1; i<argc; i++){
    if (0 == strcmp("--", argv[i])){
      /* We are done reading the application arguments;
         all the remaining arguments will be passed to the PaRSEC engine */
      for (int j=i+1; j<argc; j++){
        parsec_argv[parsec_argc++]=argv[j];
      }
      break;
    }
  }

  g_parsec = parsec_init(nthreads, &parsec_argc, &parsec_argv);

  free(parsec_argv);

#endif

  // Create the computational domain
  if (nullptr == CartesianMesh3D::getInstance(nx, ny, nz, ds)){
    std::cerr << "Unable to create the mesh. Exiting...\n";
    return -1;
  }
  auto pCartesianMesh=CartesianMesh3D::getInstance();

  // Create a finite difference operator
  CentralFDOperator fdo(pCartesianMesh->dx(), pCartesianMesh->dy(), pCartesianMesh->dz());

  // Create the domain partitioning
  if (nullptr == Mesh3DPartitioning::getInstance(cx, cy, cz,
                                                 fdo.hnx(), fdo.hny(), fdo.hnz(),
                                                 P, Q, R)){
    std::cerr << "Unable to partition the mesh. Exiting...\n";
    return -1;
  }
  auto pMeshPartitioning=Mesh3DPartitioning::getInstance();

  // Create a dataset:
  /* - Lamé parameters (lambda, mu)
     - Material density (rho)
     - velocities of primary and seconday waves
  */
  if (nullptr == DataSet::getInstance(pm)){
    std::cerr << "Unable to create the dataset. Exiting...\n";
    return -1;
  }
  auto pDataSet=DataSet::getInstance();

  // Create the seismic wave model
  if (nullptr == LinearSeismicWaveModel::getInstance(fdo, pm, nt, tmax)){
    std::cerr << "Unable to create the linear seismic wave model. Exiting...\n";
    return -1;
  }
  auto pSEWAS=LinearSeismicWaveModel::getInstance();

  int rank=0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

#ifdef VISUALIZE_EXECUTION
  if (0 == rank){
    // Create and start the visualization engine : only the master process
    if (nullptr == VisualizationManager::getInstance(pSEWAS->nt(),
                                                     pMeshPartitioning->nxx(), pMeshPartitioning->nyy(), pMeshPartitioning->nzz(),
                                                     nthreads)){
      std::cerr << "Unable to create the visualization engine. Exiting...\n";
      return -1;
    }
    tRenderer = std::thread(render);
  }
#endif

  MPI_Barrier(MPI_COMM_WORLD);

  mm->stop("Initialization");

  mm->start("Core simulation");

  // Let's rock!
  pSEWAS->propagate();

  mm->stop("Core simulation");

#ifdef VISUALIZE_EXECUTION
  if (0 == rank){
    bContinueRendering = false;
    tRenderer.join();
  }
#endif

  std::cout << "I am " << rank << " and I completed the simulation...\n";

  MPI_Barrier(MPI_COMM_WORLD);

#if 0
  const auto & Sxx=pSEWAS->sigma(SWS::XX);
  const auto & Syy=pSEWAS->sigma(SWS::YY);
  const auto & Szz=pSEWAS->sigma(SWS::ZZ);
  const auto & Sxy=pSEWAS->sigma(SWS::XY);
  const auto & Sxz=pSEWAS->sigma(SWS::XZ);
  const auto & Syz=pSEWAS->sigma(SWS::YZ);

  const auto & Vx=pSEWAS->v(SWS::X);
  const auto & Vy=pSEWAS->v(SWS::Y);
  const auto & Vz=pSEWAS->v(SWS::Z);

  // int j0=Vx(0,0,0).js()[0]; // output plane
  // int k0=Vx(0,0,0).ks()[0]; // output plane

  // Vx(0,0,0).plot2D(k0, std::string("Vx-0.0.0") + "-" + std::to_string(rank));

  // Vx(0,0,0).plot1D(j0, k0, std::string("Vx-0.0.0") + "-" + std::to_string(rank));
  // Vy(0,0,0).plot1D(j0, k0, std::string("Vy-0.0.0") + "-" + std::to_string(rank));
  // Vz(0,0,0).plot1D(j0, k0, std::string("Vz-0.0.0") + "-" + std::to_string(rank));

  // Sxx(0,0,0).plot1D(j0, k0, std::string("Sxx-0.0.0") + "-" + std::to_string(rank));
  // Syy(0,0,0).plot1D(j0, k0, std::string("Syy-0.0.0") + "-" + std::to_string(rank));
  // Szz(0,0,0).plot1D(j0, k0, std::string("Szz-0.0.0") + "-" + std::to_string(rank));
  // Sxy(0,0,0).plot1D(j0, k0, std::string("Sxy-0.0.0") + "-" + std::to_string(rank));
  // Sxz(0,0,0).plot1D(j0, k0, std::string("Sxz-0.0.0") + "-" + std::to_string(rank));
  // Syz(0,0,0).plot1D(j0, k0, std::string("Syz-0.0.0") + "-" + std::to_string(rank));
#endif

  // Release memory
  CartesianMesh3D::releaseInstance();
  DataSet::releaseInstance();
  ExternalSource::releaseInstance();
  LinearSeismicWaveModel::releaseInstance();
#ifdef VISUALIZE_EXECUTION
  if (0 == rank){
    VisualizationManager::releaseInstance();
  }
#endif

#ifdef SEWAS_WITH_PARSEC
  status=parsec_fini(&g_parsec);
  PARSEC_CHECK_ERROR(status, "parsec_fini");
#endif

  MPI_Finalize();

  mm->stop("Global");

  if (0 == rank){
    mm->show();
  }

  MetricsManager::releaseInstance();

  return status;
}
