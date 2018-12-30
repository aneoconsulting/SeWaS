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

#include "ExecutionContext.hxx"
#include "SEWASPaRSEC.hxx"
#include "LogManager.hxx"

#ifdef SEWAS_DISTRIBUTED
int ExecutionContext::world_=1;
int ExecutionContext::rank_=0;
#endif

int ExecutionContext::init(SEWASParameterManager & pm)
{
#ifdef SEWAS_DISTRIBUTED
  /* Start the MPI runtime */
  LOG(SWS::INFO, "Starting the MPI runtime");

  int provided;
  MPI_Init_thread(&pm.argc(), &pm.argv(), MPI_THREAD_SERIALIZED, &provided);
  if (MPI_THREAD_SERIALIZED != provided)
    LOG(SWS::WARN, "The thread level supported by the used MPI implementation is {}", provided);

  LOG(SWS::INFO, "MPI runtime is started");

  MPI_Comm_size(MPI_COMM_WORLD, &world_);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank_);

  LOG(SWS::INFO, "MPI size: {}", world_);
  LOG(SWS::INFO, "MPI rank: {}", rank_);
#endif

#ifdef SEWAS_WITH_PARSEC
  /* PaRSEC initialization */
  LOG(SWS::INFO, "Starting the PaRSEC runtime");

  SEWASPaRSEC::init(pm);

  LOG(SWS::INFO, "PaRSEC runtime is started");
#endif

  return 0;
}
