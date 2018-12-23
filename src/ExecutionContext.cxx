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

int ExecutionContext::init(SEWASParameterManager & pm)
{
#if SEWAS_DISTRIBUTED
  /* Start the MPI runtime */

  int provided;
  MPI_Init_thread(&pm.argc(), &pm.argv(), MPI_THREAD_SERIALIZED, &provided);
  if (MPI_THREAD_SERIALIZED != provided)
    std::cerr << "WARNING: The thread level supported by the MPI implementation is " << provided << "\n";
#endif

#ifdef SEWAS_WITH_PARSEC
  /* PaRSEC initialization */

  SEWASPaRSEC::init(pm);
#endif

  return 0;
}
