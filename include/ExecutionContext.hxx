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

#ifdef SEWAS_DISTRIBUTED
#include <mpi.h>
#endif

#ifdef SEWAS_WITH_PARSEC
#include <parsec/parsec_config.h>
#include <parsec.h>
#include <parsec/data_distribution.h>
#include <parsec/datatype.h>
#include <parsec/utils/mca_param.h>
#include <parsec/scheduling.h>
#include <parsec/arena.h>

#include "sewas.h"
#include "SEWASPaRSEC.hxx"
#endif

#include "SEWASParameterManager.hxx"

class ExecutionContext
{
public:
  ExecutionContext()
  {
  }

  ~ExecutionContext()
  {
  }

  static int init(SEWASParameterManager & pm);

  static inline auto rank()
  {
    int r=0;
#if SEWAS_DISTRIBUTED
    MPI_Comm_rank(MPI_COMM_WORLD, &r);
#endif
    return r;
  }

  static inline auto world()
  {
    int w=0;
#if SEWAS_DISTRIBUTED
    MPI_Comm_size(MPI_COMM_WORLD, &w);
#endif
    return w;
  }

  static inline void barrier()
  {
#if SEWAS_DISTRIBUTED
    MPI_Barrier(MPI_COMM_WORLD);
#endif
  }

  static inline void finalize()
  {
#ifdef SEWAS_WITH_PARSEC
    SEWASPaRSEC::finalize();
#endif

#if SEWAS_DISTRIBUTED
  MPI_Finalize();
#endif
  }

private:
};