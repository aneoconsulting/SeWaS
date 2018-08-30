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

#ifdef SEWAS_WITH_PARSEC
#include <vector>

#include <parsec/parsec_config.h>
#include <parsec.h>
#include <parsec/data_distribution.h>

#include "Config.hxx"

#include "sewas.h"

class SEWASPaRSEC{
public:
  static SEWASPaRSEC * getInstance(const int nt=1,
				  const int nxx=1, const int nyy=1, const int nzz=1);
  static void releaseInstance();

  int run();

private:
  SEWASPaRSEC(const int nt,
	     const int nxx, const int nyy, const int nzz);
  ~SEWASPaRSEC();

  void buildEngine();
  void buildDataDescriptor();
  void buildDAG();
  void enqueueDAG();

  void addArena(const short arena_idx, const SWS::Locations l=SWS::NB_LOCATIONS);


  static SEWASPaRSEC * pInstance_;

  int world_;
  int rank_;

  int nt_;

  int nxx_;
  int nyy_;
  int nzz_;

  parsec_context_t         * pEngine_;
  parsec_data_collection_t * pDDesc_;
  parsec_sewas_taskpool_t   * pDAG_;
};
#endif
