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

#include "SEWASParameterManager.hxx"

SEWASParameterManager::SEWASParameterManager(int* pargc, char*** pargv)
{
  pargc_ = pargc;
  pargv_ = pargv;

  parse();

  LOG(SWS::LOG_INFO, "Geophysical data model: {}", model_);
  LOG(SWS::LOG_INFO, "Spatial discretization step: ds = {}", ds_);
  LOG(SWS::LOG_INFO, "Temporal discretization step: dt = {}", dt_);
  LOG(SWS::LOG_INFO, "Initial global domain size: ({}, {}, {})", initial_lx_, initial_ly_, initial_lz_);
  LOG(SWS::LOG_INFO, "Adjusted global domain size: ({}, {}, {})", lx_, ly_, lz_);
  LOG(SWS::LOG_INFO, "Global spatial mesh size: (nx, ny, nz) = ({}, {}, {})", nx_, ny_, nz_);
  LOG(SWS::LOG_INFO, "Spatial blocks size: (cx, cy, cz) = ({}, {}, {})", cx(), cy(), cz());
  LOG(SWS::LOG_INFO, "Global spatial blocks count: (nxx, nyy, nzz) = ({}, {}, {})", nxx_, nyy_, nzz_);
  LOG(SWS::LOG_INFO, "Local spatial blocks count: (lnxx, lnyy, lnzz) = ({}, {}, {})", lnxx_, lnyy_, lnzz_);

  LOG(SWS::LOG_INFO, "Thread count: {}", nthreads());
}

SEWASParameterManager::~SEWASParameterManager() {}

int
SEWASParameterManager::parse()
{
  /* Parse command line arguments */
  parseArgs();

  /* Read the geophysic data file and evaluate partitioning parameters */
  parseDataFile();

  return 0;
}
