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

#include <atomic>

#include "Config.hxx"
#include "Constants.hxx"
#include "ExecutionContext.hxx"
#include "Mesh3DPartitioning.hxx"
#include "SpatialBlockField.hxx"
#include "Types.hxx"

#ifdef ENABLE_IO
#include <adios2.h>
#endif

class IOManager
{
public:
  static IOManager& getInstance(const int nt = 1, const int lnxx = 1, const int lnyy = 1, const int lnzz = 1);
  static void releaseInstance();

  template<short TaskType>
  inline void start(const int ts, const int ii, const int jj, const int kk)
  {
#ifdef ENABLE_IO

    auto performPuts = [&](const int taskCount, std::mutex& guard, adios2::Engine& writer) {
      LOG(SWS::LOG_TRACE, "processedTasks_({})({},{},{}) = {}", ts, ii, jj, kk, processedTasks_[ts]);

      std::lock_guard<std::mutex> lock(guard);

      if (taskCount == processedTasks_[ts - 2] ||
          (processedTasks_[ts - 2] >= 0 && -1 == processedTasks_[ts])) {
        /*
            We are about to start processing tasks from time-step ts while some
            tasks from time-step ts - 2 are still to be processed. Thus,
            we cannot close the time-step ts - 2. Instead of that,
            we force the writer to perform all pending "put" operations in order
            to avoid overwriting the buffers holding the data.
        */

        LOG(SWS::LOG_TRACE,
            "[PerformPuts] processedTasks_({}) = {}, processedTasks_({}) = {}",
            ts,
            processedTasks_[ts],
            ts - 2,
            processedTasks_[ts - 2]);

        writer.PerformPuts();

        if (taskCount == processedTasks_[ts - 2])
          processedTasks_[ts - 2] = -2;

        if (-1 == processedTasks_[ts])
          processedTasks_[ts] = 0;
      }
    };

    if constexpr (COMPUTE_VELOCITY == TaskType) {
      auto numberOfVelocityTasks = 3 * lnxx_ * lnyy_ * lnzz_;
      performPuts(numberOfVelocityTasks, velocityGuard_, velocityWriter_);
    } else if constexpr (COMPUTE_STRESS == TaskType) {
      auto numberOfStressTasks = 6 * lnxx_ * lnyy_ * lnzz_;
      performPuts(numberOfStressTasks, stressGuard_, stressWriter_);
    }
#endif
  }

  template<short TaskType>
  inline void stop(const int ts, const int ii, const int jj, const int kk)
  {
#ifdef ENABLE_IO
    processedTasks_[ts]++;
#endif
  }

  static int dumpVelocityWrapper(const int d, const int ts, const int ii, const int jj, const int kk);

  int dumpVelocity(const SWS::Directions d, const int ts, const int ii, const int jj, const int kk);

  static int dumpStressWrapper(const int sc, const int ts, const int ii, const int jj, const int kk);

  int dumpStress(const SWS::StressFieldComponents sc, const int ts, const int ii, const int jj, const int kk);

  int init();

  int finalize();

private:
  IOManager(const int nt, const int lnxx, const int lnyy, const int lnzz);

  ~IOManager();

#ifdef ENABLE_IO
  int dumpTile(adios2::Engine& writer,
               const SWS::SpatialBlockField<SWS::RealType>& tile3D,
               const std::string tileID);

  template<typename T>
  inline auto getUID(T arg)
  {
    return std::to_string(arg);
  }

  template<typename T, typename... Args>
  inline std::string getUID(T arg, Args... args)
  {
    return getUID(arg) + "-" + getUID(args...);
  }
#endif

  int nt_;
  int lnxx_;
  int lnyy_;
  int lnzz_;

  std::vector<std::atomic<int>> processedTasks_;

#ifdef ENABLE_IO
  std::mutex velocityGuard_;
  std::mutex stressGuard_;

  adios2::Engine velocityWriter_;
  adios2::Engine stressWriter_;

  std::unique_ptr<adios2::ADIOS> adios_;
  adios2::IO io_;
#endif
};
