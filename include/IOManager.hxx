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

#include "Types.hxx"
#include "ExecutionContext.hxx"
#include "SpatialBlockField.hxx"

#include <adios2.h>

class IOManager
{
public:
    static IOManager &getInstance();
    static void releaseInstance();

    static int dumpVelocityWrapper(const int d,
                                          const int ts,
                                          const int ii, const int jj, const int kk);

    int dumpVelocity(const SWS::Directions d,
                     const int ts,
                     const int ii, const int jj, const int kk);

    static int dumpStressWrapper(const int sc,
                                 const int ts,
                                 const int ii, const int jj, const int kk);

    int dumpStress(const SWS::StressFieldComponents sc,
                   const int ts,
                   const int ii, const int jj, const int kk);

    int init();

    int finalize();

private:
    IOManager();
    
    ~IOManager();

    int dumpTile(adios2::IO &io, adios2::Engine &writer, const SWS::SpatialBlockField<SWS::RealType> &tile3D, const std::string tileID);

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


    std::unique_ptr<adios2::ADIOS> adios_;

    std::mutex velocityGuard_;
    std::mutex stressGuard_;

    adios2::Engine velocityWriter_;
    adios2::Engine stressWriter_;

    adios2::IO vIO_;
    adios2::IO sigmaIO_;
};
