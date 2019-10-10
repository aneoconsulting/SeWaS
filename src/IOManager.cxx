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

#include <mutex>

#include "IOManager.hxx"
#include "Mesh3DPartitioning.hxx"
#include "LinearSeismicWaveModel.hxx"
#include "LogManager.hxx"

IOManager::IOManager()
{
}

IOManager::~IOManager()
{
}

IOManager &IOManager::getInstance()
{
    static IOManager instance;
    return instance;
}

void IOManager::releaseInstance()
{
}

int IOManager::dumpVelocityWrapper(const int d,
                                   const int ts,
                                   const int ii, const int jj, const int kk)
{
    auto pMesh = Mesh3DPartitioning::getInstance();

    const auto lii = pMesh->lii(ii);
    const auto ljj = pMesh->ljj(jj);
    const auto lkk = pMesh->lkk(kk);

    return getInstance().dumpVelocity((const SWS::Directions)d, ts, lii, ljj, lkk);
}

int IOManager::dumpVelocity(const SWS::Directions d, const int ts, const int ii, const int jj, const int kk)
{
    int status = 0;

#ifdef ENABLE_IO
    LOG(SWS::LOG_DEBUG, "[start] Dumping v({})({},{},{}) at ts={}", d, ii, jj, kk, ts);

    const SWS::SpatialBlockField<SWS::RealType> &tile3D = LinearSeismicWaveModel::getInstance()->v(d)(ii,jj,kk);

    auto tileID = "Velocity-" + getUID(d, ts, ii, jj, kk);

    {
        std::lock_guard<std::mutex> lock(velocityGuard_);

        status = dumpTile(vIO_, velocityWriter_, tile3D, tileID);
    }

    LOG(SWS::LOG_DEBUG, "[stop] Dumping v({})({},{},{}) at ts={}", d, ii, jj, kk, ts);
#endif

    return status;
}

int IOManager::dumpStressWrapper(const int sc,
                                 const int ts,
                                 const int ii, const int jj, const int kk)
{
    auto pMesh = Mesh3DPartitioning::getInstance();

    const auto lii = pMesh->lii(ii);
    const auto ljj = pMesh->ljj(jj);
    const auto lkk = pMesh->lkk(kk);

    return getInstance().dumpStress((const SWS::StressFieldComponents)sc, ts, lii, ljj, lkk);
}

int IOManager::dumpStress(const SWS::StressFieldComponents sc, const int ts, const int ii, const int jj, const int kk)
{
    int status = 0;

#ifdef ENABLE_IO
    LOG(SWS::LOG_DEBUG, "[start] Dumping sigma({})({},{},{}) at ts={}", sc, ii, jj, kk, ts);

    const SWS::SpatialBlockField<SWS::RealType> &tile3D = LinearSeismicWaveModel::getInstance()->sigma(sc)(ii,jj,kk);

    auto tileID = "Stress-" + getUID(sc, ts, ii, jj, kk);

    {
        std::lock_guard<std::mutex> lock(stressGuard_);
        
        status = dumpTile(sigmaIO_, stressWriter_, tile3D, tileID);
    }

    LOG(SWS::LOG_DEBUG, "[stop] Dumping sigma({})({},{},{}) at ts={}", sc, ii, jj, kk, ts);
#endif

    return status;
}

#ifdef ENABLE_IO
int IOManager::dumpTile(adios2::IO &io, adios2::Engine &writer, const SWS::SpatialBlockField<SWS::RealType> &tile3D, const std::string tileID)
{
    int status = 0;

    auto n = tile3D.n();

    auto world = ExecutionContext::world();
    auto rank = ExecutionContext::rank();

    auto lncells = Mesh3DPartitioning::getInstance()->lncells();

    auto buffer = io.DefineVariable<SWS::RealType>(tileID, {world * lncells}, {rank * lncells}, {n}, adios2::ConstantDims);

    /** Write variable for buffering */
    writer.Put<SWS::RealType>(buffer, tile3D.data().data());

    return status;
}
#endif

int IOManager::init()
{
    int status = 0;
#ifdef ENABLE_IO
#if SEWAS_DISTRIBUTED
    adios_ = std::make_unique<adios2::ADIOS>(MPI_COMM_WORLD, adios2::DebugON);
#else
    adios_ = std::make_unique<adios2::ADIOS>(adios2::DebugON);
#endif

    // std::unordered_map<short, std::string> p = {{0, "0"}};
    
    vIO_ = adios_->DeclareIO("VelocityIO");
    sigmaIO_ = adios_->DeclareIO("StressIO");

    for (auto io : {vIO_, sigmaIO_})
    {
        io.SetEngine("BP4");

        io.DefineAttribute<std::string>("adios2_schema/version_major", std::to_string(ADIOS2_VERSION_MAJOR));
        io.DefineAttribute<std::string>("adios2_schema/version_minor", std::to_string(ADIOS2_VERSION_MINOR));

        io.DefineAttribute<std::int16_t>("adios2_schema/mesh/ordering", SWS::Ordering);

        io.DefineAttribute<std::int64_t>("adios2_schema/mesh/nxx", Mesh3DPartitioning::getInstance()->nxx());
        io.DefineAttribute<std::int64_t>("adios2_schema/mesh/nyy", Mesh3DPartitioning::getInstance()->nyy());
        io.DefineAttribute<std::int64_t>("adios2_schema/mesh/nzz", Mesh3DPartitioning::getInstance()->nzz());
    }

    /** Engine derived class, spawned to start IO operations */
    velocityWriter_ = vIO_.Open("Velocity.bp", adios2::Mode::Write);
    stressWriter_ = sigmaIO_.Open("Stress.bp", adios2::Mode::Write);

    velocityWriter_.BeginStep(adios2::StepMode::Append, 0);
    stressWriter_.BeginStep(adios2::StepMode::Append, 0);
#endif
    return status;
}

int IOManager::finalize()
{
    int status = 0;

#ifdef ENABLE_IO
    stressWriter_.EndStep();
    velocityWriter_.EndStep();

    stressWriter_.Close();
    velocityWriter_.Close();
#endif

    return status;
}