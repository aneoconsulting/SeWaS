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
#include "LinearSeismicWaveModel.hxx"
#include "LogManager.hxx"

IOManager::IOManager(const int nt, const int lnxx, const int lnyy, const int lnzz)
  : nt_(nt)
  , lnxx_(lnxx)
  , lnyy_(lnyy)
  , lnzz_(lnzz)
  , processedTasks_(nt)
{
  for (auto& item : processedTasks_) {
    item = -1;
  }
}

IOManager::~IOManager() {}

IOManager&
IOManager::getInstance(const int nt, const int lnxx, const int lnyy, const int lnzz)
{
  static IOManager instance(nt, lnxx, lnyy, lnzz);
  return instance;
}

void
IOManager::releaseInstance()
{}

int
IOManager::dumpVelocityWrapper(const int d, const int ts, const int ii, const int jj, const int kk)
{
  auto pMesh = Mesh3DPartitioning::getInstance();

  const auto lii = pMesh->lii(ii);
  const auto ljj = pMesh->ljj(jj);
  const auto lkk = pMesh->lkk(kk);

  return getInstance().dumpVelocity((const SWS::Directions)d, ts, lii, ljj, lkk);
}

int
IOManager::dumpVelocity(const SWS::Directions d, const int ts, const int ii, const int jj, const int kk)
{
  int status = 0;

#ifdef ENABLE_IO
  LOG(SWS::LOG_DEBUG, "[start] Dumping v({})({},{},{}) at ts={}", d, ii, jj, kk, ts);

  const auto& tile3D = LinearSeismicWaveModel::getInstance()->v(d)(ii, jj, kk);

  auto tileID = "Velocity-" + getUID(d, ii, jj, kk);

  {
    std::lock_guard<std::mutex> lock(velocityGuard_);

    status = dumpTile(velocityWriter_, tile3D, tileID);
  }

  LOG(SWS::LOG_DEBUG, "[stop] Dumping v({})({},{},{}) at ts={}", d, ii, jj, kk, ts);
#endif

  return status;
}

int
IOManager::dumpStressWrapper(const int sc, const int ts, const int ii, const int jj, const int kk)
{
  auto pMesh = Mesh3DPartitioning::getInstance();

  const auto lii = pMesh->lii(ii);
  const auto ljj = pMesh->ljj(jj);
  const auto lkk = pMesh->lkk(kk);

  return getInstance().dumpStress((const SWS::StressFieldComponents)sc, ts, lii, ljj, lkk);
}

int
IOManager::dumpStress(const SWS::StressFieldComponents sc,
                      const int ts,
                      const int ii,
                      const int jj,
                      const int kk)
{
  int status = 0;

#ifdef ENABLE_IO
  LOG(SWS::LOG_DEBUG, "[start] Dumping sigma({})({},{},{}) at ts={}", sc, ii, jj, kk, ts);

  const auto& tile3D = LinearSeismicWaveModel::getInstance()->sigma(sc)(ii, jj, kk);

  auto tileID = "Stress-" + getUID(sc, ii, jj, kk);

  {
    std::lock_guard<std::mutex> lock(stressGuard_);

    status = dumpTile(stressWriter_, tile3D, tileID);
  }

  LOG(SWS::LOG_DEBUG, "[stop] Dumping sigma({})({},{},{}) at ts={}", sc, ii, jj, kk, ts);
#endif

  return status;
}

#ifdef ENABLE_IO
int
IOManager::dumpTile(adios2::Engine& writer,
                    const SWS::SpatialBlockField<SWS::RealType>& tile3D,
                    const std::string tileID)
{
  int status = 0;

  auto buffer = io_.InquireVariable<SWS::RealType>(tileID);

  /** Write variable for buffering */
  writer.Put<SWS::RealType>(buffer, tile3D.data().data(), adios2::Mode::Deferred);

  return status;
}
#endif

int
IOManager::init()
{
  int status = 0;
#ifdef ENABLE_IO
#if SEWAS_DISTRIBUTED
  adios_ = std::make_unique<adios2::ADIOS>(MPI_COMM_WORLD, adios2::DebugON);
#else
  adios_ = std::make_unique<adios2::ADIOS>(adios2::DebugON);
#endif

  processedTasks_[0] = 3 * lnxx_ * lnyy_ * lnzz_;
  processedTasks_[1] = 6 * lnxx_ * lnyy_ * lnzz_;

  // std::unordered_map<short, std::string> p = {{0, "0"}};

  auto world = ExecutionContext::world();
  auto rank = ExecutionContext::rank();

  auto tileSize = Mesh3DPartitioning::getInstance()->tileSize();

  auto lncells = tileSize; // Mesh3DPartitioning::getInstance()->lncells();

  io_ = adios_->DeclareIO("BP4");

  io_.SetEngine("BP4");

  io_.DefineAttribute<std::string>("adios2_schema/version_major", std::to_string(ADIOS2_VERSION_MAJOR));
  io_.DefineAttribute<std::string>("adios2_schema/version_minor", std::to_string(ADIOS2_VERSION_MINOR));

  io_.DefineAttribute<std::int16_t>("adios2_schema/mesh/ordering", SWS::Ordering);

  io_.DefineAttribute<std::int64_t>("adios2_schema/mesh/nxx", Mesh3DPartitioning::getInstance()->pm().nxx());
  io_.DefineAttribute<std::int64_t>("adios2_schema/mesh/nyy", Mesh3DPartitioning::getInstance()->pm().nyy());
  io_.DefineAttribute<std::int64_t>("adios2_schema/mesh/nzz", Mesh3DPartitioning::getInstance()->pm().nzz());

  // io_.SetParameter("StatsLevel", "0");
  // io_.SetParameter("MaxBufferSize", "2Gb");
  // io_.SetParameter("Threads", "2");

  for (int ii = 0; ii < lnxx_; ii++) {
    for (int jj = 0; jj < lnyy_; jj++) {
      for (int kk = 0; kk < lnzz_; kk++) {
        for (auto d : { SWS::X, SWS::Y, SWS::Z }) {
          auto tileID = "Velocity-" + getUID(d, ii, jj, kk);
          io_.DefineVariable<SWS::RealType>(tileID,
                                            { (unsigned int)world * lncells },
                                            { (unsigned int)rank * lncells },
                                            { (unsigned int)tileSize },
                                            adios2::ConstantDims);
        }

        for (auto sc : { SWS::XX, SWS::YY, SWS::ZZ, SWS::XY, SWS::XZ, SWS::YZ }) {
          auto tileID = "Stress-" + getUID(sc, ii, jj, kk);
          io_.DefineVariable<SWS::RealType>(tileID,
                                            { (unsigned int)world * lncells },
                                            { (unsigned int)rank * lncells },
                                            { (unsigned int)tileSize },
                                            adios2::ConstantDims);
        }
      }
    }
  }

  /** Engine derived class, spawned to start IO operations */
  velocityWriter_ = io_.Open("Velocity.bp", adios2::Mode::Write);
  stressWriter_ = io_.Open("Stress.bp", adios2::Mode::Write);

  LOG(SWS::LOG_TRACE, "Opening first time step");
  velocityWriter_.BeginStep(adios2::StepMode::Append, 0);
  stressWriter_.BeginStep(adios2::StepMode::Append, 0);
#endif
  return status;
}

int
IOManager::finalize()
{
  int status = 0;

#ifdef ENABLE_IO
  LOG(SWS::LOG_TRACE, "Closing last time step");
  stressWriter_.EndStep();
  velocityWriter_.EndStep();

  stressWriter_.Close();
  velocityWriter_.Close();
#endif

  return status;
}