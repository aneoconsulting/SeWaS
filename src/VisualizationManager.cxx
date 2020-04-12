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

#ifdef USE_VTK

#include "VisualizationManager.hxx"
#include "Mesh3DPartitioning.hxx"

#include "Config.hxx"

#include <vtkFFMPEGWriter.h>
#include <vtkSmartPointer.h>

VisualizationManager* VisualizationManager::pInstance_ = nullptr;

// For compatibility with new VTK generic data arrays
#ifdef vtkGenericDataArray_h
#define InsertNextTupleValue InsertNextTypedTuple
#endif

VisualizationManager*
VisualizationManager::getInstance(const int nt,
                                  const int nxx,
                                  const int nyy,
                                  const int nzz,
                                  const int nthreads)
{
  if (nullptr == pInstance_) {
    pInstance_ = new VisualizationManager(nt, nxx, nyy, nzz, nthreads);
    return pInstance_;
  } else {
    return pInstance_;
  }
}

void
VisualizationManager::releaseInstance()
{
  if (pInstance_) {
    delete pInstance_;
    pInstance_ = nullptr;
  }
}

int
VisualizationManager::removeVelocityTask(const SWS::Directions& d,
                                         const int& ii,
                                         const int& jj,
                                         const int& kk)
{
  /* In order to remove the task previously displayed at
     the same spatial location (during the previous time-step) */

  auto& actor = VActors_(d)(ii, jj, kk);
  auto& renderer = VRenderers_[d];

  removeActor(actor, renderer);

  return 0;
}

int
VisualizationManager::desactivateVelocityTask(const SWS::Directions& d,
                                              const int& ii,
                                              const int& jj,
                                              const int& kk)
{
  auto& actor = VActors_(d)(ii, jj, kk);
  auto& renderer = VRenderers_[d];

  desactivateActor(actor, renderer);

  return 0;
}

int
VisualizationManager::activateVelocityTask(const SWS::Directions& d,
                                           const int& ts,
                                           const int& ii,
                                           const int& jj,
                                           const int& kk,
                                           const int& tid)
{
  auto& actor = VActors_(d)(ii, jj, kk);
  auto& renderer = VRenderers_[d];

  activateActor<VELOCITY>(actor, renderer, ts, tid);

  return 0;
}

int
VisualizationManager::displayVelocity(const SWS::Directions& d,
                                      const int& ts,
                                      const int& ii,
                                      const int& jj,
                                      const int& kk,
                                      const int& tid)
{
  int status = 0;

#ifdef ENABLE_CLUSTER_RENDERING
  auto n = rank(ii, jj, kk);
  displayCore((SWS::ClusterNodes)n, tid, ts);
#endif

#ifdef ENABLE_VELOCITY_RENDERING
  status = removeVelocityTask(d, ii, jj, kk);
  status = writeOneFrame();

  status = activateVelocityTask(d, ts, ii, jj, kk, tid);
  status = writeOneFrame();

  status = desactivateVelocityTask(d, ii, jj, kk);
  status = writeOneFrame();
#endif

  return status;
}

int
VisualizationManager::removeStressTask(const SWS::StressFieldComponents& sc,
                                       const int& ii,
                                       const int& jj,
                                       const int& kk)
{
  /* In order to remove the task previously displayed at
     the same spatial location (during the previous time-step) */

#ifdef ENABLE_STRESS_RENDERING
  auto& actor = SActors_(sc)(ii, jj, kk);
  auto& renderer = SRenderers_[sc];

  removeActor(actor, renderer);
#endif

  return 0;
}

int
VisualizationManager::desactivateStressTask(const SWS::StressFieldComponents& sc,
                                            const int& ii,
                                            const int& jj,
                                            const int& kk)
{
#ifdef ENABLE_STRESS_RENDERING
  auto& actor = SActors_(sc)(ii, jj, kk);
  auto& renderer = SRenderers_[sc];

  desactivateActor(actor, renderer);
#endif

  return 0;
}

int
VisualizationManager::activateStressTask(const SWS::StressFieldComponents& sc,
                                         const int& ts,
                                         const int& ii,
                                         const int& jj,
                                         const int& kk,
                                         const int& tid)
{

#ifdef ENABLE_STRESS_RENDERING
  auto& actor = SActors_(sc)(ii, jj, kk);
  auto& renderer = SRenderers_[sc];

  activateActor<STRESS>(actor, renderer, ts, tid);
#endif

  return 0;
}

int
VisualizationManager::displayStress(const SWS::StressFieldComponents& sc,
                                    const int& ts,
                                    const int& ii,
                                    const int& jj,
                                    const int& kk,
                                    const int& tid)
{
  int status = 0;

#ifdef ENABLE_CLUSTER_RENDERING
  auto n = rank(ii, jj, kk);
  displayCore((SWS::ClusterNodes)n, tid, ts);
#endif

#ifdef ENABLE_STRESS_RENDERING
  status = removeStressTask(sc, ii, jj, kk);
  status = writeOneFrame();

  status = activateStressTask(sc, ts, ii, jj, kk, tid);
  status = writeOneFrame();

  status = desactivateStressTask(sc, ii, jj, kk);
  status = writeOneFrame();
#endif

  return status;
}

int
VisualizationManager::removeCore(const SWS::ClusterNodes& n, const int& tid)
{
  /* In order to remove the task previously displayed at
     the same spatial location (during the previous time-step) */

  auto& actor = CActors_(n, tid);
  auto& renderer = CRenderers_[n];

  removeActor(actor, renderer);

  return 0;
}

int
VisualizationManager::addCore(const SWS::ClusterNodes& n, const int& tid, const int& ts)
{
  auto& actor = CActors_(n, tid);
  auto& renderer = CRenderers_[n];

  activateActor<CLUSTER>(actor, renderer, ts, tid);

  return 0;
}

int
VisualizationManager::displayCore(const SWS::ClusterNodes& n, const int& tid, const int& ts)
{
  /* This method is called when a task has been executed
     and is intended to illustrate the activity of a single CPU core */

  int status = 0;

  status = removeCore(n, tid);
  status = writeOneFrame();

  status = addCore(n, tid, ts);
  status = writeOneFrame();

  status = removeCore(n, tid);
  status = writeOneFrame();

  return status;
}

VisualizationManager::VisualizationManager(const int nt,
                                           const int nxx,
                                           const int nyy,
                                           const int nzz,
                                           const int nthreads)
  : nt_(nt)
  , nxx_(nxx)
  , nyy_(nyy)
  , nzz_(nzz)
  , nthreads_(nthreads)
{
  cube_ = vtkCubeSource::New();
  cube_->SetBounds(0, 1, 0, 1, 0, 1);

  cubeMapper_ = vtkPolyDataMapper::New();
  cubeMapper_->SetInputConnection(cube_->GetOutputPort());
  cubeMapper_->ImmediateModeRenderingOn();

  /* Build the rendering window */
  renderingWindow_ = vtkRenderWindow::New();
  // renderingWindow_->SetFullScreen(true);
  renderingWindow_->SetSize(1920, 1080);
  // renderingWindow_->SetSize(800, 600);
  renderingWindow_->SetWindowName("SeWaS");

  VRenderers_.resize(SWS::DIM);
#ifdef ENABLE_STRESS_RENDERING
  SRenderers_.resize(SWS::NB_STRESS_FIELD_COMPONENTS);
#endif
#ifdef ENABLE_CLUSTER_RENDERING
  CRenderers_.resize(SWS::NB_CLUSTER_NODES);
#endif

  /* Define viewport ranges */
  setViewportRanges();

#ifdef ENABLE_VELOCITY_RENDERING
  /* Add Velocity renderers */
  for (auto d = 0; d < SWS::DIM; d++) {
    /* Build the renderer */
    auto& renderer = VRenderers_[d];
    renderer = vtkRenderer::New();
    setCamera<VELOCITY>(renderer);
    setViewport<VELOCITY>(d);

    renderingWindow_->AddRenderer(renderer);
  }
#endif

#ifdef ENABLE_STRESS_RENDERING
  /* Add Stress renderers */
  for (auto sc = 0; sc < SWS::NB_STRESS_FIELD_COMPONENTS; sc++) {
    /* Build the renderer */
    auto& renderer = SRenderers_[sc];
    renderer = vtkRenderer::New();
    setCamera<STRESS>(renderer);
    setViewport<STRESS>(sc);

    renderingWindow_->AddRenderer(renderer);
  }
#endif

#ifdef ENABLE_CLUSTER_RENDERING
  /* Add Cluster Node renderers */
  for (auto n = 0; n < SWS::NB_CLUSTER_NODES; n++) {
    /* Build the renderer */
    auto& renderer = CRenderers_[n];
    renderer = vtkRenderer::New();
    setCamera<CLUSTER>(renderer);
    setViewport<CLUSTER>(n);

    renderingWindow_->AddRenderer(renderer);
  }
#endif

  /* At this stage, we put all actors on the screen and switch off their visibility */
#ifdef ENABLE_VELOCITY_RENDERING
  initActors<VELOCITY>();
#endif
#ifdef ENABLE_STRESS_RENDERING
  initActors<STRESS>();
#endif
#ifdef ENABLE_CLUSTER_RENDERING
  initActors<CLUSTER>();
#endif

  /* Setup the video writer */
  initVideoWriter();
}

VisualizationManager::~VisualizationManager()
{
  cube_->Delete();
  cubeMapper_->Delete();

  renderingWindow_->Delete();

#ifdef ENABLE_VELOCITY_RENDERING
  for (auto d = 0; d < SWS::DIM; d++) {
    VRenderers_[d]->Delete();
    for (int ii = 0; ii < nxx_; ii++) {
      for (int jj = 0; jj < nyy_; jj++) {
        for (int kk = 0; kk < nzz_; kk++) {
          auto& actor = VActors_(d)(ii, jj, kk);
          if (actor) {
            actor->Delete();
          }
        }
      }
    }
  }
#endif

#ifdef ENABLE_STRESS_RENDERING
  for (auto sc = 0; sc < SWS::NB_STRESS_FIELD_COMPONENTS; sc++) {
    SRenderers_[sc]->Delete();
    for (int ii = 0; ii < nxx_; ii++) {
      for (int jj = 0; jj < nyy_; jj++) {
        for (int kk = 0; kk < nzz_; kk++) {
          auto& actor = SActors_(sc)(ii, jj, kk);
          if (actor) {
            actor->Delete();
          }
        }
      }
    }
  }
#endif

#ifdef ENABLE_CLUSTER_RENDERING
  for (auto n = 0; n < SWS::NB_CLUSTER_NODES; n++) {
    CRenderers_[n]->Delete();
    for (int tid = 0; tid < SWS::NB_CORES; tid++) {
      auto& actor = CActors_(n, tid);
      if (actor) {
        actor->Delete();
      }
    }
  }
#endif

  if (getDumpVideoFlag()) {
    writer_->End();
    writer_->Delete();
    windowToImage_->Delete();
  }
}

#endif
