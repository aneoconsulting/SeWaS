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

#ifdef USE_VTK

#include <iostream>
#include <vector>
#include <mutex>
#include <thread>

#include <vtkAutoInit.h>
#ifdef OPENGL2_SUPPORTED
VTK_MODULE_INIT(vtkRenderingOpenGL2);
#else
VTK_MODULE_INIT(vtkRenderingOpenGL);
#endif
VTK_MODULE_INIT(vtkInteractionStyle);

#include <vtkCamera.h>
#include <vtkCubeSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkProperty.h>
#include <vtkActor.h>
#include <vtkActor2D.h>
#include <vtkFFMPEGWriter.h>
#include <vtkSmartPointer.h>
#include <vtkWindowToImageFilter.h>
#include <vtkRenderWindowInteractor.h>

#include "Config.hxx"
#include "Mesh3DPartitioning.hxx"
#include "LogManager.hxx"

class VisualizationManager{
public:
  static VisualizationManager * getInstance(const int nt=1,
					    const int nxx=1, const int nyy=1, const int nzz=1,
                                            const int nthreads=1);
  static void releaseInstance();

  static inline int displayVelocityWrapper(const int d,
					   const int ts,
					   const int ii, const int jj, const int kk,
					   const int tid){
    return pInstance_->displayVelocity((const SWS::Directions) d, ts, ii, jj, kk, tid);
  }

  int displayVelocity(const SWS::Directions & d,
		      const int & ts,
		      const int & ii, const int & jj, const int & kk,
		      const int & tid);

  static inline int displayStressWrapper(const int sc,
					 const int ts,
					 const int ii, const int jj, const int kk,
					 const int tid){
    return pInstance_->displayStress((const SWS::StressFieldComponents) sc, ts, ii, jj, kk, tid);
  }

  int displayStress(const SWS::StressFieldComponents & sc,
		    const int & ts,
		    const int & ii, const int & jj, const int & kk,
		    const int & tid);

  inline void render()
  {
    renderingWindow_->Render();
  }

private:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW // To force an object of this class to be allocated as aligned

  VisualizationManager(const int nt,
		       const int nxx, const int nyy, const int nzz,
                       const int nthreads);
  ~VisualizationManager();

  enum {VELOCITY, STRESS, CLUSTER, NB_FIELDS}; // Fields to be displayed on the screen
  enum {X, Y};
  enum {MIN, MAX};

#ifdef ENABLE_CLUSTER_RENDERING
  static constexpr SWS::RealType CLUSTER_PANEL_SIZE=0.1;
#else
  static constexpr SWS::RealType CLUSTER_PANEL_SIZE=0.0;
#endif

  inline void lock()   { guard_.lock(); }
  inline void unlock() { guard_.unlock(); }

  template<int FIELD>
  inline void setCamera(vtkRenderer * renderer)
  {
    constexpr short side=3;
    const SWS::RealType diag=sqrt(nxx_*nxx_+nyy_*nyy_+nzz_*nzz_);

    renderer->SetBackground(.1, .1, .1);

    switch(FIELD){
    case VELOCITY:
    case STRESS:
      renderer->GetActiveCamera()->SetPosition(diag*side+1, diag*side+1, diag*side+1);
      renderer->GetActiveCamera()->SetFocalPoint(0, 0, 0);
      renderer->GetActiveCamera()->Zoom(0.9);
      break;
    case CLUSTER:
      renderer->GetActiveCamera()->SetPosition(0.5, 0.5, 10);
      renderer->GetActiveCamera()->SetFocalPoint(0.5, 0.5, 0);
      renderer->GetActiveCamera()->Zoom(0.5);
      break;
    default:
      break;
    }
  }

  template<int FIELD>
  inline void setViewport(const int id)
  {
    vtkRenderer * renderer=nullptr;
    const auto & vp=viewportRanges_[FIELD][id];

    switch(FIELD){
    case VELOCITY:
      renderer=VRenderers_[id];
      break;
    case STRESS:
      renderer=SRenderers_[id];
      break;
    case CLUSTER:
      renderer=CRenderers_[id];
      break;
    default:
      break;
    }

    if (renderer){
      renderer->SetViewport(vp[X][MIN], vp[Y][MIN], vp[X][MAX], vp[Y][MAX]);
    }
  }

  inline int desactivateActor(vtkActor * actor, vtkRenderer * renderer)
  {
    actor->GetProperty()->SetOpacity(0.5);

    return 0;
  }

  inline int removeActor(vtkActor * actor, vtkRenderer * renderer)
  {
    actor->VisibilityOff();

    return 0;
  }

  template <int FIELD>
  inline int activateActor(vtkActor * actor, vtkRenderer * renderer,
                           const int & ts,
                           const int & tid)
  {
    hsv in;

    SWS::ColorizationStrategies cs=getColorizationStrategy(FIELD);

    switch(cs){
    case SWS::CORE:
      in.h=360.*SWS::RealType(tid)/nthreads_;
      in.s=1.0;
      in.v=1.0;
      break;
    case SWS::TIME_STEP:
      in.h=360.*SWS::RealType(ts)/nt_;
      in.s=1.0;
      in.v=1.0;
      break;
    default:
      break;
    }

    rgb out=hsv2rgb(in);

    actor->GetProperty()->SetOpacity(1.);
    actor->GetProperty()->SetColor(out.r, out.g, out.b);

    actor->VisibilityOn();

    return 0;
  }

  template<int FIELD>
  inline int initActor(vtkActor * actor, vtkRenderer * renderer,
                       const int & ii, const int & jj, const int & kk)
  {
    /* Get my coordinates according to the process-grid */
    int lp=Mesh3DPartitioning::getInstance()->lp(ii);
    int lq=Mesh3DPartitioning::getInstance()->lq(jj);
    int lr=Mesh3DPartitioning::getInstance()->lr(kk);

    /* Evaluate the shift along the 3 dimensions relative to my rank */
    int sx=3*lp;
    int sy=3*lq;
    int sz=3*lr;

    actor->SetMapper(cubeMapper_);
    actor->SetPosition(ii+sx, jj+sy, kk+sz);
    actor->GetProperty()->SetAmbient(0.5);
    actor->GetProperty()->SetOpacity(0.5);

    renderer->AddActor(actor);

    switch(FIELD){
    case VELOCITY:
    case STRESS:
      actor->VisibilityOff();
      break;
    case CLUSTER:
      actor->VisibilityOn();
      break;
    default:
      break;
    }

    return 0;
  }

  template<int FIELD>
  inline int initActors()
  {
    switch(FIELD){
    case VELOCITY:
      for (auto d=0; d<SWS::DIM; d++){
        VActors_(d).resize(nxx_,nyy_,nzz_);
        for (int ii=0; ii<nxx_; ii++){
          for (int jj=0; jj<nyy_; jj++){
            for (int kk=0; kk<nzz_; kk++){
              auto & actor=VActors_(d)(ii,jj,kk);
              auto & renderer=VRenderers_[d];
              actor=vtkActor::New();

              initActor<VELOCITY>(actor, renderer, ii, jj, kk);
            }
          }
        }
      }
      break;
    case STRESS:
      for (auto sc=0; sc<SWS::NB_STRESS_FIELD_COMPONENTS; sc++){
        SActors_(sc).resize(nxx_,nyy_,nzz_);
        for (int ii=0; ii<nxx_; ii++){
          for (int jj=0; jj<nyy_; jj++){
            for (int kk=0; kk<nzz_; kk++){
              auto & actor=SActors_(sc)(ii,jj,kk);
              auto & renderer=SRenderers_[sc];
              actor=vtkActor::New();

              initActor<STRESS>(actor, renderer, ii, jj, kk);
            }
          }
        }
      }
      break;
    case CLUSTER:
      for (auto n=0; n<SWS::NB_CLUSTER_NODES; n++){
        for (int tid=0; tid<SWS::NB_CORES; tid++){
          auto & actor=CActors_(n,tid);
          auto & renderer=CRenderers_[n];
          actor=vtkActor::New();

          // (tid) -> (ii,jj,kk)
          auto ii=tid/2;
          auto jj=tid%2;
          auto kk=0;

          initActor<CLUSTER>(actor, renderer, ii, jj, kk);
        }
      }
      break;
    default:
      break;
    }

    return 0;
  }

  inline int initVideoWriter()
  {
    if (!getDumpVideoFlag()){
      return 0;
    }

    const auto P=Mesh3DPartitioning::getInstance()->P();
    const auto Q=Mesh3DPartitioning::getInstance()->Q();
    const auto R=Mesh3DPartitioning::getInstance()->R();

    windowToImage_ = vtkWindowToImageFilter::New();
    windowToImage_->SetInput(renderingWindow_);

    writer_=vtkFFMPEGWriter::New();
    writer_->SetInputConnection(windowToImage_->GetOutputPort());
    writer_->SetFileName(("sewas." + std::to_string(P) + "-" + std::to_string(Q) + "-" + std::to_string(R) + ".avi").c_str());
    writer_->SetQuality(20);

    writer_->Start();

    return 0;
  }

  inline int writeOneFrame()
  {
    if (!getDumpVideoFlag()){
      return 0;
    }

    lock();

    windowToImage_->Modified();
    writer_->Write();

    unlock();

    return 0;
  }

  inline void setViewportRanges()
  {
    viewportRanges_.resize(NB_FIELDS);

#ifdef ENABLE_VELOCITY_RENDERING
    // (Velocity)(X|Y|Z)(x|y)(MIN|MAX)
    viewportRanges_[VELOCITY].resize(SWS::DIM);
    for (auto d=0; d<SWS::DIM; d++){
      viewportRanges_[VELOCITY][d].resize(2);
      for (int i=0; i<2; i++){ // x,y
	viewportRanges_[VELOCITY][d][i].resize(2);
	switch(i){
	case X:
	  viewportRanges_[VELOCITY][d][X][MIN]=((SWS::RealType) d  )*(1-CLUSTER_PANEL_SIZE)/SWS::DIM;
	  viewportRanges_[VELOCITY][d][X][MAX]=((SWS::RealType) d+1)*(1-CLUSTER_PANEL_SIZE)/SWS::DIM;
	  break;
	case Y:
	  viewportRanges_[VELOCITY][d][Y][MIN]=0;
	  viewportRanges_[VELOCITY][d][Y][MAX]=1./2;
	  break;
	default:
	  break;
	}
      }
    }
#endif

#if ENABLE_STRESS_RENDERING
    // (Stress)(XX|YY|ZZ|XY|XZ|YZ)(x|y)(MIN|MAX)
    viewportRanges_[STRESS].resize(SWS::NB_STRESS_FIELD_COMPONENTS);
    for (auto sc=0; sc<SWS::NB_STRESS_FIELD_COMPONENTS; sc++){
      viewportRanges_[STRESS][sc].resize(2);
      for (int i=0; i<2; i++){ // x,y
	viewportRanges_[STRESS][sc][i].resize(2);
	switch(i){
	case X:
	  viewportRanges_[STRESS][sc][X][MIN]=((SWS::RealType) sc  )*(1-CLUSTER_PANEL_SIZE)/SWS::NB_STRESS_FIELD_COMPONENTS;
	  viewportRanges_[STRESS][sc][X][MAX]=((SWS::RealType) sc+1)*(1-CLUSTER_PANEL_SIZE)/SWS::NB_STRESS_FIELD_COMPONENTS;
	  break;
	case Y:
	  viewportRanges_[STRESS][sc][Y][MIN]=1./2;
	  viewportRanges_[STRESS][sc][Y][MAX]=1.;
	  break;
	default:
	  break;
	}
      }
    }
#endif

#ifdef ENABLE_CLUSTER_RENDERING
    // Right panel for displaying the cluster Map
    viewportRanges_[CLUSTER].resize(SWS::NB_CLUSTER_NODES);
    for (auto n=0; n<SWS::NB_CLUSTER_NODES; n++){
      viewportRanges_[CLUSTER][n].resize(2);
      for (int i=0; i<2; i++){ // x,y
	viewportRanges_[CLUSTER][n][i].resize(2);
	switch(i){
	case X:
	  viewportRanges_[CLUSTER][n][X][MIN]=(1-CLUSTER_PANEL_SIZE)+(n%2  )*(CLUSTER_PANEL_SIZE)/2;
	  viewportRanges_[CLUSTER][n][X][MAX]=(1-CLUSTER_PANEL_SIZE)+(n%2+1)*(CLUSTER_PANEL_SIZE)/2;
	  break;
	case Y:
	  viewportRanges_[CLUSTER][n][Y][MIN]=( n   /2)/(SWS::NB_CLUSTER_NODES/2.);
	  viewportRanges_[CLUSTER][n][Y][MAX]=((n/2)+1)/(SWS::NB_CLUSTER_NODES/2.);
	  break;
	default:
	  break;
	}
      }
    }
#endif
  }

  inline int rank(const int & ii, const int & jj, const int & kk)
  {
    return Mesh3DPartitioning::getInstance()->rank_of(NULL, ii, jj, kk);
  }

  inline bool getDumpVideoFlag() const
  {
    const char * dv=getenv("DUMP_VIDEO");
    if (dv){
      if (strcmp("0", dv) == 0){
        return 0;
      }
      else if (strcmp("1", dv) == 0){
        return 1;
      }
      else{
        return 0;
      }
    }
    else{
      return 0;
    }
  }

  inline SWS::ColorizationStrategies getColorizationStrategy(const int f) const
  {
    if (CLUSTER == f){
      /* Force the cores to have different colors */
      return SWS::CORE;
    }

    const char * cs=getenv("COLORIZATION_STRATEGY");
    if (cs){
      if (strcmp("CORE", cs) == 0){
        return SWS::CORE;
      }
      else if (strcmp("TIME_STEP", cs) == 0){
        return SWS::TIME_STEP;
      }
      else{
        LOG(SWS::LOG_WARN, "Unknown colorization strategy {}. Accepted values are CORE and TIME_STEP. Fallback to CORE", cs);
        return SWS::CORE;
      }
    }
    else{
      LOG(SWS::LOG_WARN, "Colorization strategy is not set. Using default strategy: CORE");
      return SWS::CORE;
    }
  }

  /*
    The code that converts hsv to rgb is taken from https://stackoverflow.com/a/6930407
  */
  typedef struct
  {
    double r;       // percent
    double g;       // percent
    double b;       // percent
  } rgb;

  typedef struct
  {
    double h;       // angle in degrees
    double s;       // percent
    double v;       // percent
  } hsv;

  inline rgb hsv2rgb(hsv in)
  {
    double      hh, p, q, t, ff;
    long        i;
    rgb         out;

    if(in.s <= 0.0){       // < is bogus, just shuts up warnings
      out.r = in.v;
      out.g = in.v;
      out.b = in.v;
      return out;
    }

    hh = in.h;
    if(hh >= 360.0) hh = 0.0;
    hh /= 60.0;
    i = (long)hh;
    ff = hh - i;
    p = in.v * (1.0 - in.s);
    q = in.v * (1.0 - (in.s * ff));
    t = in.v * (1.0 - (in.s * (1.0 - ff)));

    switch(i){
    case 0:
      out.r = in.v;
      out.g = t;
      out.b = p;
      break;

    case 1:
      out.r = q;
      out.g = in.v;
      out.b = p;
      break;

    case 2:
      out.r = p;
      out.g = in.v;
      out.b = t;
      break;

    case 3:
      out.r = p;
      out.g = q;
      out.b = in.v;
      break;

    case 4:
      out.r = t;
      out.g = p;
      out.b = in.v;
      break;
    case 5:

    default:
      out.r = in.v;
      out.g = p;
      out.b = q;
      break;
    }

    return out;
  }

  int removeVelocityTask(const SWS::Directions & d,
                         const int & ii, const int & jj, const int & kk);
  int desactivateVelocityTask(const SWS::Directions & d,
                              const int & ii, const int & jj, const int & kk);
  int activateVelocityTask(const SWS::Directions & d,
                           const int & ts,
                           const int & ii, const int & jj, const int & kk,
                           const int & tid);

  int removeStressTask(const SWS::StressFieldComponents & sc,
                       const int & ii, const int & jj, const int & kk);
  int desactivateStressTask(const SWS::StressFieldComponents & sc,
                            const int & ii, const int & jj, const int & kk);
  int activateStressTask(const SWS::StressFieldComponents & sc,
                         const int & ts,
                         const int & ii, const int & jj, const int & kk,
                         const int & tid);

  int removeCore(const SWS::ClusterNodes & n,
                 const int & tid);
  int addCore(const SWS::ClusterNodes & n,
              const int & tid,
              const int & ts);
  int displayCore(const SWS::ClusterNodes & n,
                  const int & tid,
                  const int & ts);



  static VisualizationManager * pInstance_;

  std::mutex guard_;

  const int nt_;
  const int nxx_;
  const int nyy_;
  const int nzz_;

  const int nthreads_;

  short colorizationStrategy_;

  vtkCubeSource			* cube_;
  vtkPolyDataMapper		* cubeMapper_;

  std::vector<vtkRenderer *>      VRenderers_; // Velocity
  std::vector<vtkRenderer *>      SRenderers_; // Stress
  std::vector<vtkRenderer *>      CRenderers_; // Cluster nodes

  vtkRenderWindow		* renderingWindow_;

  vtkRenderWindowInteractor	* renderingWindowInteractor_;

  // (Velocity|Stress|Cluster)(id)(x|y)(MIN|MAX)
  std::vector<std::vector<std::vector<std::vector<SWS::RealType>>>> viewportRanges_;

  vtkFFMPEGWriter		* writer_;
  vtkWindowToImageFilter	* windowToImage_;

  SWS::VTKVelocityActors          VActors_;
  SWS::VTKStressActors            SActors_;
  SWS::VTKClusterActors           CActors_;
};

#endif
