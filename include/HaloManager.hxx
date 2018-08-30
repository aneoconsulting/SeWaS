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

#include "Config.hxx"

class HaloManager{
public:
  static HaloManager * getInstance(SWS::StressField & sigma, SWS::Velocity & v,
				   const int & hnx, const int & hny, const int & hnz);
  static HaloManager * getInstance();

  static void releaseInstance();



  /* The following wrappers correspond to entry points called by the PaRSEC runtime when executing tasks' BODY */
  static int extractVelocityHaloWrapper(const int l,
					const int d,
					const int ts,
					const int ii, const int jj, const int kk,
					void * vH);

  static int updateVelocityWrapper(const int l,
				   const int d,
				   const int ts,
				   const int ii, const int jj, const int kk,
				   void * vH);

  static int extractStressHaloWrapper(const int l,
				      const int sc,
				      const int ts,
				      const int ii, const int jj, const int kk,
				      void * sigmaH);

  static int updateStressWrapper(const int l,
				 const int sc,
				 const int ts,
				 const int ii, const int jj, const int kk,
				 void * sigmaH);


  inline void setVelocityHalo(const SWS::Locations l,
			      const SWS::Directions d,
			      const int ii, const int jj, const int kk,
			      void * vH) noexcept
  {
    pInstance_->vH()(d)(ii,jj,kk)(l)=((SWS::RealType *) vH);
  }

  int extractVelocityHalo(const SWS::Locations l,
			  const SWS::Directions d,
			  const int ts,
			  const int ii, const int jj, const int kk) noexcept;

  int updateVelocity(const SWS::Locations l,
		     const SWS::Directions d,
		     const int ts,
		     const int ii, const int jj, const int kk) noexcept;


  inline void setStressHalo(const SWS::Locations l,
			    const SWS::StressFieldComponents sc,
			    const int ii, const int jj, const int kk,
			    void * sigmaH) noexcept
  {
    pInstance_->sigmaH()(sc)(ii,jj,kk)(l)=((SWS::RealType *) sigmaH);
  }

  int extractStressHalo(const SWS::Locations l,
			const SWS::StressFieldComponents sc,
			const int ts,
			const int ii, const int jj, const int kk) noexcept;

  int updateStress(const SWS::Locations l,
		   const SWS::StressFieldComponents sc,
		   const int ts,
		   const int ii, const int jj, const int kk) noexcept;


  SWS::StressField & sigma(){
    return sigma_;
  }

  SWS::Velocity & v(){
    return v_;
  }

  SWS::StressFieldHalo & sigmaH(){
    return sigmaH_;
  }

  SWS::VelocityHalo & vH(){
    return vH_;
  }

  inline const int & hnx() const { return hnx_; }
  inline const int & hny() const { return hny_; }
  inline const int & hnz() const { return hnz_; }

  size_t getHaloSize(const int ii, const int jj, const int kk) const;

  size_t getHaloSize(const SWS::Locations l,
                     const int ii, const int jj, const int kk) const;

private:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW // To force an object of this class to be allocated as aligned

  HaloManager(SWS::StressField & sigma, SWS::Velocity & v,
	      const int & hnx, const int & hny, const int & hnz);
  HaloManager();

  ~HaloManager();

  void setExtractOffsets(const SWS::Locations l,
                         const int ii, const int jj, const int kk,
                         int & iStart, int & iEnd,
                         int & jStart, int & jEnd,
                         int & kStart, int & kEnd,
                         int & iShift, int & jShift, int & kShift) noexcept;

  void setUpdateOffsets(const SWS::Locations l,
                        const int ii, const int jj, const int kk,
                        int & iStart, int & iEnd,
                        int & jStart, int & jEnd,
                        int & kStart, int & kEnd,
                        int & iShift, int & jShift, int & kShift) noexcept;

  static HaloManager * pInstance_;

  SWS::StressField & sigma_;
  SWS::Velocity & v_;

  /* These members are allocated by the PaRSEC runtime and we keep here only the pointers to reference them */
  SWS::StressFieldHalo sigmaH_;
  SWS::VelocityHalo vH_;

  const int hnx_;
  const int hny_;
  const int hnz_;
};
