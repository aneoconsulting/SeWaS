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

#include <iostream>
#include <algorithm>
#include <array>

#include "Config.hxx"
#include "LogManager.hxx"

class CentralFDOperator
{
  // This class implements the 4th order central finite difference operator
public:
  CentralFDOperator(const SWS::RealType dx, const SWS::RealType dy, const SWS::RealType dz):dx_(dx), dy_(dy), dz_(dz)
  {
    /* ci coefficients are taken in "Simulating Seismic Wave Propagation in 3D Elastic Media
       Using Staggered-Grid Finite Differences" [Graves, 1996] */
    constexpr SWS::RealType c1=-1./24;
    constexpr SWS::RealType c2=+1./24;
    constexpr SWS::RealType c3=+9./8;
    constexpr SWS::RealType c4=-9./8;

    c1_={c1/dx, c1/dy, c1/dz};
    c2_={c2/dx, c2/dy, c2/dz};
    c3_={c3/dx, c3/dy, c3/dz};
    c4_={c4/dx, c4/dy, c4/dz};
  }

  ~CentralFDOperator()
  {
  }


private:
  /*
    The variables names (c1, c2, c3, c4) and (d1, d2, d3, d4) are the same used in the paper:
    "3D Elastic Finite-Difference Modeling of Seismic Motion Using Staggered Grids with Nonuniform Spacing"
  */

  std::array<SWS::RealType, SWS::Directions::DIM> c1_;
  std::array<SWS::RealType, SWS::Directions::DIM> c2_;
  std::array<SWS::RealType, SWS::Directions::DIM> c3_;
  std::array<SWS::RealType, SWS::Directions::DIM> c4_;

  static constexpr short d1_=1;
  static constexpr short d2_=2;
  static constexpr short d3_=0;
  static constexpr short d4_=1;


  inline auto _c1(const SWS::Directions & d) const { return c1_[d]; }
  inline auto _c2(const SWS::Directions & d) const { return c2_[d]; }
  inline auto _c3(const SWS::Directions & d) const { return c3_[d]; }
  inline auto _c4(const SWS::Directions & d) const { return c4_[d]; }

  static inline constexpr auto _d1(const SWS::Directions & d) { return d1_; }
  static inline constexpr auto _d2(const SWS::Directions & d) { return d2_; }
  static inline constexpr auto _d3(const SWS::Directions & d) { return d3_; }
  static inline constexpr auto _d4(const SWS::Directions & d) { return d4_; }


public:
  static inline constexpr auto hnx() { return hnx_; }
  static inline constexpr auto hny() { return hny_; }
  static inline constexpr auto hnz() { return hnz_; }

  template<SWS::Directions d>
  inline auto apply(const SWS::SpatialBlockField<SWS::RealType> & fijk, const int i, const int j, const short DK) const noexcept
  {
    const auto cz=fijk.dimension(SWS::Z);

    const auto kStart=fijk.kStart(), kEnd=fijk.kEnd();

#ifndef EIGEN_VECTORIZATION
    SWS::SpatialBlockField1D rij(1,cz);
    rij=0.0;
#endif

#ifdef BOOST_SIMD_VECTORIZATION
    const auto nkp=rij.size()/SWS::packetSize;
    using pack_t=boost::simd::pack<SWS::RealType>;
#endif

    const auto c1=_c1(d);
    const auto c2=_c2(d);
    const auto c3=_c3(d);
    const auto c4=_c4(d);

    switch(d){
    case SWS::X:{

      constexpr auto d1=_d1(d);
      constexpr auto d2=_d2(d);
      constexpr auto d3=_d3(d);
      constexpr auto d4=_d4(d);

#ifdef EIGEN_VECTORIZATION
      return c1*fijk(i+d1,j)+c2*fijk(i-d2,j)+c3*fijk(i+d3,j)+c4*fijk(i-d4,j);
#elif  BOOST_SIMD_VECTORIZATION
      for (int kp=0; kp<nkp; kp++){
        auto rijkp=SWS::SpatialBlockField<SWS::RealType>::getPack(rij, kp);

        const auto fi1jkp=SWS::SpatialBlockField<SWS::RealType>::getPack(fijk(i+d1,j), kp);
        const auto fi2jkp=SWS::SpatialBlockField<SWS::RealType>::getPack(fijk(i-d2,j), kp);
        const auto fi3jkp=SWS::SpatialBlockField<SWS::RealType>::getPack(fijk(i+d3,j), kp);
        const auto fi4jkp=SWS::SpatialBlockField<SWS::RealType>::getPack(fijk(i-d4,j), kp);

        rijkp = c1*fi1jkp+c2*fi2jkp+c3*fi3jkp+c4*fi4jkp;

        SWS::SpatialBlockField<SWS::RealType>::setPack(rijkp, rij, kp);
      }
#else
      for (int k=kStart; k<kEnd; k++){
        rij(k) = c1*fijk(i+d1,j,k)+c2*fijk(i-d2,j,k)+c3*fijk(i+d3,j,k)+c4*fijk(i-d4,j,k); // 4 (*) + 3 (+) = 7 flops
      }
#endif

      break;
    }
    case SWS::Y:{

      constexpr auto d1=_d1(d);
      constexpr auto d2=_d2(d);
      constexpr auto d3=_d3(d);
      constexpr auto d4=_d4(d);

#ifdef EIGEN_VECTORIZATION
      return c1*fijk(i,j+d1)+c2*fijk(i,j-d2)+c3*fijk(i,j+d3)+c4*fijk(i,j-d4);
#elif  BOOST_SIMD_VECTORIZATION
      for (int kp=0; kp<nkp; kp++){
        auto rijkp=SWS::SpatialBlockField<SWS::RealType>::getPack(rij, kp);

        const auto fij1kp=SWS::SpatialBlockField<SWS::RealType>::getPack(fijk(i,j+d1), kp);
        const auto fij2kp=SWS::SpatialBlockField<SWS::RealType>::getPack(fijk(i,j-d2), kp);
        const auto fij3kp=SWS::SpatialBlockField<SWS::RealType>::getPack(fijk(i,j+d3), kp);
        const auto fij4kp=SWS::SpatialBlockField<SWS::RealType>::getPack(fijk(i,j-d4), kp);

        rijkp = c1*fij1kp+c2*fij2kp+c3*fij3kp+c4*fij4kp;

        SWS::SpatialBlockField<SWS::RealType>::setPack(rijkp, rij, kp);
      }
#else
      for (int k=kStart; k<kEnd; k++){
        rij(k) = c1*fijk(i,j+d1,k)+c2*fijk(i,j-d2,k)+c3*fijk(i,j+d3,k)+c4*fijk(i,j-d4,k);
      }
#endif

      break;
    }
    case SWS::Z:{

      const auto d1=_d1(d)+DK;
      const auto d2=_d2(d)-DK;
      const auto d3=_d3(d)+DK;
      const auto d4=_d4(d)-DK;

      auto _cz=kEnd-kStart;

#ifdef EIGEN_VECTORIZATION
      return c1*fijk.get(i,j).segment(kStart+d1,_cz)+c2*fijk.get(i,j).segment(kStart-d2,_cz)+c3*fijk.get(i,j).segment(kStart+d3,_cz)+c4*fijk.get(i,j).segment(kStart-d4,_cz);
#else
      for (int k=kStart; k<kEnd; k++){
        rij(k) = c1*fijk(i,j,k+d1)+c2*fijk(i,j,k-d2)+c3*fijk(i,j,k+d3)+c4*fijk(i,j,k-d4);
      }
#endif
      break;
    }
    default:
      LOG(SWS::ERROR, "Unknown direction {} requested within CentralFDOperator::apply()", d);
      break;
    }

#ifndef EIGEN_VECTORIZATION
    return rij;
#endif
  }

  template<SWS::Directions d, short DI, short DJ, short DK>
  inline auto apply(const SWS::SpatialBlockField<SWS::RealType> & fijk) const noexcept
  {
    const auto cx=fijk.dimension(SWS::X);
    const auto cy=fijk.dimension(SWS::Y);
    const auto cz=fijk.dimension(SWS::Z);

    /* The operator is not evaluated within the halo */
    const auto iStart=fijk.iStart(), iEnd=fijk.iEnd();
    const auto jStart=fijk.jStart(), jEnd=fijk.jEnd();

    // TODO avoid the following temporary object (can we overwrite data for the time step ts-1?)
    SWS::SpatialBlockField<SWS::RealType> r(cx,cy,cz);

    for (int i=iStart; i<iEnd; i++){
      for (int j=jStart; j<jEnd; j++){
        r(i,j)=apply<d>(fijk, i+DI, j+DJ, DK);
      }
    }

    return r;
  }


private:

  /* Halo size */
  // TODO need to add padding for vectorization efficiency
  static constexpr auto hnx_=(std::max)({d1_, d2_, d3_, d4_});
  static constexpr auto hny_=(std::max)({d1_, d2_, d3_, d4_});
  static constexpr auto hnz_=(std::max)({d1_, d2_, d3_, d4_});

  /* Spatial grid spacing */
  SWS::RealType dx_;
  SWS::RealType dy_;
  SWS::RealType dz_;
};
