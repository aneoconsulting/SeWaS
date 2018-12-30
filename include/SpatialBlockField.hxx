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
#include <fstream>
#include <tuple>

#ifdef USE_MATPLOTLIB
#include <matplotlibcpp.h>
#endif

#include <Eigen/Core>

#ifdef BOOST_SIMD_VECTORIZATION
#include <boost/simd/pack.hpp>
#include <boost/simd/function/load.hpp>
#include <boost/simd/function/store.hpp>
#include <boost/simd/function/aligned_load.hpp>
#include <boost/simd/function/aligned_store.hpp>
#endif

#include "Constants.hxx"
#include "LogManager.hxx"

/*
  This class represents a 3D block of spatial cells. It internally uses Eigen::Array to store its data.
*/

namespace SWS
{
  template<typename RealType>
  class SpatialBlockField
  {

  public:
    SpatialBlockField():nx_(1),
			ny_(1),
			nz_(1)
    {
      // Add padding along z-axis for vectorization efficiency
      px_=0;
      py_=0;
      pz_=std::ceil(RealType(nz_)/packetSize_)*packetSize_-nz_;

      nz_+=pz_;

      n_=nx_*ny_*nz_;

      data_.resize(nx_,ny_);
      for (int i=0; i<nx_; i++){
        for (int j=0; j<ny_; j++){
          auto & dij=data_(i,j);
          dij.resize(1,nz_);
        }
      }
    }

    SpatialBlockField(const int nx, const int ny, const int nz):nx_(nx),
                                                                ny_(ny),
                                                                nz_(nz)
    {
      // Add padding along z-axis for vectorization efficiency
      px_=0;
      py_=0;
      pz_=std::ceil(RealType(nz_)/packetSize_)*packetSize_-nz_;

      nz_+=pz_;

      n_=nx_*ny_*nz_;

      data_.resize(nx_,ny_);
      for (int i=0; i<nx_; i++){
        for (int j=0; j<ny_; j++){
          auto & dij=data_(i,j);
          dij.resize(1,nz_);
        }
      }
    }

    SpatialBlockField(const SpatialBlockField & o)
    {
      if (this != &o){
	nx_=o.nx_;
	ny_=o.ny_;
	nz_=o.nz_;
	n_=o.n_;

        // hnx_=o.hnx_;
        // hny_=o.hny_;
        // hnz_=o.hnz_;

        px_=o.px_;
        py_=o.py_;
        pz_=o.pz_;

        // FIXME we just need to copy the data without modifying source locations
        // hasSource_=o.hasSource_;

        // is_=o.is_;
        // js_=o.js_;
        // ks_=o.ks_;

        data_=o.data_;
      }
    }


    ~SpatialBlockField()
    {
    }

    inline auto dimension(const short d) const
    {
      size_t dim;

      switch(d){
      case X:
	dim=nx_;
	break;
      case Y:
	dim=ny_;
	break;
      case Z:
	dim=nz_;
	break;
      default:
        LogManager::getInstance()->log<SWS::CRITICAL>("Unknown dimension {} requested within SpatialBlockField::dimension()", d);
	exit(SWS::UNKNOWN_SPATIAL_DIRECTION);
      }

      return dim;
    }

    inline auto dimensions() const
    {
      return std::make_tuple<const int, const int, const int>((const int) nx_, (const int) ny_, (const int) nz_);
    }

    inline void setHaloSize(const int hnx, const int hny, const int hnz)
    {
      hnx_=hnx;
      hny_=hny;
      hnz_=hnz;
    }

    inline void resize(const int nx, const int ny, const int nz)
    {
      nx_=nx;
      ny_=ny;
      nz_=nz;

      // Add padding along z-axis for vectorization efficiency
      px_=0;
      py_=0;
      pz_=std::ceil(RealType(nz_)/packetSize_)*packetSize_-nz_;

      nz_+=pz_;

      n_=nx_*ny_*nz_;

      data_.resize(nx_,ny_);
      for (int i=0; i<nx_; i++){
        for (int j=0; j<ny_; j++){
          auto & dij=data_(i,j);
          dij.resize(1,nz_);
        }
      }
    }

    inline void addSource(const int is, const int js, const int ks)
    {
      hasSource_=true;

      is_.push_back(is);
      js_.push_back(js);
      ks_.push_back(ks);
    }


    inline const auto & nx() const { return nx_; }
    inline const auto & ny() const { return ny_; }
    inline const auto & nz() const { return nz_; }

    inline const auto & hnx() const { return hnx_; }
    inline const auto & hny() const { return hny_; }
    inline const auto & hnz() const { return hnz_; }

    inline const auto & px() const { return px_; }
    inline const auto & py() const { return py_; }
    inline const auto & pz() const { return pz_; }

    inline const auto & n() const { return n_; }

    inline auto & data() { return data_; }
    inline const auto & data() const { return data_; }

    inline const auto iStart() const { return hnx_; }
    inline const auto iEnd() const { return nx_-hnx_-px_; }

    inline const auto jStart() const { return hny_; }
    inline const auto jEnd() const { return ny_-hny_-py_; }

    inline const auto kStart() const { return hnz_; }
    inline const auto kEnd() const { return nz_-hnz_-pz_; }

    inline auto index(const int i, const int j, const int k) const { return k*nx_*ny_+j*nx_+i; }

    inline const auto & hasSource() const { return hasSource_; }

    inline const auto & is() const { return is_; }
    inline const auto & js() const { return js_; }
    inline const auto & ks() const { return ks_; }

    inline auto & operator=(const SpatialBlockField && o)
    {
      if (this != &o){
	nx_=o.nx_;
	ny_=o.ny_;
	nz_=o.nz_;

	n_=o.n_;

        hnx_=o.hnx_;
        hny_=o.hny_;
        hnz_=o.hnz_;

        px_=o.px_;
        py_=o.py_;
        pz_=o.pz_;

        // FIXME we just need to copy the data without modifying source locations
        // hasSource_=o.hasSource_;

        // is_=o.is_;
        // js_=o.js_;
        // ks_=o.ks_;

        data_=o.data_;
      }
      return *this;
    }

    inline auto & operator=(const RealType & v)
    {
      for (int i=0; i<nx_; i++){
        for (int j=0; j<ny_; j++){
          data_(i,j)=v;
        }
      }
      return *this;
    }


    inline const auto & get(const int i, const int j) const
    {
      return data_(i,j);
    }

    inline const auto & get(const int i, const int j, const int k) const
    {
      return data_(i,j)(k);
    }


    inline auto operator()(const int i, const int j)
    {
      return Eigen::VectorBlock<SpatialBlockField1D>(data_(i,j),kStart(),kEnd()-kStart());
    }

    inline const auto operator()(const int i, const int j) const
    {
      return Eigen::VectorBlock<const SpatialBlockField1D>(data_(i,j),kStart(),kEnd()-kStart());
    }


    inline auto & operator()(const int i, const int j, const int k)
    {
      return data_(i,j)(k);
    }

    inline const auto & operator()(const int i, const int j, const int k) const
    {
      return data_(i,j)(k);
    }

    inline auto operator+(const SpatialBlockField & o) const
    {
      SpatialBlockField r(nx_,ny_,nz_);
      for (int i=0; i<nx_; i++){
        for (int j=0; j<ny_; j++){
          r.data_(i,j)=data_(i,j)+o.data()(i,j);
        }
      }
      return r;
    }

    inline auto & operator+=(const SpatialBlockField & o)
    {
      for (int i=0; i<nx_; i++){
        for (int j=0; j<ny_; j++){
          data_(i,j)+=o.data()(i,j);
        }
      }
      return *this;
    }

    inline auto operator-(const SpatialBlockField & o) const
    {
      SpatialBlockField r(nx_,ny_,nz_);
      for (int i=0; i<nx_; i++){
        for (int j=0; j<ny_; j++){
          r.data_(i,j)=data_(i,j)-o.data()(i,j);
        }
      }
      return r;
    }

    inline auto & operator-=(const SpatialBlockField & o)
    {
      for (int i=0; i<nx_; i++){
        for (int j=0; j<ny_; j++){
          data_(i,j)-=o.data()(i,j);
        }
      }
      return *this;
    }

    inline auto operator*(const SpatialBlockField & o) const
    {
      SpatialBlockField r(nx_,ny_,nz_);
      for (int i=0; i<nx_; i++){
        for (int j=0; j<ny_; j++){
          r.data_(i,j)=data_(i,j)*o.data()(i,j);
        }
      }
      return r;
    }

    inline auto operator*(const RealType & v) const
    {
      SpatialBlockField r(nx_,ny_,nz_);
      for (int i=0; i<nx_; i++){
        for (int j=0; j<ny_; j++){
          r.data_(i,j)=v*data_(i,j);
        }
      }
      return r;
    }

    inline friend SpatialBlockField operator*(const RealType & v, const SpatialBlockField & o)
    {
      auto nx=o.dimension(X);
      auto ny=o.dimension(Y);
      auto nz=o.dimension(Z);

      SpatialBlockField r(nx,ny,nz);
      for (int i=0; i<nx; i++){
        for (int j=0; j<ny; j++){
          r.data_(i,j)=v*o.data_(i,j);
        }
      }
      return r;
    }

    inline auto & operator*=(const RealType & v)
    {
      for (int i=0; i<nx_; i++){
        for (int j=0; j<ny_; j++){
          data_(i,j)*=v;
        }
      }
      return *this;
    }

    inline auto & operator*=(const SpatialBlockField & o)
    {
      for (int i=0; i<nx_; i++){
        for (int j=0; j<ny_; j++){
          data_(i,j)*=o.data()(i,j);
        }
      }
      return *this;
    }

    inline friend SpatialBlockField operator/(const RealType & v, const SpatialBlockField & o)
    {
      auto nx=o.dimension(X);
      auto ny=o.dimension(Y);
      auto nz=o.dimension(Z);

      SpatialBlockField r(nx,ny,nz);
      for (int i=0; i<nx; i++){
        for (int j=0; j<ny; j++){
          r.data_(i,j)=v/o.data_(i,j);
        }
      }

      return r;
    }

    inline friend SpatialBlockField operator/(const SpatialBlockField & o, const RealType & v)
    {
      auto nx=o.dimension(X);
      auto ny=o.dimension(Y);
      auto nz=o.dimension(Z);

      SpatialBlockField r(nx,ny,nz);
      for (int i=0; i<nx; i++){
        for (int j=0; j<ny; j++){
          r.data_(i,j)=o.data_(i,j)/v;
        }
      }

      return r;
    }

    inline auto & operator/=(const RealType & v)
    {
      for (int i=0; i<nx_; i++){
        for (int j=0; j<ny_; j++){
          data_(i,j)/=v;
        }
      }
      return *this;
    }

    inline auto & operator/=(const SpatialBlockField & o)
    {
      for (int i=0; i<nx_; i++){
        for (int j=0; j<ny_; j++){
          data_(i,j)/=o.data()(i,j);
        }
      }
      return *this;
    }

#ifdef BOOST_SIMD_VECTORIZATION
    static inline auto getPack(const SpatialBlockField1D & v1D, const int kp)
    {
      const size_t offset=kp*Eigen::internal::packet_traits<RealType>::size;
      return boost::simd::aligned_load<boost::simd::pack<RealType>>((RealType *) (v1D.data()+offset));
    }

    static inline void setPack(const boost::simd::pack<RealType> & pack, SpatialBlockField1D & v1D, const int kp)
    {
      const size_t offset=kp*Eigen::internal::packet_traits<RealType>::size;
      boost::simd::aligned_store(pack, v1D.data()+offset);
    }
#endif

    inline void display() const
    {
      for (int k=kStart(); k<kEnd(); k++){
        std::cout << "k = " << k-hnz_ << " *********************************" << std::endl;
        for (int i=iStart(); i<iEnd(); i++){
          for (int j=jStart(); j<jEnd(); j++){
            std::cout << data_(i,j)(k) << " ";
          }
          std::cout << std::endl;
        }
      }
    }

    inline const auto norm2() const
    {
      const auto _iStart=iStart();
      const auto _jStart=jStart();
      const auto _kStart=kStart();

      const auto _iEnd=iEnd();
      const auto _jEnd=jEnd();
      const auto _kEnd=kEnd();

      SWS::RealType n2=0.0;

      for (int i=_iStart; i<_iEnd; i++){
        for (int j=_jStart; j<_jEnd; j++){
          for (int k=_kStart; k<_kEnd; k++){
            n2+=get(i,j,k)*get(i,j,k);
          }
        }
      }

      return n2;
    }

    inline void plot1D(const int j0, const int k0,
                       const std::string name) const
    {
#ifdef USE_MATPLOTLIB
      std::ofstream of(name + "-" + std::to_string(j0) + "-" + std::to_string(k0) + ".txt");

      namespace mpl=matplotlibcpp;

      const auto n=nx_;
      std::vector<RealType> x(n), y(n);
      for(int i=0; i<nx_; i++){
        x.at(i)=i;
        y.at(i)=data_(i,j0)(k0);

        of << int(x.at(i)) << " " << std::scientific << y.at(i) << std::endl;
      }

      of.close();


      mpl::named_plot(name, x, y);

      mpl::xlim(0, nx_);

      mpl::xlabel("i");
      mpl::ylabel(name + "(i," + std::to_string(j0) + "," + std::to_string(k0) + ")");

      mpl::legend();

      mpl::save("./" + name + ".pdf");
#endif
    }

    inline void plot2D(const int k0,
                       const std::string name) const
    {
#ifdef USE_MATPLOTLIB
      std::ofstream of(name + "-" + std::to_string(k0) + ".txt");

      for(int i=0; i<nx_; i++){
        for(int j=0; j<ny_; j++){
          of << i << " " << j << " " << std::scientific << data_(i,j)(k0) << std::endl;
        }
      }

      of.close();
#endif
    }

    inline void plot3D() const
    {
    }

    inline void plotHalo(const SWS::Locations l, const std::string name) const
    {
#ifdef USE_MATPLOTLIB

      switch(l){
      case LEFT:{
        for (int i=0; i<hnx_; i++){
          std::ofstream of(name + "-" + std::to_string(i) + ".txt");

          for (int j=0; j<ny_; j++){
            for (int k=0; k<nz_; k++){
              of << j << " " << k << " " << std::scientific << data_(i,j)(k) << std::endl;
            }
          }

          of.close();
        }
        break;
      }
      case RIGHT:{
        for (int i=nx_-hnx_; i<nx_; i++){
          std::ofstream of(name + "-" + std::to_string(i) + ".txt");

          for (int j=0; j<ny_; j++){
            for (int k=0; k<nz_; k++){
              of << j << " " << k << " " << std::scientific << data_(i,j)(k) << std::endl;
            }
          }

          of.close();
        }
        break;
      }
      default:
        break;
      }

#endif
    }

  private:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW // To force an object of this class to be allocated as aligned

    enum{
      X=0,
      Y,
      Z,
      DIM
    };

    const int packetSize_=Eigen::internal::packet_traits<RealType>::size;

    int nx_;
    int ny_;
    int nz_;

    /* Halo size along x, y and z */
    // FIXME find a better way to define the halo size!!!
    int hnx_=2;
    int hny_=2;
    int hnz_=2;

    /* Padding along x, y and z */
    int px_=0;
    int py_=0;
    int pz_=0;

    int n_;

    /* Indicates whether this spatial block contains a source */
    bool hasSource_=false;

    /* Coordinates of the sources located within this spatial block */
    std::vector<int> is_;
    std::vector<int> js_;
    std::vector<int> ks_;

    // (i,j)(k)
    Eigen::Array<SpatialBlockField1D, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor|Eigen::AutoAlign> data_;
  };

}
