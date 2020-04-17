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

#include <cmath>
#include <fstream>
#include <iomanip>
#include <map>
#include <string>
#include <vector>

#include <boost/program_options.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include "Config.hxx"
#include "LogManager.hxx"

class SEWASParameterManager
{
public:
  SEWASParameterManager(int* pargc, char*** pargv);
  ~SEWASParameterManager();

  inline auto& argc() { return *pargc_; }
  inline auto& argv() { return *pargv_; }

  inline const auto& cx() const { return getOption<int>("cx"); }
  inline const auto& cy() const { return getOption<int>("cy"); }
  inline const auto& cz() const { return getOption<int>("cz"); }

  inline const auto& nthreads() const { return getOption<int>("nthreads"); }

  inline const auto& P() const { return getOption<int>("P"); }
  inline const auto& Q() const { return getOption<int>("Q"); }
  inline const auto& R() const { return getOption<int>("R"); }

  inline const auto& dfile() const { return getOption<std::string>("dfile"); }

  inline const auto& tmax() const { return tmax_; }
  inline const auto& dt() const { return dt_; }

  inline const auto& nt() const { return nt_; }

  inline const auto& lx() const { return lx_; }
  inline const auto& ly() const { return ly_; }
  inline const auto& lz() const { return lz_; }
  inline const auto& ds() const { return ds_; }

  inline const auto& nx() const { return nx_; }
  inline const auto& ny() const { return ny_; }
  inline const auto& nz() const { return nz_; }

  inline const auto& nxx() const { return nxx_; }
  inline const auto& nyy() const { return nyy_; }
  inline const auto& nzz() const { return nzz_; }

  inline const auto& lnxx() const { return lnxx_; }
  inline const auto& lnyy() const { return lnyy_; }
  inline const auto& lnzz() const { return lnzz_; }

  inline const auto& layers() const { return layers_; }

  inline const auto& start() const { return start_; }
  inline const auto& end() const { return end_; }

  inline const auto& Vp() const { return Vp_; }
  inline const auto& Vs() const { return Vs_; }
  inline const auto& rho() const { return rho_; }
  // inline const auto & Q() const{ return Q_; }

  inline const auto& sources() const { return sources_; }

  inline const auto& xs() const { return xs_; }
  inline const auto& ys() const { return ys_; }
  inline const auto& zs() const { return zs_; }

private:
  template<typename T>
  inline const T& getOption(const std::string& o) const
  {
    if (vm_.count(o)) {
      return vm_[o].as<T>();
    } else {
      std::cerr << std::quoted(o) << " is not set. Exiting...\n";
      exit(-1);
    }
  }

  int parse();

  inline int parseArgs()
  {
    namespace po = boost::program_options;

    auto& argc = this->argc();
    auto& argv = this->argv();

    constexpr char caption[] = "Seismic Wave Simulator (SeWaS) \n\nAllowed options";

    po::options_description desc(caption);

    constexpr int nb_options = 8; // TODO use an enumeration to track all the following options
    desc.add_options()("help,h", "Produce help message")("cx,x", po::value<int>(), "Block size along x-axis")(
      "cy,y", po::value<int>(), "Block size along y-axis")(
      "cz,z", po::value<int>(), "Block size along z-axis")(
      "P,P", po::value<int>(), "Number of MPI processes along x-axis")(
      "Q,Q", po::value<int>(), "Number of MPI processes along y-axis")(
      "R,R", po::value<int>(), "Number of MPI processes along z-axis")(
      "nthreads,T", po::value<int>(), "Number of threads per MPI process")(
      "dfile,D",
      po::value<std::string>(),
      "Path to the Json file containing mechanical properties and kinematic source parameter")(
      "config,C", po::value<std::string>(), "Configuration file");

    try {

      po::store(po::parse_command_line(argc, argv, desc), vm_);

      if (vm_.count("config")) {
        std::ifstream config(getOption<std::string>("config"));
        if (config) {
          po::store(po::parse_config_file(config, desc), vm_);
        }
      }

      po::notify(vm_);

      if (vm_.size() < nb_options) {
        std::cout << desc << "\n";
        exit(-1);
      }

      if (vm_.count("help")) {
        std::cout << desc << "\n";
        return 0;
      }
    } catch (const po::error& e) {
      std::cerr << e.what() << "\n";
      return -1;
    }

    return 0;
  }

  inline void evaluateLayerBoundaries(const std::string& layerName,
                                      const boost::property_tree::ptree& s,
                                      const SWS::Directions d,
                                      const SWS::RealType initial_l,
                                      const SWS::RealType l)
  {
    start_[layerName][d] = s.get<SWS::RealType>("start");
    end_[layerName][d] = s.get<SWS::RealType>("end");

    // Add padding on the last layer
    if (end_[layerName][d] == initial_l) {
      end_[layerName][d] = l;
    }
  }

  inline void parseDataFile()
  {
    namespace pt = boost::property_tree;

    pt::ptree root;
    pt::read_json(dfile(), root);

    model_ = root.get<std::string>("Model");

    tmax_ = root.get<SWS::RealType>("tmax");
    dt_ = root.get<SWS::RealType>("dt");

    nt_ = std::ceil(tmax_ / dt_);

    lx_ = initial_lx_ = root.get<SWS::RealType>("lx");
    ly_ = initial_ly_ = root.get<SWS::RealType>("ly");
    lz_ = initial_lz_ = root.get<SWS::RealType>("lz");

    ds_ = root.get<SWS::RealType>("ds");

    /* Evaluate and adjust discretization/partitioning parameters */
    nx_ = std::ceil(lx_ / ds_);
    ny_ = std::ceil(ly_ / ds_);
    nz_ = std::ceil(lz_ / ds_);

    nxx_ = (nx_ + cx() - 1) / cx();
    nyy_ = (ny_ + cy() - 1) / cy();
    nzz_ = (nz_ + cz() - 1) / cz();

    lnxx_ = (nxx_ + P() - 1) / P();
    lnyy_ = (nyy_ + Q() - 1) / Q();
    lnzz_ = (nzz_ + R() - 1) / R();

    /* Recompute parameters for taking into account all padding */
    nxx_ = lnxx_ * P();
    nyy_ = lnyy_ * Q();
    nzz_ = lnzz_ * R();

    nx_ = nxx_ * cx();
    ny_ = nyy_ * cy();
    nz_ = nzz_ * cz();

    lx_ = SWS::RealType(nx_) * ds_;
    ly_ = SWS::RealType(ny_) * ds_;
    lz_ = SWS::RealType(nz_) * ds_;

    /* Layers */
    for (auto& layer : root.get_child("Layers")) {
      const auto& layerName = layer.second.get<std::string>("Name");

      layers_.push_back(layerName);

      start_[layerName].resize(SWS::DIM);
      end_[layerName].resize(SWS::DIM);

      const auto& size = layer.second.get_child("Size");
      pt::ptree s;
      for (auto d : { SWS::X, SWS::Y, SWS::Z }) {
        switch (d) {
          case SWS::X:
            s = size.get_child("X");

            evaluateLayerBoundaries(layerName, s, d, initial_lx_, lx_);

            break;
          case SWS::Y:
            s = size.get_child("Y");

            evaluateLayerBoundaries(layerName, s, d, initial_ly_, ly_);

            break;
          case SWS::Z:
            s = size.get_child("Z");

            evaluateLayerBoundaries(layerName, s, d, initial_lz_, lz_);

            break;
          default:
            break;
        }

        Vp_[layerName] = layer.second.get<SWS::RealType>("Vp");
        Vs_[layerName] = layer.second.get<SWS::RealType>("Vs");
        rho_[layerName] = layer.second.get<SWS::RealType>("rho");
        // Q_[name]=layer.second.get<SWS::RealType>("Q");
      } // d

    } // layer

    /* Sources */
    for (auto& source : root.get_child("Sources")) {
      const auto& sourceName = source.second.get<std::string>("Name");

      sources_.push_back(sourceName);

      const auto& xs = source.second.get<SWS::RealType>("xs");
      const auto& ys = source.second.get<SWS::RealType>("ys");
      const auto& zs = source.second.get<SWS::RealType>("zs");

      xs_[sourceName] = xs;
      ys_[sourceName] = ys;
      zs_[sourceName] = zs;
    }
  }

  int* pargc_;
  char*** pargv_;

  boost::program_options::variables_map vm_;

  std::string model_;

  /* Time discretization */
  SWS::RealType tmax_;
  SWS::RealType dt_;

  int nt_;

  /* Spatial discretization */
  SWS::RealType initial_lx_;
  SWS::RealType initial_ly_;
  SWS::RealType initial_lz_;

  SWS::RealType lx_;
  SWS::RealType ly_;
  SWS::RealType lz_;

  SWS::RealType ds_;

  int nx_;
  int ny_;
  int nz_;

  int nxx_;
  int nyy_;
  int nzz_;

  /* Spatial blocks count per local MPI process */
  int lnxx_;
  int lnyy_;
  int lnzz_;

  /* Layers structure */
  std::vector<std::string> layers_;

  std::map<std::string, std::vector<SWS::RealType>> start_;
  std::map<std::string, std::vector<SWS::RealType>> end_;

  std::map<std::string, SWS::RealType> Vp_;
  std::map<std::string, SWS::RealType> Vs_;
  std::map<std::string, SWS::RealType> rho_;
  std::map<std::string, SWS::RealType> Q_;

  /* Number of velocity source points */
  std::vector<std::string> sources_;

  /* Velocity source coordinates */
  std::map<std::string, SWS::RealType> xs_;
  std::map<std::string, SWS::RealType> ys_;
  std::map<std::string, SWS::RealType> zs_;
};
