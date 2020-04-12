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

  inline void parseDataFile()
  {
    namespace pt = boost::property_tree;

    pt::ptree root;
    pt::read_json(dfile(), root);

    model_ = root.get<std::string>("Model");

    tmax_ = root.get<SWS::RealType>("tmax");
    dt_ = root.get<SWS::RealType>("dt");

    nt_ = ceil(tmax_ / dt_);

    lx_ = root.get<SWS::RealType>("lx");
    ly_ = root.get<SWS::RealType>("ly");
    lz_ = root.get<SWS::RealType>("lz");

    ds_ = root.get<SWS::RealType>("ds");

    // TODO handle the case where ds l{x,y,z}_ is not a multiple of ds_
    nx_ = lx_ / ds_;
    ny_ = ly_ / ds_;
    nz_ = lz_ / ds_;

    /* Layers */
    for (auto& layer : root.get_child("Layers")) {
      const auto& name = layer.second.get<std::string>("Name");

      layers_.push_back(name);

      start_[name].resize(SWS::DIM);
      end_[name].resize(SWS::DIM);

      const auto& size = layer.second.get_child("Size");
      pt::ptree s;
      for (auto d : { SWS::X, SWS::Y, SWS::Z }) {
        switch (d) {
          case SWS::X:
            s = size.get_child("X");
            break;
          case SWS::Y:
            s = size.get_child("Y");
            break;
          case SWS::Z:
            s = size.get_child("Z");
            break;
          default:
            break;
        }

        start_[name][d] = s.get<SWS::RealType>("start");
        end_[name][d] = s.get<SWS::RealType>("end");

        Vp_[name] = layer.second.get<SWS::RealType>("Vp");
        Vs_[name] = layer.second.get<SWS::RealType>("Vs");
        rho_[name] = layer.second.get<SWS::RealType>("rho");
        // Q_[name]=layer.second.get<SWS::RealType>("Q");
      } // d

    } // layer

    /* Sources */
    for (auto& source : root.get_child("Sources")) {
      const auto& name = source.second.get<std::string>("Name");

      sources_.push_back(name);

      const auto& xs = source.second.get<SWS::RealType>("xs");
      const auto& ys = source.second.get<SWS::RealType>("ys");
      const auto& zs = source.second.get<SWS::RealType>("zs");

      xs_[name] = xs;
      ys_[name] = ys;
      zs_[name] = zs;
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
  SWS::RealType lx_;
  SWS::RealType ly_;
  SWS::RealType lz_;
  SWS::RealType ds_;

  int nx_;
  int ny_;
  int nz_;

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
