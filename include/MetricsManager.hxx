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

#include <chrono>
#include <ctime>
#include <map>
#include <string>
#include <vector>
#include <tuple>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <utility>
#include <queue>

#ifdef COLLECT_STATS
#include <tbb/concurrent_queue.h>
#include <tbb/concurrent_vector.h>
#include <tbb/concurrent_unordered_map.h>

#ifdef USE_PAPI
#include <papi.h>
#endif
#endif // COLLECT_STATS

class MetricsManager{
public:
  static MetricsManager * getInstance(const std::string applicationName);
  static MetricsManager * getInstance();

  static void releaseInstance();

private:
  enum Metric{ELAPSED, CPU, FLOPS, NB_METRICS};
  enum Measurement{START, STOP, DURATION};
  enum MetricSummary{AVG, MIN, MAX, TOT, SD, NB_STATISTICS};

#ifdef COLLECT_STATS
  // (ELAPSED|CPU)(START, STOP, DURATION)(measurements)
  typedef tbb::concurrent_vector< tbb::concurrent_vector< tbb::concurrent_queue<double> > > CollectedMeasurements;

  typedef std::array<double, NB_STATISTICS> ComputedStatistics;
  typedef std::array<ComputedStatistics, NB_METRICS> ComputedMetrics;

  const std::chrono::high_resolution_clock::time_point oElapsedTime_=std::chrono::high_resolution_clock::now();
  const std::clock_t oCPUTime_=std::clock();
#endif

#ifdef COLLECT_STATS
  template<Metric t, Measurement m>
  inline auto get(const std::string & id)
  {
    auto & events=events_[id][t][m];

    double e=0.0;
    events.try_pop(e);

    return e;
  }

  template<Metric m>
  inline auto now() const
  {
    double n;

    switch(m){
    case Metric::ELAPSED:{
      // Returns the elapsed time
      auto now=std::chrono::high_resolution_clock::now();
      n=std::chrono::duration<double>(now-oElapsedTime_).count();
      break;
    }
    case Metric::CPU:{
      // Returns the CPU time
      auto now=std::clock();
      n=double(now-oCPUTime_)/CLOCKS_PER_SEC;
      break;
    }
    case Metric::FLOPS:{
      // Returns the flops
#ifdef USE_PAPI
      /* Setup PAPI library and begin collecting data from the counters */
      if(PAPI_flops(&papiElapsedTime_, &papiCPUTime_, &papiFlpIns_, &papiMFlops_) < PAPI_OK){
        printf("Error starting the counters, aborting.\n");
        exit(-1);
      }

      auto now=papiMFlops_;
      n=double(now-oPAPIMFlops_);
#endif
      break;
    }
    default:
      break;
    }

    return n;
  }

  // inline auto evalStats(std::deque<double> & q)
  inline auto evalStats(tbb::concurrent_queue<double> & q)
  {
    double avg=0.0, min=1.e99, max=0.0, tot=0.0, sd=0.0;

    auto ncalls=q.unsafe_size();

    while (!q.empty()){
      double measure=0.0;
      q.try_pop(measure);

      min=(measure < min) ? measure : min;
      max=(measure > max) ? measure : max;
      tot+=measure;
      sd+=measure*measure;
    }
    sd/=ncalls;

    avg=tot/ncalls;

    sd=std::sqrt(sd-avg*avg);

    return ComputedStatistics{avg, min, max, tot, sd};
  }
#endif

public:
  inline void setSampledEvents(const std::initializer_list<std::string> events)
  {
#ifdef COLLECT_STATS
    for (auto & id : events){
      if (events_.find(id) == events_.end()){
        events_[id].resize(NB_METRICS);
        for (auto & event : events_[id]){
          event.resize(NB_STATISTICS);
        }
      }
    }
#endif
  }

  inline void start(const std::string id, const double flops=0)
  {
#ifdef COLLECT_STATS
    if (events_.find(id) != events_.end()){
      events_[id][Metric::ELAPSED][START].push(now<Metric::ELAPSED>());
      events_[id][Metric::CPU][START].push(now<Metric::CPU>());
    }
#endif
  }

  inline void stop(const std::string id)
  {
#ifdef COLLECT_STATS
    if (events_.find(id) != events_.end()){
      events_[id][Metric::ELAPSED][STOP].push(now<Metric::ELAPSED>());
      events_[id][Metric::CPU][STOP].push(now<Metric::CPU>());

      auto cpu    =get<Metric::CPU,     Measurement::STOP>(id)-get<Metric::CPU,     Measurement::START>(id);
      auto elapsed=get<Metric::ELAPSED, Measurement::STOP>(id)-get<Metric::ELAPSED, Measurement::START>(id);

      events_[id][Metric::CPU][DURATION].push(cpu);
      events_[id][Metric::ELAPSED][DURATION].push(elapsed);
    }
#endif
  }

  inline void show()
  {
#ifdef COLLECT_STATS
    std::cout << "++++++++++++++++++++++++++++++++++++++++++\n";
    std::cout << " Statistics \n";
    std::cout << "------------------------------------------\n";
    std::cout << " " << applicationName_ << " \n";
    std::cout << "++++++++++++++++++++++++++++++++++++++++++\n";

    std::cout << std::setw(32) << " ";
    std::cout << std::left << std::setw(12) << "AVG\t";
    std::cout << std::left << std::setw(12) << "MIN\t";
    std::cout << std::left << std::setw(12) << "MAX\t";
    std::cout << std::left << std::setw(12) << "TOT\t";
    std::cout << std::left << std::setw(12) << "SD";
    std::cout << std::endl;

    /* Compute statistics from collected data */
    for (auto & it : events_){
      auto sname=it.first;

      ComputedMetrics cm;
      for (auto & metric : {ELAPSED, CPU}){
        auto qmetric=it.second[metric][DURATION];
        cm[metric]=evalStats(qmetric);
      }

      cmetrics_.push_back(std::make_pair<std::string, ComputedMetrics>(std::move(sname), std::move(cm)));
    }

    /* Sort the statistics according to the avarage elapsed time */
    sort<ELAPSED, AVG>();

    /* Print the computed statistics */
    for (auto & it : cmetrics_){
      std::cout << it.first << "\n";

      for (auto metric : {ELAPSED, CPU}){
        const auto ncalls=events_[it.first][metric][DURATION].unsafe_size();

        switch(metric){
        case ELAPSED:
          std::cout << "\t#Calls           : " << "\t" << ncalls << std::endl;
          std::cout << "\tElapsed Time (s) : ";
          break;
        case CPU:
          std::cout << "\tCPU Time (s)     : ";
          break;
        default:
          break;
        }

        std::cout << std::scientific;

        const auto & stats=it.second[metric];
        std::cout << "\t" << stats[AVG];
        std::cout << "\t" << stats[MIN];
        std::cout << "\t" << stats[MAX];
        std::cout << "\t" << stats[TOT];
        std::cout << "\t" << stats[SD];
        std::cout << std::endl;
      }
    }
#endif
  }

private:
  MetricsManager(const std::string applicationName);
  ~MetricsManager();

#ifdef COLLECT_STATS
  template<Metric m, MetricSummary ms>
  inline void sort()
  {
    std::sort(cmetrics_.begin(), cmetrics_.end(),
              [](const auto & cm1, const auto & cm2)
              {
                return cm1.second[m][ms] > cm2.second[m][ms];
              });
  }
#endif

  static MetricsManager * pInstance_;

  const std::string applicationName_;

#ifdef COLLECT_STATS
  // {"event id" : [ELAPSED, CPU][START, STOP]}
  tbb::concurrent_unordered_map<std::string, CollectedMeasurements> events_;

  std::vector<std::pair<std::string, ComputedMetrics>> cmetrics_;
#endif
};
