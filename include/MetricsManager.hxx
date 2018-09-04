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
#include <unordered_map>
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
#include <thread>
#include <atomic>
#include <mutex>

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
  static MetricsManager * getInstance(const std::string applicationName,
                                      const std::initializer_list<std::string> sampledEvents);
  static MetricsManager * getInstance();

  static void releaseInstance();

private:
  enum Metric{ELAPSED, CPU, FLOPS, NB_METRICS};
  enum Measurement{START, STOP, DURATION};
  enum MetricSummary{NCALLS, AVG, MIN, MAX, TOT, SD, NB_STATISTICS};

#ifdef COLLECT_STATS
  // (ELAPSED|CPU)(START, STOP, DURATION)(measurements)
  typedef tbb::concurrent_vector< tbb::concurrent_vector< tbb::concurrent_queue<double> > > CollectedMeasurements;

  typedef std::array<double, NB_STATISTICS> ComputedStatistics;
  typedef std::array<ComputedStatistics, NB_METRICS> ComputedMetrics;

  /* Origin of time */
  const std::chrono::high_resolution_clock::time_point oElapsedTime_=std::chrono::high_resolution_clock::now();
  const std::clock_t oCPUTime_=std::clock();
#endif

#ifdef COLLECT_STATS
  template<Metric t, Measurement m>
  inline auto get(const int tn, const std::string & eid)
  {
    auto & events=events_[tn][eid][t][m];

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

    return ComputedStatistics{double(ncalls), avg, min, max, tot, sd};
  }

  inline auto getThreadNumber()
  {
    static std::atomic<int> tn{0};

    auto tid=std::this_thread::get_id();

#ifdef COLLECT_STATS
    if (threadId2threadNumber_.find(tid) == threadId2threadNumber_.end()){
      std::lock_guard<std::mutex> lock(guard_);

      threadId2threadNumber_[tid]=tn++;

      /* Allocate space for storing collected metrics per thread and for all sampled events */
      buildEvents(threadId2threadNumber_[tid]);
    }
#endif

    return threadId2threadNumber_[tid];
  }
#endif

  inline void buildEvents(const int tn)
  {
#ifdef COLLECT_STATS
    for (auto eid : sampledEvents_){
      events_[tn][eid].resize(NB_METRICS);
      for (auto & event : events_[tn][eid]){
        event.resize(NB_STATISTICS);
      }
    }
#endif
  }

public:
  inline void start(const std::string eid, const double flops=0)
  {
#ifdef COLLECT_STATS
    auto tn=getThreadNumber();

    if (events_[tn].find(eid) != events_[tn].end()){
      events_[tn][eid][Metric::ELAPSED][START].push(now<Metric::ELAPSED>());
      events_[tn][eid][Metric::CPU][START].push(now<Metric::CPU>());
    }
#endif
  }

  inline void stop(const std::string eid)
  {
#ifdef COLLECT_STATS
    auto tn=getThreadNumber();

    if (events_[tn].find(eid) != events_[tn].end()){
      events_[tn][eid][Metric::ELAPSED][STOP].push(now<Metric::ELAPSED>());
      events_[tn][eid][Metric::CPU][STOP].push(now<Metric::CPU>());

      auto cpu    =get<Metric::CPU,     Measurement::STOP>(tn, eid)-get<Metric::CPU,     Measurement::START>(tn, eid);
      auto elapsed=get<Metric::ELAPSED, Measurement::STOP>(tn, eid)-get<Metric::ELAPSED, Measurement::START>(tn, eid);

      events_[tn][eid][Metric::CPU][DURATION].push(cpu);
      events_[tn][eid][Metric::ELAPSED][DURATION].push(elapsed);
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

    /* Compute per-thread statistics from collected data */
    auto nthreads=events_.size();
    cmetrics_.resize(nthreads+1);
    for (auto tn=0; tn<nthreads; tn++){
      for (auto & it : events_[tn]){
        auto sevent=it.first;

        ComputedMetrics cm;
        for (auto & metric : {ELAPSED, CPU}){
          auto qmetric=it.second[metric][DURATION];
          cm[metric]=evalStats(qmetric);
        }

        cmetrics_[tn].push_back(std::make_pair<std::string, ComputedMetrics>(std::move(sevent), std::move(cm)));
      }
    }

    /* Aggregate statistics for all threads */
    for (auto sevent : sampledEvents_){

      ComputedMetrics cm;
      for (auto metric : {ELAPSED, CPU}){
        cm[metric][NCALLS]=0.0;
        cm[metric][AVG]=0.0;
        cm[metric][MIN]=1.e99;
        cm[metric][MAX]=0.0;
        cm[metric][TOT]=0.0;
        cm[metric][SD]=0.0;
      }

      // FIXME optimize the following reduction
      for (auto tn=0; tn<nthreads; tn++){
        for (auto & it : cmetrics_[tn]){
          const auto & cm_n=it.second;

          if (sevent != it.first)
            continue;

          for (auto metric : {ELAPSED, CPU}){
            cm[metric][NCALLS]+=cm_n[metric][NCALLS];

            if (0 == cm_n[metric][NCALLS]){
              continue;
            }

            cm[metric][MIN]=std::min(cm[metric][MIN], cm_n[metric][MIN]);
            cm[metric][MAX]=std::max(cm[metric][MAX], cm_n[metric][MAX]);
            cm[metric][TOT]+=cm_n[metric][TOT];

            cm[metric][SD]+=(cm_n[metric][SD]*cm_n[metric][SD] + cm_n[metric][AVG]*cm_n[metric][AVG])*cm_n[metric][NCALLS];
          }
        }
      } // tn

      for (auto metric : {ELAPSED, CPU}){
        cm[metric][AVG]=cm[metric][TOT]/cm[metric][NCALLS];
        cm[metric][SD]=std::sqrt(cm[metric][SD]/cm[metric][NCALLS] - cm[metric][AVG]*cm[metric][AVG]);
      }

      cmetrics_[nthreads].push_back(std::make_pair<std::string, ComputedMetrics>(std::move(sevent), std::move(cm)));
    } // sevent


    /* Sort the statistics according to the average elapsed time */
    sort<ELAPSED, AVG>();


    /* Print the computed statistics */
#ifdef SHOW_PER_THREAD_STATS
    for (auto tn=0; tn<nthreads+1; tn++){
      if (nthreads == tn)
        std::cout << "Aggregated statistics\n";
      else
        std::cout << "Thread ID: " << tn << "\n";
      std::cout << "------------------------------------------\n";
#else
    for (auto tn=nthreads; tn<nthreads+1; tn++){
#endif
      for (auto & it : cmetrics_[tn]){
        std::cout << it.first << "\n";

        for (auto metric : {ELAPSED, CPU}){
          const auto & stats=it.second[metric];

          if (0 == stats[NCALLS]){
            std::cout << "\tNONE" << std::endl;
            continue;
          }

          switch(metric){
          case ELAPSED:
            std::cout << "\t#Calls           : " << "\t" << int(stats[NCALLS]) << std::endl;
            std::cout << "\tElapsed Time (s) : ";
            break;
          case CPU:
            std::cout << "\tCPU Time (s)     : ";
            break;
          default:
            break;
          }

          std::cout << std::scientific;
          std::cout << "\t" << stats[AVG];
          std::cout << "\t" << stats[MIN];
          std::cout << "\t" << stats[MAX];
          std::cout << "\t" << stats[TOT];
          std::cout << "\t" << stats[SD];
          std::cout << std::endl;
        }
      }
      std::cout << "------------------------------------------\n";
    } // tn
#endif
  }

private:
  MetricsManager(const std::string applicationName,
                 const std::initializer_list<std::string> sampledEvents);
  ~MetricsManager();

#ifdef COLLECT_STATS
  template<Metric m, MetricSummary ms>
  inline void sort()
  {
    for (auto tn=0; tn<cmetrics_.size(); tn++){
      std::sort(cmetrics_[tn].begin(), cmetrics_[tn].end(),
                [](const auto & cm1, const auto & cm2)
                {
                  return cm1.second[m][ms] > cm2.second[m][ms];
                });
    }
  }
#endif

  static MetricsManager * pInstance_;

  const std::string applicationName_;

  const std::vector<std::string> sampledEvents_;

#ifdef COLLECT_STATS
  std::mutex guard_;

  // {"thread num" : {"event id" : [ELAPSED, CPU][START, STOP]}}
  std::unordered_map<int, tbb::concurrent_unordered_map<std::string, CollectedMeasurements>> events_;

  /* The list of computed metrics from collected data */
  std::vector<std::vector<std::pair<std::string, ComputedMetrics>>> cmetrics_;

  // Map a thread id to a unique number
  std::unordered_map<std::thread::id, int> threadId2threadNumber_;
#endif
};
