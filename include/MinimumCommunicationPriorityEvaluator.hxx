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

class MinimumCommunicationPriorityEvaluator{
public:
  MinimumCommunicationPriorityEvaluator(const int nxx, const int nyy, const int nzz,
                                        const int lnxx, const int lnyy, const int lnzz) : nxx_(nxx), nyy_(nyy), nzz_(nzz),
                                                                                          lnxx_(lnxx), lnyy_(lnyy), lnzz_(lnzz)
  {
    taskGroups_(COMPUTE_VELOCITY)=4;
    taskGroups_(DISPLAY_VELOCITY)=4;
    taskGroups_(COMPUTE_STRESS)=4;
    taskGroups_(DISPLAY_STRESS)=4;
    taskGroups_(UPDATE_VELOCITY)=3;
    taskGroups_(UPDATE_STRESS)=3;
    taskGroups_(EXTRACT_VELOCITY_HALO)=2;
    taskGroups_(EXTRACT_STRESS_HALO)=2;
    taskGroups_(INITIALIZE_FIELDS)=1;
  }

  ~MinimumCommunicationPriorityEvaluator()
  {
  }

  inline void buildPriorityField(SWS::TaskPriorities & taskPriorities)
  {
    for (auto taskType : {COMPUTE_VELOCITY, EXTRACT_VELOCITY_HALO, UPDATE_VELOCITY, DISPLAY_VELOCITY, COMPUTE_STRESS, EXTRACT_STRESS_HALO, UPDATE_STRESS, DISPLAY_STRESS, INITIALIZE_FIELDS}){
      taskPriorities(taskType).resize(nxx_, nyy_, nzz_);
      taskPriorities(taskType).setZero();
    }
  }

  template<typename MESH3D_PARTITIONING>
  inline void evaluate(SWS::TaskPriorities & taskPriorities)
  {
    buildPriorityField(taskPriorities);

    // Only one MPI node
    if (nxx_==lnxx_ && nyy_==lnyy_ && nzz_==lnzz_){
      return;
    }

    // Compute Maximum priority
    int nCoucheMax = (std::max)({lnxx_, lnyy_, lnzz_});
    int offsetRankZero = 0;
    if ( nxx_!=lnxx_ ){
      nCoucheMax = (std::min)(nCoucheMax, lnxx_);
      offsetRankZero = (std::max)(offsetRankZero, lnxx_);
    }

    if ( nyy_!=lnyy_ ){
      nCoucheMax = (std::min)(nCoucheMax, lnyy_);
      offsetRankZero = (std::max)(offsetRankZero, lnyy_);
    }

    if ( nzz_!=lnzz_ ){
      nCoucheMax = (std::min)(nCoucheMax, lnzz_);
      offsetRankZero = (std::max)(offsetRankZero, lnzz_);
    }

    const int priorityMax = 4*( nCoucheMax ) * ( ( offsetRankZero ) + ( offsetRankZero - nCoucheMax + 1  ) ) / 2 - 1;


    //We compute the maximum distance to a communication border for each ii

    std::vector<int> distanceToCommunicationBorderXArray(nxx_);
    for(int ii=0; ii<nxx_; ii++){
      const int lii=MESH3D_PARTITIONING::getInstance()->lii(ii);

      if(nxx_!=lnxx_){
        if( ii < lnxx_){ // Node 0
          distanceToCommunicationBorderXArray[ii]=lnxx_-lii-1;
        } else if (ii > nxx_-lnxx_) { // Last node
          distanceToCommunicationBorderXArray[ii]=lii;
        } else{
          distanceToCommunicationBorderXArray[ii]=(std::min)(lnxx_-lii-1, lii);
        }
      }
    }


    //We compute the maximum distance to a communication border for each ii

    std::vector<int> distanceToCommunicationBorderYArray(nyy_);
    for(int jj=0; jj<nyy_;jj++){
      const int ljj=MESH3D_PARTITIONING::getInstance()->ljj(jj);

      if(nyy_!=lnyy_) {
        if( jj < lnyy_) { // Node 0
          distanceToCommunicationBorderYArray[jj]=lnyy_-ljj-1;
        } else if ( jj > nyy_-lnyy_) { // Last node
          distanceToCommunicationBorderYArray[jj]=ljj;
        } else {
          distanceToCommunicationBorderYArray[jj]=(std::min)(lnyy_-ljj-1, ljj);
        }
      }
    }

    //We compute the maximum distance to a communication border for each ii

    std::vector<int> distanceToCommunicationBorderZArray(nzz_);
    for(int kk=0; kk<nzz_;kk++){
      const int lkk=MESH3D_PARTITIONING::getInstance()->lkk(kk);

      if(nzz_!=lnzz_){
        if( kk < lnzz_){ // Node 0
          distanceToCommunicationBorderZArray[kk]=lnzz_-lkk-1;
        } else if ( kk > nzz_-lnzz_) { // Last node
          distanceToCommunicationBorderZArray[kk]=lkk;
        } else {
          distanceToCommunicationBorderZArray[kk]=(std::min)(lnzz_-lkk-1, lkk);
        }
      }
    }


    // We compute priorities
    int minDistanceToCommunicationBorderX;
    int maxDistanceToCommunicationBorderX;
    int minDistanceToCommunicationBorderXY;
    int maxDistanceToCommunicationBorderXY;
    int minDistanceToCommunicationBorder;
    int maxDistanceToCommunicationBorder;

    for(int ii=0; ii<nxx_; ii++){
      if(nxx_!=lnxx_){
        minDistanceToCommunicationBorderX = distanceToCommunicationBorderXArray[ii];
        maxDistanceToCommunicationBorderX = distanceToCommunicationBorderXArray[ii];
      } else{
        minDistanceToCommunicationBorderX = (std::max)({lnxx_, lnyy_, lnzz_});
        maxDistanceToCommunicationBorderX = 0;
      }

      for(int jj=0; jj<nyy_; jj++){
        if(nyy_!=lnyy_) {
          minDistanceToCommunicationBorderXY = (std::min)(minDistanceToCommunicationBorderX, distanceToCommunicationBorderYArray[jj]);
          maxDistanceToCommunicationBorderXY = (std::max)(maxDistanceToCommunicationBorderX, distanceToCommunicationBorderYArray[jj]);
        } else {
          minDistanceToCommunicationBorderXY = minDistanceToCommunicationBorderX;
          maxDistanceToCommunicationBorderXY = maxDistanceToCommunicationBorderX;
        }

        for(int kk=0; kk<nzz_; kk++){
          if(nzz_!=lnzz_){
            minDistanceToCommunicationBorder = (std::min)(minDistanceToCommunicationBorderXY, distanceToCommunicationBorderZArray[kk]);
            maxDistanceToCommunicationBorder = (std::max)(maxDistanceToCommunicationBorderXY, distanceToCommunicationBorderZArray[kk]);
          }else{
            minDistanceToCommunicationBorder = minDistanceToCommunicationBorderXY;
            maxDistanceToCommunicationBorder = maxDistanceToCommunicationBorderXY;
          }

          const int nCouche = minDistanceToCommunicationBorder;

          for (auto taskType : {COMPUTE_VELOCITY, EXTRACT_VELOCITY_HALO, UPDATE_VELOCITY, DISPLAY_VELOCITY, COMPUTE_STRESS, EXTRACT_STRESS_HALO, UPDATE_STRESS, DISPLAY_STRESS, INITIALIZE_FIELDS}){
            auto offsetFactor=taskGroups_(taskType);

            const int inversePriority =
              nTaskGroups_*nCouche*( ( offsetRankZero ) + ( offsetRankZero - nCouche + 1 ) ) / 2
              + (offsetRankZero - nCouche)*(offsetFactor - 1)
              + (maxDistanceToCommunicationBorder - minDistanceToCommunicationBorder);

            taskPriorities(taskType)(ii,jj,kk) = priorityMax - inversePriority;
          }

        } // kk
      } // jj
    } //ii



    // Version de Wil
    // const int a=(std::max)({lnxx_, lnyy_, lnzz_});
    // const int b=(std::min)(ii+1, lnxx_-lii) + (std::min)(jj+1, lnyy_-ljj) + (std::min)(kk+1, lnzz_-lkk);
    // const int c=(lnxx_ + lnyy_ + lnzz_)/2;

    // p=3*a-b-c;

    // Version de Salli
    // p+=(lii - lnxx_/2)*(lii - lnxx_/2);
    // p+=(ljj - lnyy_/2)*(ljj - lnyy_/2);
    // p+=(lkk - lnzz_/2)*(lkk - lnzz_/2);

    // fprintf(stdout, "Priority(%d,%d,%d)=%d\n", ii, jj, kk, p);

    // Version Cedric par fonction
    // La priorité correspond à la valeur maximale du tableau moins la valeur de la case correspondante
    // Ici p26-valeur
    // Extract
    // +----+----+----+----+----+----+----+
    // |  0 |  1 |  2 |  3 |  2 |  1 |  0 |
    // +----+----+----+----+----+----+----+
    // |  1 | 12 | 13 | 14 | 13 | 12 |  1 |
    // +----+----+----+----+----+----+----+
    // |  2 | 13 | 21 | 22 | 21 | 13 |  2 |
    // +----+----+----+----+----+----+----+
    // |  1 | 12 | 13 | 14 | 13 | 12 |  1 |
    // +----+----+----+----+----+----+----+
    // |  0 |  1 |  2 |  3 |  2 |  1 |  0 |
    // +----+----+----+----+----+----+----+
    // Update
    // +----+----+----+----+----+----+----+
    // |  4 |  5 |  6 |  7 |  6 |  5 |  4 |
    // +----+----+----+----+----+----+----+
    // |  5 | 15 | 16 | 17 | 16 | 15 |  5 |
    // +----+----+----+----+----+----+----+
    // |  6 | 16 | 23 | 24 | 23 | 16 |  6 |
    // +----+----+----+----+----+----+----+
    // |  5 | 15 | 16 | 17 | 16 | 15 |  5 |
    // +----+----+----+----+----+----+----+
    // |  4 |  5 |  6 |  7 |  6 |  5 |  4 |
    // +----+----+----+----+----+----+----+
    // Compute
    // +----+----+----+----+----+----+----+
    // |  8 |  9 | 10 | 11 | 10 |  9 |  8 |
    // +----+----+----+----+----+----+----+
    // |  9 | 18 | 19 | 20 | 19 | 18 |  9 |
    // +----+----+----+----+----+----+----+
    // | 10 | 19 | 25 | 26 | 25 | 19 | 10 |
    // +----+----+----+----+----+----+----+
    // |  9 | 18 | 19 | 20 | 19 | 18 |  9 |
    // +----+----+----+----+----+----+----+
    // |  8 |  9 | 10 | 11 | 10 |  9 |  8 |
    // +----+----+----+----+----+----+----+
  }

private:
  SWS::TaskGroups taskGroups_;

  const int nTaskGroups_=4;

  /* Number of sub-blocks on each of the three axes */
  int nxx_;
  int nyy_;
  int nzz_;

  /* Local number of sub-blocks on each of the three axes */
  int lnxx_;
  int lnyy_;
  int lnzz_;
};
