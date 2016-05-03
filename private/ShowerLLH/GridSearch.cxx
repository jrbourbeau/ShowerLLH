#include "ShowerLLH/GridSearch.h"

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <limits>

using namespace std;

template <typename Archive>
void GridSearch::serialize(Archive &ar, unsigned version) {
  ar & make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
}


//vector<double> GridSearch::csearch(const vector<double> &gridx,
//        const vector<double> &gridy, const vector<double> &llhTable, 
//        const vector<int> &binDims, const vector<string> &binNames,
//        int Zbin, const vector<int> &Sbin, const vector<int> &Cbin,
//        double theta, double phi, const vector<double> &tankx,
//        const vector<double> &tanky, const vector<double> &tankz,
//        const vector<double> &Dbins) const {

vector<double> GridSearch::csearch(const vector<double> &gridx,
        const vector<double> &gridy, const vector<double> &llhTable, 
        const vector<double> &binDims, const vector<string> &binNames,
        vector<double> binnedVals,
        double theta, double phi, const vector<double> &tankx,
        const vector<double> &tanky, const vector<double> &tankz,
        const vector<double> &Dbins) const {

  // Convert double vectors that should be ints (ugh)
  vector<int> binDims_i, binnedVals_i;
  for (unsigned int i=0; i < binDims.size(); ++i)
    binDims_i.push_back(int(binDims[i]));
  for (unsigned int i=0; i < binnedVals.size(); ++i)
    binnedVals_i.push_back(int(binnedVals[i]));

  int nDims  = binNames.size();
  int nG     = gridx.size();
  int eIdx   = find(binNames.begin(), binNames.end(), "E") - binNames.begin();
  int nE     = binDims_i[eIdx];
  int ntanks = tankx.size();
  double inf = numeric_limits<double>::infinity();

  int e_idx = 0;
  int d_idx = 3;

  double RecoMaxLLH = -inf;
  double grid_idx = 0;
  double argmax = 0;
  vector<int> Dbin(ntanks, 0);

  double Z = 1947;
  double ex = -sin(theta) * cos(phi);
  double ey = -sin(theta) * sin(phi);
  double ez = -cos(theta);

  // Loop over every grid position
  for (int i=0; i < nG; ++i) {

    double X = gridx[i];
    double Y = gridy[i];

    // Calculate and bin closest approach distances
    for (int j=0; j < ntanks; ++j) {
      double hx = X - tankx[j];
      double hy = Y - tanky[j];
      double hz = Z - tankz[j];
      double s = ex*hx + ey*hy + ez*hz;
      double x1 = tankx[j] + s*ex;
      double y1 = tanky[j] + s*ey;
      double z1 = tankz[j] + s*ez;
      double dist = sqrt(pow(x1-X,2) + pow(y1-Y,2) + pow(z1-Z,2));
      int bin = 0;
      while (dist > Dbins[bin+1])
        bin += 1;
      binnedVals_i[j*nDims + d_idx] = bin;
    }

    // Calculate llh for every energy bin, return most likely
    int idx, c;
    for (int e=0; e < nE; ++e) {
      double llh = 0.0;
      for (int p=0; p < ntanks; ++p) {
        // Calculate index
        binnedVals_i[p*nDims + e_idx] = e;
        idx = 0;
        c = 1;
        for (int k=nDims; k-- > 0; ) {
          idx += c * binnedVals_i[p*nDims + k];
          c *= binDims_i[k];
        }
        llh += llhTable[idx];
      }
      if (llh > RecoMaxLLH) {
        RecoMaxLLH = llh;
        grid_idx = i;
        argmax = e;
      }
    }
  }

  vector<double> output(3);
  output[0] = RecoMaxLLH;
  output[1] = grid_idx;
  output[1] = argmax;

  return output;
}

I3_SERIALIZABLE(GridSearch);

