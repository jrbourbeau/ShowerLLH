#ifndef GRIDSEARCH_H_INCLUDED
#define GRIDSEARCH_H_INCLUDED

#include <icetray/I3FrameObject.h>
#include <icetray/serialization.h>
#include <cstdio>
#include <string>
#include <vector>

using namespace std;

class GridSearch : public I3FrameObject {
public:

  /*
  vector<double> csearch(const vector<double> &gridx, 
                  const vector<double> &gridy, const vector<double> &llhTable, 
                  const vector<int> &binDims, const vector<string> &binNames,
                  int Zbin, const vector<int> &Sbin, const vector<int> &Cbin,
                  double theta, double phi, const vector<double> &tankx,
                  const vector<double> &tanky, const vector<double> &tankz,
                  const vector<double> &Dbins) const;
  */

   vector<double> csearch(const vector<double> &gridx, 
                  const vector<double> &gridy, const vector<double> &llhTable, 
                  const vector<double> &binDims, const vector<string> &binNames,
                  vector<double> binnedVals, 
                  double theta, double phi, const vector<double> &tankx,
                  const vector<double> &tanky, const vector<double> &tankz,
                  const vector<double> &Dbins) const;

private:
  friend class boost::serialization::access;
  template <typename Archive>
  void serialize(Archive &ar, unsigned version);

};

I3_POINTER_TYPEDEFS(GridSearch);

#endif // GRIDSEARCH_H_INCLUDED
