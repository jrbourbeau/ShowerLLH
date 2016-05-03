#ifndef TIMESCRAMBLE_H_INCLUDED
#define TIMESCRAMBLE_H_INCLUDED

#include <icetray/I3FrameObject.h>
#include <icetray/serialization.h>
#include <cstdio>
#include <string>
#include <vector>

using namespace std;

class TimeScramble : public I3FrameObject {
public:

   vector<double> cscramble(const vector<double> mjdList, const vector<double> thetaList, 
                  const vector<double> phiList, int nSide, int nBGResample, 
                  string method) const;
   //cscramble(const vector<double> mjdList, const vector<double> thetaList, 
   //               const vector<double> phiList, int nSide, int nBGResample, 
   //               vector<double> BGMap, string method) const;

private:
  friend class boost::serialization::access;
  template <typename Archive>
  void serialize(Archive &ar, unsigned version);

};

I3_POINTER_TYPEDEFS(TimeScramble);

#endif // TIMESCRAMBLE_H_INCLUDED
