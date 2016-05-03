#include "ShowerLLH/TimeScramble.h"

#include <slalib/slalib.h>
#include <gsl/gsl_rng.h>

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <stdexcept>

#include <healpix-cxx/cxxsupport/constants.h>
//#include <healpix-cxx/cxxsupport/fitshandle.h>
#include <healpix-cxx/healpix/healpix_map.h>
#include <healpix-cxx/healpix/healpix_map_fitsio.h>
#include <healpix-cxx/cxxsupport/pointing.h>

#include <astro/astro.h>
#include <astro/Astro_Time.h>
#include <astro/Astro_Coords.h>
#include <astro/Astro_Detector.h>

using namespace std;

template <typename Archive>
void TimeScramble::serialize(Archive &ar, unsigned version)
{
  ar & make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
}


//void TimeScramble::cscramble(const vector<double> mjdList, 
vector<double> TimeScramble::cscramble(const vector<double> mjdList, 
                const vector<double> thetaList, const vector<double> phiList, 
                int nSide, int nBGResample, 
                //vector<double> BGMap, string method) const {
                string method) const
{
  Astro::IceCubeDetector ice;
  Astro::LocalCoord local;
  Astro::Time t;
  Astro::EqCoord eqApparent;
  Astro::EqCoord eqSolar;

  Healpix_Map<double> BGInt;
  BGInt.SetNside(nSide, RING);
  BGInt.fill(0.);
  int npix = BGInt.Npix();
  vector<double> BGMap (npix);

  pointing sphereDir;
  int pixelId;
  double theta, phi, rndMJD;

  const gsl_rng_type *T = gsl_rng_default;
  gsl_rng *r = gsl_rng_alloc (T);

  // Should find some way to set the random seed
  //gRandom->SetSeed(0);
  int nEntries = mjdList.size();
  for (int iEntry = 0; iEntry < nEntries; ++iEntry) {

    theta = thetaList[iEntry];
    phi = phiList[iEntry];
    for (int k = 0; k < nBGResample; ++k) {

      // Scramble the time
      int rndIndex = gsl_rng_uniform_int(r, nEntries);
      double rndMJD = mjdList[rndIndex];
      t.SetTime(rndMJD);

      //Generate the new equatorial coordinates
      local.SetLocalRad(theta, phi);
      if (method == "anti_sid") {
        eqApparent = ice.LocalToEquatorial_FromAntiSid(local, t);
        sphereDir.phi = eqApparent.GetRaRad();
      }
      if (method == "ext_sid") {
        eqApparent = ice.LocalToEquatorial_FromExtendedSid(local, t);
        sphereDir.phi = eqApparent.GetRaRad();
      }
      if (method == "solar") {
        eqApparent = ice.LocalToEquatorial(local, t);
        eqSolar = ice.PlanetToEquatorial(0, t);
        sphereDir.phi = eqApparent.GetRaRad() - eqSolar.GetRaRad();
      }
      else {
        eqApparent = ice.LocalToEquatorial(local, t);
        sphereDir.phi = eqApparent.GetRaRad();
      }

      sphereDir.theta = M_PI/2. - eqApparent.GetDecRad();
      if (sphereDir.phi < 0)
        sphereDir.phi += 2.*M_PI;

      pixelId = BGInt.ang2pix(sphereDir);
      //BGInt[pixelId] += 1.0;
      BGMap[pixelId] += 1.0;
    }
  }
  return BGMap;
}

I3_SERIALIZABLE(TimeScramble);

