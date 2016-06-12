#include <icetray/serialization.h>
#include <ShowerLLH/ShowerLLH.h>

#include <string>
#include <vector>
#include <cmath>
#include <limits>

#include "icetray/I3Context.h"
#include "icetray/I3Frame.h"

#include <dataclasses/I3Double.h>
#include <dataclasses/physics/I3Particle.h>
#include <dataclasses/physics/I3RecoPulse.h>
#include <dataclasses/geometry/I3Geometry.h>
#include <dataclasses/status/I3DetectorStatus.h>
#include <dataclasses/calibration/I3Calibration.h>
#include <dataclasses/status/I3DOMStatus.h>


// Constants for saturation
const double SAT_LG = 90000.;   // Upper bound for LG in PE
const double SAT_HG = 3000.;    // ...for HG (same as topeventcleaning)

I3_MODULE(ShowerLLH);

//----------------------------------------------------------------------------
ShowerLLH::ShowerLLH(const I3Context& ctx) : 
    I3ConditionalModule(ctx) {

  AddOutBox("OutBox");

  // Default values
  directionFit_ = "ShowerPlane";
  recoPulses_ = "IceTopHLCVEMPulses";
  composition_ = "proton";

  // Default values
  vector<int> gridSteps_, nside_;
  gridSteps_.push_back(20);
  gridSteps_.push_back(5);
  nside_.push_back(6);
  nside_.push_back(6);

  // Options
  AddParameter("DirectionFit", "I3Particle track to get direction from", 
      directionFit_);
  AddParameter("RecoPulses", "Pulse series used for reconstruction",
      recoPulses_);
  AddParameter("LLHTable", "Input LLH table", 0);
  AddParameter("eBins", "Bin edges in energy", 0);
  //AddParameter("zBins", "Bin edges in zenith", 0);
  AddParameter("sBins", "Bin edges in snow depth", 0);
  AddParameter("dBins", "Bin edges in distance", 0);
  AddParameter("cBins", "Bin edges in charge", 0);
  AddParameter("gridSteps", "Step sizes (m) to use for iterative search", 
      gridSteps_);
  AddParameter("nside", "Number of spots per side to use for iterative search", 
      nside_);
  AddParameter("gridx", "Input grid for initial ShowerLLH search", 0);
  AddParameter("gridy", "Input grid for initial ShowerLLH search", 0);
  AddParameter("comp", "Composition used to make likelihood tables",
      composition_);
}

void ShowerLLH::Configure() {

  GetParameter("DirectionFit", directionFit_);
  GetParameter("RecoPulses", recoPulses_);
  GetParameter("LLHTable", llhTable_);
  GetParameter("eBins", ebins_);
  //GetParameter("zBins", zbins_);
  GetParameter("sBins", sbins_);
  GetParameter("dBins", dbins_);
  GetParameter("cBins", cbins_);
  GetParameter("gridSteps", gridSteps_);
  GetParameter("nside", nside_);
  GetParameter("gridx", gridx_);
  GetParameter("gridy", gridy_);
  GetParameter("comp", composition_);
}

void ShowerLLH::Geometry(I3FramePtr frame) {
  if (!frame->Has("I3Geometry"))
    log_fatal("No Geometry in Frame!");
  PushFrame(frame, "OutBox");
}  

void ShowerLLH::Calibration(I3FramePtr frame) {
  if (!frame->Has("I3Calibration"))
    log_fatal("No Calibration in Frame!");
  PushFrame(frame, "OutBox");
}

void ShowerLLH::DetectorStatus(I3FramePtr frame) {
  if (!frame->Has("I3DetectorStatus"))
    log_fatal("No DetectorStatus in Frame!");
  PushFrame(frame, "OutBox");
}

//----------------------------------------------------------------------------
void ShowerLLH::Physics(I3FramePtr frame) {

  // Initial cuts
  // NOT SHOWN - check for minimum number of tanks?
  // Is this a good place for cuts???

  if (!frame->Has(directionFit_) || !frame->Has(recoPulses_)){
    PushFrame(frame, "OutBox");
    return;
  }

  const I3DetectorStatus &status = frame->Get<I3DetectorStatus>();
  const I3Geometry &geometry = frame->Get<I3Geometry>();
  const I3Calibration &calibration = frame->Get<I3Calibration>();
  const I3OMGeoMap &om_map = geometry.omgeo;
  const std::map<OMKey, I3DOMStatus> &status_map = status.domStatus;
  const std::map<OMKey, I3VEMCalibration> &vemcal_map = calibration.vemCal;
  I3StationGeoMap smap = geometry.stationgeo;

  I3ParticleConstPtr directionFit = 
      frame->Get<I3ParticleConstPtr>(directionFit_);
  I3RecoPulseSeriesMapConstPtr recoPulses = 
      frame->Get<I3RecoPulseSeriesMapConstPtr>(recoPulses_);

  if (directionFit->GetFitStatus() != I3Particle::OK) {
    PushFrame(frame, "OutBox");
    return;
  }

  double inf = numeric_limits<double>::infinity();
  double pi = M_PI;

  int nE = ebins_.size() - 1;
  //int nZ = zbins_.size() - 1;
  int nS = sbins_.size() - 1;
  int nD = dbins_.size() - 1;
  int nC = cbins_.size() - 1;
  double cmax = 1.01 * cbins_[nC-1];

  // Setup vectors for storage of tank parameters
  std::set<vector<int> > str_tnk;
  vector<int> om_vec(2);
  int tankID;
  vector<double> tankx, tanky, tankz, snowheight, vemcharges;
  double vem;

  // Loop over the pulse series to get list of (string,tank) pairs
  std::map<OMKey, I3RecoPulseSeries>::const_iterator p_om;
  for (p_om = recoPulses->begin(); p_om != recoPulses->end(); p_om++) {
    OMKey dom_key = p_om->first;
    TankKey tk(dom_key);
    tankID = tk.tank==TankKey::TankA?0:1;
    om_vec[0] = dom_key.GetString();
    om_vec[1] = tankID;
    str_tnk.insert(om_vec);
  }

  // Loop over the DOMs in the current detector setup
  std::map<OMKey, I3DOMStatus>::const_iterator i_om;
  for (i_om = status_map.begin(); i_om != status_map.end(); i_om++) {

    // Get and check the OMGeo
    OMKey dom_key = i_om->first;
    I3OMGeo om = om_map.find(dom_key)->second;
    if (om.omtype != I3OMGeo::IceTop)
      continue;

    // Setup vector for storing string and tankID
    TankKey tk(dom_key);
    tankID = tk.tank==TankKey::TankA?0:1;
    om_vec[0] = dom_key.GetString();
    om_vec[1] = tankID;

    // BAD STATIONS CHECK?  CHECK DOMS ARE ON??

    // Get and check its pulse series
    I3Map<OMKey, I3RecoPulseSeries>::const_iterator iter = 
        recoPulses->find(dom_key);

    // Ignore untriggered doms in taken (string,tank)
    if ((str_tnk.find(om_vec)!=str_tnk.end()) && (iter==recoPulses->end()))
      continue;

    // If dom_key not found in pulse series, fill with 0
    I3RecoPulseSeries pulses;
    vem = -1.;
    if (iter == recoPulses->end())
      vem = 0.;
    // If in pulse map but no pulses found, no information(??)
    else {
      pulses = iter->second;
      if (pulses.size() == 0) {
        //vem = 0.;
        continue;
      }
    }

    // Extract charge from pulse information
    if (vem == -1) {

      // Could write as a loop over pulses, only look at first for now...
      I3RecoPulseSeries::const_iterator it_pulse;
      it_pulse = pulses.begin();
      vem = it_pulse->GetCharge();

      // Determine saturation values
      I3VEMCalibration vemCalib = vemcal_map.find(dom_key)->second;
      double pe_per_vem = vemCalib.pePerVEM / vemCalib.corrFactor;
      double lg_sat = SAT_LG / pe_per_vem;
      I3DOMStatus::DOMGain gain = i_om->second.domGainType;

      // If saturated and HG, ignore
      if ((vem != vem) && (gain == I3DOMStatus::High))
        continue;
      // If saturated and LG, push back max charge
      if ((vem == vem) && (vem > lg_sat) && (gain == I3DOMStatus::Low))
        vem = cmax;
    }

    // Fill vectors
    vemcharges.push_back(vem);

    I3StationGeoMap::const_iterator siter = smap.find(dom_key.GetString());
    if (siter==smap.end())
      log_fatal("Station %d not in StationGeoMap!", dom_key.GetString());
    const I3TankGeo &tankGeo = siter->second.at(tankID);
    double h = tankGeo.snowheight;
    double x = (om.position).GetX();
    double y = (om.position).GetY();
    double z = (om.position).GetZ();
    tankx.push_back(x);
    tanky.push_back(y);
    tankz.push_back(z);
    snowheight.push_back(h);

    // Add om to "claimed" set (automatically erases duplicates)
    str_tnk.insert(om_vec);

  }

  int ntanks = tankx.size();
  const double theta = directionFit->GetZenith();
  const double phi = directionFit->GetAzimuth();

  // Bin elements
  //int bzenith = 0;
  int j = 0;
  //while (theta > zbins_[j+1])
  //  j += 1;
  //bzenith = j;
  vector<int> bsnow(ntanks);
  vector<int> bcharge(ntanks);
  for (int i = 0; i < ntanks; ++i) {
    j = 0;
    while (snowheight[i] > sbins_[j+1])
      j += 1;
    bsnow[i] = j;
    j = 0;
    while (vemcharges[i] > cbins_[j+1])
      j += 1;
    bcharge[i] = j;
  }

  vector<double> temp_gridx, temp_gridy;
  vector<int> bdist(ntanks);
  double xmax, ymax;
  double x, y;
  temp_gridx = gridx_;
  temp_gridy = gridy_;
  vector< vector<double> > maxLLHs, argmax;
  //double inf = numeric_limits<double>::infinity();
  int emax, gmax;
  double maxLLH;

  // Starting parameters
  double Z = 1947;
  double ex = -sin(theta) * cos(phi);
  double ey = -sin(theta) * sin(phi);
  double ez = -cos(theta);

  // Iterate through each set of grids
  for (unsigned int g = 0; g <= gridSteps_.size(); ++g) {
    int ngrid = temp_gridx.size();
    maxLLH = -inf;
    gmax = 0;
    emax = 0;

    // Iterate over all grid positions
    for (int g_idx = 0; g_idx < ngrid; ++g_idx) {
      double X = temp_gridx[g_idx];
      double Y = temp_gridy[g_idx];

      // Calculate closest approach distances
      for (int t_idx=0; t_idx < ntanks; ++t_idx) {
        double hx = X - tankx[t_idx];
        double hy = Y - tanky[t_idx];
        double hz = Z - tankz[t_idx];
        double s = ex*hx + ey*hy + ez*hz;
        double x1 = tankx[t_idx] + s*ex;
        double y1 = tanky[t_idx] + s*ey;
        double z1 = tankz[t_idx] + s*ez;
        double dist = sqrt(pow(x1-X,2) + pow(y1-Y,2) + pow(z1-Z,2));
        int dbin = 0;
        while (dist > dbins_[dbin+1])
          dbin += 1;
        bdist[t_idx] = dbin;
      }

      int llh_idx;
      // Calculate llh for every energy bin, return most likely
      for (int e_idx = 0; e_idx < nE; ++e_idx) {
        double llh = 0.0;
        for (int t = 0; t < ntanks; ++t) {
          llh_idx = 0;
          llh_idx += bcharge[t];
          llh_idx += nC * bdist[t];
          llh_idx += nC * nD * bsnow[t];
          llh_idx += nC * nD * nS * e_idx;
          //llh_idx += nC * nD * nS * bzenith;
          //llh_idx += nC * nD * nS * nZ * e_idx;
          llh += llhTable_[llh_idx];
        }
        if (llh > maxLLH) {
          maxLLH = llh;
          gmax = g_idx;
          emax = e_idx;
        }
      }
    }

    xmax = temp_gridx[gmax];
    ymax = temp_gridy[gmax];

    // Create new hexagon of search points centered on most likely position
    if (g < gridSteps_.size()) {
      temp_gridx.clear();
      temp_gridy.clear();
      for (int i = 0; i < nside_[g]; ++i) {
        for (int j = 0; j < (2*nside_[g]-i-1); ++j) {
          x = (-(2*nside_[g]-i-2)*gridSteps_[g])/2.0 + j*gridSteps_[g];
          y = (sqrt(3)*i*gridSteps_[g])/2.0;
          temp_gridx.push_back(x + xmax);
          temp_gridy.push_back(y + ymax);
          if (y != 0) {
            temp_gridx.push_back(x + xmax);
            temp_gridy.push_back(ymax - y);
          }
        }
      }
    }
  }


  I3ParticlePtr track (new I3Particle());
  track->SetShape(I3Particle::InfiniteTrack);
  track->SetFitStatus(I3Particle::OK);
  track->SetPos(xmax, ymax, Z);
  track->SetDir(theta, phi);

  // Calculate energy from bin
  double ereco = pow(10, (ebins_[emax] + ebins_[emax+1]) / 2.0);
  track->SetEnergy(ereco);
  I3DoublePtr maxLLH_out(new I3Double(maxLLH));

  frame->Put("ShowerLLH_"+composition_, track);
  std::cout << "Running..." << std::endl;
  ShowerLLHFitParamsPtr params(new ShowerLLHFitParams());
  params->maxLLH = double(maxLLH);
  frame->Put("ShowerLLHParams_"+composition_, params);
  PushFrame(frame, "OutBox");
}


void ShowerLLH::Finish() {
  log_info("ShowerLLH finished");
}













