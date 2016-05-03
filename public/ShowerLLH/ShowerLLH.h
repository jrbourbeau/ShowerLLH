#ifndef SHOWERLLH_H_INCLUDED
#define SHOWERLLH_H_INCLUDED

#include <icetray/I3ConditionalModule.h>
#include <ShowerLLH/ShowerLLHFitParams.h>

#include <string>
#include <vector>

using namespace::std;

class ShowerLLH : public I3ConditionalModule {

	public:

    // Constructor
    ShowerLLH(const I3Context& ctx);
    // Destructor
    ~ShowerLLH(){}

    void Configure();
    void Geometry(I3FramePtr frame);
    void Calibration(I3FramePtr frame);
    void DetectorStatus(I3FramePtr frame);
    void Physics(I3FramePtr frame);
    void Finish();

  private:

    ShowerLLH();
    ShowerLLH(const ShowerLLH &source);
    ShowerLLH& operator=(const ShowerLLH &source);

    string directionFit_;
    string recoPulses_;
    vector<double> llhTable_;
    vector<double> ebins_;
    //vector<double> zbins_;
    vector<double> sbins_;
    vector<double> dbins_;
    vector<double> cbins_;
    vector<int> gridSteps_;
    vector<int> nside_;
    vector<double> gridx_;
    vector<double> gridy_;
    string composition_;

};


#endif //SHOWERLLH_H_INCLUDED

