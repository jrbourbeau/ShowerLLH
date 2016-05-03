#include <icetray/serialization.h>
#include "icetray/I3FrameObject.h"
#include <ShowerLLH/ShowerLLHFitParams.h>

using namespace std;

template <class Archive>
void ShowerLLHFitParams::serialize(Archive & ar, unsigned version){

	ar & make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
	//ar & make_nvp("RecoRanCut", RecoRanCut);
	ar & make_nvp("maxLLH", maxLLH);

  /* Removed due to large file size */
  /*
  ar & make_nvp("coarse_llhs", coarse_llhs);
  ar & make_nvp("med_llhs", med_llhs);
  ar & make_nvp("fine_llhs", fine_llhs);
  */
}

BOOST_CLASS_VERSION(ShowerLLHFitParams, 0);
I3_SERIALIZABLE(ShowerLLHFitParams);



