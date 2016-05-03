#ifndef SHOWERLLHFITPARAMS_H
#define SHOWERLLHFITPARAMS_H

#include <icetray/I3ConditionalModule.h>
#include "icetray/I3FrameObject.h"

#include <string>
#include <vector>

#include <boost/shared_ptr.hpp>
#include <boost/python.hpp>

//class I3FrameObject;

class ShowerLLHFitParams : public I3FrameObject {
	public:
		//bool RecoRanCut;
		double maxLLH;
    /* Removed due to large file size */
    /*
    std::vector<double> coarse_llhs;
    std::vector<double> med_llhs;
    std::vector<double> fine_llhs;
    */

		ShowerLLHFitParams() :
			//RecoRanCut(0),
			maxLLH(0.0)
      /*
      coarse_llhs(0.0),
      med_llhs(0.0),
      fine_llhs(0.0),
      */
			{  };

		virtual ~ShowerLLHFitParams() {  };
	protected:
		friend class boost::serialization::access;
		template <class Archive> void serialize(Archive & ar,
		    unsigned version);
};

I3_POINTER_TYPEDEFS(ShowerLLHFitParams);


#endif //SHOWERLLHFITPARAMS_H

