
//#include <icetray/load_project.h>
#include <ShowerLLH/ShowerLLHFitParams.h>
#include <ShowerLLH/ShowerLLH.h>
#include <ShowerLLH/GridSearch.h>
#include <ShowerLLH/converter/ShowerLLHFitParamsConverter.h>
#include <tableio/converter/pybindings.h>

namespace bp = boost::python;

void
register_GridSearch()
{
  bp::class_<GridSearch, GridSearchPtr, 
                  bp::bases<I3FrameObject> >("GridSearch")
    .def("csearch", &GridSearch::csearch,
        (bp::arg("self"), bp::arg("gridx"), bp::arg("gridy"), 
        bp::arg("llhTable"),
        bp::arg("binDims"), bp::arg("binNames"), bp::arg("binnedVals"), 
        bp::arg("theta"), bp::arg("phi"),
        bp::arg("tankx"), bp::arg("tanky"), bp::arg("tankz"), 
        bp::arg("Dbins") ),
        "Perform grid search, returning most likely position and energy")
  ;
}

/*
void
register_ShowerLLH()
{
  bp::class_<ShowerLLH, bp::bases<I3FrameObject>, 
                  boost::shared_ptr<ShowerLLH> >("ShowerLLH")
    .def("ShowerLLH", &ShowerLLH::ShowerLLH,
      (bp::arg("self"),
      bp::arg("DirectionFit"),
      bp::arg("RecoPulses"),
      bp::arg("LLHTable"),
      bp::arg("eBins"),
      bp::arg("sBins"),
      bp::arg("dBins"),
      bp::arg("cBins"),
      bp::arg("gridSteps"),
      bp::arg("nside"),
      bp::arg("gridx"),
      bp::arg("gridy"),
      bp::arg("comp") ),
      "C++ implementation of ShowerLLH grid search")
  ;
}
*/

void
register_ShowerLLHFitParams()
{
	bp::class_<ShowerLLHFitParams, 
	    bp::bases<I3FrameObject>, ShowerLLHFitParamsPtr >("ShowerLLHFitParams")
	    //.def_readwrite("RecoRanCut", &ShowerLLHFitParams::RecoRanCut)
	    .def_readwrite("maxLLH", &ShowerLLHFitParams::maxLLH)
      /* Removed due to large file size */
      /*
	    .def_readwrite("coarse_llhs", &ShowerLLHFitParams::coarse_llhs)
	    .def_readwrite("med_llhs", &ShowerLLHFitParams::med_llhs)
	    .def_readwrite("fine_llhs", &ShowerLLHFitParams::fine_llhs)
      */
  ;

  register_pointer_conversions<ShowerLLHFitParams>();

}

void
register_tableio_converters()
{
	I3CONVERTER_NAMESPACE(ShowerLLH);
	I3CONVERTER_EXPORT_DEFAULT(ShowerLLHFitParamsConverter,
	    "Dumps the fit statistics from ShowerLLH reconstructions");
}


