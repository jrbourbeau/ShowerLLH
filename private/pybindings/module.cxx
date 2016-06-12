
#include <icetray/load_project.h>
#include <ShowerLLH/ShowerLLHFitParams.h>
#include <ShowerLLH/GridSearch.h>
#include <ShowerLLH/ShowerLLH.h>
#include <ShowerLLH/converter/ShowerLLHFitParamsConverter.h>
#include <tableio/converter/pybindings.h>

namespace bp = boost::python;

void register_GridSearch();
//void register_ShowerLLH();
void register_ShowerLLHFitParams();
void register_tableio_converters();


I3_PYTHON_MODULE(ShowerLLH) {

  load_project("ShowerLLH", false);

  register_GridSearch();
  //register_ShowerLLH();
	register_ShowerLLHFitParams();
	register_tableio_converters();
}


