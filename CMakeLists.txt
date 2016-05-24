i3_project(ShowerLLH
	DOCS_DIR resources/docs
	PYTHON_DIR python
)

i3_add_library(ShowerLLH
	private/ShowerLLH/ShowerLLH.cxx
	private/ShowerLLH/GridSearch.cxx
  private/ShowerLLH/ShowerLLHFitParams.cxx
	private/ShowerLLH/converter/ShowerLLHFitParamsConverter.cxx

	USE_TOOLS boost python sla gsl
	USE_PROJECTS icetray dataclasses tableio healpix-cxx
)

i3_add_pybindings(ShowerLLH
  private/pybindings/module.cxx
  private/pybindings/ShowerLLHFitParams.cxx
  USE_TOOLS boost python sla gsl
  USE_PROJECTS icetray dataclasses tableio healpix-cxx ShowerLLH
)

