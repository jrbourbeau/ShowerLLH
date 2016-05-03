import numpy
import pickle
from icecube import icetray
from icecube import ShowerLLH
LLHTables = ShowerLLH.SupportFunctions.LLHTables

ShowerLLH_keys = ['pML_llhs', 'fML_llhs', 'pMLtrack', 'fMLtrack']
ShowerLLH_keys += [dict(key='passed', converter=ShowerLLH.converters.ShowerLLHFitParamsConverter())]

@icetray.traysegment
def GridLLHSegment(tray, name,
		pulses='IceTopVEMPulses_0',
		DirectionFit='ShowerPlane',
		If = lambda frame: True):
	resourcedir = '/home/fmcnally/ShowerLLH/resources/'
	gridList = [[],[],[]]
	gridList[0] = pickle.load(open(resourcedir + 'BinSpots_Coarse.pkl', 'rb'))
	gridList[1] = pickle.load(open(resourcedir + 'BinSpots_Middle.pkl', 'rb'))
	gridList[2] = pickle.load(open(resourcedir + 'BinSpots_Fine.pkl', 'rb'))
	LLHlist = LLHTables(resourcedir + 'FunctionHist.hdf5')
	EDGEbins = pickle.load(open(resourcedir + 'EDGEbins.pkl', 'rb'))
	EDGEbins2 = numpy.array(range(8)+[14,15,23,24,33,34,44,45,56,57,69,70,81,82,92,93,102,103,111,112]+range(119,127))
	EDGElist = [EDGEbins, EDGEbins2]
	tray.AddModule(ShowerLLH.GridLLH, 'ShowerLLHGrid', DirectionFit=DirectionFit, LLHTables=LLHlist, gridList=gridList, edgeList=EDGElist)


