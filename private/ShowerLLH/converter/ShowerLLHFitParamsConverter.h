/**
 * copyright  (C) 2012
 * The Icecube Collaboration
 *
 * $Id: MillipedeFitParamsConverter.h 76611 2011-06-14 19:30:18Z jvansanten $
 *
 * @version $Revision: 76611 $
 * @date $LastChangedDate: 2011-06-14 14:30:18 -0500 (Tue, 14 Jun 2011) $
 * @author Frank McNally <fmcnally@wisc.edu> $LastChangedBy: jvansanten $
 */

#include <tableio/I3Converter.h>
#include <ShowerLLH/ShowerLLHFitParams.h>

class ShowerLLHFitParamsConverter : public I3ConverterImplementation<ShowerLLHFitParams > {
private:
    I3TableRowDescriptionPtr CreateDescription(const ShowerLLHFitParams & params); 
    size_t FillRows(const ShowerLLHFitParams& params, I3TableRowPtr rows);
};
    
