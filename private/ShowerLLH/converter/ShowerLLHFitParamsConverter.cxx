/**
 * copyright  (C) 2012
 * The Icecube Collaboration
 *
 * $Id: ShowerLLHFitParamsConverter.cxx 84012 2012-01-17 17:58:12Z nwhitehorn $
 *
 * @version $Revision: 84012 $
 * @date $LastChangedDate: 2012-01-17 11:58:12 -0600 (Tue, 17 Jan 2012) $
 * @author Frank McNally  <fmcnally@wisc.edu> $LastChangedBy: nwhitehorn $
 */

#include "ShowerLLHFitParamsConverter.h"

I3TableRowDescriptionPtr ShowerLLHFitParamsConverter::CreateDescription(const ShowerLLHFitParams& params) {
    I3TableRowDescriptionPtr desc(new I3TableRowDescription());

    //desc->AddField<bool>("RecoRanCut", "", "successful reconstruction");
    desc->AddField<double>("maxLLH", "", "max likelihood for composition");

    /* Likelihood descriptors: removed due to large file size */
    /*
    desc->AddField<double>("coarse_llhs", "", "LLH space for coarse grid",145);
    desc->AddField<double>("med_llhs", "", "LLH space for med grid", 91);
    desc->AddField<double>("fine_llhs", "", "LLH space for fine grid", 91);
    */

    return desc;
}
    
size_t ShowerLLHFitParamsConverter::FillRows(const ShowerLLHFitParams& params,
        I3TableRowPtr rows) {

    //rows->Set<bool>("RecoRanCut",  params.RecoRanCut);
    rows->Set<double>("maxLLH",  params.maxLLH);

    /* Fill likelihood descriptors */
    /*
    double *coarse_llhs = rows->GetPointer<double>("coarse_llhs");
    double *med_llhs = rows->GetPointer<double>("med_llhs");
    double *fine_llhs = rows->GetPointer<double>("fine_llhs");
    std::copy(params.coarse_llhs.begin(),params.coarse_llhs.end(),coarse_llhs);
    std::copy(params.med_llhs.begin(), params.med_llhs.end(), med_llhs);
    std::copy(params.fine_llhs.begin(), params.fine_llhs.end(), fine_llhs);
    */
    return 1;
}

