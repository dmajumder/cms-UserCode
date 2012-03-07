#define REWEIGHTMC_CXX
#include "reweightMC.h"
#include "HEADERDIR/format.h" 
#include "BTagWeight.h" 

#include <iostream>
#include <cmath>
#include <TLorentzVector.h> 

using namespace std;

reweightMC :: reweightMC (const std::string& infiles,const bool& isData,const double& crosssection, const double& weight,const double& lumiInt,const double& events,const std::string&outfile) 
	: bprimeAnalyzer (infiles,isData,crosssection,weight,lumiInt,events,outfile)
{
  initHists_(histBprimeCandWithBTagSF_,"BprimeCandWithBTagSF");  
}

reweightMC :: ~reweightMC () { 
}

void reweightMC :: evtLoop_ () {

  GenInfoBranches    GenInfo ; 
  EvtInfoBranches    EvtInfo ; 
  VertexInfoBranches VtxInfo ; 
  LepInfoBranches    LepInfo ; 
  JetInfoBranches    JetInfo ; 

  GenInfo.Register(thisChain_) ; 
  EvtInfo.Register(thisChain_) ; 
  VtxInfo.Register(thisChain_) ; 
  LepInfo.Register(thisChain_,"PFLepInfo") ; 
  JetInfo.Register(thisChain_,"PFJetInfo") ; 

  /**\ Looping over events */ 
  for (int iEntry=0;iEntry<thisChain_->GetEntries();++iEntry) { 

    thisChain_->GetEntry(iEntry) ; 

    /**\ Calculate b-tagging weight */ 

  }

  normalizeHists_(histBprimeCandWithBTagSF_) ; 

}
