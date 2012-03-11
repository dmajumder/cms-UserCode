#ifndef REWEIGHTMC_H
#define REWEIGHTMC_H
#include "bprimeAnalyzer.h"

#include "HEADERDIR/format.h" 

class LepInfoBranches ; 
class JetInfoBranches ; 

class reweightMC : public bprimeAnalyzer { 

	public:
  reweightMC  (const std::string& infiles=0,const bool& isData=0,const double& crosssection=0,const double& weight=0,const double& lumiInt=0,const double& events=0,const std::string&outfile=0) ; 
  ~reweightMC  () ;  
  void process_() ; 

	private:
  TFile outfile_ ;

  /**\ Event loop function */
  void initHists_() ; 
  void evtLoop_ () ;
  inline bool jetAndBprimeSel_(LepInfoBranches*,JetInfoBranches*,int&,int&,int&,TString&) ; 

  std::map<TString,TH1F*> hmeeBjet_; 

}; 

#endif
