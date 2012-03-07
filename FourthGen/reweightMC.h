#ifndef REWEIGHTMC_H
#define REWEIGHTMC_H
#include "bprimeAnalyzer.h"

class reweightMC : public bprimeAnalyzer { 

	public:
  reweightMC  (const std::string& infiles=0,const bool& isData=0,const double& crosssection=0,const double& weight=0,const double& lumiInt=0,const double& events=0,const std::string&outfile=0) ; 
  ~reweightMC  () ;  

	private:

  /**\ Event loop function */
  virtual void evtLoop_ () ;

  TChain* thisChain_ ; 

  Histograms* histBprimeCandWithBTagSF_ ;

}; 

#endif
