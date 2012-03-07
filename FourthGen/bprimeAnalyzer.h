#ifndef DATASET
#define DATASET
struct Dataset {
  Dataset(
    const bool&   isData=false,		  
    const double& csIn=0,
    const double& wtIn=0,
    const double& lumiIn=0,
    const long&   evtIn=0 
  ) {;} 
  bool   isData ;  
  double cs     ; 
  double wt     ; 
  double lumi   ; 
  long   evts   ; 
} ; 
#endif 

#ifndef HISTS
#define HISTS

#include <TH1.h>
#include <TH2.h>
#include <TH3.h>

struct Hists{
  TH1F* hPUWeight_      ; 
  TH1F* hnPU_           ; 
  TH1F* hnGoodVtxs_     ; 
  TH1F* hlepId_         ; 
  TH1F* hnRecoEl_       ; 
  TH1F* hnPFEl_         ; 
  TH1F* h1stElPt_       ; 
  TH1F* h2ndElPt_       ; 
  TH1F* h1stElEta_      ; 
  TH1F* h2ndElEta_      ; 
  TH1F* hnZcands_       ; 
  TH1F* hZPt_           ; 
  TH1F* hZEta_          ; 
  TH1F* hZPhi_          ; 
  TH1F* hmee_           ; 
  TH1F* hnJets_         ; 
  TH1F* hallJetsPt_     ; 
  TH1F* hallJetsEta_    ; 
  TH1F* hnBjets_        ; 
  TH1F* hallBJetsPt_    ; 
  TH1F* hallBJetsEta_   ; 
  TH1F* hallBJetsPhi_   ; 
  TH1F* hleadingBJetPt_ ; 
  TH1F* hleadingBJetEta_; 
  TH1F* hleadingBJetPhi_; 
  TH1F* hnBprimes_      ; 
  TH1F* hmeeBjet_       ; 
  TH1F* heeBjetPt_      ; 
  TH1F* heeBjetEta_     ; 
  TH1F* heeBjetPhi_     ; 
};
#endif

#ifndef HISTOGRAMS
#define HISTOGRAMS
struct Histograms : public Hists { 

  Histograms() ;    
  virtual ~Histograms() {} ;
  virtual void fill() ;
	       
};
#endif 

#ifndef BTAGWEIGHTER
#define BTAGWEIGHTER

#include <vector>

class BTagWeighter { 
	private:
 static const double ptmin[14] ; 
 static const double ptmax[14] ;
 static const double SFb_error[14] ;

	public:
  BTagWeighter () {} ;
  ~BTagWeighter () {} ;

  /**\ B- C-jet SFs */ 
  double sfbc_TCHPM (TString meanminmax,double pt, double eta) { 
    if (pt >=30. && pt <= 670. ) return  0.616456*((1.+(0.145816*pt))/(1.+(0.0904067*pt))) ; 
    else if (pt > 670.) { 
      pt = 670. ;
      return 0.616456*((1.+(0.145816*pt))/(1.+(0.0904067*pt))) ; 
    } else return 0; 
  }

 /**\ B-jet SF errors */  
 double sferrorb_TCHPM (TString meanminmax,double pt, double eta) {
   int ptbin(-1) ;
   double error(0.) ; 
   for (int ii=0;ii<14;++ii) { 
     if (pt>=ptmin[ii] && pt<ptmax[ii]) { ptbin = ii ; break ; }
     else if (pt>ptmax[13]) { ptbin = 13 ; break ; }
   } 
   if (ptbin>=0) error = SFb_error[ptbin] ; 
   else error = 0 ; 
   return error;
 }

 /**\ C-jet SF errors */  
 double sferrorc_TCHPM (TString meanminmax,double pt, double eta) {
   int ptbin(-1) ;
   double error(0.) ; 
   for (int ii=0;ii<14;++ii) { 
     if (pt>=ptmin[ii] && pt<ptmax[ii]) { ptbin = ii ; break ; }
     else if (pt>ptmax[13]) { ptbin = 13 ; break ; }
   } 
   if (ptbin>=0) error = 2.*SFb_error[ptbin] ; 
   else error = 0 ; 
   return error;
 }

 /**\ Light-jet SF and SF errors */  
 double sflight_TCHPM (TString meanminmax,Float_t x, Float_t eta) { 
 
   if( eta >= 0.0 && eta < 0.8)   {
   	if( meanminmax == "mean" ) return ((1.27011+(-0.000869141*x))+(2.49796e-06*(x*x)))+(-2.62962e-09*(x*(x*x)));
   	if( meanminmax == "min" )  return ((1.12949+(-0.000678492*x))+(2.02219e-06*(x*x)))+(-2.21675e-09*(x*(x*x)));
   	if( meanminmax == "max" )  return ((1.41077+(-0.00105992*x))+(2.97373e-06*(x*x)))+(-3.0425e-09*(x*(x*x)));
   } else if( eta >= 0.8 && eta < 1.6)   {
   	if( meanminmax == "mean" ) return ((1.36167+(-0.00153237*x))+(4.54567e-06*(x*x)))+(-4.38874e-09*(x*(x*x)));
   	if( meanminmax == "min" )  return ((1.21289+(-0.00126411*x))+(3.81676e-06*(x*x)))+(-3.75847e-09*(x*(x*x)));
   	if( meanminmax == "max" )  return ((1.51053+(-0.00180085*x))+(5.27457e-06*(x*x)))+(-5.01901e-09*(x*(x*x)));
   }  else if( eta >= 1.6 && eta < 2.4)   {
   	if( meanminmax == "mean" ) return ((1.22696+(0.000249231*x))+(9.55279e-08*(x*x)))+(-1.04034e-09*(x*(x*x)));
   	if( meanminmax == "min" )  return ((1.07572+(0.00055366*x))+(-9.55796e-07*(x*x)))+(-3.73943e-11*(x*(x*x)));
   	if( meanminmax == "max" )  return ((1.3782+(-5.52498e-05*x))+(1.14685e-06*(x*x)))+(-2.04329e-09*(x*(x*x)));
   } else return 0; 
 
   return 0; 
 
 } ; 

};

const double BTagWeighter::ptmin[14] = {30., 40., 50., 60., 70., 80., 100., 120., 160., 210., 260., 320., 400., 500.}; 
const double BTagWeighter::ptmax[14] = {40., 50., 60., 70., 80.,100., 120., 160., 210., 260., 320., 400., 500., 670.};

const double BTagWeighter::SFb_error[14] = { 
 0.0365776,
 0.036307,
 0.0261062,
 0.0270308,
 0.0276016,
 0.0175067,
 0.0179022,
 0.0198104,
 0.0197836,
 0.024912,
 0.0273767,
 0.0398119,
 0.0418751,
 0.0605975 
}; 

#endif

#ifndef BPRIMEANALYZER_00_00_00_H
#define BPRIMEANALYZER_00_00_00_H

#include <TFile.h>
#include <TTree.h> 
#include <TChain.h>

#include <string>

/**\ \class bprimeAalyser 
 *
 * Main class for implementing event selection in bprime analysis
 * Default constructor: inputs
 *                    - name of input text file specifying names of rootfiles containing ntuples
 *                    - whether isData
 *                    - input cross section for porcesses (MC only)
 *                    - integrated lumi for dataset (MC only)
 *                    - name of outfile 
 * 
 * */

class bprimeAnalyzer {

	public:

  /**\ Default constructor */
  bprimeAnalyzer (const std::string& infiles=0,const bool& isData=0,const double& crosssection=0,const double& weight=0,const double& lumiInt=0,const double& events=0,const std::string&outfile=0); 
  /**\ Default destructor */
  virtual ~bprimeAnalyzer (); 
  /** Process TChain */ 
  virtual void process_() ; 

	protected: 

  /**\ Initialize histograms */
  virtual void initHists_(Histograms*&, std::string="") ; 
  /**\ Call to fill histograms */
  virtual void fillHists_(Histograms*)  ;
  /**\ Call to normalize histograms */
  virtual void normalizeHists_(Histograms*) ; 
  /**\ Set HLT bits */
  virtual void setHLTBits_(); 
  /**\ Event loop function */
  virtual void evtLoop_ () ;
  /**\ Returns true for electrons pasing electron ID */
  virtual bool isGoodElectron_ (int&) ; 
  /**\ Returns true for muons passing muon ID */
  virtual bool isGoodMuon_ (int&) ; 

	private:

  std::string inputfiles_ ;
  TFile outfile_ ;
  bool isData_ ; 
  Dataset dataset_ ;
  TChain* thisChain_ ; 

  std::vector<int>hltBits_ ; 

  Histograms* histBeforeSelection_ ; 
  Histograms* histElSelection_ ; 
  Histograms* histZCand_ ; 
  Histograms* histZPtGt95_ ; 
  Histograms* histBprimeCand_ ; 

  TTree* bprimeEvtsTree_ ; 

  /**\ Event variables  */
  unsigned int evtProc_ ; 
  unsigned int evtAccept_ ; 
  unsigned int nzcands_ ; 
  unsigned int nzPtGt95_; 
  unsigned int nbjets_ ; 
  unsigned int nBprimeCands_ ; 
  std::vector<int> el1Index_ ; 
  std::vector<int> el2Index_ ; 
  std::vector<int> jetIndex_ ; 
  std::vector<int> bjetIndex_ ; 
  
};

#endif 
