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
  TH1F* hweightPU_      ; 
  TH1F* hweightBtagging_; 
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

class LepInfoBranches ; 

static const int nmaxBprimes =  15 ; 

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
  virtual bool isGoodElectron_ (LepInfoBranches*,int&,int&) ; 
  /**\ Returns true for muons passing muon ID */
  virtual bool isGoodMuon_ (LepInfoBranches*,int&,int&) ; 

  std::string inputfiles_ ;
  TString outfileName_ ; 
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
  Histograms* histBprimeCandBtagScaled_ ; 

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

  struct BPRIMECANDS {
    int _ele1Index[nmaxBprimes] ; 
    int _ele2Index[nmaxBprimes] ; 
    int _bjetIndex[nmaxBprimes] ; 
    double _bprimeMass[nmaxBprimes] ; 
    double _bprimePt[nmaxBprimes] ; 
    double _bprimeEta[nmaxBprimes] ; 
    double _bprimePhi[nmaxBprimes] ; 
    double _zMass[nmaxBprimes] ; 
    double _zPt[nmaxBprimes]  ; 
    double _zEta[nmaxBprimes] ; 
    double _zPhi[nmaxBprimes] ; 
    double _bjetPt[nmaxBprimes] ; 
    double _bjetEta[nmaxBprimes] ; 
    double _bjetPhi[nmaxBprimes] ; 
  }; 
  
};

#endif 
