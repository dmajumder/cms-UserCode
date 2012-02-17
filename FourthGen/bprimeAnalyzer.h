#ifndef bprimeAnalyzer_00_00_00_h
#define bprimeAnalyzer_00_00_00_h

#include <TFolder.h>
#include <TFile.h>
#include <TString.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TCut.h>
#include <TChain.h>

#include <string>

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

struct Hists{
  TH1F* hPUWeight_      ; 
  TH1F* hnPU_           ; 
  TH1F* hnGoodVtxs_     ; 
  TH1F* hlepId_         ; 
  TH1F* hnRecoEl_       ; 
  TH1F* hnPFEl_         ; 
  TH2F* hptPFElRecoEl_  ; 
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

struct Histograms : public Hists { 

  Histograms() ;    
  virtual ~Histograms() {} ;
  virtual void fill() ;
	       
};

////class EvtInfoBranches ; 
////class LepInfoBranches ; 
////class JetInfoBranches ; 
////class PairInfoBranches ; 
////class PhotonInfoBranches ; 
////class VertexInfoBranches ; 
////class GenInfoBranches ; 

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
  ~bprimeAnalyzer ();
  /** Process TChain */ 
  virtual void process_() ; 

	protected: 

  /**\ Initialize histograms */
  virtual void initHists_(Histograms*&, std::string="") ; 
  /**\ Call to fill histograms */
  virtual void fillHists_(Histograms*)  ;
  /**\ Call to normalize histograms */
  virtual void normalizeHists_(Histograms*) ; 
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

  Histograms* histBeforeSelection_ ; 
  Histograms* histElSelection_ ; 
  Histograms* histZCand_ ; 
  Histograms* histZPtGt95_ ; 
  Histograms* histBprimeCand_ ; 

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
