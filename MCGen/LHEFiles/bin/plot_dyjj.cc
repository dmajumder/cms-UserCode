#include <cmath>
#include <TFile.h>
#include <TH1D.h>
#include <TLorentzVector.h>

#include "LHEF.h"

#define pow2(a) pow(a,2.)

using namespace std;

void readerExample() {

  //// Open a stream connected to an event file:
  std::ifstream ifs1("/afs/cern.ch/user/d/devdatta/Madgraph/madgraph/proc_DYToJJ_2Partons/proc_DYToJJ_2Partons/Events/run_01/unweighted_events.lhe"); 

  //// Create the Reader object:
  LHEF::Reader reader1(ifs1);

  //// Print out the header information:
  std::cout << reader1.headerBlock;

  //// Print out the addinional comments in the init block:
  //std::cout << reader1.initComments;

  //// Print out the beam energies:
  std::cout << "Beam A: " << reader1.heprup.EBMUP.first << " GeV, Beam B: "
    << reader1.heprup.EBMUP.second << " GeV." << std::endl;

  // Now loop over all events:
  long ieve = 0;

  TFile *rootoutput = new TFile("dytojj.root","RECREATE");
  rootoutput->cd(); 
  TH1D* hmz   = new TH1D("hmz",  "M(Z)",         100, 0., 200.)  ; 
  TH1D* hptz  = new TH1D("hptz", "pT(Z)",        500, 0., 1000.) ; 
  TH1D* hmjj  = new TH1D("hmjj", "M(jj), Z>jj",  100, 0., 200.)  ; 
  TH1D* hptjj = new TH1D("hptjj","pT(jj), Z>jj", 500, 0., 1000.) ; 

  while ( reader1.readEvent() ) {

    ++ieve;

    // Some events may be specially tagged "# SPECIAL" in the comment
    // lines, in that case write out the number of particles in that
    // event:
    if ( reader1.eventComments.find("# SPECIAL") != std::string::npos ) 
      std::cout << "Event " << ieve << " contained a special event with "
        << reader1.hepeup.NUP << " particles." << std::endl;

    double pxz, pyz, pzz, ez;
    double pxjj, pyjj, pzjj, ejj;  
    for (int ii=0;ii< reader1.hepeup.NUP;ii++) {

      if(reader1.hepeup.IDUP[ii]==23){
        pxz = reader1.hepeup.PUP[ii][0];
        pyz = reader1.hepeup.PUP[ii][1];
        pzz = reader1.hepeup.PUP[ii][2];
        ez  = reader1.hepeup.PUP[ii][3];
      }
      if(abs(float(reader1.hepeup.IDUP[ii])) < 6 || abs(float(reader1.hepeup.IDUP[ii])) == 21) { 
        if ( abs(float(reader1.hepeup.MOTHUP[ii].first)) == 23 && abs(float(reader1.hepeup.MOTHUP[ii].second)) == 23) {   
          pxjj += reader1.hepeup.PUP[ii][0];
          pyjj += reader1.hepeup.PUP[ii][1];
          pzjj += reader1.hepeup.PUP[ii][2];
          ejj  += reader1.hepeup.PUP[ii][3];
        } 
      }
    } //// Particle loop 
    TLorentzVector p4z(pxz, pyz, pzz, ez) ; 
    TLorentzVector p4jj(pxjj, pyjj, pzjj, ejj) ; 
    hmz  -> Fill(p4z.Mag())  ; 
    hmjj -> Fill(p4jj.Mag()) ; 
    hptz -> Fill(p4z.Pt())   ;  
    hptjj-> Fill(p4jj.Pt())  ;  
  } //// Event loop 

  rootoutput->Write();
  rootoutput->Close(); 

  // Now we are done.

}

int main(){

  readerExample();

  return 0;

}

