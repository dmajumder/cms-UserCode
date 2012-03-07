#define BTAGWEIGHT_CXX 
#include "BTagWeight.h"

#include <iostream>

int fact(int n) {
  if(n < 1) return 1;
  int r=1;
  for(int i=n;i > 1; i--) r*=i;
  
  return r;
}

int comb(int n, int k) {
  return fact(n)/fact(k)/fact(n-k);
}

BTagWeight::BTagWeight(const std::vector<JetVariables>& JetsIn, const bool& isDataIn,const double& discIn,const int& minTagsIn,const int& maxTagsIn,const int& nbjetsIn) :
  isData(isDataIn), 
  disc(discIn), 
  Jets(JetsIn),
  maxTags(minTagsIn), 
  minTags(maxTagsIn), 
  nbjets(nbjetsIn) 
{
}

const double BTagWeight::ptmin[14] = {30., 40., 50., 60., 70., 80., 100., 120., 160., 210., 260., 320., 400., 500.}; 
const double BTagWeight::ptmax[14] = {40., 50., 60., 70., 80.,100., 120., 160., 210., 260., 320., 400., 500., 670.};

const double BTagWeight::SFb_error[14] = { 
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

bool BTagWeight::filter(const int& t) {
  return (t >= minTags && t <= maxTags);
}

double BTagWeight::getSF(int& flavour, double& pt, double& eta, TString& string) {
  double sf(0.)	; 
  if (flavour==5 || flavour==4) {
    sf = sfbc_TCHPM (string,pt,eta) ; 
  } else if (flavour<4) {
    sf = sflight_TCHPM(string,pt,eta) ; 
  } else sf = 0. ; 
  return sf ; 
}

double BTagWeight::getEff(bool& isData,int& flavour, double& pt, double& eta, TString& string,double& disc) {
  double eff(0.) ; 
  if (flavour==5 || flavour==4) {
    eff = btagEff_bc(isData,string,flavour,disc) ; 
  } else if (flavour<4) {
    eff = btagEff_light(isData,string,pt,eta) ; 
  } else eff = 0 ; 
  return eff ; 
} 

#include <TString.h> 

double BTagWeight::weight(const int& tags=1) { 

  //if(!filter(tags)) {
  //    return 0;
  //}

  for (unsigned iJet=0;iJet<Jets.size();++iJet) {
    TString mean = "mean" ; 
    Jets[iJet].sf = getSF(Jets[iJet].partonFlavour,Jets[iJet].ptj,Jets[iJet].etaj,mean); 
    Jets[iJet].eff = getEff(isData,Jets[iJet].partonFlavour,Jets[iJet].ptj,Jets[iJet].etaj,mean,disc);  
  }

  /**\ Event weight for b-taggig */
  double weight(0.) ; 
  if (isData) { weight = 1. ; } else {
    if (Jets.size()==1 && nbjets==1) { weight = Jets[0].sf ; std::cout << "Comment c'est possible!!\n" ; }
    if (Jets.size()==2 && nbjets==1) {
      int nbtagged(0) ; 
      double sfb, sfNotb ;
      double eb, eNotb ; 
      for(unsigned iJet=0;iJet<Jets.size();++iJet) {
        if (Jets[iJet].isBtag==true) { 
          sfb = Jets[iJet].sf ; 
	  ++nbtagged; 
	} else if (Jets[iJet].isBtag==false) { 
	  sfNotb = Jets[iJet].sf ;  
	  eNotb = Jets[iJet].eff ;
	}
      }
      weight = sfb*(1-eNotb*sfNotb)/sfNotb ; 
    }
  }
  return weight ; 

  ////int njets=b+c+l;

  ////double pMC=0;
  ////double pData=0;
  ////for(int ib=0;ib<=b;ib++)
  ////  for(int ic=0;ic<=c;ic++)
  ////    for(int il=0;il<=l;il++) {
  ////        int t=ib+ic+il;
  ////        if(!filter(t)) continue;
  ////     
  ////        // how many equivalent ways 
  ////        int totalCombinations=comb(b,ib)*comb(c,ic)*comb(l,il);
  ////    
  ////        //       std::cout <<  totalCombinations << std::endl;
  ////        pMC+=1.*totalCombinations * pow(meanBeff,ib)*pow(1.-meanBeff,b-ib)*pow(meanCeff,ic)*pow(1.-meanCeff,c-ic)*pow(meanUDSGeff,il)*pow(1.-meanUDSGeff,l-il);
  ////        pData+=1.*totalCombinations *  pow(meanBeff*sfb,ib)*pow(1.-meanBeff*sfb,b-ib)*pow(meanCeff*sfc,ic)*pow(1.-meanCeff*sfc,c-ic)*pow(meanUDSGeff*sfl,il)*pow(1.-meanUDSGeff*sfl,l-il);
  ////        std::cout << totalCombinations << " " << pow(meanBeff,ib)*pow(1.-meanBeff,b-ib)*pow(meanCeff,ic)*pow(1.-meanCeff,c-ic)*pow(meanUDSGeff,il)*pow(1.-meanUDSGeff,l-il) << " ";   
  ////        std::cout <<  pow(meanBeff*sfb,ib)*pow(1.-meanBeff*sfb,b-ib)*pow(meanCeff*sfc,ic)*pow(1.-meanCeff*sfc,c-ic)*pow(meanUDSGeff*sfl,il)*pow(1.-meanUDSGeff*sfl,l-il) << std::endl ;
  ////        /* std::cout << pow(meanBeff*sfb,ib)*pow(1.-meanBeff*sfb,b-ib)*pow(meanCeff*sfc,ic)*pow(1.-meanCeff*sfc,c-ic)*pow(meanUDSGeff*sfl,il)*pow(1.-meanUDSGeff*sfl,l-il)/(pow(meanBeff,ib)*pow(1.-meanBeff,b-ib)*pow(meanCeff,ic)*pow(1.-meanCeff,c-ic)*pow(meanUDSGeff,il)*pow(1.-meanUDSGeff,l-il)) << std::endl; 
  ////    std::cout << pMC << " " << pData <<  " " <<  pData/pMC << std::endl; 
  ////        */   
  ////   } 
  ////if(pMC==0) return 0; 
  ////return pData/pMC;

  return 0 ;

}


//int main() {
//
//  BTagWeight b(1,2,0.6,0.08,0.01);
//
//  for(int ib=2;ib<3;ib++)
//    for(int ic=0;ic<1;ic++)
//      for(int il=0;il<1;il++) {
//	  /*      
//      double w= b.weight(ib,ic,il,0.9,0.9,1.1,1);
//      double wbcp= b.weight(ib,ic,il,1,1,1.1,1);
//      double wbcm= b.weight(ib,ic,il,0.8,0.8,1.1,1);
//      double wlp= b.weight(ib,ic,il,0.9,0.9,1.2,1);
//      double wlm= b.weight(ib,ic,il,0.9,0.9,1.0,1);
//      double err=sqrt((wbcp-wbcm)*(wbcp-wbcm)+(wlp-wlm)*(wlp-wlm));
//      std::cout << ib << "   " << ic << "  " << "  " << il << "  " << b.weight(ib,ic,il,0.9,1,1,1) << " +- " <<  err << std::endl; 
//	  */ 
//	  std::cout << ib << "   " << ic << "  " << "  " << il << "  " << b.weight(ib,ic,il,0.9,1,1,2)  << std::endl; 
//
//      } 
//
//  BTagWeight bloose(2,100,0.8,0.15,0.08);   // at least 2 loose
//  BTagWeight btigth(1,100,0.5,0.07,0.005); // at least 1 tight
//
//  for(int ib=2;ib<3;ib++)
//    for(int ic=0;ic<1;ic++)
//      for(int il=0;il<2;il++) 	{
//	  double w1= btigth.weight(ib,ic,il,0.9,0.9,1.1,1);
//	  double w2= bloose.weight(ib,ic,il,0.95,0.95,1.2,2);
//	  std::cout << ib << "   " << ic << "  " << "  " << il << "  " << w1*w2 << std::endl;
//	}
//
//  return 0; 
//
//}

