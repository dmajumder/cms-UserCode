#ifndef BTAGWEIGHT_H
#define BTAGWEIGHT_H

#include <cmath> 
#include <vector> 
#include <TString.h> 

struct JetVariables { 
  JetVariables(bool& isBtagIn,int& partonFlavourIn, double& ptjIn, double& etajIn): 
    isBtag(isBtagIn),
    partonFlavour(partonFlavourIn),
    ptj(ptjIn),
    etaj(etajIn) {} ; 
  bool isBtag; 
  int partonFlavour;
  double ptj;
  double etaj; 
  double eff;
  double effErrorP;
  double effErrorM;
  double sf; 
  double sfErrorP; 
  double sfErrorM; 
};

class BTagWeight {

 public:

  explicit BTagWeight( const std::vector<JetVariables>& JetsIn,const bool& isData=false, const double& discIn=1.93,const int& minTagsIn=0, const int& maxTagsIn=1,const int& nbjetsIn=1) ; 
  ~BTagWeight() {}; 
  double weight(const TString&); 

 private:

  bool isData ; 
  double disc ; 
  std::vector<JetVariables> Jets; 
  int maxTags;
  int minTags;
  int nbjets ; 

  bool filter(const int& t);
  double getSF(int& flavour, double& pt, double& eta, TString& string) ; 
  double getEff(int& flavour, double& pt, double& eta, TString& string,double& disc) ; 

  static const double ptmin[14] ; 
  static const double ptmax[14] ;
  static const double SFb_error[14] ;

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
  double sflight_TCHPM (TString meanminmax,double x, double eta) { 
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
  }  

  /**\ Btagging efficiency for b and c quarks */
  double btagEff_bc (TString& mode,int& flavour,double& disc) {
    double eff(0) ; 
    if (isData) {
      if (mode=="mean") eff = 9.83842428415e-06*disc*disc*disc*disc +  -0.000556835427293*disc*disc*disc +  0.0123891144567*disc*disc +  -0.141658673059*disc +  0.804455651041 ; 
      if (mode=="errorP" || mode=="errorM") eff = 1.23607953977e-05*disc*disc*disc*disc + -0.000639251438742*disc*disc*disc + 0.0131227817518*disc*disc + -0.14535620858*disc + 0.854260405647 ; 
    } 
    if (!isData) {
      if (/*mode=="mean" &&*/ TMath::Abs(flavour)==5) eff = 1.26119661124e-05*disc*disc*disc*disc +  -0.000683198597977*disc*disc*disc +  0.0145106168149*disc*disc +  -0.159575511553*disc +  0.887707865272 ; 
      if(/*mode=="mean" &&*/ TMath::Abs(flavour)==4) eff = 0.451288118581*exp(-0.0213290505241*disc*disc*disc + 0.356020789904*disc*disc + -2.20158883207*disc + 1.84838018633 ) ; 
    } 
    return eff ; 
  }
  
  /**\ Btagging efficiency for light quarks */ 
  double btagEff_light (TString& mode,double& pt,double& eta) {
    double eff(0.) ; 
    if (pt < 20.) {
      eff = 0;
      return eff ; 
    }
    else if (pt > 670.) pt = 670 ; 
    if(eta>=0.0 && eta<0.8) {
      /*if( mode == "mean" )*/ eff = (((-0.00464673+(0.000247485*pt))+(9.13236e-07*(pt*pt)))+(-2.49994e-09*(pt*(pt*pt))))+(1.65678e-12*(pt*(pt*(pt*pt))));
    }
    if(eta>=0.8 && eta<1.6) {
      /*if( mode == "mean" )*/ eff = (((-0.0060878+(0.000297422*pt))+(1.13369e-06*(pt*pt)))+(-2.84945e-09*(pt*(pt*pt))))+(1.64721e-12*(pt*(pt*(pt*pt))));
    }
    if(eta>=1.6 && eta<2.4) { 
      /*if( mode == "mean" )*/ eff = (((-0.00836219+(0.000391889*pt))+(2.78156e-07*(pt*pt)))+(-6.14017e-10*(pt*(pt*pt))))+(-1.30592e-13*(pt*(pt*(pt*pt))));
    }
    //if(eta>=0.0 && eta<2.4) {
    //  if( mode == "mean" ) eff = (((-0.00609327+(0.00029725*pt))+(9.46617e-07*(pt*pt)))+(-2.71065e-09*(pt*(pt*pt))))+(1.77968e-12*(pt*(pt*(pt*pt))));
    //}
    return eff ; 
  }

};

#endif

