#ifndef BOOKHISTOGRAMS_H
#define BOOKHISTOGRAMS_H

#include <TH1D.h>
#include <iostream>

class BookHistograms {

	public: 

  BookHistograms() ;
  ~BookHistograms() {} ;  
  void fill(TH1D* hist, double bin, double wt=1.0) ; 

  void setSumw2() ; 

  TH1D* m_hnvtx_ ; 

  TH1D* m_hpte0_ ; 
  TH1D* m_hpte1_ ; 
  TH1D* m_hptm0_ ; 
  TH1D* m_hptm1_ ; 

  TH1D* m_hpte0_dy_ ; 
  TH1D* m_hpte1_dy_ ; 
  TH1D* m_hptm0_dy_ ; 
  TH1D* m_hptm1_dy_ ; 

  TH1D* m_hpte0_notdy_ ; 
  TH1D* m_hpte1_notdy_ ; 
  TH1D* m_hptm0_notdy_ ; 
  TH1D* m_hptm1_notdy_ ; 

  TH1D* m_hetae0_ ; 
  TH1D* m_hetae1_ ; 
  TH1D* m_hetam0_ ; 
  TH1D* m_hetam1_ ; 

  TH1D* m_hetae0_dy_ ; 
  TH1D* m_hetae1_dy_ ; 
  TH1D* m_hetam0_dy_ ; 
  TH1D* m_hetam1_dy_ ; 

  TH1D* m_hetae0_notdy_ ; 
  TH1D* m_hetae1_notdy_ ; 
  TH1D* m_hetam0_notdy_ ; 
  TH1D* m_hetam1_notdy_ ; 

  TH1D* m_hmee_ ; 
  TH1D* m_hmmumu_ ; 
  TH1D* m_hmemu_ ; 

  TH1D* m_hnjets_ ; 
  TH1D* m_hnbjets_ ; 

  TH1D* m_hnjets_dy_ ; 
  TH1D* m_hnbjets_dy_ ; 

  TH1D* m_hnjets_notdy_ ; 
  TH1D* m_hnbjets_notdy_ ; 

  TH1D* m_hmet_ ; 
  TH1D* m_hst_ ; 

  TH1D* m_hmet_dy_ ; 
  TH1D* m_hst_dy_ ; 

  TH1D* m_hmet_notdy_ ; 
  TH1D* m_hst_notdy_ ; 

}; 


BookHistograms::BookHistograms () {

  int bins ;
  double x0, xn ; 

  std::cout << " Histograms will be booked\n" ; 

  bins = 1000 ; 
  x0 = 0. ; 
  xn = 1000. ; 

  m_hpte0_ = new TH1D("hpte0","",bins,x0,xn) ;  
  m_hpte1_ = new TH1D("hpte1","",bins,x0,xn) ;  
  m_hptm0_ = new TH1D("hptm0","",bins,x0,xn) ;  
  m_hptm1_ = new TH1D("hptm1","",bins,x0,xn) ;  

  m_hpte0_dy_ = new TH1D("hpte0_dy","",bins,x0,xn) ;  
  m_hpte1_dy_ = new TH1D("hpte1_dy","",bins,x0,xn) ;  
  m_hptm0_dy_ = new TH1D("hptm0_dy","",bins,x0,xn) ;  
  m_hptm1_dy_ = new TH1D("hptm1_dy","",bins,x0,xn) ;  

  m_hpte0_notdy_ = new TH1D("hpte0_notdy","",bins,x0,xn) ;  
  m_hpte1_notdy_ = new TH1D("hpte1_notdy","",bins,x0,xn) ;  
  m_hptm0_notdy_ = new TH1D("hptm0_notdy","",bins,x0,xn) ;  
  m_hptm1_notdy_ = new TH1D("hptm1_notdy","",bins,x0,xn) ;  

  bins = 400 ; 
  x0 = -4. ; 
  xn = 4. ; 

  m_hetae0_ = new TH1D("hetae0","",bins,x0,xn) ;  
  m_hetae1_ = new TH1D("hetae1","",bins,x0,xn) ;  
  m_hetam0_ = new TH1D("hetam0","",bins,x0,xn) ;  
  m_hetam1_ = new TH1D("hetam1","",bins,x0,xn) ;  

  m_hetae0_dy_ = new TH1D("hetae0_dy","",bins,x0,xn) ;  
  m_hetae1_dy_ = new TH1D("hetae1_dy","",bins,x0,xn) ;  
  m_hetam0_dy_ = new TH1D("hetam0_dy","",bins,x0,xn) ;  
  m_hetam1_dy_ = new TH1D("hetam1_dy","",bins,x0,xn) ;  

  m_hetae0_notdy_ = new TH1D("hetae0_notdy","",bins,x0,xn) ;  
  m_hetae1_notdy_ = new TH1D("hetae1_notdy","",bins,x0,xn) ;  
  m_hetam0_notdy_ = new TH1D("hetam0_notdy","",bins,x0,xn) ;  
  m_hetam1_notdy_ = new TH1D("hetam1_notdy","",bins,x0,xn) ;  

  bins = 1000 ; 
  x0 = 0. ; 
  xn = 1000. ; 

  m_hmee_ = new TH1D("hmee","",bins,x0,xn);  
  m_hmmumu_ = new TH1D("hmmumu","",bins,x0,xn) ;  
  m_hmemu_ = new TH1D("hmemu","",bins,x0,xn) ;  

  bins = 20 ; 
  x0 = 0. ; 
  xn = 20. ; 

  m_hnvtx_ = new TH1D("hnvtx","",bins,x0,xn) ;  
  m_hnjets_ = new TH1D("hnjets","",bins,x0,xn) ;  
  m_hnbjets_ = new TH1D("hnbjets","",bins,x0,xn) ;  

  m_hnjets_dy_ = new TH1D("hnjets_dy","",bins,x0,xn) ;  
  m_hnbjets_dy_ = new TH1D("hnbjets_dy","",bins,x0,xn) ;  

  m_hnjets_notdy_ = new TH1D("hnjets_notdy","",bins,x0,xn) ;  
  m_hnbjets_notdy_ = new TH1D("hnbjets_notdy","",bins,x0,xn) ;  

  bins = 1000 ; 
  x0 = 0. ; 
  xn = 1000. ; 

  m_hmet_ = new TH1D("hmet","",bins,x0,xn) ;  
  m_hst_ = new TH1D("hst","",bins,x0,xn) ;  

  m_hmet_dy_ = new TH1D("hmet_dy","",bins,x0,xn) ;  
  m_hst_dy_ = new TH1D("hst_dy","",bins,x0,xn) ;  

  m_hmet_notdy_ = new TH1D("hmet_notdy","",bins,x0,xn) ;  
  m_hst_notdy_ = new TH1D("hst_notdy","",bins,x0,xn) ;  

  std::cout << " Histograms have been booked: going to set Sumw2()\n" ; 

  setSumw2() ; 	

}

void BookHistograms::setSumw2() {

  m_hnvtx_ ->Sumw2() ;  

  m_hpte0_ ->Sumw2() ;  
  m_hpte1_ ->Sumw2() ;  
  m_hptm0_ ->Sumw2() ;  
  m_hptm1_ ->Sumw2() ;  

  m_hpte0_dy_ ->Sumw2() ;  
  m_hpte1_dy_ ->Sumw2() ;  
  m_hptm0_dy_ ->Sumw2() ;  
  m_hptm1_dy_ ->Sumw2() ;  

  m_hpte0_notdy_ ->Sumw2() ;  
  m_hpte1_notdy_ ->Sumw2() ;  
  m_hptm0_notdy_ ->Sumw2() ;  
  m_hptm1_notdy_ ->Sumw2() ;  

  m_hetae0_ ->Sumw2() ;  
  m_hetae1_ ->Sumw2() ;  
  m_hetam0_ ->Sumw2() ;  
  m_hetam1_ ->Sumw2() ;  

  m_hetae0_dy_ ->Sumw2() ;  
  m_hetae1_dy_ ->Sumw2() ;  
  m_hetam0_dy_ ->Sumw2() ;  
  m_hetam1_dy_ ->Sumw2() ;  

  m_hetae0_notdy_ ->Sumw2() ;  
  m_hetae1_notdy_ ->Sumw2() ;  
  m_hetam0_notdy_ ->Sumw2() ;  
  m_hetam1_notdy_ ->Sumw2() ;  

  m_hmee_ ->Sumw2() ; 
  m_hmmumu_ ->Sumw2() ;  
  m_hmemu_ ->Sumw2() ;  

  m_hnjets_ ->Sumw2() ;  
  m_hnbjets_ ->Sumw2() ;  

  m_hnjets_dy_ ->Sumw2() ;  
  m_hnbjets_dy_ ->Sumw2() ;  

  m_hnjets_notdy_ ->Sumw2() ;  
  m_hnbjets_notdy_ ->Sumw2() ;  

  m_hmet_ ->Sumw2() ;  
  m_hst_ ->Sumw2() ;  

  m_hmet_dy_ ->Sumw2() ;  
  m_hst_dy_ ->Sumw2() ;  

  m_hmet_notdy_ ->Sumw2() ;  
  m_hst_notdy_ ->Sumw2() ;  

  return ; 
	
}

void BookHistograms::fill(TH1D* hist, double bin, double wt) {

  hist->Fill(bin,wt) ; 
   
  return ; 

}

#endif 
