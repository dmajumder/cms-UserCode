#ifndef CommHistProducer_h
#define CommHistProducer_h

/**\ 
 *
 * Class to fill in histograms of quantities for different jet flavours (MC) and data
 * Histograms booked in file opened by TFileService 
 *
 */ 

#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TString.h>
#include <string>
#include <map>
#include <vector> 
#include <cmath>
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

class CommHistProducer {

	public: 

		CommHistProducer(edm::Service<TFileService>&, const bool isDataIn) ; 
		virtual ~CommHistProducer() ; 

		virtual void     AddHisto(TString name, TString title,  int nbins, float min, float max);
		virtual void     AddHisto2D(TString name,TString title,int nbins,float min,float max,int nbins2,float min2,float max2);   
		virtual void     FillHisto_int(int flavour, bool isGS, int number, int value, float weight);
		virtual void     FillHisto_float(int flavour, bool isGS, int number, float value, float weight);
		virtual void     FillHisto_floatFromMap(TString name, int flavour, bool isGS, float value, float weight);
		virtual void     FillHisto_intFromMap(TString name, int flavour, bool isGS, int value, float weight);
		virtual void     FillHisto2D_int_floatFromMap(TString name, int flavour, bool isGS, int value, float value2, float weight);
		virtual void     FillHisto2D_float_floatFromMap(TString name, int flavour, bool isGS, float value, float value2, float weight); 

		edm::Service<TFileService> fs ; 
		std::vector<TH1D*>   HistoBtag;  
		std::vector<TH2D*>   HistoBtag2D;  

		std::map<TString, int>   HistoBtag_map; 
		std::map<TString, int>   HistoBtag2D_map; 

		const bool isData; 
		double x_section[10]; 
		double nmc_evt_vect[10];
		float WeightXS;
		float sum_xs;
		int choice;

		int numb_histo;
		int numb_histo2D;   
		int numb_histo2;

}; 

#endif 

#ifdef CommHistProducer_h

using namespace std ; 

CommHistProducer::CommHistProducer (edm::Service<TFileService>& fsIn, const bool isDataIn=false) : 

	        fs(fsIn), 
		isData(isDataIn),  
		numb_histo  (0),
		numb_histo2D(0),   
		numb_histo2 (0) 

{


}

CommHistProducer::~CommHistProducer() {
}

void CommHistProducer::AddHisto(TString name, TString title, int nbins, float min, float max)  {

	TH1D* h_b      = fs->make<TH1D>(name+"_b",title+"_b",nbins,min,max);
	TH1D* h_bfromg = fs->make<TH1D>(name+"_bfromg",title+"_bfromg",nbins,min,max);  
	TH1D* h_c      = fs->make<TH1D>(name+"_c",title+"_c",nbins,min,max);  
	TH1D* h_l      = fs->make<TH1D>(name+"_l",title+"_l",nbins,min,max);
	TH1D* h_data   = fs->make<TH1D>(name+"_data",title+"_data",nbins,min,max);

	h_b        ->Sumw2();
	h_bfromg   ->Sumw2();  
	h_c        ->Sumw2();  
	h_l        ->Sumw2(); 
	h_data     ->Sumw2();

	HistoBtag.push_back(h_b);
	HistoBtag.push_back(h_bfromg);  
	HistoBtag.push_back(h_c);  
	HistoBtag.push_back(h_l);  
	HistoBtag.push_back(h_data);  
	HistoBtag_map[name] = numb_histo;

	numb_histo++;

}

void CommHistProducer::FillHisto_float(int flavour, bool isGS, int number, float value, float weight)  {

	if (!isData){
		if (fabs(flavour)==5 && !isGS)                  HistoBtag[number*5 +0]->Fill(value,weight);
		else if (fabs(flavour)==5 && isGS)              HistoBtag[number*5 +1]->Fill(value,weight);
		else if (fabs(flavour)==4)                      HistoBtag[number*5 +2]->Fill(value,weight); 
		else if (fabs(flavour)< 4 || fabs(flavour)==21) HistoBtag[number*5 +3]->Fill(value,weight);

	}  
	else                                              HistoBtag[number*5 +4]->Fill(value);


}
void CommHistProducer::FillHisto_int(int flavour, bool isGS, int number, int value, float weight)  {

	if (!isData){
		if (fabs(flavour)==5 && !isGS)             HistoBtag[number*5 +0]->Fill(value,weight);
		else if (fabs(flavour)==5 && isGS)              HistoBtag[number*5 +1]->Fill(value,weight);
		else if (fabs(flavour)==4)                      HistoBtag[number*5 +2]->Fill(value,weight); 
		else if (fabs(flavour)< 4 || fabs(flavour)==21) HistoBtag[number*5 +3]->Fill(value,weight);
	}  
	else                                         HistoBtag[number*5 +4]->Fill(value);

}


void CommHistProducer::FillHisto_floatFromMap(TString name, int flavour, bool isGS, float value, float weight)  {


	int number = HistoBtag_map[name] ;
	if (!isData){
		if (fabs(flavour)==5 && !isGS)                  HistoBtag[number*5 +0]->Fill(value,weight);
		else if (fabs(flavour)==5 && isGS)              HistoBtag[number*5 +1]->Fill(value,weight);
		else if (fabs(flavour)==4)                      HistoBtag[number*5 +2]->Fill(value,weight); 
		else if (fabs(flavour)< 4 || fabs(flavour)==21) HistoBtag[number*5 +3]->Fill(value,weight);

	}  
	else                                              HistoBtag[number*5 +4]->Fill(value);


}

void CommHistProducer::FillHisto_intFromMap(TString name, int flavour, bool isGS, int value, float weight)  {

	int number = HistoBtag_map[name] ;
	if (!isData){
		if (fabs(flavour)==5 && !isGS)                  HistoBtag[number*5 +0]->Fill(value,weight);
		else if (fabs(flavour)==5 && isGS)              HistoBtag[number*5 +1]->Fill(value,weight);
		else if (fabs(flavour)==4)                      HistoBtag[number*5 +2]->Fill(value,weight); 
		else if (fabs(flavour)< 4 || fabs(flavour)==21) HistoBtag[number*5 +3]->Fill(value,weight);
	}  
	else                                         HistoBtag[number*5 +4]->Fill(value);

}

void CommHistProducer::AddHisto2D(TString name, TString title, int nbins, float min, float max, int nbins2, float
		min2, float max2)  {

	TH2D* h_b      = fs->make<TH2D>(name+"_b",title+"_b",nbins,min,max,nbins2,min2,max2);
	TH2D* h_bfromg = fs->make<TH2D>(name+"_bfromg",title+"_bfromg",nbins,min,max,nbins2,min2,max2);  
	TH2D* h_c      = fs->make<TH2D>(name+"_c",title+"_c",nbins,min,max,nbins2,min2,max2);  
	TH2D* h_l      = fs->make<TH2D>(name+"_l",title+"_l",nbins,min,max,nbins2,min2,max2);
	TH2D* h_data   = fs->make<TH2D>(name+"_data",title+"_data",nbins,min,max,nbins2,min2,max2);


	h_b        ->Sumw2();
	h_bfromg   ->Sumw2();  
	h_c        ->Sumw2();  
	h_l        ->Sumw2(); 
	h_data     ->Sumw2();

	HistoBtag2D.push_back(h_b);
	HistoBtag2D.push_back(h_bfromg);  
	HistoBtag2D.push_back(h_c);  
	HistoBtag2D.push_back(h_l);  
	HistoBtag2D.push_back(h_data);  
	HistoBtag2D_map[name] = numb_histo2D;
	numb_histo2D++;

}

void CommHistProducer::FillHisto2D_int_floatFromMap(TString name, int flavour, bool isGS, int value, float value2, float weight)  {


	int number = HistoBtag2D_map[name] ;
	if (!isData){
		if (fabs(flavour)==5 && !isGS)                  HistoBtag2D[number*5 +0]->Fill(value,value2,weight);
		else if (fabs(flavour)==5 && isGS)              HistoBtag2D[number*5 +1]->Fill(value,value2,weight);
		else if (fabs(flavour)==4)                      HistoBtag2D[number*5 +2]->Fill(value,value2,weight); 
		else if (fabs(flavour)< 4 || fabs(flavour)==21) HistoBtag2D[number*5 +3]->Fill(value,value2,weight);

	}  
	else                                              HistoBtag2D[number*5 +4]->Fill(value,value2);


}

void CommHistProducer::FillHisto2D_float_floatFromMap(TString name, int flavour, bool isGS, float value, float value2, float weight)  {


	int number = HistoBtag2D_map[name] ;
	if (!isData){
		if (fabs(flavour)==5 && !isGS)                  HistoBtag2D[number*5 +0]->Fill(value,value2,weight);
		else if (fabs(flavour)==5 && isGS)              HistoBtag2D[number*5 +1]->Fill(value,value2,weight);
		else if (fabs(flavour)==4)                      HistoBtag2D[number*5 +2]->Fill(value,value2,weight); 
		else if (fabs(flavour)< 4 || fabs(flavour)==21) HistoBtag2D[number*5 +3]->Fill(value,value2,weight);

	}  
	else                                              HistoBtag2D[number*5 +4]->Fill(value,value2);


}

#endif 
