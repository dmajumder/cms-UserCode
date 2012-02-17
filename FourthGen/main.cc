#include <iostream>
#include "bprimeAnalyzer.h"

int main () {

  time_t start, stop ; 
  double time_elapsed ; 
  time(&start) ; 

  std::string ip ; 
  std::string op ; 
  double cs(0.) ;
  double wt(0.) ; 
  bool isData(0) ;
  double lumi(0) ; 
  double events(0);

  std::cout << " enter input\n" ;
  std::cin >> ip ;
  std::cout << " enter output\n" ; 
  std::cin >> op ; 
  std::cout << " enter isData (=1 for data =0 for MC)\n" ;
  std::cin >> isData ;
  if (!isData) {
    std::cout << " enter sample cross-section\n" ;
    std::cin >> cs ; 
    std::cout << " enter sample weight\n" ;
    std::cin >> wt ; 
    std::cout << " enter integrated luminosity \n" ;
    std::cin >> lumi ; 
    std::cout << " enter total number of events\n" ;
    std::cin >> events;
  }

  //bprimeAnalyzer* bprimeanalyzer = new bprimeAnalyzer(ip,isData,cs,op) ; 
  bprimeAnalyzer* bprimeanalyzer = new bprimeAnalyzer(ip,isData,cs,wt,lumi,events,op) ; 
  bprimeanalyzer->process_() ; 
  delete bprimeanalyzer; 
  //bprimeanalyzer->~bprimeAnalyzer() ; 
  
  time(&stop) ; 
  time_elapsed = difftime(stop, start) ; 
  std::cout << "\n Time taken for program " << time_elapsed << " seconds \n" ;

  return 0 ; 

}


