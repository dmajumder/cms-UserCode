#include <iostream>
#include <string>
#include "MultileptonBprime/OS2LSelector/interface/OS2LSelector.h" 

using namespace std ;

int main (int argc, char *argv[]) {

	time_t start, stop ;
	double time_elapsed ;
	time(&start) ;

	std::string ip ;
	std::string op ;
	bool isData(0) ;
	std::string jsonfile;

	std::cout << " enter input\n" ;
	std::cin >> ip ;
	std::cout << " enter output\n" ;
	std::cin >> op ;
	std::cout << " enter isData (=1 for data =0 for MC)\n" ;
	std::cin >> isData ;
	std::cout << " enter json file if data\n" ;
	std::cin >> jsonfile ;

	OS2LSelector* selector = new OS2LSelector(ip,isData,jsonfile,op) ;
	selector->process() ;
	delete selector;

	time(&stop) ;
	time_elapsed = difftime(stop, start) ;
	std::cout << "\n Time taken for program " << time_elapsed << " seconds \n" ;

	return 0 ;

}

