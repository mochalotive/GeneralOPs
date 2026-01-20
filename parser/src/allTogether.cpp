#include "parser/include/basictype.h"
//#include "parser/include/basicconverter.h"
//#include "parser/include/converters.h"
#include "parser/include/parseinput.h"
//#include "parser/include/primativeconverters.h"
#include <iostream>


//int main(){
//int a=5;
//std::string b="hello!";
//double c = 10.01;
//BasicType d(a,b,c);
//d.Print();
////basic type declaration
////just to test
//
////void pointer init
//std::string temp[] = {"5","hello!","10.01"};//data from cmd
//std::vector<std::string> data(temp, temp + sizeof(temp)/sizeof(temp[0]));
/////////
//
//fp func(Converters::ReturnConverter("BASICTYPE"));//function
//
////output for data
//std::vector<BasicType> myParsedData;
//void * sendToConverter = &myParsedData;
//(*fp(Converters::ReturnConverter("BASICTYPE")))(data, sendToConverter);
//
////BasicConverter::converter(data, sendToConverter);
//myParsedData[0].Print();
//return(0);
//}


//int main() {
//std::vector<fp> v;
//v.push_back(&testfunc);
//int a=10;
//void *b = &a;
//std::string name="name";
//v[0](name, b);
//std::cout << a << std::endl;
//return(0);
//
//}

int main(int argc, const char* argv[]){
  if (argc==1) {std::cout << "Nothin to do, no arguments given." << std::endl; throw std::exception();}
  std::vector<ArgumentPair> argp; //an argument pair is a typedef, for an std::pair containing a string, and a vector of strings. This holds your cmd data, something like --numbers 1 2 3
  if (argc>1) argp = ParsingHelperNS::parseArgvToPair(argc, argv);//simplified to just take argc and argv}

  //now define all expected arguments
  //ExpectedArgument - a class which contains --term as a string, expected number of inputs as an int, a bool to determine if mandatory, a string to indicate data type, and a string to contain error  
  const ExpectedInput eargpt[] =        {ExpectedInput("files", 3, true, "STRING","Error: Expecting Three File Names in --files"),
                                         ExpectedInput("voronoi", 1, true, "BOOL", "Error: on or off,  expected for --voronoi"), 
                                         ExpectedInput("side-length", 1, false, "REAL", "Error: Expecting one Real number for side length of cube."), 
                                         ExpectedInput("numcells", 3, true, "INT", "Error: Expecting 3 integers for number of cells for order parameters."),
                                         ExpectedInput("test", 3, true, "BASICTYPE", "Error: Expecting an int, string, and double in that order!")};
  //will be replacing this soon
  std::vector<ExpectedInput> eargp(eargpt, eargpt+sizeof(eargpt)/sizeof(eargpt[0]));
 
  //finally input defaults for non mandatory arguments
  //just input a pseudo argv, and do the same as argp
//  std::cout << "break!!" << std::endl;
  std::vector<ArgumentPair> dargp;
  const char * dargpt[] = {"--side-length", "10"};
  //std::vector<std::string> dargpt2(dargpt1, dargpt1 + sizeof(dargpt1)/sizeof(dargpt1[0]));
  dargp = ParsingHelperNS::parseArgvToPair(2, dargpt, false );//not my favorite, not sure how c++ actually makes argc, just pass the real number of args, the array, and false
//  std::cout << dargp[0].second.size() << std::endl;
//  ParsedInput parsedInputs(&argp, &dargp, &eargp); 
  ParsedInput inputs(&argp, &dargp, &eargp);
  std::vector<std::string> myParsedData;
  inputs.WriteOut("files", myParsedData, "STRING");
//    std::cout << Converters::IsRegistered("INT") << std::endl;
  for (size_t i=0; i<myParsedData.size(); i++){std::cout << myParsedData[i] << std::endl;}
 return(0); 
}
