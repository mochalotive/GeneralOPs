#include <fstream>
#include <sstream>
#include <string>

#include "parser/include/parseinput.h" 
#include "converters/include/primativeconverters.h"

int main(int argc, const char* argv[]){
{

//EXPECTED ARGUMENTS
//quick reference
//ExpectedInput(std::string parameter_, int arguments_, bool mandatory_, std::string type_, std::string errorMessage_, bool strict_=true);
  const ExpectedInput eargpt[] ={//define what you'd like to find in the input, and whether or not it is mandatory
                                 ExpectedInput("plumedout", 1, true, "PLUMED_FORCE", "Plumed's output on a restraint force for --plumedout.", true),
  std::vector<ExpectedInput> eargp(eargpt, eargpt+sizeof(eargpt)/sizeof(eargpt[0]));//convert to vector
  
//PRINT HELP MESSAGE
  //print instructions
  if (argc==1) 
  {
  std::cout << "plumed2opd converts a plumed output to a binary opd file. \n"
            << "just feed it the output of plumed, printing out the restraint force as often as you want \n"
            << "plumed2opd will parse through it, and convert it to an opd file." << std::endl;
  //print all arguments
  for (size_t i=0; i<eargp.size(); i++)
    {
    eargp[i].PrintExpected();
    }
    throw std::exception();
  }
//PARSE CMD INPUT FOR COMPARISION TO EXPECTED ARGUMENTS
  std::vector<ArgumentPair> argp; //an argument pair is a typedef, for an std::pair containing a string, and a vector of strings. This holds your cmd data, something like --numbers 1 2 3
  if (argc>1) argp = ParsingHelperNS::parseArgvToPair(argc, argv);//simplified to just take argc and argv directly

//DEFINE DEFAULT ARGUMENTS FOR NON-MANDATORY INPUTS
//just pretend the dargpt is the command line
  std::vector<ArgumentPair> dargp;
  const char * dargpt[] = {"--plumedout", "restaint.dat"};
  size_t defSize = sizeof(dargpt)/sizeof(dargpt[0]);
  dargp = ParsingHelperNS::parseArgvToPair(defSize, dargpt, false );//not my favorite, not sure how c++ actually makes argc, just pass the real number of args, the array, and false

//finally, the parsed input class hold all this info, and contains all the machinery to be passed a vector of the type you're after
//and then populate it  
  ParsedInput inputs(&argp, &dargp, &eargp);



return 0;
}
