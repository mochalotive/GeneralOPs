#include <fstream>
#include <vector>
#include <string>

#include "plumedinput/include/plumedfile.h"
#include "plumedinput/include/kernelfile.h"
#include "plumedinput/include/plumedparameters.h"
#include "plumedinput/include/kernelparameters.h"
#include "plumedinput/include/fileconverters.h"

#include "parser/include/converters.h" 
#include "parser/include/parseinput.h"
#include "parser/include/primativeconverters.h"

#include "crystdist/include/crystalops.h"
#include "crystdist/include/moleculemap.h"
#include "crystdist/include/statparameters.h"
#include "crystdist/include/mmpfile.h"
#include "crystdist/include/xtpfile.h"
#include "crystdist/include/copfile.h"

#include "mymol/include/geometry.h"
#include "mymol/include/file_formats/fileformats.h"
#include "mymol/include/system.h"
int main(int argc, const char* argv[]){
  
//std::vector<std::string> arguments;

//read in cmd inputs
//you can copy paste this to other files
  if (argc==1) {std::cout << "Nothin to do, no arguments given." << std::endl; throw std::exception();}
  std::vector<ArgumentPair> argp; //an argument pair is a typedef, for an std::pair containing a string, and a vector of strings. This holds your cmd data, something like --numbers 1 2 3
  if (argc>1) argp = ParsingHelperNS::parseArgvToPair(argc, argv);//simplified to just take argc and argv}

//define expected inputs
//ExpectedArgument - a class which contains --term as a string, expected number of inputs as an int, a bool to determine if mandatory, a string to indicate data type, and a string to contain error  
  const ExpectedInput eargpt[] ={ExpectedInput("xtp", 1, true, "CRYSTAL_DISTRIBUTION_PARAMETERS","Error: error, expecting one filepath for --xtp")};
  //will be replacing this soon to just pass [] array directly to constructor
  std::vector<ExpectedInput> eargp(eargpt, eargpt+sizeof(eargpt)/sizeof(eargpt[0]));

//no default arguments this time


std::string pdbFileName, mmpFileName, xtpFileName, copFileName, plumedFileName;


for(int i = 0; i < argc; ++i)
//  arguments.push_back(std::string(argv[i]));
xtpFileName=argv[1];
mmpFileName=argv[2];
copFileName=argv[3];
pdbFileName=argv[4];
plumedFileName=argv[5];
//part out names


//for now, just do xtp,mmp,cop,pdb in that order
//no input checking living on the wild side 
//xtp
CrystalDistributionParameters myCrystalParameters;
XTPFile myXTPFile(xtpFileName,IN); 
//mmp
std::vector<MoleculeMap> myMoleculeMaps;
MMPFile myMMPFile(mmpFileName); //not readwrite mode, only in
//cop
CrystalOrderParameters myParameters;
COPFile myCOPFile(copFileName, IN);
//pdb
System<Geometry> mySystem;
PDBFile myPDBFile(pdbFileName, IN);

//now read info into data structures

//xtp
myXTPFile >> myCrystalParameters;
  if((!myXTPFile && !myXTPFile.eof()) ||
     myCrystalParameters.concentrationParameters.size() < 1 ||
     myCrystalParameters.means.size() < 1 ||
     myCrystalParameters.normalizationFactors.size() < 1)
  {
    std::cerr << "Error reading crystal distribution parameters file " << xtpFileName << std::endl
              << "Use copcalc -help for help." << std::endl;
    return -1;
  }
//mmp
myMMPFile >> myMoleculeMaps;
  if((!myMMPFile && !myMMPFile.eof()) ||
     myMoleculeMaps.size() < 1 ||
     (myMoleculeMaps.size() > 0 && myMoleculeMaps[0].framePoints.size() < 1))
  {
    std::cerr << "Error reading molecule map file " << mmpFileName << std::endl
              << "Use copcalc -help for help." << std::endl;
    return -1;
  }

//cop
  myCOPFile >> myParameters;
//didn't check this one 8)
//pdb
myPDBFile >> mySystem;

  if((!myPDBFile && !myPDBFile.eof()) ||
     mySystem.size() < 1 ||
     (mySystem.size() > 0 && mySystem[0].numAtoms() == 0))
  {
    std::cerr << "Error reading geometry file " << pdbFileName << std::endl
              << "Use copcalc -help for help." << std::endl;
    return -1;
  }
//all read in, now << to write
//shoddy at best 8-)


PlumedFile myPlumedFile(plumedFileName); //no IN for plumed 
PlumedParameters myPlumedParameters(&myCrystalParameters,&myMoleculeMaps,&myParameters,&mySystem); //internal data types for xtp, mmp, cop, pdb in that order


//KernelParameters *myKernelParameters[3];
KernelFile kernel1("dopsKernel.dat"), kernel2("bopsKernel.dat"),kernel3("ropsKernel.dat");
KernelParameters myKernelParameters(&myCrystalParameters, DOPS_KERNEL);
KernelParameters myKernelParameters2(&myCrystalParameters, BOPS_KERNEL);
KernelParameters myKernelParameters3(&myCrystalParameters, ROPS_KERNEL);




myPlumedFile << (myPlumedFile, myPlumedParameters);
kernel1 << (kernel1, myKernelParameters); 
kernel2 << (kernel2, myKernelParameters2); 
kernel3 << (kernel3, myKernelParameters3);

return 0;}

