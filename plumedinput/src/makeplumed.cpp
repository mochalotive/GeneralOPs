#include <fstream>
#include <vector>
#include <string>

#include "mymol/include/geometry.h"
#include "mymol/include/file_formats/fileformats.h"
#include "mymol/include/system.h"
#include "mymol/include/lattice.h"

#include "crystdist/include/crystalops.h"
#include "crystdist/include/moleculemap.h"
#include "crystdist/include/statparameters.h"
#include "crystdist/include/mmpfile.h"
#include "crystdist/include/xtpfile.h"
#include "crystdist/include/copfile.h"
#include "crystdist/include/relativeconfiguration.h"

//#include "converters/include/converters.h" 
#include "parser/include/parseinput.h" //comes bundled with both now
#include "converters/include/primativeconverters.h"
//#include "parser/include/primativeconverters.h"

#include "plumedinput/include/plumedfile.h"
#include "plumedinput/include/forcefile.h"
#include "plumedinput/include/kernelfile.h"
#include "plumedinput/include/plumedparameters.h"
#include "plumedinput/include/kernelparameters.h"
#include "plumedinput/include/fileconverters.h"
#include "plumedinput/include/equationcontainer.h"

#include <algorithm>

//TODO
//grab residues from mmp files instead
//make the help message a bit more formal, more automated
//a non-mandatory argument with a non-strict number of inputs could potentially be an issue if the number of arguments depends on other things
//for example COM should have the same number of elements as MMP and COP, so just be careful when relying on default arguments
//as it is right now I just make it the max number


  int main(int argc, const char* argv[]){
  //define expected arguments
  //if i were a smarty, I could just grab the residues from the mmp files, but I'll leave it manual for now.
  const ExpectedInput eargpt[] ={
                                 ExpectedInput("mmp", 5, true, "MOLECULE_MAPS", "Expecting Once MMP filepath for --mmp.", false),
                                 ExpectedInput("cop", 5, true, "CRYSTAL_ORDER_PARAMETERS", "Expecting one COP filepath for --cop", false),
                                 ExpectedInput("cop-old", 5, true, "CRYSTAL_ORDER_PARAMETERS", "Expecting one old COP filepath for --cop", false),
                                 ExpectedInput("xtp", 5, true, "CRYSTAL_DISTRIBUTION_PARAMETERS","Error: error, expecting one filepath for --xtp", false),
                                 ExpectedInput("pdb", 1 , true, "SYSTEM_GEOMETRY", "Expecting one pbb filepath for --pdb"),
                                 ExpectedInput("resname", 5 , true, "STRING", "Expecting a residue to be restrained for --resname", false),
                                 ExpectedInput("box-size", 3, false, "REAL", "Expecting the three side lengths for a 90 90 90 rectangle for --box-size"),
                                 ExpectedInput("COM" , 5, false, "BOOL", "Expecting on or off for each residue type", false),
                                 ExpectedInput("OP", 5, true, "STRING", "Expecting ROPS, BOPS, DOPS, LOPS, or NOPS for a given residue", false),
                                 ExpectedInput("cutoff", 5, true, "REAL", "Expecting a cutoff for the OPs", false),
                                 ExpectedInput("contact",5, false, "BOOL", "Speed up with MASK?", false),
                                 ExpectedInput("contact-cutoff", 5, false, "REAL", "Cutoff for MASK vector", false),
                                 ExpectedInput("eqf", 5, true, "EQUATION_CONTAINER", "Output of reweight.py", false),
                                 ExpectedInput("numstep", 1, true, "INT", "Number of steps", true),
                                 ExpectedInput("numdragstep", 1, true, "INT", "Number of drag steps", true),
                                 ExpectedInput("kappa", 5, true, "REAL", "Force Constants for restraints", false)
};
  std::vector<ExpectedInput> eargp(eargpt, eargpt+sizeof(eargpt)/sizeof(eargpt[0]));//convert to vector
  //if no arguments given, print output
  if (argc==1) 
  {
   //pre-message
    std::cout 
     << "\n\nmakeplumed produces kernel files (files containing means, concentration parameters) for the 3 main kinds of order parameters defined by Santiso and Trout.\n\n" 
     <<  "It also produces a PLUMED input file for the given system for use in restrained dynamics\n\n" 
     << "It needs exactly one pdb file, and then cop, old_cop, mmp, and xtp files for all residues in the pdb which whill be restrained.\n\n" 
     << "MAKE SURE the --resname flags are in the same order as the --cop, --mmp, and --xtp files!\n\n" 
     << "These files don't have any linking to the resnames, so make sure they're in the same order. eg --resname SFD ACN, make sure the rest are in the order of SFD and ACN.\n\n" 
     <<" Kernel Files - they will be produced according to the resnames given eg SFDDOPS.dat, SFDBOPS.dat, etc. It just makes all of them, so 3 per residue. They're not big.\n\n"
     << "ALSO ensure that the number of ops in the cop files match, and also match the numbers given here. I could add a check here, or just have it grab the box size from the cop file direcrly. I haven't done this yet.\n\n" << std::endl;
    //print all possible arguments, and whether they are mandatory and strict or not 
    for (size_t i=0; i<eargp.size(); i++)
    {
    eargp[i].PrintExpected();
    
    }
    //post messgae
    throw std::exception();

  }



  std::vector<ArgumentPair> argp; //an argument pair is a typedef, for an std::pair containing a string, and a vector of strings. This holds your cmd data, something like --numbers 1 2 3
  if (argc>1) argp = ParsingHelperNS::parseArgvToPair(argc, argv);//simplified to just take argc and argv}

//define expected inputs
//ExpectedArgument - a class which contains --term as a string, expected number of inputs as an int, a bool to determine if mandatory, a string to indicate data type, and a string to contain error  
  std::vector<ArgumentPair> dargp;
  const char * dargpt[] = {/*"--dops-kernel", "DOPS.dat", "--rops-kernel", "ROPS.dat", "--bops-kernel", "BOPS.dat",*/ "--box-size", "0", "0", "0", "--COM", "on", "on", "on", "on", "on", "--contact", "off", "off", "off", "off", "off", "--contact-cutoff", "0.0", "0.0", "0.0", "0.0", "0.0"};
  //std::vector<std::string> dargpt2(dargpt1, dargpt1 + sizeof(dargpt1)/sizeof(dargpt1[0]));
  size_t defSize = sizeof(dargpt)/sizeof(dargpt[0]);
  dargp = ParsingHelperNS::parseArgvToPair(defSize, dargpt, false );//not my favorite, not sure how c++ actually makes argc, just pass the real number of args, the array, and false

  ParsedInput inputs(&argp, &dargp, &eargp);

  //now we can write everything out via WriteOut

  std::vector<CrystalDistributionParameters> myCrystalDistributionParameterVector;//xtp data
  inputs.WriteOut("xtp", myCrystalDistributionParameterVector, "CRYSTAL_DISTRIBUTION_PARAMETERS");
  
  std::vector<std::vector<MoleculeMap> > myMoleculeMapsVector;//mmp data
  inputs.WriteOut("mmp", myMoleculeMapsVector, "MOLECULE_MAPS");

//std::cout << myMoleculeMapsVector[0][0].framePoints[0].atoms[0] << " " << myMoleculeMapsVector[0][0].framePoints[1].atoms[0] << " " << myMoleculeMapsVector[0][0].framePoints[2].atoms[0] << std::endl;

  
  std::vector<CrystalOrderParameters> myCrystalOrderParameterVector;//cop data
  inputs.WriteOut("cop", myCrystalOrderParameterVector, "CRYSTAL_ORDER_PARAMETERS");

  std::vector<CrystalOrderParameters> myOldCrystalOrderParameterVector;//cop data
  inputs.WriteOut("cop-old", myOldCrystalOrderParameterVector, "CRYSTAL_ORDER_PARAMETERS");
  
  std::vector<System<Geometry > > mySystemGeometryVector;//pdb data
  inputs.WriteOut("pdb", mySystemGeometryVector, "SYSTEM_GEOMETRY");


  std::vector<EquationContainer> myEquations;
  inputs.WriteOut("eqf", myEquations, "EQUATION_CONTAINER");

  std::vector<int> numstep;
  std::vector<int> numdragstep;
  inputs.WriteOut("numstep", numstep, "INT");
  inputs.WriteOut("numdragstep", numdragstep, "INT");
  if (numstep[0] < 1)
  {
    std::cout << "Need a positive number of steps!" << std::endl;
    throw std::exception();
  }
  if (numdragstep[0] > numstep[0])
  {
  std::cout << "numdragsteps must be less than or equal to numstep" << std::endl;
  throw std::exception();
  }
  if (numdragstep[0] > 0)
  if (numstep[0] % numdragstep[0] != 0 || numdragstep[0] % 5000 !=0 || numstep[0] % 5000 != 0)
  {
  std::cout << "try to be smart about drag steps, make sure they are cleanly divisible by 5000" << std::endl;
  throw std::exception();
  }

  std::pair<int,int> steps(numstep[0], numdragstep[0]);
  SystemDataPackage mySystemDataPackage(mySystemGeometryVector[0],
                    myCrystalDistributionParameterVector,
                    myMoleculeMapsVector,
                    myCrystalOrderParameterVector,
                    myOldCrystalOrderParameterVector,
                    myEquations,
                    steps);//not too happy with 5 arguments, but its ok for now.



  //now plumed parameters not from other files
  std::vector<Real> box_size;//how many boxes to separate
  inputs.WriteOut("box-size", box_size, "REAL");

  std::vector<bool> COM;//boolean, using center of masses or not?
  inputs.WriteOut("COM", COM, "BOOL");


  std::vector<std::string> residueNames;
  inputs.WriteOut("resname", residueNames, "STRING");
 

  //Jake 2025 added new params
  std::vector<Real> cutoff;
  std::vector<bool> contact;
  std::vector<Real> contact_cutoff;
  std::vector<Real> kappa; 
 
  inputs.WriteOut("cutoff", cutoff, "REAL");
  inputs.WriteOut("contact", contact, "BOOL");
  inputs.WriteOut("contact-cutoff", contact_cutoff, "REAL");
  inputs.WriteOut("kappa", kappa, "REAL");

  PlumedSystemParameters myPlumedSystemParameters(COM, box_size, residueNames, cutoff, contact, contact_cutoff, kappa);  


  //write out kernel files (gotta do it somewhere!)
  for (size_t i=0; i<residueNames.size(); i++)
  {
  KernelFile kernelDops(residueNames[i] + "DOPS");
  KernelFile kernelBops(residueNames[i] + "BOPS");
  KernelFile kernelRops(residueNames[i] + "ROPS");

  KernelParameters myKernelParametersDops(&myCrystalDistributionParameterVector[i], DOPS_KERNEL);
  KernelParameters myKernelParametersBops(&myCrystalDistributionParameterVector[i], BOPS_KERNEL);
  KernelParameters myKernelParametersRops(&myCrystalDistributionParameterVector[i], ROPS_KERNEL);
  
  kernelDops << (kernelDops, myKernelParametersDops); 
  kernelBops << (kernelBops, myKernelParametersBops); 
  kernelRops << (kernelRops, myKernelParametersRops); 
  }


PlumedFile myPlumedFile("myPlumedFile"); //no IN for plumed
std::vector<std::string> OPS;
inputs.WriteOut("OP", OPS, "STRING");
PlumedParameters myPlumedParameters(OPS, mySystemDataPackage, myPlumedSystemParameters );


myPlumedFile << (myPlumedFile, myPlumedParameters);

//print all the force files
for (size_t res=0; res < residueNames.size(); res++)
  for (size_t i=0; i < size_t(box_size[0]*box_size[1]*box_size[2]); ++i)
    {
    std::string myNumString = std::to_string(i);
    std::string forceFileName = residueNames[res] + myNumString;
    ForceFile temp(forceFileName);
    PlumedParameters myPlumedParameters(OPS, mySystemDataPackage, myPlumedSystemParameters, i, res );
    temp << (temp, myPlumedParameters);
    }
return 0;}

