#include "smcv/include/opdfile.h"
#include "parser/include/parseinput.h" //comes bundled with both now
#include "converters/include/primativeconverters.h"







int main(int argc, const char* argv[]){
//define expected arguments
//if i were a smarty, I could just grab the residues from the mmp files, but I'll leave it manual for now.
const ExpectedInput eargpt[] ={
                                 ExpectedInput("kappa", 5, false, "REAL", "Force constant used for MD", false),
                                 ExpectedInput("gradF", 1, true, "STRING", "formatted gradients of potential" ),
                                 ExpectedInput("gradPhi", 1, true, "STRING", "formatted gradients of OPs wrt atoms"),
//                                 ExpectedInput("num-atoms", 1 , true, "INT", "number of atoms"),
//                                 ExpectedInput("steps", 1, true, "INT", "number of steps to calculate"),
                                 ExpectedInput("num-skip", 1, true, "INT", "number of steps to skip")
//                                 ExpectedInput("numops", 1, true, "INT", "number of order parameters"),
 //                                ExpectedInput("numres", 1, true, "INT", "number of residues")

};
  std::vector<ExpectedInput> eargp(eargpt, eargpt+sizeof(eargpt)/sizeof(eargpt[0]));//convert to vector
  //if no arguments given, print output
  if (argc==1) 
  {
   //pre-message
   std::cout << "total number of lines (commas) in matrix file should be steps * residues * ops * numatoms" << std::endl;
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
  const char * dargpt[] = {"--kappa", "0"};  //crashes if this is empty haha  
  size_t defSize = sizeof(dargpt)/sizeof(dargpt[0]);
  dargp = ParsingHelperNS::parseArgvToPair(defSize, dargpt, false );//not my favorite, not sure how c++ actually makes argc, just pass the real number of args, the array, and false

  ParsedInput inputs(&argp, &dargp, &eargp);

//std::vector<int> numstep;
std::vector<int> numskip;
//std::vector<int> numops;
//std::vector<int> numres;
//std::vector<int> numatoms;

//inputs.WriteOut("num-atoms", numstep ,"INT");
//inputs.WriteOut("steps",numskip,"INT");
inputs.WriteOut("num-skip",numskip,"INT");
//inputs.WriteOut("numops",numres,"INT");
//inputs.WriteOut("numres",numatoms,"INT");

std::vector<std::string> gradF;
std::vector<std::string> gradPhi;
inputs.WriteOut("gradF", gradF, "STRING");
inputs.WriteOut("gradPhi", gradPhi, "STRING");


std::vector<Real> kappa;
inputs.WriteOut("kappa", kappa, "REAL");


//read in files now
   std::ifstream gradFFile(gradF[0]);
    if (!gradFFile.is_open()) {
        std::cerr << "Error: Could not open gradF file.\n";
        return 1;
    }

    std::string gradFline;
    std::vector<std::string> gradFData;

    while (std::getline(gradFFile, gradFline, ',')) {
        gradFData.push_back(gradFline);
        }
    gradFFile.close();

  std::ifstream gradPhiFile(gradPhi[0]);
    if (!gradPhiFile.is_open()) {
        std::cerr << "Error: Could not open gradPhi file.\n";
        return 1;
    }
    std::string gradPhiline;
    std::vector<std::string> gradPhiData;

    while (std::getline(gradPhiFile, gradPhiline, ',')) {
        gradPhiData.push_back(gradPhiline);
        }
    gradFFile.close();

int numRes, numStep, numAtom, numCell;
numRes = std::stoi(gradFData[1]);
numStep = std::stoi(gradFData[0]);
numCell = std::stoi(gradFData[2]);
numAtom = std::stoi(gradPhiData[3]);




if (numRes != std::stoi(gradPhiData[1])) 
  {
  std::cerr << "inconsistent residue numbers" << std::endl;
  throw std::exception();
  }
if (numStep != std::stoi(gradPhiData[0])) 
  {
  std::cerr << "inconsistent step numbers" << std::endl;
  throw std::exception();
  }
if (numCell != std::stoi(gradPhiData[2])) 
  {
  std::cerr << "inconsistent number of OPs" << std::endl;
  throw std::exception();
  }
if (numStep  < numskip[0] || numskip[0] < 0)
  {
  std::cerr << "number of steps to skip needs to be less than overall number of steps, and also a positive integer. Zero means no skipping.\n" <<
               "The number of steps should be in terms of how many times the file was written to (usually every 5000 steps or so, not the overall number of steps." << std::endl; 
  throw std::exception();
  }
  

std::cout << "number of steps:        " << numStep << "\n" <<
             "number of residues:     " << numRes << "\n" << 
             "number of order params: " << numCell << "\n" << 
             "number of atoms:        " << numAtom << "\n" <<
             "skipping "<< numskip[0] << " steps\n" <<  std::endl;
//split this data up into workable files
std::vector<std::vector<std::vector<Real> > > gradFFinal;
  for (size_t step=0; step < numStep; ++step)
  {
    std::vector< std::vector < Real> > tempStep;
    for (size_t res=0; res < numRes; ++res)
    {
      std::vector < Real > tempRes;
      for (size_t cell=0; cell < numCell ; ++cell)
      {
      tempRes.push_back(stod(gradFData[3 + step*(numRes*numCell) + res*numCell + cell]));
      }
      tempStep.push_back(tempRes);
    }
    gradFFinal.push_back(tempStep);
  }

std::vector<std::vector< std::vector< std::vector <Vector3D> > > > gradPhiFinal;



for (size_t step=0; step < numStep; ++step)
{
  std::vector<std::vector< std::vector< Vector3D> > > tempStep;
  for (size_t res=0; res < numRes; ++res)
  {
  std::vector<std::vector< Vector3D > > tempRes; 
    for (size_t cell=0; cell < numCell ; ++cell)
    {
      std::vector<Vector3D> tempCell;
      for (size_t atom=0; atom<numAtom; ++atom)
      {
      Vector3D tempVec;
      Real x,y,z;
      std::stringstream ts(gradPhiData[4 + step*(numRes*numCell*numAtom) + res*(numCell*numAtom) + cell*numAtom + atom]);
      std::string sx;
      std::string sy;
      std::string sz;
      std::getline(ts, sx, ' ');
      std::getline(ts, sy, ' ');
      std::getline(ts, sz, ' ');
//std::cout << sx << " " << sy << " " << sz << std::endl;
      x = stod(sx);
      y = stod(sy);
      z = stod(sz);
      tempVec.x = x; tempVec.y = y; tempVec.z = z;
      tempCell.push_back(tempVec);
      }
    tempRes.push_back(tempCell);
    }
    tempStep.push_back(tempRes);
  }
  gradPhiFinal.push_back(tempStep);
}

std::cout << gradPhiFinal.size() << std::endl;
std::cout << gradPhiFinal[0].size() << std::endl;
std::cout << gradPhiFinal[0][0].size() << std::endl;
std::cout << gradPhiFinal[0][0][0].size() << std::endl;

std::string opdFileName = "myOPDFile";
OPDFile myOPDFile(opdFileName, OUT);
#define opindex 1
for (size_t step=numskip[0]; step < numStep; ++step)
  {
  ForceData myForceData;
  myForceData.clear();
  for (size_t numKappa=0; numKappa < 4; ++numKappa)
    myForceData.forceConstants.push_back(kappa[numKappa]);

  for (size_t i=0; i< 4; ++i)
    {
    if (i==opindex)
    {
    for (size_t res=0; res<numRes; ++res)
      {
      for (size_t cell=0; cell<numCell; ++cell)
        {
          myForceData.restraintForces.push_back(gradFFinal[step][res][cell]); 
        }
      }
    }
    else
    {
    for (size_t res=0; res<numRes; ++res)
      {
      for (size_t cell=0; cell<numCell; ++cell)
        {
          myForceData.restraintForces.push_back(0); 
        }
      }
      }

    } 
// now do metric tensor 
std::cout << "beginning metric tensor calculations" << std::endl;
myForceData.mMatrix.resize(numRes*numCell*4);
for (size_t i=0; i<4; ++i)
{
for (size_t resi=0; resi < numRes; ++resi)
  {
    for (size_t celli=0; celli < numCell; ++celli)
      {
       for (size_t resj=resi; resj < numRes; ++resj)
       {
      if (resi == resj)
        for (size_t cellj=celli; cellj < numCell; ++cellj)
          {
          Real Mij = 0.0;
          if (i==opindex){
          for (size_t atom=0; atom < numAtom; ++atom)
            Mij += gradPhiFinal[step][resi][celli][atom]*gradPhiFinal[step][resj][cellj][atom]; 
          myForceData.mMatrix(resi*numCell + celli,resj*numCell + cellj) = Mij;
           }
          else 
          {
          myForceData.mMatrix(resi*numCell + celli,resj*numCell + cellj) = Mij;

          }         
           }
      else
        for (size_t cellj=0; cellj < numCell; ++cellj)
          {
          Real Mij=0.0;
          if (i==opindex)
           {
          for (size_t atom=0; atom < numAtom; ++atom)
             Mij += gradPhiFinal[step][resi][celli][atom]*gradPhiFinal[step][resj][cellj][atom];
          myForceData.mMatrix(resi*numCell + celli,resj*numCell + cellj) = Mij;
         }
         else myForceData.mMatrix(resi*numCell + celli,resj*numCell + cellj) = Mij;

          }
        }
      }

    }
}
  std::cout << myForceData << std::endl;
  myOPDFile << myForceData;
  }

return 0;


}

