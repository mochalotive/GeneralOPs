#include "plumedinput/include/plumedfile.h"
#include "plumedinput/include/plumedparameters.h"
#include "mymath/include/vector3D.h"
//static int compareInternal(std::string& s, std::string t)//notice t is a copy :)
//  {
//   std::string::iterator end_pos = std::remove(t.begin(), t.end(), ' ');
//   t.erase(end_pos, t.end()); 
//   return (s.compare(t));
//  }
#define R_CUTOFF 0.01
PlumedFile::PlumedFile():IOFile(){};
PlumedFile::PlumedFile(std::string const &fileName):IOFile(fileName, OUT){};
//operator

PlumedFile &operator<<(PlumedFile &plumedFile, PlumedParameters &parameters) 
{

//pretty much always operating in angstroms because we're coming from pdbs
//and all our code uses pdbs 
plumedFile << "UNITS LENGTH=A" << std::endl;//always need that
//could change in the future if doing a different origin
//used in the arounds
plumedFile << "o: FIXEDATOM AT=0,0,0" << std::endl;
//if any of the parameters use COM, include coms.dat only once
for (size_t i=0; i<parameters.myPlumedSystemParameters.residueNames.size(); i++)
  if (parameters.myPlumedSystemParameters.COM[i])
    {
    plumedFile << "INCLUDE FILE=coms.dat" << std::endl;
    break;
    }

//print a vector of ones for making neighbor vectors
for (size_t i=0; i< parameters.myPlumedSystemParameters.residueNames.size(); ++i)
    plumedFile << "num" << parameters.myPlumedSystemParameters.residueNames[i] << ": ONES SIZE=" << parameters.residueNumbers[i] << "\n";


//put COMS in stdout, make sure to pipe this file.
for (size_t i=0; i<parameters.myPlumedSystemParameters.residueNames.size(); i++)//very very very bad, but I made the indices of myPlumedSystemParameters.residueNames match the first indices of allRes, so its ok
  {//looping over residues to be restrained eg ACN,SFD. Order within pdb doesn't matter
  if (parameters.myPlumedSystemParameters.COM[i])
    {
    //2025 hotfix - cout this, better to have this on a separate file for readability, just > coms.dat. 
    //in the future add a proper file writing
    int COMcounter=0;//when writing COM [number] statements
    int atomIndex=0; //when giving the COM statement a list of atoms
    for (size_t j=0; j<parameters.mySystemDataPackage.mySystem.size(); j++)
      { //loopinf over each molecule within the pdb file, and check if it matches the current residue
      if (parameters.compareInternal(parameters.myPlumedSystemParameters.residueNames[i], parameters.mySystemDataPackage.mySystem[j].name()) == 0)//check for residue matching, takes 2 strings, 2nd string is a copy
        {//if matching res, add a new COM in plumed.dat
        //plumedFile << "com" << parameters.myPlumedSystemParameters.residueNames[i] << COMcounter << ": COM ATOMS=";
        std::cout << "com" << parameters.myPlumedSystemParameters.residueNames[i] << COMcounter << ": COM ATOMS=";
        for (size_t k=0; k<(int)parameters.mySystemDataPackage.mySystem[j].numAtoms(); k++)
          {//loop over number of atoms within molecule
          if (k==(int)parameters.mySystemDataPackage.mySystem[j].numAtoms()-1) {/*plumedFile << atomIndex+1;*/std::cout << atomIndex+1; }// 
          else {/*plumedFile << atomIndex+1 << ",";*/std::cout << atomIndex+1 << ",";}
          atomIndex++;//increment atom index
          }
          //plumedFile << "\n";
          std::cout << "\n";
          COMcounter++;   
        }  
      else {atomIndex=atomIndex+(int)parameters.mySystemDataPackage.mySystem[j].numAtoms();}//if residues don't match just add the number of atoms for current residue
      }
    }
  }//all coms, if needed

//define shortcuts for all types
for (size_t i=0; i< parameters.myPlumedSystemParameters.residueNames.size(); i++)
  {
  if (parameters.myPlumedSystemParameters.COM[i])
    {
    plumedFile << parameters.myPlumedSystemParameters.residueNames[i] << ": GROUP ATOMS="; 
    for (size_t j=0; j<parameters.residueNumbers[i]; j++)
      {
      if (j==parameters.residueNumbers[i]-1){plumedFile << "com" << parameters.myPlumedSystemParameters.residueNames[i]<< j << "\n";}
      else plumedFile << "com" << parameters.myPlumedSystemParameters.residueNames[i]<< j << ",";
      }
    }
  else //i.e if we're using framepoints (faster in plumed) 
    {
    int resCounter=0;
    int moleculeIndex=0;
    plumedFile << parameters.myPlumedSystemParameters.residueNames[i] << ": GROUP ATOMS=";
    for (size_t j=0; j< parameters.mySystemDataPackage.mySystem.size(); j++)
      {
      if (parameters.compareInternal(parameters.myPlumedSystemParameters.residueNames[i], parameters.mySystemDataPackage.mySystem[j].name()) == 0 &&
          resCounter==parameters.residueNumbers[i]-1)
          {
          plumedFile <<moleculeIndex +1+ parameters.mySystemDataPackage.myMoleculeMaps[i][0].framePoints[0].atoms[0] << "\n";
          resCounter++;
          moleculeIndex = moleculeIndex + (int)parameters.mySystemDataPackage.mySystem[j].numAtoms();
          }
      else if (parameters.compareInternal(parameters.myPlumedSystemParameters.residueNames[i], parameters.mySystemDataPackage.mySystem[j].name())
                == 0)
         {
           plumedFile <<moleculeIndex +1+ parameters.mySystemDataPackage.myMoleculeMaps[i][0].framePoints[0].atoms[0] << ",";
           resCounter++;
           moleculeIndex = moleculeIndex + (int)parameters.mySystemDataPackage.mySystem[j].numAtoms();
         }
      else moleculeIndex = moleculeIndex + (int)parameters.mySystemDataPackage.mySystem[j].numAtoms();
      }
    }
  }

Real coreCutoff = 0.0;
coreCutoff = parameters.myPlumedSystemParameters.contact_cutoff[0];
for (size_t i=0; i< parameters.myPlumedSystemParameters.residueNames.size(); ++i)
  plumedFile << "core" << parameters.myPlumedSystemParameters.residueNames[i] << ": CONTACT_MATRIX GROUPA=" <<
  parameters.myPlumedSystemParameters.residueNames[i] << " GROUPB=o COMPONENTS SWITCH={RATIONAL R_0=" << R_CUTOFF << " D_0=" << coreCutoff*1.414 << " D_MAX=" << 
  coreCutoff*1.414+1 << "}\n";
for (size_t i=0; i< parameters.myPlumedSystemParameters.residueNames.size(); ++i)
  plumedFile << "f" <<  parameters.myPlumedSystemParameters.residueNames[i] << ":  FLATTEN ARG=core" <<  parameters.myPlumedSystemParameters.residueNames[i] <<
                ".w\n";


////if using MASKS, use first res as contact, if you need more you need to do it yourself 
////itll make the input more of a mess than it already is
//for (size_t i=1; i< parameters.myPlumedSystemParameters.residueNames.size(); ++i)
//  if (parameters.myPlumedSystemParameters.contact[i]){
//    plumedFile << parameters.myPlumedSystemParameters.residueNames[0] << parameters.myPlumedSystemParameters.residueNames[i] <<
//                  ": CONTACT_MATRIX GROUPA=" <<  parameters.myPlumedSystemParameters.residueNames[i] <<
//                  " GROUPB=" <<  parameters.myPlumedSystemParameters.residueNames[0] << " SWITCH={RATIONAL R_0="<< R_CUTOFF<< " D_0=";
//    if (parameters.myPlumedSystemParameters.contact_cutoff[i] > 0.0) plumedFile <<  parameters.myPlumedSystemParameters.contact_cutoff[i] <<                        
//    " D_MAX=" << parameters.myPlumedSystemParameters.contact_cutoff[i]+2.0  << "}\n";
//    else  plumedFile << parameters.myPlumedSystemParameters.cutoff[i] <<                        
//    " D_MAX=" << parameters.myPlumedSystemParameters.cutoff[i]+2.0  << "}\n";
//}
//
////make no norm neighbor vector
//for (size_t i=1; i< parameters.myPlumedSystemParameters.residueNames.size(); ++i)
//  if (parameters.myPlumedSystemParameters.contact[i]){
//      
//      plumedFile << "no_norm_" <<  parameters.myPlumedSystemParameters.residueNames[i] << 
//      ":  MATRIX_VECTOR_PRODUCT ARG=" <<parameters.myPlumedSystemParameters.residueNames[0] << parameters.myPlumedSystemParameters.residueNames[i] << 
//      "," << "num" << parameters.myPlumedSystemParameters.residueNames[0] << "\n"
//      << "f" <<  parameters.myPlumedSystemParameters.residueNames[i]     << ": MORE_THAN ARG=" << "no_norm_"<< 
//          parameters.myPlumedSystemParameters.residueNames[i] << " SWITCH={EXP R_0=" << R_CUTOFF << "D_0=0.5 D_MAX=3.0}\n";
//  }

//contact matrices
//dont mask the solute

  plumedFile << parameters.myPlumedSystemParameters.residueNames[0] << "c: CONTACT_MATRIX GROUP=" << parameters.myPlumedSystemParameters.residueNames[0] << 
                " SWITCH={RATIONAL R_0=" << R_CUTOFF <<" D_0=" << parameters.myPlumedSystemParameters.cutoff[0] << " D_MAX=" << 
                parameters.myPlumedSystemParameters.cutoff[0] + 1 << "} COMPONENTS ";
  if (parameters.myPlumedSystemParameters.contact[0]) plumedFile << " MASK=f" <<  parameters.myPlumedSystemParameters.residueNames[0] << "\n";
  else plumedFile << "\n";

for (size_t i=1; i< parameters.myPlumedSystemParameters.residueNames.size(); ++i){
  plumedFile << parameters.myPlumedSystemParameters.residueNames[i] << "c: CONTACT_MATRIX GROUP=" << parameters.myPlumedSystemParameters.residueNames[i] << 
                " SWITCH={RATIONAL R_0=" << R_CUTOFF <<" D_0=" << parameters.myPlumedSystemParameters.cutoff[i] << " D_MAX=" << 
                parameters.myPlumedSystemParameters.cutoff[i] + 1 << "} COMPONENTS ";
  if (parameters.myPlumedSystemParameters.contact[i]) plumedFile << "NL_CUTOFF=" << parameters.myPlumedSystemParameters.cutoff[i]+1.5<< " NL_STRIDE=10\n";
  else plumedFile << "\n";
}

for (size_t i=0; i< parameters.myPlumedSystemParameters.residueNames.size(); ++i)
  {
  plumedFile << parameters.myPlumedSystemParameters.residueNames[i] << "n: MATRIX_VECTOR_PRODUCT ARG=" <<
                parameters.myPlumedSystemParameters.residueNames[i] << "c.w,num" << parameters.myPlumedSystemParameters.residueNames[i] << "\n";
  }

//Do Quaternions, if needed. 
for (size_t i=0; i<parameters.myPlumedSystemParameters.residueNames.size(); i++)
  {
  switch(parameters.fileType[i])
    {
    case DOPS://only case where we don't need Quaterions
      {break;}
       
    case NOPS:
      {break;}
    case LOPS:
      {break;}
    default:
     {//all other cases need quaternions
     plumedFile << "quat" << parameters.myPlumedSystemParameters.residueNames[i] << ": QUATERNION ";
     int resCounter=0;
     int moleculeIndex=0;
     for (size_t j=0; j< parameters.mySystemDataPackage.mySystem.size(); j++)
       {
       if (parameters.compareInternal(parameters.myPlumedSystemParameters.residueNames[i], parameters.mySystemDataPackage.mySystem[j].name()) == 0 && 
           resCounter==parameters.residueNumbers[i]-1)
           {  plumedFile <<"ATOMS" << resCounter+1 << "=" <<  moleculeIndex +1+ parameters.mySystemDataPackage.myMoleculeMaps[i][0].framePoints[0].atoms[0] << 
                           "," << moleculeIndex +1+ parameters.mySystemDataPackage.myMoleculeMaps[i][0].framePoints[1].atoms[0]<<"," << 
                           moleculeIndex +1+ parameters.mySystemDataPackage.myMoleculeMaps[i][0].framePoints[2].atoms[0];
           resCounter++;
           moleculeIndex = moleculeIndex + (int)parameters.mySystemDataPackage.mySystem[j].numAtoms();
           }
       else if (parameters.compareInternal(parameters.myPlumedSystemParameters.residueNames[i], parameters.mySystemDataPackage.mySystem[j].name())
                       == 0)
         {  plumedFile << "ATOMS" << resCounter+1 << "=" << moleculeIndex +1+ parameters.mySystemDataPackage.myMoleculeMaps[i][0].framePoints[0].atoms[0] << 
                          "," << moleculeIndex +1+ parameters.mySystemDataPackage.myMoleculeMaps[i][0].framePoints[1].atoms[0]<<"," << 
                           moleculeIndex +1+ parameters.mySystemDataPackage.myMoleculeMaps[i][0].framePoints[2].atoms[0] << " ";

            resCounter++;
            moleculeIndex = moleculeIndex + (int)parameters.mySystemDataPackage.mySystem[j].numAtoms();
         }
         else moleculeIndex = moleculeIndex + (int)parameters.mySystemDataPackage.mySystem[j].numAtoms();
        }
       if (parameters.myPlumedSystemParameters.contact[i]) plumedFile << " " << "MASK=" << "f" <<  
                                                                                  parameters.myPlumedSystemParameters.residueNames[i] << "\n";
       else plumedFile << "\n";
     break; 
     } 
  }
}

//now do all quaternion and quaternion bond products
//SFDqbp: QUATERNION_BOND_PRODUCT_MATRIX ARG=quatSFDROPS.*,sc.* 
//MEOHqbp: QUATERNION_BOND_PRODUCT_MATRIX ARG=quatMEOHROPS.*,mc.* MASK=neighbors

for (size_t i=0; i<parameters.myPlumedSystemParameters.residueNames.size(); i++)
  {
  switch(parameters.fileType[i])
    {
    case DOPS://only case where we don't need Quaterions
      {break;}

    case NOPS:
      {break;}
    case LOPS:
      {break;}
    case TOPS: //total ops, need both kinds of product 
      {
      plumedFile << "qbp" << parameters.myPlumedSystemParameters.residueNames[i] << ": QUATERNION_BOND_PRODUCT_MATRIX ARG=" <<
                    "quat" << parameters.myPlumedSystemParameters.residueNames[i] << ".*," << 
                     parameters.myPlumedSystemParameters.residueNames[i] << "c.*";
      if (parameters.myPlumedSystemParameters.contact[i]) 
        plumedFile << " MASK=f" << parameters.myPlumedSystemParameters.residueNames[i] << "\n";
      else plumedFile << "\n"; 

      //this always gets a mask
      plumedFile << "qp" <<  parameters.myPlumedSystemParameters.residueNames[i] << ": QUATERNION_PRODUCT_MATRIX ARG=" <<
      "quat" << parameters.myPlumedSystemParameters.residueNames[i] << ".*,quat" << parameters.myPlumedSystemParameters.residueNames[i] << ".* " <<
       "MASK=" << parameters.myPlumedSystemParameters.residueNames[i] << "c.w\n";
      break;
      }

    case BOPS:
{
    plumedFile << "qbp" << parameters.myPlumedSystemParameters.residueNames[i] << ": QUATERNION_BOND_PRODUCT_MATRIX ARG=" <<
                    "quat" << parameters.myPlumedSystemParameters.residueNames[i] << ".*," <<
                     parameters.myPlumedSystemParameters.residueNames[i] << "c.*";
      if (parameters.myPlumedSystemParameters.contact[i])
        plumedFile << " MASK=f" << parameters.myPlumedSystemParameters.residueNames[i] << "\n";
      else plumedFile << "\n";
break;


}
    case ROPS:
{
      plumedFile << "qp" <<  parameters.myPlumedSystemParameters.residueNames[i] << ": QUATERNION_PRODUCT_MATRIX ARG=" <<
      "quat" << parameters.myPlumedSystemParameters.residueNames[i] << ".*,quat" << parameters.myPlumedSystemParameters.residueNames[i] << ".* " <<
       "MASK=" << parameters.myPlumedSystemParameters.residueNames[i] << "c.w\n";

break;
}
    default:
      {
      break;      


      }
    }
  }

//make op eqn
for (size_t i=0; i<parameters.myPlumedSystemParameters.residueNames.size(); i++)
  {
  plumedFile << "ops" << parameters.myPlumedSystemParameters.residueNames[i] << "p: CUSTOM ARG=";
  switch(parameters.fileType[i])
  {
  case DOPS:
  {
    plumedFile << parameters.myPlumedSystemParameters.residueNames[i]  << "c.* VAR=d,x,y,z FUNC=";
    size_t eqSize = parameters.mySystemDataPackage.myEq[i].getSize();
    for (size_t eq=0; eq<eqSize ; eq++)
      {
        if (eq == eqSize - 1) plumedFile << parameters.mySystemDataPackage.myEq[i].getString(eq) << " PERIODIC=NO\n";
        else plumedFile << parameters.mySystemDataPackage.myEq[i].getString(eq) << "+";
      }
  break;
  }
  case BOPS:
  {
    plumedFile << parameters.myPlumedSystemParameters.residueNames[i]  << "c.*,qbp" <<
                  parameters.myPlumedSystemParameters.residueNames[i] << ".* VAR=d,x,y,z,o,l,m,n FUNC=";
    size_t eqSize = parameters.mySystemDataPackage.myEq[i].getSize();
    for (size_t eq=0; eq<eqSize ; eq++)
      {
        if (eq == eqSize - 1) plumedFile << parameters.mySystemDataPackage.myEq[i].getString(eq) << " PERIODIC=NO\n";
        else plumedFile << parameters.mySystemDataPackage.myEq[i].getString(eq) << "+";
      }
  break;
  }

  case ROPS:
  {
  plumedFile << parameters.myPlumedSystemParameters.residueNames[i]  << "c.*,qp" <<
                  parameters.myPlumedSystemParameters.residueNames[i] << ".* VAR=d,x,y,z,w,i,j,k FUNC=";
    size_t eqSize = parameters.mySystemDataPackage.myEq[i].getSize();
    for (size_t eq=0; eq<eqSize ; eq++)
      {
        if (eq == eqSize - 1) plumedFile << parameters.mySystemDataPackage.myEq[i].getString(eq) << " PERIODIC=NO\n";
        else plumedFile << parameters.mySystemDataPackage.myEq[i].getString(eq) << "+";
      }
  break;

  }
  case TOPS:
  {
    plumedFile << parameters.myPlumedSystemParameters.residueNames[i]  << "c.*,qbp" <<
                  parameters.myPlumedSystemParameters.residueNames[i] << ".*,qp" <<
                  parameters.myPlumedSystemParameters.residueNames[i] << ".* VAR=d,x,y,z,o,l,m,n,w,i,j,k FUNC=";
    size_t eqSize = parameters.mySystemDataPackage.myEq[i].getSize();
    for (size_t eq=0; eq<eqSize ; eq++)
      {
        if (eq == eqSize - 1) plumedFile << parameters.mySystemDataPackage.myEq[i].getString(eq) << " PERIODIC=NO\n";
        else plumedFile << parameters.mySystemDataPackage.myEq[i].getString(eq) << "+";
      }
  break;
 
  }
  default: {break;}
  }
  }

//now matrix vector product 

for (size_t i=0; i< parameters.myPlumedSystemParameters.residueNames.size(); ++i)
  plumedFile << "ops" << parameters.myPlumedSystemParameters.residueNames[i] << "t: MATRIX_VECTOR_PRODUCT ARG=ops" <<
                parameters.myPlumedSystemParameters.residueNames[i] << "p,num" << parameters.myPlumedSystemParameters.residueNames[i] << "\n";

for (size_t i=0; i< parameters.myPlumedSystemParameters.residueNames.size(); ++i)
  plumedFile << "ops" << parameters.myPlumedSystemParameters.residueNames[i] <<  ": CUSTOM ARG=ops" <<
  parameters.myPlumedSystemParameters.residueNames[i] << "t," << parameters.myPlumedSystemParameters.residueNames[i] << "n FUNC=x/(y+1) PERIODIC=NO\n";

for (size_t i=0; i< parameters.myPlumedSystemParameters.residueNames.size(); ++i)
  plumedFile << parameters.myPlumedSystemParameters.residueNames[i] << "x: FLATTEN ARG=core" <<
                parameters.myPlumedSystemParameters.residueNames[i] << ".x\n" << 
                parameters.myPlumedSystemParameters.residueNames[i] << "y: FLATTEN ARG=core" <<
                parameters.myPlumedSystemParameters.residueNames[i] << ".y\n" << 
                parameters.myPlumedSystemParameters.residueNames[i] << "z: FLATTEN ARG=core" <<
                parameters.myPlumedSystemParameters.residueNames[i] << ".z\n"; 


  //this is really unbelievable pollution of my code I managed to eat up every common varaible name
  //and boy do i not care
  //this just makes AROUNDs
  Lattice l = parameters.mySystemDataPackage.mySystem.lattice();

  Vector3D x = l.latticeVector(0);
  Vector3D y = l.latticeVector(1);
  Vector3D z = l.latticeVector(2);
  x[0] = 2*coreCutoff;
  y[1] = 2*coreCutoff;
  z[2] = 2*coreCutoff;
  int a = parameters.myPlumedSystemParameters.box_size[0];
  int b = parameters.myPlumedSystemParameters.box_size[1];
  int c = parameters.myPlumedSystemParameters.box_size[2];
  

  int numResidue = 0;
  for (size_t i=0; i<parameters.residueNumbers.size();i++) numResidue = numResidue+parameters.residueNumbers[i];
  Real m,n,o;
  m = x[0]/(Real)a;
  n = y[1]/(Real)b;
  o = z[2]/(Real)c;
  
  //a,b,c is number of boxes
  //x,y,z are box vectors
  //i,j,k are indices ranging over a,b,c
  
  //gromacs supports all kinds of box types - https://manual.gromacs.org/current/reference-manual/algorithms/periodic-boundary-conditions.html
  //and the default center is x/2, y/2, z/2 https://manual.gromacs.org/current/onlinehelp/gmx-editconf.html check -center flag


for (size_t res=0; res < parameters.myPlumedSystemParameters.residueNames.size(); res++){
int aroundCounter = 0;
for (size_t k=0; k<c; k++)
  {
  for (size_t j=0; j<b; j++)
    {
    for (size_t i=0; i<a; i++)
      {
      plumedFile << "a" << parameters.myPlumedSystemParameters.residueNames[res] << aroundCounter << ": CUSTOM ARG=" <<
                           parameters.myPlumedSystemParameters.residueNames[res] << "x," <<
                           parameters.myPlumedSystemParameters.residueNames[res] << "y," <<
                           parameters.myPlumedSystemParameters.residueNames[res] << "z," << 
                           "ops" <<parameters.myPlumedSystemParameters.residueNames[res] << " VAR=x,y,z,o FUNC=";
        //here goes xlower, xupper etc
        plumedFile << "(1-tanh((" << -x[0]/2.0 + i*m << "-x)/" << R_CUTOFF << "))*(1-tanh((x-" << -x[0]/2.0 + (i+1)*m << ")/" << R_CUTOFF << "))*";
        plumedFile << "(1-tanh((" << -y[1]/2.0 + j*n << "-y)/" << R_CUTOFF << "))*(1-tanh((y-" << -y[1]/2.0 + (j+1)*n << ")/" << R_CUTOFF << "))*";
        plumedFile << "(1-tanh((" << -z[2]/2.0 + k*o << "-z)/" << R_CUTOFF << "))*(1-tanh((z-" << -z[2]/2.0 + (k+1)*o << ")/" << R_CUTOFF << "))*o/64";
        plumedFile << " PERIODIC=NO MASK=f" << parameters.myPlumedSystemParameters.residueNames[res] << "\n";
        aroundCounter++;
      }
      }
    } 
}
//now need restraints
//skipping moving restraints for now?

//loop thru number of boxes, and apply restraints from COP data.
//

//sum up arounds times ops 





//ra (ops times around)
//for (size_t res=0; res < parameters.myPlumedSystemParameters.residueNames.size(); res++){
//size_t raCounter=0;
//for (size_t k=0; k<c; k++)
//  {
//  for (size_t j=0; j<b; j++)
//    {
//    for (size_t i=0; i<a; i++)
//      {
//      plumedFile << "ra" << parameters.myPlumedSystemParameters.residueNames[res] 
//      << raCounter << ": CUSTOM ARG=ops" 
//      <<  parameters.myPlumedSystemParameters.residueNames[res] << ",a" 
//      << parameters.myPlumedSystemParameters.residueNames[res] << raCounter 
//      << " FUNC=x*y PERIODIC=NO\n";
//      raCounter++;
//}}}}

//summing arounds
for (size_t res=0; res < parameters.myPlumedSystemParameters.residueNames.size(); res++){
size_t saCounter=0;
for (size_t k=0; k<c; k++)
  {
  for (size_t j=0; j<b; j++)
    {
    for (size_t i=0; i<a; i++)
      {
      plumedFile << "sa" << parameters.myPlumedSystemParameters.residueNames[res] 
      << saCounter << ": SUM ARG=a"
      << parameters.myPlumedSystemParameters.residueNames[res] << saCounter 
      << " PERIODIC=NO\n";
      saCounter++;
}}}}
////summing ops within box
//for (size_t res=0; res < parameters.myPlumedSystemParameters.residueNames.size(); res++){
//size_t sraCounter=0;
//for (size_t k=0; k<c; k++)
//  {
//  for (size_t j=0; j<b; j++)
//    {
//    for (size_t i=0; i<a; i++)
//      {
//      plumedFile << "sra" << parameters.myPlumedSystemParameters.residueNames[res] 
//      << sraCounter << ": SUM ARG=ra"
//      << parameters.myPlumedSystemParameters.residueNames[res] << sraCounter 
//      << " PERIODIC=NO\n";
//      sraCounter++;
//}}}}
////finally normalizing summed ops with summed arounds
//for (size_t res=0; res < parameters.myPlumedSystemParameters.residueNames.size(); res++){
//size_t opfCounter=0;
//for (size_t k=0; k<c; k++)
//  {
//  for (size_t j=0; j<b; j++)
//    {
//    for (size_t i=0; i<a; i++)
//      {
//      plumedFile << "opf" << parameters.myPlumedSystemParameters.residueNames[res] 
//      << opfCounter << ": CUSTOM ARG=" 
//      << "sa" << parameters.myPlumedSystemParameters.residueNames[res] 
//      << opfCounter << ","
//      << "sra" 
//      << parameters.myPlumedSystemParameters.residueNames[res] << opfCounter
//      << " FUNC=(y*(1/x)) PERIODIC=NO\n";
//      opfCounter++;
//}}}}

int numstep = parameters.mySystemDataPackage.steps.first;
int numdragstep = parameters.mySystemDataPackage.steps.second;
//check if we need to drag
for (size_t res=0; res < parameters.myPlumedSystemParameters.residueNames.size() ; ++res)
{
if (numdragstep == 0 ) 
  {
  plumedFile << "r" << parameters.myPlumedSystemParameters.residueNames[res] << ": RESTRAINT ARG=";
  for (size_t i=0; i < a*b*c; ++i)
    {
    if (i == a*b*c-1)
      plumedFile << "sa" << parameters.myPlumedSystemParameters.residueNames[res] << i << " "; 
      
    else plumedFile << "sa" << parameters.myPlumedSystemParameters.residueNames[res] << i << ","; 
    }
  //now print OPs
  plumedFile << "AT=";
  switch (parameters.fileType[res])
    {
    case (DOPS): 
      {
      for (size_t i=0; i< a*b*c ; ++i)
        {
        if (i == a*b*c -1 ) plumedFile << parameters.mySystemDataPackage.myCrystalOrderParameters[res].distance[i] << " "; 
        else plumedFile << parameters.mySystemDataPackage.myCrystalOrderParameters[res].distance[i] << ",";
        }//loop over ops
        break;
      }//end dops
    case (BOPS): 
      {
      for (size_t i=0; i< a*b*c ; ++i)
        {
        if (i == a*b*c -1 ) plumedFile << parameters.mySystemDataPackage.myCrystalOrderParameters[res].bondOrientation[i] << " "; 
        else plumedFile << parameters.mySystemDataPackage.myCrystalOrderParameters[res].bondOrientation[i] << ",";
        }//loop over ops
      break;
      }//end bops
    case (ROPS): 
      {
      for (size_t i=0; i< a*b*c ; ++i)
        {
        if (i == a*b*c -1 ) plumedFile << parameters.mySystemDataPackage.myCrystalOrderParameters[res].relativeOrientation[i] << " "; 
        else plumedFile << parameters.mySystemDataPackage.myCrystalOrderParameters[res].relativeOrientation[i] << ",";
        }//loop over ops
        break;
      }//end rops

    case (TOPS): 
      {
      for (size_t i=0; i< a*b*c ; ++i)
        {
        if (i == a*b*c -1 ) plumedFile << parameters.mySystemDataPackage.myCrystalOrderParameters[res].total[i] << " "; 
        else plumedFile << parameters.mySystemDataPackage.myCrystalOrderParameters[res].total[i] << ",";
        }//loop over ops
        break;
      }//end tops
    default: {break;}
    }//end switch
  plumedFile << "KAPPA=";
  for (size_t i=0; i< a*b*c; ++i)
    {
    if (i== a*b*c-1) plumedFile << parameters.myPlumedSystemParameters.kappa[res] << "\n";
    else plumedFile << parameters.myPlumedSystemParameters.kappa[res] << ",";
    }
  }//end case with no dragging
//complicated one
//MOVINGRESTRAINT ...
//   ARG=d 
//   STEP0=0 AT0=1.0 KAPPA0=100.0 
//   STEP1=1000 AT1=2.0 
//   STEP2=2000 AT2=1.0 
//   STEP3=2500 KAPPA3=0.0 
//...
else
  {
  //loop over residues
  plumedFile << "MOVINGRESTRAINT ...\n";
  plumedFile << "  ARG=";
  for (size_t j=0; j< a*b*c ; ++j)
    {
    if (j == a*b*c - 1) plumedFile << "sa" << parameters.myPlumedSystemParameters.residueNames[res] << j << "\n";
    else plumedFile << "sa" << parameters.myPlumedSystemParameters.residueNames[res] << j << ",";
    }
  int stepCounter = 0;
  int numDragIncr = numdragstep/ 5000; 
  for (size_t i=0; i< numDragIncr+1; ++i) 
  {
   plumedFile << "   STEP" <<  i << "=" << 5000*i << " AT" << i << "="; //now loop over ops
 

  switch (parameters.fileType[res])
  {
  case DOPS:
  {  
    for (size_t j=0; j< a*b*c ; ++j)
    {
    if (i == a*b*c -1 ) plumedFile << parameters.mySystemDataPackage.myOldCrystalOrderParameters[res].distance[j]  + 
    (parameters.mySystemDataPackage.myCrystalOrderParameters[res].distance[j]-parameters.mySystemDataPackage.myOldCrystalOrderParameters[res].distance[j])*i/Real(numDragIncr)                                        << " "; 
    else plumedFile << parameters.mySystemDataPackage.myOldCrystalOrderParameters[res].distance[j]  +
    (parameters.mySystemDataPackage.myCrystalOrderParameters[res].distance[j]-parameters.mySystemDataPackage.myOldCrystalOrderParameters[res].distance[j])*i/Real(numDragIncr)                                        << ",";
    }
break;
  }
  case BOPS:
  {  
    for (size_t j=0; j< a*b*c ; ++j)
    {
    if (i == a*b*c -1 ) plumedFile << parameters.mySystemDataPackage.myOldCrystalOrderParameters[res].bondOrientation[j]  +
    (parameters.mySystemDataPackage.myCrystalOrderParameters[res].bondOrientation[j]-parameters.mySystemDataPackage.myOldCrystalOrderParameters[res].bondOrientation[j])*i/Real(numDragIncr) << " "; 
    else plumedFile << parameters.mySystemDataPackage.myOldCrystalOrderParameters[res].bondOrientation[j]  +
    (parameters.mySystemDataPackage.myCrystalOrderParameters[res].bondOrientation[j]-parameters.mySystemDataPackage.myOldCrystalOrderParameters[res].bondOrientation[j])*i/Real(numDragIncr) << ",";
    }

break;
  }
  case ROPS:
  {  
    for (size_t j=0; j< a*b*c ; ++j)
    {
    if (i == a*b*c -1 ) plumedFile << parameters.mySystemDataPackage.myOldCrystalOrderParameters[res].relativeOrientation[j]  +
    (parameters.mySystemDataPackage.myCrystalOrderParameters[res].relativeOrientation[j]-parameters.mySystemDataPackage.myOldCrystalOrderParameters[res].relativeOrientation[j])*i/Real(numDragIncr) << " "; 
    else plumedFile << parameters.mySystemDataPackage.myOldCrystalOrderParameters[res].relativeOrientation[j]  +
    (parameters.mySystemDataPackage.myCrystalOrderParameters[res].relativeOrientation[j]-parameters.mySystemDataPackage.myOldCrystalOrderParameters[res].relativeOrientation[j])*i/Real(numDragIncr) << ",";
    }
break;
  }
  case TOPS:
  {  
    for (size_t j=0; j< a*b*c ; ++j)
    {
    if (i == a*b*c -1 ) plumedFile << parameters.mySystemDataPackage.myOldCrystalOrderParameters[res].total[j]  +
    (parameters.mySystemDataPackage.myCrystalOrderParameters[res].total[j]-parameters.mySystemDataPackage.myOldCrystalOrderParameters[res].total[j])*i/Real(numDragIncr) << " "; 
    else plumedFile << parameters.mySystemDataPackage.myOldCrystalOrderParameters[res].total[j]  +
    (parameters.mySystemDataPackage.myCrystalOrderParameters[res].total[j]-parameters.mySystemDataPackage.myOldCrystalOrderParameters[res].total[j])*i/Real(numDragIncr) << ",";

//    parameters.mySystemDataPackage.myCrystalOrderParameters[i].total[j] 
    }
   break;
  }
  default: {break;}
  }
  
  plumedFile << " KAPPA" <<  i << "="; //now loop over ops
    for (size_t j=0; j< a*b*c ; ++j)
    {
    if (j == a*b*c - 1) plumedFile << parameters.myPlumedSystemParameters.kappa[res] << "\n";
    else plumedFile << parameters.myPlumedSystemParameters.kappa[res] << ",";
    }
}
//print final step and finish 
plumedFile << "...\n";

  }
}

//now print the individual OP's so I can calculate dF/dphi
for (size_t res=0; res < parameters.myPlumedSystemParameters.residueNames.size(); ++res)
for (size_t i=0; i<a*b*c; ++i) 
  {
  switch(parameters.fileType[res])
    {
    case (DOPS) :
    {
      plumedFile << "#csa" << parameters.myPlumedSystemParameters.residueNames[res] << i << ":  CUSTOM ARG=sa" <<
                  parameters.myPlumedSystemParameters.residueNames[res] << i << " FUNC=-" << parameters.myPlumedSystemParameters.kappa[res] <<
                  "*(x-" <<       parameters.mySystemDataPackage.myCrystalOrderParameters[res].distance[i] << ") PERIODIC=NO\n"; 
    break;
    }
    case (BOPS) :
    {
      plumedFile << "#csa" << parameters.myPlumedSystemParameters.residueNames[res] << i << ":  CUSTOM ARG=sa" <<
                  parameters.myPlumedSystemParameters.residueNames[res] << i << " FUNC=-" << parameters.myPlumedSystemParameters.kappa[res] <<
                  "*(x-" <<       parameters.mySystemDataPackage.myCrystalOrderParameters[res].bondOrientation[i] << ") PERIODIC=NO\n"; 
    break;
    }

    case (ROPS) :
    {
      plumedFile << "#csa" << parameters.myPlumedSystemParameters.residueNames[res] << i << ":  CUSTOM ARG=sa" <<
                  parameters.myPlumedSystemParameters.residueNames[res] << i << " FUNC=-" << parameters.myPlumedSystemParameters.kappa[res] <<
                  "*(x-" <<       parameters.mySystemDataPackage.myCrystalOrderParameters[res].relativeOrientation[i] << ") PERIODIC=NO\n"; 
    break;
    }
    case (TOPS) :
    {
      plumedFile << "#csa" << parameters.myPlumedSystemParameters.residueNames[res] << i << ":  CUSTOM ARG=sa" <<
                  parameters.myPlumedSystemParameters.residueNames[res] << i << " FUNC=-" << parameters.myPlumedSystemParameters.kappa[res] <<
                  "*(x-" <<       parameters.mySystemDataPackage.myCrystalOrderParameters[res].total[i] << ") PERIODIC=NO\n"; 
    break;
    default: {break;}
    }
    }
  }
//print

for (size_t res=0; res < parameters.myPlumedSystemParameters.residueNames.size(); ++res)
  {
  plumedFile << "psa" << parameters.myPlumedSystemParameters.residueNames[res] << ": PRINT ARG=";
  for (size_t i=0; i < a*b*c; ++i) 
    {
    if (i== a*b*c-1)
      plumedFile << "sa" << parameters.myPlumedSystemParameters.residueNames[res] << (a*b*c-1-i) << " STRIDE=10000\n";
    else
      plumedFile << "sa" << parameters.myPlumedSystemParameters.residueNames[res] << (a*b*c-1-i) << ",";
    }  
  }



for (size_t res=0; res < parameters.myPlumedSystemParameters.residueNames.size(); ++res)
  {
  plumedFile << "#csa" << parameters.myPlumedSystemParameters.residueNames[res] << ": PRINT ARG=";
  for (size_t i=0; i<a*b*c; ++i) 
    {
    if (i== a*b*c -1 )
      plumedFile << "csa" << parameters.myPlumedSystemParameters.residueNames[res] << i << "\n";
    else
      plumedFile << "csa" << parameters.myPlumedSystemParameters.residueNames[res] << i << ",";
    }  
  }
return plumedFile; 
}













