#include "plumedinput/include/fileconverters.h"
namespace FileINConverters {
void XTPConverterIN(std::vector<std::string>&s, void* data)
{
//shortcut to read in all the xtp files and make crystalop objects without nonsense
  std::vector<CrystalDistributionParameters>* myCrystalParameters = static_cast<std::vector<CrystalDistributionParameters>*>(data);
  for (size_t i=0; i<s.size(); i++) {
  XTPFile tempFile(s[i],IN);
  //error checking
  CrystalDistributionParameters tempCrystalParameters;
  tempFile >> tempCrystalParameters;
  if((!tempFile && !tempFile.eof()) ||
     tempCrystalParameters.concentrationParameters.size() < 1 ||
     tempCrystalParameters.means.size() < 1 ||
     tempCrystalParameters.normalizationFactors.size() < 1)
  {
    std::cerr << "Error reading crystal distribution parameters file " << s[i] << std::endl;
 //             << "Use copcalc -help for help." << std::endl;
    throw std::exception();
  }
  myCrystalParameters->push_back(tempCrystalParameters);
  tempFile.close();
  }
  return;
}

void MMPConverterIN(std::vector<std::string>&s, void* data)
{
  std::vector <std::vector<MoleculeMap> >* myMoleculeMaps = static_cast<std::vector< std::vector<MoleculeMap > > *>(data);
  for (size_t i=0; i<s.size(); i++)
  {
   std::vector<MoleculeMap> tempMaps;
   MMPFile myTempMMPFile(s[i]); //not readwrite mode, only in
   myTempMMPFile >> tempMaps;

   if((!myTempMMPFile && !myTempMMPFile.eof()) ||
     tempMaps.size() < 1 ||
     (tempMaps.size() > 0 && tempMaps[0].framePoints.size() < 1))
  {
    std::cerr << "Error reading molecule map file " << s[i] << std::endl;
 //             << "Use copcalc -help for help." << std::endl``;
    throw std::exception();
  }
  myMoleculeMaps->push_back(tempMaps);
  myTempMMPFile.close();

 
 
  }
  return;
}
void COPConverterIN(std::vector<std::string>&s, void* data)
{
  //hopefully somewhere else checks this somewhere :)
  std::vector <CrystalOrderParameters >* myCrystalOrderParameters = static_cast<std::vector<CrystalOrderParameters>*  >(data);
  for (size_t i=0; i<s.size(); i++)
  {
  CrystalOrderParameters myTempParameters;
  COPFile myTempCOPFile(s[i], IN);
  myTempCOPFile >> myTempParameters;
  myCrystalOrderParameters->push_back(myTempParameters);
  myTempCOPFile.close(); 
  }
  return;

}

void PDBConverterIN(std::vector<std::string>&s, void* data)
{
std::vector<System<Geometry> >* mySystem = static_cast<std::vector<System<Geometry> > * >(data);
  for (size_t i=0; i<s.size(); i++)
  {
  PDBFile myPDBFile(s[i], IN);
  System<Geometry> temp;
  myPDBFile >> temp;
  if((!myPDBFile && !myPDBFile.eof()) ||
     temp.size() < 1 ||
     (temp.size() > 0 && temp[0].numAtoms() == 0))
  {
    std::cerr << "Error reading geometry file " << s[i] << std::endl;
    throw std::exception();
  }
  mySystem->push_back(temp);
  myPDBFile.close();
  }
  return;
}

void EQFConverterIN(std::vector<std::string>&s, void* data)
{
std::vector<EquationContainer>* myEq = static_cast<std::vector<EquationContainer> * >(data);
for (size_t i=0; i< s.size(); ++i)
  {
  EquationFile myEqFile(s[i]);
  EquationContainer temp;
  myEqFile >> temp;
  if (temp.getSize() < 1) std::cerr << "zero size in equation container. Error Reading file " << s[i] << std::endl;
  myEq->push_back(temp);
  myEqFile.close();
  }
  return;
}



}

REGISTER_CLASS(CrystalDistrubtionParameters, "CRYSTAL_DISTRIBUTION_PARAMETERS",&FileINConverters::XTPConverterIN);
REGISTER_CLASS(MoleculeMaps, "MOLECULE_MAPS",&FileINConverters::MMPConverterIN);
REGISTER_CLASS(CrystalOrderParameters, "CRYSTAL_ORDER_PARAMETERS",&FileINConverters::COPConverterIN);
REGISTER_CLASS(SystemGeometry, "SYSTEM_GEOMETRY",&FileINConverters::PDBConverterIN);
REGISTER_CLASS(EquationContainer, "EQUATION_CONTAINER",&FileINConverters::EQFConverterIN);

