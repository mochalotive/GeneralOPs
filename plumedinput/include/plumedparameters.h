#ifndef H_PLUMEDPARAMETERS
#define H_PLUMEDPARAMETERS


#include "mymol/include/system.h"
#include "mymol/include/geometry.h"
#include "crystdist/include/crystalops.h"
#include <vector>
#include <string>
#include "crystdist/include/moleculemap.h"
#include "crystdist/include/statparameters.h"
#include <algorithm>
#include <iostream>
#include "plumedinput/include/equationcontainer.h"
//a class to store lots of parameters for a plumed input file
//could be a struct possibly
//added by Jake McKibbin march 2024



//not sure how im gonna do this later on
//enum class UnitName {LENGTH, ENERGY, TIME, CHARGE, MASS};
//struct Units {
//std::string LENGTH = "nm";
//std::string ENERGY = "kj/mol";
//std::string TIME = "ps";
//std::string CHARGE = "e";
//std::string MASS = "amu";
//};
enum PlumedFileType {DOPS=0, ROPS=1, BOPS=2, LOPS=3, NOPS=4, TOPS=5}; 
std::ostream &operator << ( std::ostream& strm, PlumedFileType tt );

struct SystemDataPackage
  {
  System<Geometry>& mySystem;//pdb data, only one passed
  std::vector<CrystalDistributionParameters>& myCrystalDistributionParameters;//xtp datas
  std::vector<std::vector<MoleculeMap> >& myMoleculeMaps;//mmp datas
  std::vector<CrystalOrderParameters>& myCrystalOrderParameters;//cop files
  std::vector<CrystalOrderParameters>& myOldCrystalOrderParameters;//old files for moving restraint
  std::vector<EquationContainer>& myEq;
  std::pair<int,int> steps;
  SystemDataPackage(System<Geometry>& mySystem_,
                    std::vector<CrystalDistributionParameters>& myCrystalDistributionParameters_,
                    std::vector<std::vector<MoleculeMap> >& myMoleculeMaps_,
                    std::vector<CrystalOrderParameters>& myCrystalOrderParameters_,
                    std::vector<CrystalOrderParameters>& myOldCrystalOrderParameters_i,
                    std::vector<EquationContainer>& myEq_,
                    std::pair<int,int> steps_);
  };
//all these classes exist elsewhere, so I'm confortable doing references
//they should all outlive instances of this class


struct PlumedSystemParameters
  {
  std::vector<Real> box_size; //box size
  std::vector<bool> COM; //are we using COM or frame
  std::vector<std::string> residueNames;
  std::vector<Real> cutoff;
  std::vector<bool> contact;
  std::vector<Real> contact_cutoff;
  std::vector<Real> kappa;
  PlumedSystemParameters(std::vector<bool> COM_,
                         std::vector<Real> box_size_, 
                         std::vector<std::string> residueNames_,
                         std::vector<Real> cutoff_,
                         std::vector<bool> contact_,
                         std::vector<Real> contact_cutoff_,
                         std::vector<Real> kappa_);
  };
//will pass most of this as arguments, so not making them references



class PlumedParameters {//no inheritence needed I think
private:
public:
//calculated within constructor
std::vector<size_t> residueNumbers; //made in constructor
std::vector<std::string> allRes; //made in constructor

//taken as arguments
std::vector<PlumedFileType> fileType;
SystemDataPackage& mySystemDataPackage; 
PlumedSystemParameters& myPlumedSystemParameters;
size_t op; //for usage in computing forces with ForceFile
size_t resid;


int compareInternal(std::string s, std::string t);//just make two damn copies
//PlumedParameters(bool COM_, std::vector<int> resNum_, std::vector<int> atomPerMolecule_, std::vector<size_t> frame_, std::vector<int> centerOfBox_, std::string kernelFile_)
PlumedParameters(std::vector<PlumedFileType> fileType_, SystemDataPackage& mySystemDataPackage_, PlumedSystemParameters& myPlumedSystemParameters_, size_t op_= 0, size_t resid_ = 0);
PlumedParameters(std::vector<std::string> fileType_, SystemDataPackage& mySystemDataPackage_, PlumedSystemParameters& myPlumedSystemParameters_ , size_t op_ = 0, size_t resid_= 0);

};
#endif
