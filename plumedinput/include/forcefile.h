#include "common/include/iofile.h"
#include <string>
#include "plumedinput/include/plumedparameters.h"


//#include "mymol/include/file_formats/fileformats.h"
//#include "mymol/include/system.h"
//#include "crystdist/include/pointmolecule.h"
//#include "crystdist/include/mmpfile.h"
//#include "crystdist/include/xtpfile.h"
//#include "crystdist/include/copfile.h"
//#include "crystdist/include/moleculemap.h"
//#include "mymol/include/geometry.h"


#ifndef H_PLUMEDFORCEFILE
#define H_PLUMEDFORCEFILE

//nearly a direct copy of the PlumedFile class, because I'm printing nearly the same thing, just replacing 
//the actual restraints with printing instead!

//struct DataPackage{
//CrystalDistributionParameters const &distparameters;
//std::vector<MoleculeMap> &moleculeMaps;
//CrystalOrderParameters &parameters; 
//System<Geometry> &system;
//};


class ForceFile : public IOFile {
//constructors 
public:
ForceFile();  //empty plumed input file
ForceFile(std::string const &fileName); //plumed input with given filename

//interface 
void setFile(std::string const &fileName );

friend ForceFile &operator<<(ForceFile &plumedFile, PlumedParameters &parameters); //plumed, xtp, mmp, cop. pdb
// Write parameters to file

};
//end plumed input class
//inline for setFile, only mode is out because I'm not taking any in 
inline void ForceFile::setFile(std::string const &fileName)
{
  IOFile::setFile(fileName, OUT);
}
#endif
