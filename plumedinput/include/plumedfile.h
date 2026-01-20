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


#ifndef H_PLUMEDFILE
#define H_PLUMEDFILE

//produces an input file for use in plumed. Requires a PDB file for molecule information, and mmp,xtp, and cop file for setting parameters in input
//output only for the moment, no real information to gain from a plumed input at the moment 


//struct DataPackage{
//CrystalDistributionParameters const &distparameters;
//std::vector<MoleculeMap> &moleculeMaps;
//CrystalOrderParameters &parameters; 
//System<Geometry> &system;
//};


class PlumedFile : public IOFile {
//constructors 
public:
PlumedFile();  //empty plumed input file
PlumedFile(std::string const &fileName); //plumed input with given filename

//interface 
void setFile(std::string const &fileName );

friend PlumedFile &operator<<(PlumedFile &plumedFile, PlumedParameters &parameters); //plumed, xtp, mmp, cop. pdb
// Write parameters to file

};
//end plumed input class
//inline for setFile, only mode is out because I'm not taking any in 
inline void PlumedFile::setFile(std::string const &fileName)
{
  IOFile::setFile(fileName, OUT);
}
#endif
