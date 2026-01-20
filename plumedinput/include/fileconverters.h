#ifndef H_FILE_CONVERTERS
#define H_FILE_CONVERTERS
#include <vector>
#include <string>
#include "mymol/include/geometry.h"
#include "crystdist/include/crystalops.h"
#include "crystdist/include/moleculemap.h"
#include "crystdist/include/statparameters.h"
//
#include "mymol/include/lattice.h"
#include "mymol/include/file_formats/fileformats.h"
#include "mymol/include/system.h"
#include "crystdist/include/mmpfile.h"
#include "crystdist/include/xtpfile.h"
#include "crystdist/include/copfile.h"
//#include "converters/include/convertermacro.h"
#include "plumedinput/include/equationfile.h"

#include "converters/include/converters.h"




namespace FileINConverters {//shortcut for reading in files
////xtp
//CrystalDistributionParameters myCrystalParameters;
//XTPFile myXTPFile(xtpFileName,IN); 
////mmp
//std::vector<MoleculeMap> myMoleculeMaps;
//MMPFile myMMPFile(mmpFileName); //not readwrite mode, only in
////cop
//CrystalOrderParameters myParameters;
//COPFile myCOPFile(copFileName, IN);
////pdb
//System<Geometry> mySystem;
//PDBFile myPDBFile(pdbFileName, IN);

void XTPConverterIN(std::vector<std::string>&s, void* data);
void MMPConverterIN(std::vector<std::string>&s, void* data);
void COPConverterIN(std::vector<std::string>&s, void* data);
void PDBConverterIN(std::vector<std::string>&s, void* data);
void PLUMEDConverterIN(std::vector<std::string>&s, void* data);
void KERNELConverterIN(std::vector<std::string>&s, void* data);
void EQFConverterIN(std::vector<std::string>&s, void* data);
}

#endif
