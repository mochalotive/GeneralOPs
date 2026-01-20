#include "common/include/iofile.h"
#include <string>
#include "plumedinput/include/kernelparameters.h"


//#include "mymol/include/file_formats/fileformats.h"
//#include "mymol/include/system.h"
//#include "crystdist/include/pointmolecule.h"
//#include "crystdist/include/mmpfile.h"
//#include "crystdist/include/xtpfile.h"
//#include "crystdist/include/copfile.h"
//#include "crystdist/include/moleculemap.h"
//#include "mymol/include/geometry.h"





#ifndef H_KERNELFILE
#define H_KERNELFILE 


class KernelFile : public IOFile {
private:
public:
KernelFile();
KernelFile(std::string const &fileName);

void setFile(std::string const &fileName);
friend KernelFile &operator<<(KernelFile &kernelFile, KernelParameters &parameters);
};

//inline setfile


inline void KernelFile::setFile(std::string const &fileName) 
{
IOFile::setFile(fileName, OUT);
}


#endif
