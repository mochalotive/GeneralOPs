#include "common/include/iofile.h"
#include <string>
#include "plumedinput/include/equationcontainer.h"
#include <cassert>

#ifndef H_EQUATIONFILE
#define H_EQUATIONFILE

//in only class, maybe just include a debug option for equationcontainer if you need it
class EquationFile : public IOFile {
private:
public:
EquationFile();
EquationFile(std::string const &fileName);

void setFile(std::string const &fileName);
friend EquationFile &operator>>(EquationFile &equationFile, EquationContainer &container);
};

//inline setfile

//IN only filetype
inline void EquationFile::setFile(std::string const &fileName) 
{
IOFile::setFile(fileName, IN);
}


#endif
