#include "plumedinput/include/equationfile.h"
EquationFile::EquationFile():IOFile(){};
EquationFile::EquationFile(std::string const &fileName):IOFile(fileName, IN){};

EquationFile &operator>>(EquationFile &equationFile, EquationContainer &container)
{
  assert(equationFile.is_open());
  std::string s;
  while(std::getline(equationFile, s))
    container.addString(s);

  return (equationFile);
}
