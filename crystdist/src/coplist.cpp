#include <cmath>
#include "common/include/assert.h"
#include "crystdist/include/coplist.h"




COPList::COPList()
:
IOFile()
{}

COPList::COPList(std::string const &fileName, IOMode const &mode)
:
IOFile(fileName, mode)
{}


COPList &operator>>(COPList &copList, CrystalOrderParameters &parameters)
{
  // Read order parameters from file
  parameters.clear();

  std::string r_label;  // Label read from file
  //grab number of ops
  do
  {
    std::getline(copList, r_label);
    if(r_label.find("!NUM") != r_label.npos) break;
  }
  while(!copList.eof());
  if(copList.eof())
  {
    std::cerr << "Error reading file " << copList.fileName() << ": NUM record not found." << std::endl;
    return copList;
  }
int numCells;
copList >> numCells;
  if(!copFile)
  {
    std::cerr << "Error reading number of OPs from file " << copFile.fileName() << std::endl;
    return copFile;
  }
if (numCells < 1)
{
std::cerr << "Number of OPs should be a positive integer!" << std::endl;
return copFile;
}
  copFile.ignore(LINE_LENGTH, '\n');  // Skip to next line

  do
  {
    std::getline(copFile, r_label);
    if(r_label.find("!OP") != r_label.npos) break;
  }
  while(!copFile.eof());
  if(copFile.eof())
  {
    std::cerr << "Error reading file " << copFile.fileName() << ": CUTOFF record not found." << std::endl;
    return copFile;
  }
  

