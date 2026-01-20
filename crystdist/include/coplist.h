#ifndef H_COPLIST
#define H_COPLIST

#include "common/include/iofile.h"
#include "crystdist/include/crystalops.h"


class COPList 
{

public:
  enum OP_TYPE = {DOPS=0, BOPS=1, ROPS=2, IOPS=3, LOPS=4, TOPS=5};
  COPList();
  COPList(std::string const &fileName, IOMode const &mode);

  friend COPList &operator<<(COPList &copFile, CrystalOrderParameters const &parameters); // Write ops to file
  friend COPList &operator>>(COPList &copFile, CrystalOrderParameters &parameters);       // Read ops from file


};




#endif
