/*
** Copyright 2008-2011 Erik Santiso.
** This file is part of crystdist.
** crystdist is free software: you can redistribute it and/or modify
** it under the terms of the GNU Lesser General Public License
** version 2.1 as published by the Free Software Foundation.
** 
**
** crystdist is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU Lesser General Public License for more details.
**
** You should have received a copy of the GNU Lesser General Public License
** along with crystdist. If not, see <http://www.gnu.org/licenses/>.
*/

/*
** smcv v. 0.1
**
** smcv contains utilities to use the string method in collective variables
** using crystal order parameters.
*/

#include <fstream>
#include "smcv/include/forcedata.h"
#include "smcv/include/opdfile.h"

//#define FIRST_PASS
#define SECOND_PASS

int main(int argc, char* argv[])
{
/*
  ForceData myData;
#ifdef FIRST_PASS
  OPDFile myOPDInputFile("benzene_cry.opd", IN);
  OPDFile myOPDOutputFile("benzene_copy.opd", OUT);
  myOPDInputFile >> myData;
  myData.clear();
#else
  OPDFile myOPDInputFile("benzene_copy.opd", IN);
  OPDFile myOPDOutputFile("benzene_another_copy.opd", OUT);
#endif
  myOPDInputFile >> myData;
  std::ofstream constantsFile("constants.out", std::ios::out);
  std::ofstream forceFile("forces.out", std::ios::out);
  std::ofstream matrixFile("matrices.out", std::ios::out);
  for(size_t i = 0; i < myData.forceConstants.size(); ++i)
    constantsFile << myData.forceConstants[i] << std::endl;
  for(size_t i = 0; i < myData.numParameters(); ++i)
    forceFile << myData.restraintForces[i] << std::endl;
  for(size_t i = 0; i < myData.numParameters(); ++i)
  for(size_t j = i; j < myData.numParameters(); ++j)
    if(myData.mMatrix(i, j) != 0.0) matrixFile << i << " " << j << " " << myData.mMatrix(i, j) << std::endl;
  myOPDOutputFile << myData;
*/
  // TEST
  ForceData myData;
  OPDFile myOPDFile("test.opd", IN);
  myOPDFile >> myData;
  myOPDFile >> myData;
  std::ofstream constantsFile("constants.out", std::ios::out);
  std::ofstream forceFile("forces.out", std::ios::out);
  std::ofstream matrixFile("matrices.out", std::ios::out);
  for(size_t i = 0; i < myData.forceConstants.size(); ++i)
    constantsFile << myData.forceConstants[i] << std::endl;
  for(size_t i = 0; i < myData.numParameters(); ++i)
    forceFile << myData.restraintForces[i] << std::endl;
  for(size_t i = 0; i < myData.numParameters(); ++i)
  for(size_t j = i; j < myData.numParameters(); ++j)
    if(myData.mMatrix(i, j) != 0.0) matrixFile << i << " " << j << " " << myData.mMatrix(i, j) << std::endl;
  return 0;
}
