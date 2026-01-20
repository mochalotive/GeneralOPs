/*
** Copyright 2009-2011 Erik Santiso.
** This file is part of smcv.
** smcv is free software: you can redistribute it and/or modify
** it under the terms of the GNU Lesser General Public License
** version 2.1 as published by the Free Software Foundation.
** 
**
** smcv is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU Lesser General Public License for more details.
**
** You should have received a copy of the GNU Lesser General Public License
** along with smcv. If not, see <http://www.gnu.org/licenses/>.
*/

/*
** Force data - stores force constants, restraint forces and m-matrix
*/

#ifndef H_FORCE_DATA
#define H_FORCE_DATA

#include <vector>
#include "common/include/types.h"
#include "common/include/assert.h"
#include "mymath/include/vectorND.h"
#include "mymath/include/matrixND.h"

struct ForceData
{
  std::vector<Real> forceConstants; // Force constants 
  VectorND restraintForces;         // Restraint forces
  MatrixND mMatrix;                 // M-matrix

  inline void clear() // Clear all data
  {
    forceConstants.clear();
    restraintForces.clear();
    mMatrix.clear();
  }

  inline void resize(size_t const numParameters)  // Change number of parameters
  {
    forceConstants.resize(numParameters, 0);
    restraintForces.resize(numParameters, 0.0);
    mMatrix.resize(numParameters);
  }

  inline size_t const size() const // Return total number of parameters
  {
    return restraintForces.size();
  }

  inline Real const &operator[](size_t const i) const // Access force components by index
  {
    assert(i < restraintForces.size());
    return restraintForces[i]; 
  }

  inline friend std::ostream& operator<<(std::ostream &str, ForceData const &data) // Output
  {
    size_t const numParameters = data.size();
    str << std::endl << "Number of parameters: " << numParameters;
    str << std::endl << "Force constants: ";
    for(size_t i = 0; i < data.forceConstants.size(); ++i)
      str << std::endl << data.forceConstants[i];
    str << std::endl << "Restraint forces: ";
    for(size_t i = 0; i < numParameters; ++i)
      str << std::endl << data.restraintForces[i];
    str << std::endl << "M-matrix: ";
    for(size_t i = 0; i < numParameters; ++i)
      for(size_t j = 0; j < numParameters; ++j)
        if(data.mMatrix(i, j) != 0.0)
          str << std::endl << i << " " << j << " " << data.mMatrix(i, j);
    return str;
  }
};

/*
** End of class ForceData
*/

#endif

