/*
** Copyright 2009 Erik Santiso.
** This file is part of smcv.
** smcv is free software: you can redistribute it and/or modify
** it under the terms of the GNU Lesser General Public License as published by
** the Free Software Foundation, either version 3 of the License, or
** (at your option) any later version.
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
** .opd file format
**
** Note: Reading of internal OPs is not implemented yet
*/

#include <sstream>
#include <vector>
#include "common/include/assert.h"
#include "smcv/include/opdfile.h"

// Constructors

OPDFile::OPDFile()
:
IOFile()
{}

OPDFile::OPDFile(std::string const &fileName, IOMode const &mode)
:
IOFile(fileName, mode, BINARY)
{
  if(mode == IN) getForceConstants();
  else firstWrite_ = true; 
}

// Interface

void OPDFile::setFile(std::string const &fileName, IOMode const &mode)
{
  IOFile::setFile(fileName, mode);
  if(!*this) return; // An error ocurred
  if(mode == IN) getForceConstants();
  else firstWrite_ = true;
}

// Operators

OPDFile &operator<<(OPDFile &opdFile, ForceData const &data)
{
  assert(opdFile.is_open() && opdFile.mode() == OUT);
  if(opdFile.firstWrite_) opdFile.writeForceConstants(data.forceConstants);

  size_t numParameters = data.size(); 
  opdFile.write(reinterpret_cast<char *>(&numParameters), sizeof(int));

  // Write restraint forces 
  for(size_t i = 0; i < numParameters; ++i)
  {
    Real force = data.restraintForces[i];
    opdFile.write(reinterpret_cast<char *>(&force), sizeof(Real));
  }

  for(size_t i = 0; i < numParameters; ++i)
  for(size_t j = i; j < numParameters; ++j)
  {
    Real Mij = data.mMatrix(i, j);
    if(Mij != 0.0)
    {
      opdFile.write(reinterpret_cast<char *>(&i), sizeof(size_t));
      opdFile.write(reinterpret_cast<char *>(&j), sizeof(size_t));
      opdFile.write(reinterpret_cast<char *>(&Mij), sizeof(Real));
    }
  }
  int iend = -999;
  opdFile.write(reinterpret_cast<char *>(&iend), sizeof(int));
  opdFile.firstWrite_ = false;
  return opdFile;
}

OPDFile &operator>>(OPDFile &opdFile, ForceData &data)
{
  assert(opdFile.is_open() && opdFile.mode() == IN);
  data.clear();

  // Copy force constants
  data.forceConstants = opdFile.forceConstants_;

  size_t r_numParameters; // Number of parameters read from file
  opdFile.read(reinterpret_cast<char *>(&r_numParameters), sizeof(size_t));
  if(opdFile.eof()) return opdFile; // No more data to read
  if(!r_numParameters || !opdFile)
  {
    std::cerr << "Error reading number of parameters from file " << opdFile.fileName() << std::endl;
    return opdFile;
  }
  // Read restraint forces
  for(size_t i = 0; i < r_numParameters; ++i)
  {
    Real r_force; // Force read from file
    opdFile.read(reinterpret_cast<char *>(&r_force), sizeof(Real));
    if(opdFile.eof())
    {
      std::cerr << "Unexpected end of file while reading file " << opdFile.fileName() << std::endl;
      data.clear();
      return opdFile;
    }
    if(!opdFile)
    {
      std::cerr << "Error reading restraint forces from file " << opdFile.fileName() << std::endl;
      data.clear();
      return opdFile;
    }
    data.restraintForces.push_back(r_force);
  }
  // Read m-matrix
  data.mMatrix.resize(r_numParameters);
  for(;;)
  {
    int r_row, r_col; // Row and column read from file
    Real r_Mij;       // Matrix element read from file
    opdFile.read(reinterpret_cast<char *>(&r_row), sizeof(int));
    if(r_row < 0) break;
    opdFile.read(reinterpret_cast<char *>(&r_col), sizeof(int));
    opdFile.read(reinterpret_cast<char *>(&r_Mij), sizeof(Real));
    if(opdFile.eof()) 
    {
      std::cerr << "Unexpected end of file while reading file " << opdFile.fileName() << std::endl;
      data.clear();
      return opdFile;
    }
    if(!opdFile)
    {
      std::cerr << "Error reading M-matrix from file " << opdFile.fileName() << std::endl;
      data.clear();
      return opdFile;
    }
    data.mMatrix(r_row, r_col) = r_Mij;
    data.mMatrix(r_col, r_row) = r_Mij; // Symmetry
  }
  return opdFile;
}

// Write force constants

void OPDFile::writeForceConstants(std::vector<Real> const &forceConstants)
{
  size_t numConstants = forceConstants.size();
  write(reinterpret_cast<char *>(&numConstants), sizeof(size_t)); 
  for(size_t i = 0; i < numConstants; ++i)
  {
    Real constant = forceConstants[i];
    write(reinterpret_cast<char *>(&constant), sizeof(Real));
  }
  if(!*this)
    std::cerr << "Error writing force constants to file " << fileName() << std::endl;
}

// Read force constants

void OPDFile::getForceConstants()
{
  size_t r_numConstants; // Number of force constants read from file
  read(reinterpret_cast<char *>(&r_numConstants), sizeof(size_t));
  if(!r_numConstants || !*this)
  {
    std::cerr << "Error reading number of force constants from file " << fileName() << std::endl;
    return;
  }
  for(size_t i = 0; i < r_numConstants; ++i)
  {
    Real r_forceConstant; // Force constant read from file
    read(reinterpret_cast<char *>(&r_forceConstant), sizeof(Real));
    if(!*this)
    {
      std::cerr << "Error reading force constants from file " << fileName() << std::endl;
      forceConstants_.clear();
    }
    forceConstants_.push_back(r_forceConstant);
  }
}

