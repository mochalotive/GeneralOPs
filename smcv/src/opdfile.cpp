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
IOFile(fileName, mode)
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

  int numParameters = (int)data.size();
//  std::cout << "opd numParameters " << numParameters << std::endl; 
  opdFile.write(reinterpret_cast<char *>(&numParameters), sizeof(int));

  // Write restraint forces 
  for(int i = 0; i < numParameters; ++i)
  {
    Real force = data.restraintForces[(size_t)i];
    opdFile.write(reinterpret_cast<char *>(&force), sizeof(Real));
  }

  for(int i = 0; i < numParameters; ++i)
  for(int j = i; j < numParameters; ++j)
  {
    Real Mij = data.mMatrix((size_t)i, (size_t)j);
    if(Mij != 0.0)
    {
      opdFile.write(reinterpret_cast<char *>(&i), sizeof(int));
      opdFile.write(reinterpret_cast<char *>(&j), sizeof(int));
      opdFile.write(reinterpret_cast<char *>(&Mij), sizeof(Real));
    }
  }
  int iend = -999; 
  opdFile.write(reinterpret_cast<char *>(&iend), sizeof(int));
  opdFile.firstWrite_ = false;
 // std::cout << "opd iend " <<iend <<std::endl;
  return opdFile;
}

OPDFile &operator>>(OPDFile &opdFile, ForceData &data)
{
  assert(opdFile.is_open() && opdFile.mode() == IN);
  data.clear();

  // Copy force constants
  data.forceConstants = opdFile.forceConstants_;

  int r_numParameters; // Number of parameters read from file

  opdFile.read(reinterpret_cast<char *>(&r_numParameters), sizeof(int));
//  std::cout << "opd r_numParameters " << r_numParameters <<std::endl; 

if(opdFile.eof()) return opdFile; // No more data to read
  if(!r_numParameters || !opdFile)
  {
    std::cerr << "Error reading number of parameters from file " << opdFile.fileName() << std::endl;
    return opdFile;
  }
  // Read restraint forces
  for(int i = 0; i < r_numParameters; ++i)
  {
    Real r_force; // Force read from file
    opdFile.read(reinterpret_cast<char *>(&r_force), sizeof(Real));
//    std::cout << "opd r_force " << r_force <<std::endl; 
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
//JAKE
//    std::cout << "opd r_cow r_col r_Mij " << r_row << " " <<r_col << " " << r_Mij << std::endl;	 
   if(r_row > r_numParameters || r_col > r_numParameters || r_col < 0)
    { 
      std::cerr << "Unexpected input in file " << opdFile.fileName() << std::endl;
      data.clear();
      return opdFile;
    }
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
  int numConstants = (int)forceConstants.size();
  
  write(reinterpret_cast<char *>(&numConstants), sizeof(int)); 
  for(int i = 0; i < numConstants; ++i)
  {
    Real constant = forceConstants[(size_t)i];
    write(reinterpret_cast<char *>(&constant), sizeof(Real));
  }
  if(!*this)
    std::cerr << "Error writing force constants to file " << fileName() << std::endl;
}

// Read force constants

void OPDFile::getForceConstants()
{
  int r_numConstants; // Number of force constants read from file
  read(reinterpret_cast<char *>(&r_numConstants), sizeof(int));
// std::cout << "opd num Force Constants " << r_numConstants << std::endl; 
  if(!r_numConstants || !*this)
  {
    std::cerr << "Error reading number of force constants from file " << fileName() << std::endl;
    return;
  }
  for(int i = 0; i < r_numConstants; ++i)
  {
    Real r_forceConstant; // Force constant read from file
    read(reinterpret_cast<char *>(&r_forceConstant), sizeof(Real));
//    std::cout << "opd r_forceConstant " << r_forceConstant <<std::endl;
    if(!*this)
    {
      std::cerr << "Error reading force constants from file " << fileName() << std::endl;
      forceConstants_.clear();
    }
    forceConstants_.push_back(r_forceConstant);
  }
}

//JAKE 2025 -> adding read IN for formatted OPD files (easier to make compatible with plumed
//added a polymorph, if there's no format given, nothing is different, so just call the original again. 
//otherwise, read in the formatted one. 
//void OPDFile::getForceConstants(IOFormat &format)
//{
//  if (format == BINARY){ getForceConstants(); return;}
//  //in the case where it's FORMATTED
//  int r_numConstants; // Number of force constants read from file
//  read(reinterpret_cast<char *>(&r_numConstants), sizeof(int));
//// std::cout << "opd num Force Constants " << r_numConstants << std::endl; 
//  if(!r_numConstants || !*this)
//  {
//    std::cerr << "Error reading number of force constants from file " << fileName() << std::endl;
//    return;
//  }
//  for(int i = 0; i < r_numConstants; ++i)
//  {
//    Real r_forceConstant; // Force constant read from file
//    read(reinterpret_cast<char *>(&r_forceConstant), sizeof(Real));
////    std::cout << "opd r_forceConstant " << r_forceConstant <<std::endl;
//    if(!*this)
//    {
//      std::cerr << "Error reading force constants from file " << fileName() << std::endl;
//      forceConstants_.clear();
//    }
//    forceConstants_.push_back(r_forceConstant);
//  }
//}
