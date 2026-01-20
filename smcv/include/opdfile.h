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
** .opd file format - used to read and write force constants, restraint
** forces, and the metric tensor.
*/

#ifndef H_OPDFILE
#define H_OPDFILE

#include "common/include/types.h"
#include "common/include/iofile.h"
#include "smcv/include/forcedata.h"

class OPDFile: public IOFile
{
public:

// Constructors

  OPDFile();                                                // Defines an empty OPDFile
  OPDFile(std::string const &fileName, IOMode const &mode); // Defines an OPDFile with a given file name and I/O mode 
                                                                                     // default value is binary to preserve compatibility with old code
                                                                                     // in 99% of cases it'll be a binary opd file we are reading
                                                                                     // JAKE 2025
                                                                                    

// Interface

  void setFile(std::string const &fileName, IOMode const &mode); // Set the file name and I/O mode

// Read from and write to ForceData objects

  friend OPDFile &operator<<(OPDFile &opdFile, ForceData const &data); // Write data to file
  friend OPDFile &operator>>(OPDFile &opdFile, ForceData &data);       // Read data from file

private:

  bool firstWrite_;                  // Whether the file is being written for the first time
  std::vector<Real> forceConstants_; // Force constants from the OPD file header

// Private functions

  void getForceConstants();                                          // Get the force constants from the OPD file header
  void writeForceConstants(std::vector<Real> const &forceConstants); // Write force constants to the OPD file header
};

/*
** End of class OPDFile
*/

#endif

