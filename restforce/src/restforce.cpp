/*
 * p
** Copyright 2009-2011 Erik Santiso.
** This file is part of restforce.
** restforce is free software: you can redistribute it and/or modify
** it under the terms of the GNU Lesser General Public License
** version 2.1 as published by the Free Software Foundation.
** 
**
** restforce is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU Lesser General Public License for more details.
**
** You should have received a copy of the GNU Lesser General Public License
** along with restforce. If not, see <http://www.gnu.org/licenses/>.
*/

/*
** restforce version 0.1
**
** Functions to implement the string method in collective variables
** using crystal order parameters as defined in the crystdist library.
*/

#include <iostream>
#include <string>
#include <vector>
#include "crystdist/include/crystalops.h"
#include "crystdist/include/copfile.h"
#include "smcv/include/forcedata.h"
#include "smcv/include/opdfile.h"

void printHelpMessage()
{
  std::cout << std::endl
            << "restforce - Calculate average restraint forces and metric tensors" << std::endl
            << "            from a restrained MD simulation using crystal OPs" << std::endl << std::endl; 
  std::cout << "Usage: restforce [-help] " << std::endl
            << "                 [-avg numSkip] " << std::endl
            << "                 [-blk numPoints] " << std::endl
            << "                 [-opd opd_filename] " << std::endl
            << "                 [-cop cop_filename] " << std::endl
            << "                 [-out output_prefix] " << std::endl << std::endl
            << "Description: restforce can be used to compute average restraint forces" << std::endl
            << "             and the average M-matrix from a restrained dynamics" << std::endl
            << "             simulation of a molecular crystal." << std::endl << std::endl
            << "Options: " << std::endl
            << "         -help              - Print this help message" << std::endl
            << "         -avg numSkip       - Compute average restraint forces" << std::endl
            << "                              and average M-matrix. Skip the first" << std::endl
            << "                              numSkip frames" << std::endl
            << "         -blk numPoints     - Calculate block averages of the" << std::endl
            << "                              restraint forces and M-matrix. Each" << std::endl
            << "                              block contains numPoints frames" << std::endl
            << "         -opd opd_filename  - Output file from the restrained dynamics" << std::endl
            << "                              simulation" << std::endl
            << "         -cop cop_filename  - File containing target values of" << std::endl
            << "                              the crystal order parameters" << std::endl
            << "         -out output_prefix - Prefix for the output file names" << std::endl << std::endl
            << "         Exactly one of the -help, -avg or -blk options must be given." << std::endl
            << "         If the -out option is omitted, output is written directly" << std::endl
            << "         to standard output." << std::endl
            << "         The output file cannot have the same name as the input opd file." << std::endl
            << "         The -avg and -blk options require specifying an opd file using" << std::endl
            << "         the -opd option. If the -cop option is also given, average" << std::endl
            << "         values of the order parameters are computed from the" << std::endl
            << "         restraint forces." << std::endl << std::endl;
}

void doAverages(size_t numSkip,
                std::string opdFileName,
                std::string copFileName,
                std::string outPrefix)
{
  // Calculate averages, skipping numSkip frames
  OPDFile myOPDInputFile(opdFileName, IN);
  COPFile myCOPFile;
  CrystalOrderParameters myCrystalOPs, avgOPs;
  size_t numOPsFromCOPFile;
  if(copFileName != "") 
  {
    myCOPFile.setFile(copFileName, IN);
    if(!myCOPFile)
    {
      std::cerr << "Error opening file " << copFileName << std::endl
                << "Exiting..." << std::endl;
      return;
    }
    myCOPFile >> myCrystalOPs;
    size_t numInternalDOFs = myCrystalOPs.internal[0].size();
    numOPsFromCOPFile = myCrystalOPs.grid.numCells()*(4 + numInternalDOFs);
    avgOPs.resize(myCrystalOPs.grid.numCells(), numInternalDOFs);
    avgOPs.grid = myCrystalOPs.grid;
    avgOPs.cutoffSq = myCrystalOPs.cutoffSq;
    avgOPs.switchWidth = myCrystalOPs.switchWidth;
  }
 
  // Skip numSkip frames
  for(size_t iFrame = 0; iFrame < numSkip; ++iFrame)
  {
    ForceData dummyData;
    myOPDInputFile >> dummyData;
    if(myOPDInputFile.eof())
    {
      std::cerr << "Error while skipping " << numSkip << " frames from file " << opdFileName << ": " << std::endl
                << "Not enough frames in file." << std::endl;
      return;
    }
  }
  // Calculate averages
  ForceData avgData;
  //std::cout << "DEBUG PRINTING:   " << std::endl;
  //std::cout << avgData << std::endl;
  myOPDInputFile >> avgData;
  size_t numParameters = avgData.size();
//  std::cout << "restforce numParameters " << numParameters << std::endl;
  size_t numPointsRead = 1;
  while(!myOPDInputFile.eof())
  {
    ForceData currentData;
    myOPDInputFile >> currentData;
    if(myOPDInputFile.eof()) break; // Done reading
    if(currentData.size() != avgData.size())
    {
      std::cerr << "Inconsistent data in file " << opdFileName << std::endl;
      std::cerr << "Frames read: " << numPointsRead << std::endl;
      break;
    }
//    std::vector<Real> &forceConstants = currentData.forceConstants;
    std::vector<Real> forceConstants = currentData.forceConstants;
    forceConstants[2] = 10;
    forceConstants[3] = 10;

    
    if(copFileName != "") // Calculate OPs from restraint forces
    {
      if(numOPsFromCOPFile != avgData.size())
      {
        std::cerr << "Data in COP file " << myCOPFile.fileName() 
                  << " is not consistent with opd file " << myOPDInputFile.fileName() << std::endl
                  << numOPsFromCOPFile << " != " << avgData.size() << std::endl;
        return;
      }
      //JAKE adding fix for floating point negatives
      size_t numCells = myCrystalOPs.grid.numCells();
//      std::cout << numCells << std::endl; 
      for(size_t i = 0; i < numCells; ++i)
      {
        double temp = currentData.restraintForces[i + numCells]/forceConstants[1]
                                     + myCrystalOPs.bondOrientation[i];
        avgOPs.distance[i] += -currentData.restraintForces[i]/forceConstants[0]
                              + myCrystalOPs.distance[i];
        avgOPs.bondOrientation[i] += -currentData.restraintForces[i + numCells]/forceConstants[1]
                                     + myCrystalOPs.bondOrientation[i];
        avgOPs.relativeOrientation[i] += -currentData.restraintForces[i + 2*numCells]/forceConstants[2]
                                         + myCrystalOPs.relativeOrientation[i];
        //JAKE mar 2024
        avgOPs.localDensity[i] += -currentData.restraintForces[i+3*numCells]/forceConstants[3] + myCrystalOPs.localDensity[i];

        for(size_t j = 0; j < myCrystalOPs.internal[0].size(); ++j)
        {
          avgOPs.internal[i][j] += -currentData.restraintForces[i + (j + 4)*numCells]/forceConstants[j]
                                   + myCrystalOPs.internal[i][j];
        }
      }
    }
    for(size_t i = 0; i < numParameters; ++i)
    {
      avgData.restraintForces[i] += currentData.restraintForces[i];
      for(size_t j = 0; j < numParameters; ++j)
        avgData.mMatrix(i, j) += currentData.mMatrix(i, j);
    }
    ++numPointsRead;
  }
  for(size_t i = 0; i < numParameters; ++i)
  {
    avgData.restraintForces[i] /= (Real)numPointsRead;
    for(size_t j = 0; j < numParameters; ++j)
      avgData.mMatrix(i, j) /= (Real)numPointsRead;
  }
  if(copFileName != "")
  {
    size_t numCells = avgOPs.grid.numCells();
    for(size_t i = 0; i < numCells; ++i)
    {
      avgOPs.distance[i] /= (Real)numPointsRead;
      avgOPs.bondOrientation[i] /= (Real)numPointsRead;
      avgOPs.relativeOrientation[i] /= (Real)numPointsRead;
      avgOPs.localDensity[i] /=(Real)numPointsRead;
      for(size_t j = 0; j < avgOPs.internal[0].size(); ++j)
        avgOPs.internal[i][j] /= (Real)numPointsRead;
    }
  }
  if(outPrefix != "")
  {
    // Write to output file
    OPDFile myOPDOutputFile(outPrefix+".opd", OUT);
    myOPDOutputFile << avgData;
    if(copFileName != "")
    {
      COPFile copOutFile(outPrefix+"_OPs.cop", OUT);
      copOutFile << avgOPs;
      copOutFile.close();
    }
//       std::cout << "Average restraint forces: " << std::endl;
//    for(size_t i = 0; i < avgData.restraintForces.size(); ++i)
//      std::cout << avgData.restraintForces[i] << std::endl;
//    for(size_t i = 0; i < avgData.mMatrix.size(); ++i)
//      for(size_t j = i; j < avgData.mMatrix.size(); ++j)
//        if(avgData.mMatrix[i][j] != 0.0)
//          std::cout << i << " " << j << " " << avgData.mMatrix[i][j] << std::endl;
   }
  else
  {
    // Write to standard output (useful for debugging)
    std::cout << "Average restraint forces: " << std::endl;
    for(size_t i = 0; i < avgData.restraintForces.size(); ++i)
      std::cout << avgData.restraintForces[i] << std::endl;
    for(size_t i = 0; i < avgData.mMatrix.size(); ++i)
      for(size_t j = i; j < avgData.mMatrix.size(); ++j)
        if(avgData.mMatrix[i][j] != 0.0)
          std::cout << i << " " << j << " " << avgData.mMatrix[i][j] << std::endl;
    if(copFileName != "")
      std::cout << avgOPs;
  }
}

void doBlocks(size_t numBlock,
              std::string opdFileName,
              std::string copFileName,
              std::string outPrefix)
{
  // Do block averages with size numBlock
  if(!numBlock) // Stupidity check
  {
    std::cerr << "Number of blocks must be > 0." << std::endl;
    return;
  }
  OPDFile myOPDInputFile(opdFileName, IN);
  COPFile myCOPFile, copOutFile;
  OPDFile myOPDOutputFile;

  CrystalOrderParameters myCrystalOPs;
  size_t numOPsFromCOPFile;
  if(copFileName != "") 
  {
    myCOPFile.setFile(copFileName, IN);
    if(!myCOPFile)
    {
      std::cerr << "Error opening file " << copFileName << std::endl
                << "Exiting..." << std::endl;
      return;
    }
    myCOPFile >> myCrystalOPs;
  }

  if(outPrefix != "")
  {
    myOPDOutputFile.setFile(outPrefix+".opd", OUT);  
    if(copFileName != "")
      copOutFile.setFile(outPrefix+"_OPs.cop", OUT);
  }
  while(!myOPDInputFile.eof())
  {
    ForceData avgData;
    CrystalOrderParameters avgOPs;
    if(copFileName != "")
    {
      size_t numInternalDOFs = myCrystalOPs.internal[0].size();
      numOPsFromCOPFile = myCrystalOPs.grid.numCells()*(4 + numInternalDOFs);
      avgOPs.resize(myCrystalOPs.grid.numCells(), numInternalDOFs);
      avgOPs.grid = myCrystalOPs.grid;
      avgOPs.cutoffSq = myCrystalOPs.cutoffSq;
      avgOPs.switchWidth = myCrystalOPs.switchWidth;
    }
    size_t numPointsRead = 0;
    size_t numParameters;

    for(size_t i = 0; i < numBlock; ++i)
    {
      ForceData currentData;
      myOPDInputFile >> currentData;
      std::vector<Real> &forceConstants = currentData.forceConstants; 
      numParameters = currentData.size();
      if(i == 0) // Initialize accumulators
      { 
        avgData.clear();
        avgData.resize(numParameters); 
      }
      if(myOPDInputFile.eof()) break;
      if(copFileName != "") // Calculate OPs from restraint forces
      {
        if(numOPsFromCOPFile != numParameters)
        {
          std::cerr << "Data in COP file " << myCOPFile.fileName() 
                    << " is not consistent with opd file " << myOPDInputFile.fileName() << std::endl;
          return;
        }
        size_t numCells = myCrystalOPs.grid.numCells();
        for(size_t i = 0; i < numCells; ++i)
        {
          avgOPs.distance[i] += -currentData.restraintForces[i]/forceConstants[0]
                                + myCrystalOPs.distance[i];
          avgOPs.bondOrientation[i] += -currentData.restraintForces[i + numCells]/forceConstants[1]
                                      + myCrystalOPs.bondOrientation[i];
          avgOPs.relativeOrientation[i] += -currentData.restraintForces[i + 2*numCells]/forceConstants[2]
                                          + myCrystalOPs.relativeOrientation[i];
          for(size_t j = 0; j < myCrystalOPs.internal[0].size(); ++j)
          {
            avgOPs.internal[i][j] += -currentData.restraintForces[i + (j + 3)*numCells]/forceConstants[j]
                                    + myCrystalOPs.internal[i][j];
          }
        }
      }
      for(size_t i = 0; i < numParameters; ++i)
      {
        avgData.restraintForces[i] += currentData.restraintForces[i];
        for(size_t j = 0; j < numParameters; ++j)
          avgData.mMatrix(i, j) += currentData.mMatrix(i, j);
      }
      ++numPointsRead;
    } // End of for loop for current block
    if(myOPDInputFile.eof()) break;

    for(size_t i = 0; i < numParameters; ++i)
    {
      avgData.restraintForces[i] /= (Real)numPointsRead;
      for(size_t j = 0; j < numParameters; ++j)
        avgData.mMatrix(i, j) /= (Real)numPointsRead;
    }
    if(copFileName != "")
    {
      size_t numCells = avgOPs.grid.numCells();
      for(size_t i = 0; i < numCells; ++i)
      {
        avgOPs.distance[i] /= (Real)numPointsRead;
        avgOPs.bondOrientation[i] /= (Real)numPointsRead;
        avgOPs.relativeOrientation[i] /= (Real)numPointsRead;
        for(size_t j = 0; j < avgOPs.internal[0].size(); ++j)
          avgOPs.internal[i][j] /= (Real)numPointsRead;
      }
    }
    if(outPrefix != "")
    {
      // Write to output files
      myOPDOutputFile << avgData;
      if(copFileName != "")
        copOutFile << avgOPs;
      // Write to standard output
      std::cout << "Average restraint forces: " << std::endl;
      for(size_t i = 0; i < avgData.restraintForces.size(); ++i)
        std::cout << avgData.restraintForces[i] << std::endl;
      for(size_t i = 0; i < avgData.mMatrix.size(); ++i)
        for(size_t j = i; j < avgData.mMatrix.size(); ++j)
          if(avgData.mMatrix[i][j] != 0.0)
            std::cout << i << " " << j << " " << avgData.mMatrix[i][j] << std::endl;
      if(copFileName != "")
        std::cout << avgOPs;



    }
    else
    {
      // Write to standard output
      std::cout << "Average restraint forces: " << std::endl;
      for(size_t i = 0; i < avgData.restraintForces.size(); ++i)
        std::cout << avgData.restraintForces[i] << std::endl;
      for(size_t i = 0; i < avgData.mMatrix.size(); ++i)
        for(size_t j = i; j < avgData.mMatrix.size(); ++j)
          if(avgData.mMatrix[i][j] != 0.0)
            std::cout << i << " " << j << " " << avgData.mMatrix[i][j] << std::endl;
      if(copFileName != "")
        std::cout << avgOPs;
    }
  } // End of main loop to read opd file
}

int main(int argc, char* argv[])
{
  enum CommandLineOption { HELP, AVERAGE, BLOCK };

  if(argc < 2)
  {
    // No input options, print help message
    printHelpMessage();
    return 0;
  }
  // Copy argument list to a vector of strings
  std::vector<std::string> arguments;
  for(int i = 0; i < argc; ++i)
    arguments.push_back(std::string(argv[i]));
  // Search for command line options
  bool foundInputOption = false;
  bool foundOutputOption = false;
  bool foundOPDFileOption = false;
  bool foundCopFileOption = false;
  size_t avgOptionIndex = 0;
  size_t blkOptionIndex = 0;
  size_t outputOptionIndex = 0;
  size_t opdOptionIndex = 0;
  size_t copOptionIndex = 0;
  CommandLineOption inputOption = HELP;
  for(size_t i = 0; i < arguments.size(); ++i)
  {
    std::string &argument = arguments[i];
    size_t switchIndex = argument.find("-");
    if(switchIndex != argument.npos)
    {
      // Found a command-line switch
      if(foundInputOption && (argument == "-help" || argument == "avg" || argument == "blk"))
      {
        std::cerr << "Invalid option combination." << std::endl
                  << "Use scmv -help for help." << std::endl;
        return -1;
      }
      if(argument == "-help")
      {
        foundInputOption = true;
        inputOption = HELP;
      }
      else if(argument == "-avg")
      {
        if(foundInputOption)
        {
          std::cerr << "Only one of the -avg and -blk options can be specified." << std::endl
                    << "Use restforce -help for help." << std::endl;
          return -1;
        }
        foundInputOption = true;
        avgOptionIndex = i;
        inputOption = AVERAGE;
      }
      else if(argument == "-blk")
      {
        if(foundInputOption)
        {
          std::cerr << "Only one of the -avg and -blk options can be specified." << std::endl
                    << "Use restforce -help for help." << std::endl;
          return -1;
        }
        foundInputOption = true;
        blkOptionIndex = i;
        inputOption = BLOCK;
      }
      else if(argument == "-opd")
      {
        foundOPDFileOption = true;
        opdOptionIndex = i;
      }
      else if(argument == "-cop")
      {
        foundCopFileOption = true;
        copOptionIndex = i;
      }
      else if(argument == "-out")
      {
        foundOutputOption = true;
        outputOptionIndex = i;
      }
      else
      {
        std::cerr << "Invalid option: " << argument << "." << std::endl
                  << "Use restforce -help for help." << std::endl;
        return -1;
      }
    }
  }
  switch(inputOption)
  {
  case HELP:
    printHelpMessage();
    return 0;
  case AVERAGE:
    {
      std::string numSkipString = arguments[avgOptionIndex + 1];
      if(numSkipString.find("-") != numSkipString.npos)
      {
        std::cerr << "Number of frames to skip must be given after -avg option." << std::endl
                  << "Use restforce -help for help." << std::endl;
        return -1;
      }
      size_t numSkip = atoi(numSkipString.c_str());
      std::string opdFileName, copFileName, outPrefix;
      if(!foundOPDFileOption)
      {
        std::cerr << "opd file name is required for -avg option." << std::endl
                  << "Use restforce -help for help." << std::endl;
        return -1;
      }
      if(opdOptionIndex + 1 >= arguments.size())
      {
        std::cerr << "opd file name should be given after -opd option." << std::endl
                  << "Use restforce -help for help." << std::endl;
        return -1;
      }
      opdFileName = arguments[opdOptionIndex + 1];
      if(opdFileName.find("-") != opdFileName.npos)
      {
        std::cerr << "opd file name should be given after -opd option." << std::endl
                  << "Use restforce -help for help." << std::endl;
        return -1;
      }
      if(foundCopFileOption)
      {
        if(copOptionIndex + 1 >= arguments.size())
        {
          std::cerr << "Target op file name should be given after -cop option." << std::endl
                    << "Use restforce -help for help." << std::endl;
          return -1;
        }
        copFileName = arguments[copOptionIndex + 1];
        if(copFileName.find("-") != copFileName.npos)
        {
          std::cerr << "Target op file name should be given after -cop option." << std::endl
                    << "Use restforce -help for help." << std::endl;
          return -1;
        }
      }
      else copFileName.clear();
      if(foundOutputOption)
      {
        if(outputOptionIndex + 1 >= arguments.size())
        {
          std::cerr << "Output file prefix should be given after -out option." << std::endl
                    << "Use restforce -help for help." << std::endl;
          return -1;
        }
        outPrefix = arguments[outputOptionIndex + 1];
        if(outPrefix.find("-") != outPrefix.npos)
        {
          std::cerr << "Output file prefix should be given after -out option." << std::endl
                    << "Use restforce -help for help." << std::endl;
          return -1;
        }
        if(opdFileName.find(outPrefix+".out") != opdFileName.npos)
        {
          std::cerr << "Output file cannot have the same name as the input opd file." << std::endl
                    << "Use restforce -help for help." << std::endl;
          return -1;
        }
      }
      else outPrefix.clear();
      doAverages(numSkip, opdFileName, copFileName, outPrefix);
    }
    break;
  case BLOCK:
    {
      std::string numBlockString = arguments[blkOptionIndex + 1];
      if(numBlockString.find("-") != numBlockString.npos)
      {
        std::cerr << "Block size must be given after -blk option." << std::endl
                  << "Use restforce -help for help." << std::endl;
        return -1;
      }
      size_t numBlock = atoi(numBlockString.c_str());
      std::string opdFileName, copFileName, outPrefix;
      if(!foundOPDFileOption)
      {
        std::cerr << "opd file name is required for -blk option." << std::endl
                  << "Use restforce -help for help." << std::endl;
        return -1;
      }
      if(opdOptionIndex + 1 >= arguments.size())
      {
        std::cerr << "opd file name should be given after -opd option." << std::endl
                  << "Use restforce -help for help." << std::endl;
        return -1;
      }
      opdFileName = arguments[opdOptionIndex + 1];
      if(opdFileName.find("-") != opdFileName.npos)
      {
        std::cerr << "opd file name should be given after -opd option." << std::endl
                  << "Use restforce -help for help." << std::endl;
        return -1;
      }
      if(foundCopFileOption)
      {
        if(copOptionIndex + 1 >= arguments.size())
        {
          std::cerr << "Target op file name should be given after -cop option." << std::endl
                    << "Use restforce -help for help." << std::endl;
          return -1;
        }
        copFileName = arguments[copOptionIndex + 1];
        if(copFileName.find("-") != copFileName.npos)
        {
          std::cerr << "Target op file name should be given after -cop option." << std::endl
                    << "Use restforce -help for help." << std::endl;
          return -1;
        }
      }
      else copFileName.clear();
      if(foundOutputOption)
      {
        if(outputOptionIndex + 1 >= arguments.size())
        {
          std::cerr << "Output file prefix should be given after -out option." << std::endl
                    << "Use restforce -help for help." << std::endl;
          return -1;
        }
        outPrefix = arguments[outputOptionIndex + 1];
        if(outPrefix.find("-") != outPrefix.npos)
        {
          std::cerr << "Output file prefix should be given after -out option." << std::endl
                    << "Use restforce -help for help." << std::endl;
          return -1;
        }
        if(opdFileName.find(outPrefix+".out") != opdFileName.npos)
        {
          std::cerr << "Output file cannot have the same name as the input opd file." << std::endl
                    << "Use restforce -help for help." << std::endl;
          return -1;
        }
      }
      else outPrefix.clear();
      doBlocks(numBlock, opdFileName, copFileName, outPrefix);
    }
    break;
  }
  return 0;
}

