/*
** Copyright 2009-2011 Erik Santiso.
** This file is part of smcvstep.
** smcvstep is free software: you can redistribute it and/or modify
** it under the terms of the GNU Lesser General Public License
** version 2.1 as published by the Free Software Foundation.
** 
**
** smcvstep is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU Lesser General Public License for more details.
**
** You should have received a copy of the GNU Lesser General Public License
** along with smcvstep. If not, see <http://www.gnu.org/licenses/>.
*/

/*
** smcvstep version 0.1
**
** Carry out a step of the (simplified) string method in collective
** variables using crystal order parameters.
**
** See Maragliano et al., J. Chem. Phys. 125, 024106 (2006) and
** E et al., J. Chem. Phys. 126, 164103 (2007) for details
*/

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include "common/include/types.h"
#include "common/include/iofile.h"
#include "mymath/include/integration.h"
#include "mymath/include/interpolation.h"
#include "mymath/include/metrics.h"
#include "crystdist/include/crystalops.h"
#include "crystdist/include/copfile.h"
#include "smcv/include/forcedata.h"
#include "smcv/include/opdfile.h"

void printHelpMessage()
{
  std::cout << std::endl
            << "smcvstep - Carry out a step of the string method in collective variables" << std::endl
            << "           using crystal order parameters." << std::endl << std::endl; 
  std::cout << "Usage: smcvstep [-help] " << std::endl
            << "                [-cop list_filename] " << std::endl
            << "                [-opd list_filename] " << std::endl
            << "                [-bspline]" << std::endl
            << "                [-linear]" << std::endl
            << "                [-tstep time_step] " << std::endl
            << "                [-out output_prefix] " << std::endl
            << "                [-nimg num_out_images] " << std::endl 
            << "                [-fzero zero_force_tol] " << std::endl << std::endl
            << "Description: smcvstep carries out a step of the string method in collective" << std::endl
            << "             variables given a set of files containing the original values" << std::endl
            << "             of the crystal order parameters and a set of files containing" << std::endl
            << "             the average restraing forces and metric tensors from a set of" << std::endl 
            << "             restrained dynamics simulations at those values of order" << std::endl 
            << "             parameters. The output is a set of files contaning the new" << std::endl
            << "             values of the order parameters and a file containing the" << std::endl
            << "             potential of mean force along the original string." << std::endl
            << "             smcvstep also writes to standard output the Frechet distance" << std::endl
            << "             between the old and new curves, which is useful to test for" << std::endl
            << "             convergence of the string method." << std::endl << std::endl
            << "Options: " << std::endl
            << "         -help                   - Print this help message" << std::endl
            << "         -cop cop_list_filename  - File containing a list of the files" << std::endl
            << "                                   containing the original order parameters," << std::endl
            << "                                   one file name per line." << std::endl
            << "         -opd opd_list_filename  - File containing a list of the files" << std::endl
            << "                                   containing the average restraint forces," << std::endl
            << "                                   and metric tensors in opd format, one " << std::endl
            << "                                   file name per line." << std::endl
            << "         -bspline                - Reparameterize using cubic B-spline" << std::endl
            << "                                   interpolation." << std::endl
            << "         -linear                 - Reparameterize using piecewise linear" << std::endl
            << "                                   interpolation." << std::endl
            << "         -tstep time_step        - Time step for the integration of the " << std::endl
            << "                                   string equation of motion." << std::endl
            << "         -out output_prefix      - Prefix used to construct the name of" << std::endl
            << "                                   the output files containing the order" << std::endl
            << "                                   parameters." << std::endl
            << "         -nimg num_out_images    - Number of replicas of the system (i.e." << std::endl 
            << "                                   of order parameter files) to output." << std::endl
            << "         -fzero zero_force_tol   - Zero out any forces and OPs for which" << std::endl
            << "                                   the restraint force is smaller than the" << std::endl
            << "                                   given tolerance. This is useful when" << std::endl
            << "                                   using a reduced set of parameters." << std::endl
            << "         -chords                 - Reparameterize using constant distance" << std::endl
            << "                                   from the initial point instead of" << std::endl
            << "                                   constant arc length." << std::endl << std::endl
            << "         The -cop, -opd, -tstep and -out options are mandatory." << std::endl
            << "         The -bspline and -linear option are mutually exclusive, if none" << std::endl
            << "         is given a B-spline is used by default." << std::endl
            << "         If the -nimg option is not given, the number of output images" << std::endl
            << "         will be the same as the number of input images." << std::endl
            << "         If the number of output images is zero, the program only finds" << std::endl
            << "         the potential of mean force." << std::endl
            << "         If the force tolerance is not given, it is assumed to be zero," << std::endl
            << "         i.e. all order parameters will be used. " << std::endl << std::endl;
}

int main(int argc, char* argv[])
{
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
  bool found_cop = false;
  bool found_opd = false;
  bool found_tstep = false;
  bool found_out = false;
  bool found_nimg = false;  
  bool found_bspline = false;
  bool found_linear = false;
  bool found_zerotol = false;
  bool found_chords = false;
  size_t numOutputImgs = 0; 
  std::string copListFileName, opdListFileName, outPrefix; 
  Real timeStep = 0.0;
  Real zeroTolerance = 0.0;

  for(size_t i = 0; i < arguments.size(); ++i)
  {
    if(arguments[i].find("-") != arguments[i].npos)
    {
      if(arguments[i] == "-help")
      {
        printHelpMessage();
        return 0;
      }
      else if(arguments[i] == "-cop")
      {
        found_cop = true;
        if(i + 1 >= arguments.size() || 
           arguments[i + 1].find("-") != arguments[i + 1].npos)
        {
          std::cerr << "Name of file containing the cop file names must be given after -cop option." << std::endl
                    << "Use smcvstep -help for help." << std::endl;
          return -1;
        }
        else copListFileName = arguments[i + 1];
      } // End of "-cop" option
      else if(arguments[i] == "-opd")
      {
        found_opd = true;
        if(i + 1 >= arguments.size() ||
           arguments[i + 1].find("-") != arguments[i + 1].npos)
        {
          std::cerr << "Name of file containing the opd file names must be given" << std::endl
                    << "after the -opd option." << std::endl
                    << "Use smcvstep -help for help." << std::endl;
          return -1;
        }
        else opdListFileName = arguments[i + 1];
      } // End of "-opd" option
      else if(arguments[i] == "-bspline")
      {
        found_bspline = true;
      } // End of "-bspline" option
      else if(arguments[i] == "-linear")
      {
        found_linear = true;
      } // End of "-linear" option
      else if(arguments[i] == "-tstep")
      {
        found_tstep = true;
        if(i + 1 >= arguments.size() ||
           arguments[i + 1].find("-") != arguments[i + 1].npos)
        {
          std::cerr << "Time step for string integration must be given after" << std::endl
                    << "the -tstep option." << std::endl
                    << "Use smcvstep -help for help." << std::endl;
          return -1;
        }
        // Read time step
        std::stringstream timeStepStream(arguments[i + 1]);
        timeStepStream >> timeStep;
        if((!timeStepStream && !timeStepStream.eof()) ||
           timeStep <= 0.0)
        {
          std::cerr << "Integration time step must be a positive real number." << std::endl
                    << "Use smcvstep -help for help." << std::endl;
          return -1;
        }     
      } // End of "-tstep" option
      else if(arguments[i] == "-out")
      {
        found_out = true;
        if(i + 1 >= arguments.size() ||
           arguments[i + 1].find("-") != arguments[i + 1].npos)
        {
          std::cerr << "Output file prefix must be given after -out option." << std::endl
                    << "Use smcvstep -help for help." << std::endl;
          return -1;
        }
        else outPrefix = arguments[i + 1];
      } // End of "-out" option
      else if(arguments[i] == "-nimg")
      {
        found_nimg = true;
        if(i + 1 >= arguments.size() ||
           arguments[i + 1].find("-") != arguments[i + 1].npos)
        {
          std::cerr << "Number of output replicas must be given after the -nimg option." << std::endl
                    << "Use smcvstep -help for help." << std::endl;
          return -1;
        }
        // Read number of images
        std::stringstream numImgStream(arguments[i + 1]);
        numImgStream >> numOutputImgs;
        if(!numImgStream && !numImgStream.eof())
        {
          std::cerr << "Number of output replicas must be an integer." << std::endl
                    << "Use smcvstep -help for help." << std::endl;
          return -1;
        }
      } // End of "-nimg" option
      else if(arguments[i] == "-fzero")
      {
        found_zerotol = true;
        // Read force tolerance
        std::stringstream zeroTolStream(arguments[i + 1]);
        zeroTolStream >> zeroTolerance;
        if((!zeroTolStream && !zeroTolStream.eof()) ||
           zeroTolerance <= 0.0)
        { 
          std::cerr << "Force tolerance must be a nonnegative real number." << std::endl
                    << "Use smcvstep -help for help." << std::endl;
          return -1;
        }
      } // End of "-fzero" option
      else if(arguments[i] == "-chords")
      {
        found_chords = true;
      }
      else
      {
        std::cerr << "Unknown option: " << arguments[i] << std::endl
                  << "Use smcvstep -help for help." << std::endl;
        return -1;
      }
    }
  } // End of loop over arguments

  if(!found_cop || !found_opd || !found_tstep || !found_out)
  {
    std::cerr << "Error: Need to specify files containing the lists of order parameter files," << std::endl
              << "restraint forces and metric tensors, the time step, and the output file prefix" << std::endl 
              << "Use smcvstep -help for help." << std::endl;
    return -1;
  }
  
  if(found_bspline && found_linear)
  {
    std::cerr << "Error: The -bspline and -linear option are mutually exclusive" << std::endl
              << "Use smcvstep -help for help." << std::endl;
    return -1;
  }

  if(!found_bspline && !found_linear)
    found_bspline = true; // Revert to default

  // Read lists of files

  std::vector<std::string> copFiles;
  std::vector<std::string> opdFiles;

  // COP files
  IOFile copListFile(copListFileName, IN);
  if(!copListFile)
  {
    std::cerr << "Error opening order parameter file list " << copListFileName << std::endl
              << "Exiting..." << std::endl;
    return -1;
  }
  size_t numInputImgs = 0;
  while(!copListFile.eof())
  {
    std::string currentFile;
    copListFile >> currentFile;
    if(currentFile.length() < 1 && copListFile.eof()) break;
    copFiles.push_back(currentFile);
    ++numInputImgs;
  }
  if(numInputImgs < 4)
  {
    std::cerr << "Error: Too few images in order parameter list - need at least 4" << std::endl
              << "Exiting..." << std::endl;
    return -1;
  }
  if(!found_nimg) numOutputImgs = numInputImgs;

  // OPD files 
  IOFile opdListFile(opdListFileName, IN);  
  if(!opdListFile)
  {
    std::cerr << "Error opening opd file list " << opdListFileName << std::endl
              << "Exiting..." << std::endl;
    return -1;
  }
  size_t numOPDImgs = 0;
  while(!opdListFile.eof())
  {
    std::string currentFile;
    opdListFile >> currentFile;
    if(currentFile.length() < 1 && opdListFile.eof()) break;
    opdFiles.push_back(currentFile);
    ++numOPDImgs;
  }
  if(numOPDImgs != numInputImgs)
  {
    std::cerr << "Error: Inconsistent numbers of images in order parameter and opd lists" << std::endl
              << "Exiting..." << std::endl;
    return -1;
  }

  // Read order parameter and force data files

  typedef std::vector<CrystalOrderParameters> COPString;
  typedef std::vector<ForceData> ForceString; 
 
  COPString imageOPs; 
  ForceString imageForceData;
  size_t numParameters;

  for(size_t i = 0; i < numInputImgs; ++i)
  {
    COPFile currentCOPFile(copFiles[i], IN);
    CrystalOrderParameters currentOPs;
    currentCOPFile >> currentOPs;
    if(i == 0) numParameters = currentOPs.size();

    if((!currentCOPFile && !currentCOPFile.eof()) ||
       currentOPs.size() < 1)
    {
      std::cerr << "Error reading order parameter file " << copFiles[i] << std::endl
                << "Exiting..." << std::endl;
      return -1;
    }
    if(currentOPs.size() != numParameters)
    {
      std::cerr << "Inconsistent data in order parameters file" << std::endl
                << "Exiting..." << std::endl;
      return -1;
    }
    currentCOPFile.close();
    imageOPs.push_back(currentOPs);   
 
    OPDFile currentOPDFile(opdFiles[i], IN);
    ForceData currentForceData;
    currentOPDFile >> currentForceData;

    if((!currentOPDFile && !currentOPDFile.eof()) ||
       currentForceData.size() < 1)
    {
      std::cerr << "Error reading force data file " << opdFiles[i] << std::endl
                << "Exiting..." << std::endl;
      return -1;
    }
    if(currentForceData.size() != numParameters)
    {
      std::cerr << "Inconsistent data in force data file" << std::endl
                << "Exiting..." << std::endl;
      return -1;
    }
    currentOPDFile.close(); 
    imageForceData.push_back(currentForceData);
  }

  if(found_zerotol && zeroTolerance > 0.0)  // Apply tolerance rules
    for(size_t i = 0; i < numInputImgs; ++i)
      for(size_t j = 0; j < numParameters; ++j)
        if(fabs(imageForceData[i].restraintForces[j]) < zeroTolerance)
        {
          imageForceData[i].restraintForces[j] = 0.0;
          imageOPs[i][j] = 0.0;
        }

  // Calculate potential of mean force
  std::vector<Real> parameters;
  if(found_chords)
  {
    parameters.push_back(0);
    for(size_t i = 1; i < numInputImgs; ++i)
    {
      Real distSq = 0.0;
      for(size_t j = 0; j < numParameters; ++j)
      {
        Real delta = imageOPs[i][j] - imageOPs[0][j];
        distSq += delta*delta;
      }
      parameters.push_back(sqrt(distSq));
    }
  }
  else
    parameters = trapezoidCumulativeArcLength<COPString, Real>(imageOPs);
  std::vector<Real> cumulativePMF =
    trapezoidCumulativeIntegral<COPString, ForceString, Real>(imageOPs, imageForceData);

  // Write potential of mean force to file

  IOFile myPMFFile(outPrefix+".pmf", OUT);
  for(size_t i = 0; i < numInputImgs; ++i)
    myPMFFile << parameters[i] << " " << cumulativePMF[i] << std::endl;

  // Exit if nimg = 0, give error message if nimg < 4
  if(numOutputImgs == 0) return 0;
  if(numOutputImgs < 4)
  {
    std::cerr << "Number of output images must be greater than 3 " << std::endl
              << "in order to generate new images" << std::endl
              << "Exiting..." << std::endl;
    return -1;
  }

  // Do a forward Euler step to generate new string

  COPString imageOPsAfterEuler = imageOPs;
  for(size_t i = 0; i <  numInputImgs; ++i)
  {
    VectorND EulerStep = timeStep*imageForceData[i].mMatrix*imageForceData[i].restraintForces;
    for(size_t j = 0; j < numParameters; ++j)
      imageOPsAfterEuler[i][j] -= EulerStep[j];
  }

  // Reparameterize to constant arc length or distance to initial point

  std::vector<Real> parametersAfterEuler;

  if(found_chords) // Use constant distance to initial point
  {
    parametersAfterEuler.push_back(0.0);
    for(size_t i = 1; i < numInputImgs; ++i)
    {
      Real distSq = 0.0;
      for(size_t j = 0; j < numParameters; ++j)
      {
        Real delta = imageOPsAfterEuler[i][j] - imageOPsAfterEuler[0][j];
        distSq += delta*delta;
      }
      parametersAfterEuler.push_back(sqrt(distSq));
    }
  }
  else // Use constant arc length
    parametersAfterEuler = trapezoidCumulativeArcLength<COPString, Real>(imageOPsAfterEuler); 

  Real totalLength = parametersAfterEuler[numInputImgs - 1];

  COPString newImageOPs;
  for(size_t i = 0; i < numOutputImgs; ++i)
  {
    if(found_bspline)
      newImageOPs.push_back(evaluateBSpline<CrystalOrderParameters>(
                             parametersAfterEuler, 
                             imageOPsAfterEuler, 
                             totalLength*(Real)i/(Real)(numOutputImgs - 1) ) );
    else if(found_linear)
      newImageOPs.push_back(evaluateLinearSpline<CrystalOrderParameters>(
                             parametersAfterEuler, 
                             imageOPsAfterEuler, 
                             totalLength*(Real)i/(Real)(numOutputImgs - 1) ) );
    else // Should never get here! 
    {
      std::cerr << "Internal error - bad interpolation rules" << std::endl
                << "Exiting..." << std::endl;
      return -1;
    }
  }

  // Find the old image closest to each new image and write to file

  std::stringstream ImgFileNameStream;
  ImgFileNameStream << outPrefix << "imgs.out";
  std::ofstream imgFile(ImgFileNameStream.str().c_str(), std::ios::out);
  for(size_t i = 0; i < numOutputImgs; ++i)
  {
    size_t nearest = -1;
    Real minDistSq = REAL_VERYBIG;
    for(size_t j = 0; j < numInputImgs; ++j)
    {
      Real distSq = 0.0;
      for(size_t k = 0; k < numParameters; ++k)
      {
        Real delta = newImageOPs[i][k] - imageOPs[j][k];
        distSq += delta*delta;
      }
      if(distSq < minDistSq)
      {
        nearest = j;
        minDistSq = distSq;
      }
    }
    imgFile << nearest << std::endl;
  }
  imgFile.close();

  // Write images to file
  
  size_t numWidth = (size_t)(log10((Real)numOutputImgs)) + 1;
  for(size_t i = 0; i < numOutputImgs; ++i)
  {
    std::stringstream COPFileNameStream;
    COPFileNameStream << outPrefix << std::setfill('0') << std::setw(numWidth) << i << ".cop";
    COPFile currentCOPFile(COPFileNameStream.str(), OUT);
    currentCOPFile << newImageOPs[i];
    currentCOPFile.close();
  }
 
  // Calculate the Frechet distance between the old and new strings and write to standard output

  std::cout << "Frechet distance between old and new strings: "
            << frechetDistance<COPString, CrystalOrderParameters, Real>(imageOPs, newImageOPs)
            << std::endl;

  return 0;
}

