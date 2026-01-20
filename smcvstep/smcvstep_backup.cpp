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
**
** Notes:
**
** - For the reparameterization, it may be a good idea to iterate to
**   make sure that the parameter constraint is imposed more accurately.
**   Especially with chords, now it's still a bit off.
*/

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <cstring>
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
            << "                                   constant arc length." << std::endl
            << "         -frcweight              - Weigh the arc (or chord) length by" << std::endl
            << "                                   the norm of the restraint force vector." << std::endl
            << "         -pmfweight              - Weigh the arc (or chord) length by" << std::endl
            << "                                   the potental of mean force." << std::endl
            << "         -extweight filename     - Weigh the arc (or chord) length using" << std::endl
            << "                                   weights read from and external file" << std::endl
            << "         -extpars filename       - Read the parameter values to interpolate" << std::endl
            << "                                   (as a fraction from zero to one) from" << std::endl
            << "                                   the given file." << std::endl << std::endl
            << "         The -cop, -opd, -tstep and -out options are mandatory." << std::endl
            << "         The -bspline and -linear option are mutually exclusive, if none" << std::endl
            << "         is given a B-spline is used by default." << std::endl
            << "         The -frcweight, -pmfweight and -extweight options are mutually" << std::endl
            << "         exclusive, if none is given, no weights are applied." << std::endl
            << "         It is best to not use force or pmf weights until the string is" << std::endl
            << "         close to convergence." << std::endl
            << "         Using the -extpars option in combination with weights means that" << std::endl
            << "         the parameter values are for the weighted arclengths/chords. In" << std::endl
            << "         general it is simpler to not use weights in this case." << std::endl
            << "         If the -nimg option is not given, the number of output images" << std::endl
            << "         will be the same as the number of input images." << std::endl
            << "         If the number of output images is zero, the program only finds" << std::endl
            << "         the potential of mean force." << std::endl
            << "         Weighing by the restraint force vector produces more images in" << std::endl
            << "         areas where the free energy changes the most, weighing by" << std::endl
            << "         the potential of mean force produces more images where the" << std::endl
            << "         free energy is higher (e.g. the transition region)." << std::endl
            << "         If the force tolerance is not given, it is assumed to be zero," << std::endl
            << "         i.e. all order parameters will be used. " << std::endl
            //JAKE climbing image added Apr 2024
            << "         -climbing_image is optional. It expects a number to be given." << std::endl
            << "          Only use this image when you are reasonably sure you've found the saddle point." << std::endl
            << "          Give the number of the image close to the saddle point." << std::endl
            //JAKE fixing images June 2024
            << "          -fixed_climbing fixes the climbing image to let the path evolve better." << std::endl 
            << "          -fixed_others does the opposite and fixes everything else." << std::endl;
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
  bool found_frcweight = false;
  bool found_pmfweight = false;
  bool found_extweight = false;
  bool found_extpars = false;

//JAKE - CLIMBING IMAGE apr 2024
  bool found_climbing = false;
//  bool begin_climbing = false;
  size_t climbing_image = 0;
//JAKE - FIXING (as in, making them not move, not fix as in repair) THINGS JUNE 2024
  bool found_fixed_climbing = false;
  bool found_fixed_others = false;

  size_t numOutputImgs = 0; 
  std::string copListFileName, opdListFileName, outPrefix, weightsFileName, parsFileName; 
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
           timeStep < 0.0)
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
        found_chords = true;
      else if(arguments[i] == "-frcweight")
        found_frcweight = true;
      else if(arguments[i] == "-pmfweight")
        found_pmfweight = true;
      else if(arguments[i] == "-extweight")
      {
        found_extweight = true;
        if(i + 1 >= arguments.size() ||
           arguments[i + 1].find("-") != arguments[i + 1].npos)
        {
          std::cerr << "Name of file containing the weights must be given" << std::endl
                    << "after the -extweight option." << std::endl
                    << "Use smcvstep -help for help." << std::endl;
          return -1;
        }
        else weightsFileName = arguments[i + 1];
      } // End of "-extweight" option
      else if (arguments[i] == "-extpars")
      {
        found_extpars = true;
        if(i + 1 >= arguments.size() ||
           arguments[i + 1].find("-") != arguments[i + 1].npos)
        {
          std::cerr << "Name of file containing the relative parameter values" << std::endl
                    << "must be given after the -extpars option." << std::endl
                    << "Use smcvstep -help for help." << std::endl;
          return -1;
        }
        else parsFileName = arguments[i + 1];
      } // End of "-extpars" option
      //JAKE - CLIMBING images option apr 2024
      else if (arguments[i] == "-climbing_image") 
        {
        found_climbing=true;
        if(i + 1 >= arguments.size() ||
          arguments[i + 1].find("-") != arguments[i + 1].npos)
          {
          std::cerr << "Must give an image number for climbing option" << std::endl
                    << "Use smcvstep -help for help." << std::endl;
          return -1;
          }

          climbing_image = strtol(arguments[i+1].c_str(), NULL ,10);
          if (climbing_image == 0) {std::cerr << "Climbing Image is zero, this is probably a mistake, meaning strtol returned 0 due to a letter or something." << std::endl; return -1; }
        }// End of "-climbing_images" option //JAKE
        //JAKE fixed_climbing option, this means fix the position of the climbing image, and not the others 
        else if (arguments[i] == "-fixed_climbing")
         {
         found_fixed_climbing = true;
          }//end fixed_climbing
        //JAKE fixed_others - meaning fix everything except the climbing image
        else if (arguments[i] == "-fixed_others")
          {
          found_fixed_others=true; 
          }//end fixed_others
          //JAKE JUNE 2024, adding options to fix images
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

  if((found_frcweight && found_pmfweight) || (found_frcweight && found_extweight) ||
     (found_pmfweight && found_extweight) || //JAKE apr 2024
     (found_frcweight && found_climbing) || 
     (found_extweight && found_climbing) || 
     (found_pmfweight && found_climbing))//JAKE
  {
    std::cerr << "Error: The -frcweight, -pmfweight and -extweight and -climbing_image options are mutually exclusive" << std::endl
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
  // JAKE june `24, reversing order so we integrate from the liquid side. I could also switch the images in the files I guess
  COPString reversedImageOPs;
  ForceString reversedImageForceData;
  
  size_t imageOPsSize;
  size_t imageForceSize;
  imageOPsSize = imageOPs.size();
  imageForceSize = imageForceData.size();
  for (size_t i=0; i<imageOPs.size(); i++) reversedImageOPs.push_back(imageOPs[imageOPsSize-1-i]);//aka reverse order
  for (size_t i=0; i< imageForceData.size(); i++) reversedImageForceData.push_back(imageForceData[imageForceSize-1-i]);
  //JAKE
  //this bit is the original calculation 
  std::vector<Real> parameters = 
    trapezoidCumulativeArcLength<COPString, Real>(imageOPs);
  std::vector<Real> cumulativePMF =
    trapezoidCumulativeIntegral<COPString, ForceString, Real>(imageOPs, imageForceData);


  //and this is the reversed calculation
  std::vector<Real> reversedParameters = 
    trapezoidCumulativeArcLength<COPString, Real>(reversedImageOPs);
  std::vector<Real> reversedCumulativePMF =
    trapezoidCumulativeIntegral<COPString, ForceString, Real>(reversedImageOPs, reversedImageForceData);
  // Write potential of mean force to file

  IOFile myPMFFile(outPrefix+".pmf", OUT);
  for(size_t i = 0; i < numInputImgs; ++i)
    myPMFFile << parameters[i] << " " << cumulativePMF[i] << std::endl;
  //JAKE write reversed file too
  IOFile myReversedPMFFile(outPrefix+"_reversed.pmf", OUT);
  for(size_t i = 0; i < numInputImgs; ++i)
    myReversedPMFFile << reversedParameters[numInputImgs-1 - i] << " " << reversedCumulativePMF[numInputImgs-1-i] << std::endl;

  // Exit if nimg = 0, give error message if nimg < 4
  if(numOutputImgs == 0) return 0;
  if(numOutputImgs < 4)
  {
    std::cerr << "Number of output images must be greater than 3 " << std::endl
              << "in order to generate new images" << std::endl
              << "Exiting..." << std::endl;
    return -1;
  }
  //JAKE
  //check -climbing_image option for validity
  if (climbing_image > numInputImgs || climbing_image < 0) {std::cerr << "Climbing image greater than number of images! or less than 0" << std::endl; return -1;}


  //JAKE, do a spline thru current images
  COPString imageOPsBeforeEuler = imageOPs;//copy in case I need to modify anything
  // make parameters before euler
  std::vector<Real> parametersBeforeEuler;

  if(found_chords) // Use constant distance to initial point
  {
    parametersBeforeEuler.push_back(0.0);
    for(size_t i = 1; i < numInputImgs; ++i)
    {
      Real distSq = 0.0;
      for(size_t j = 0; j < numParameters; ++j)
      {
        Real delta = imageOPsBeforeEuler[i][j] - 
                     imageOPsBeforeEuler[0][j];
        distSq += delta*delta;
      }
      parametersBeforeEuler.push_back(sqrt(distSq));
    }
    //JAKE
    }
  else // Use constant arc length
  {
    parametersBeforeEuler.push_back(0.0);
    for(size_t i = 1; i < numInputImgs; ++i)
    {
      Real distSq = 0.0;
      for(size_t j = 0; j < numParameters; ++j)
      {
        Real delta = imageOPsBeforeEuler[i][j] - 
                     imageOPsBeforeEuler[i - 1][j];
        distSq += delta*delta;
      }
      parametersBeforeEuler.push_back(parametersBeforeEuler[i - 1] + sqrt(distSq));
    }
  }
  if (found_extpars && found_climbing){ std::cerr << "Climbing Images and extpars are mutually exclusive!!!" << std::endl; return -1;}
//tangent vector to putative saddle point calcuated via derivative of BSpline

//JAKE derivatives
CrystalOrderParameters tangentCOP;
CrystalOrderParameters tangentCOPBackward;
CrystalOrderParameters tangentCOPBackward2;
if (found_climbing){
if (found_bspline){
 
tangentCOPBackward2 = evaluateBSplineDerivative<CrystalOrderParameters>(parametersBeforeEuler,
                               imageOPsBeforeEuler,
                               parametersBeforeEuler[climbing_image]);
}
else {//linear option
//tangentCOP = evaluateLinearSplineDerivative<CrystalOrderParameters>(
//                             parametersBeforeEuler, 
//                             imageOPsBeforeEuler, 
//                             parametersBeforeEuler[climbing_image]);
//try backward difference to see if it looks the same?
tangentCOPBackward2 = (imageOPsBeforeEuler[climbing_image] - imageOPsBeforeEuler[climbing_image-1]);
//tangentCOPBackward2 = (imageOPsBeforeEuler[climbing_image+1] - imageOPsBeforeEuler[climbing_image-1]);
Real tangentNorm2 = 0;
//Real tangentNorm = 0;
for (size_t i=0; i<tangentCOPBackward2.size(); i++) {
tangentNorm2+=(imageOPsBeforeEuler[climbing_image][i] - imageOPsBeforeEuler[climbing_image-1][i])*(imageOPsBeforeEuler[climbing_image][i] - imageOPsBeforeEuler[climbing_image-1][i]);
//tangentNorm2+=(imageOPsBeforeEuler[climbing_image+1][i] - imageOPsBeforeEuler[climbing_image-1][i])*(imageOPsBeforeEuler[climbing_image+1][i] - imageOPsBeforeEuler[climbing_image-1][i]);

}

for (size_t i=0; i<numParameters; ++i) {
//tangentCOPBackward[i]/=sqrt(tangentNorm);
tangentCOPBackward2[i]/=sqrt(tangentNorm2);}
//for (int i=0;i<numParameters; i++){
//std::cerr  << tangentCOPBackward[i] << " " << tangentCOPBackward2[i] << std::endl;
//}
}
}
if ((!found_climbing && found_fixed_climbing) || (!found_climbing && found_fixed_others))
  {
  std::cerr << "need to specify a climbing image if fixing things!" << std::endl; return -1;
  }

#define TimeStepReductionClimbing 4.0
#define FrictionParameter 1.0
  // Do a forward Euler step to generate new string
  COPString imageOPsAfterEuler = imageOPs;
  if (!found_climbing)  
    {
    for(size_t i = 0; i <  numInputImgs; ++i)
      {
        VectorND EulerStep = timeStep*imageForceData[i].mMatrix*imageForceData[i].restraintForces;
        for(size_t j = 0; j < numParameters; ++j)
        imageOPsAfterEuler[i][j] -= EulerStep[j];
      }
    }
   else {
    //need spline for this one, through original images
    for(size_t i = 0; i <  numInputImgs; ++i)
      {
      if (i!= climbing_image && !found_fixed_others)
        {
        VectorND EulerStep = timeStep*imageForceData[i].mMatrix*imageForceData[i].restraintForces;//nabla V
        for(size_t j = 0; j < numParameters; ++j)
        imageOPsAfterEuler[i][j] -= EulerStep[j];
        }
      else if (i!=climbing_image && found_fixed_others)
        {
         for(size_t j = 0; j < numParameters; ++j)
         imageOPsAfterEuler[i][j] = imageOPsBeforeEuler[i][j];
        }
      else if (found_fixed_climbing)//only happens if i==climbing_image, because it'll enter the other if elses before this one otherwise 
        {
         for(size_t j = 0; j < numParameters; ++j)
           imageOPsAfterEuler[i][j] = imageOPsBeforeEuler[i][j];
        }
      else //climbing image
        {
         VectorND ClimbingStep = imageForceData[i].mMatrix*imageForceData[i].restraintForces;//nabla V
         Real dotSum = 0;
         Real orthoSum = 0;
         Real TSum = 0;
          for (size_t j=0; j< numParameters; j++) 
          {
          dotSum+=ClimbingStep[j]*tangentCOPBackward2[j];
          TSum += tangentCOPBackward2[j]*tangentCOPBackward2[j];
          }
          for (size_t j=0; j< numParameters; j++) orthoSum += (ClimbingStep[j] - dotSum*tangentCOPBackward2[j])*(ClimbingStep[j] - dotSum*tangentCOPBackward2[j]);
        
         std::cerr << dotSum << " " << norm(ClimbingStep) << " " << sqrt(orthoSum) << std::endl;
     
         dotSum*=2;
         for (size_t j=0; j< numParameters; j++) tangentCOPBackward2[j]*=dotSum;
         ClimbingStep*=-1;
         for (size_t j=0; j< numParameters; j++)ClimbingStep[j]+=tangentCOPBackward2[j];
         for (size_t j=0; j< numParameters; j++)ClimbingStep[j]*=timeStep*FrictionParameter;//TimeStepReductionClimbing;
         for (size_t j=0; j< numParameters; j++)imageOPsAfterEuler[i][j]+=ClimbingStep[j];
        }
      }
        }

  //for (size_t i=0; i<numParameters; i++) std::cerr << imageOPsAfterEuler[climbing_image][i] << std::endl;
  // Interpolate to calculate weights (this step uses arclenghts)
//calculate a different sort of metric for convergence.  
IOFile myConvergenceFile("C.out", OUT);
Real Bops_C;
Real Lops_C;
for (size_t j=0; j<numInputImgs; j++){
Real Bops_temp=0;
Real Lops_temp=0;
for (size_t i=0; i<numParameters/4; i++){
  Bops_temp += (imageOPsAfterEuler[j][i+64] - imageOPsBeforeEuler[j][i+64])*(imageOPsAfterEuler[j][i+64] - imageOPsBeforeEuler[j][i+64]);
  Lops_temp += (imageOPsAfterEuler[j][i+192] - imageOPsBeforeEuler[j][i+192])*(imageOPsAfterEuler[j][i+192] - imageOPsBeforeEuler[j][i+192]);
}
Bops_C += Bops_temp;
Lops_C += Lops_temp;
}
myConvergenceFile << Bops_C/20.0 << " " << Lops_C/20.0  << std::endl;


std::vector<Real> arcLengthsAfterEuler =
    trapezoidCumulativeArcLength<COPString, Real>(imageOPsAfterEuler); 

  std::vector<Real> weights;

  if(found_frcweight)
  {
    // Calculate norms of restraint forces
    std::vector<Real> forceNorms;
    for(size_t i = 0; i < numInputImgs; ++i)
      forceNorms.push_back(norm(imageForceData[i].restraintForces));

    // Interpolate to calculate weights using scaled arclengths
    std::vector<Real> scaledBeforeEuler, scaledAfterEuler;
    for(size_t i = 0; i < numInputImgs; ++i){
      scaledBeforeEuler.push_back(parameters[i]/parameters[numInputImgs - 1]);
      scaledAfterEuler.push_back(arcLengthsAfterEuler[i]/
                                  arcLengthsAfterEuler[numInputImgs - 1]);
    }
    for(size_t i = 0; i < numInputImgs; ++i)
      weights.push_back(evaluateLinearSpline<Real>(
                         scaledBeforeEuler,
                         forceNorms,
                         scaledAfterEuler[i]));
  }
  else if(found_pmfweight)
  {
    // Make PMF non-negative
    Real minPMF = REAL_VERYBIG;
    for(size_t i = 0; i < numInputImgs; ++i)
      if(cumulativePMF[i] < minPMF) minPMF = cumulativePMF[i];
    for(size_t i = 0; i < numInputImgs; ++i)
      cumulativePMF[i] += minPMF;
    // Interpolate to calculate weights using scaled arclengths
    std::vector<Real> scaledBeforeEuler, scaledAfterEuler;
    for(size_t i = 0; i < numInputImgs; ++i)
    {
      scaledBeforeEuler.push_back(parameters[i]/parameters[numInputImgs - 1]);
      scaledAfterEuler.push_back(arcLengthsAfterEuler[i]/
                                  arcLengthsAfterEuler[numInputImgs - 1]);
    }
    for(size_t i = 0; i < numInputImgs; ++i)
      weights.push_back(evaluateLinearSpline<Real>(
                         scaledBeforeEuler,
                         cumulativePMF,
                         scaledAfterEuler[i]));
  }
  else if(found_extweight)
  {
    IOFile weightsFile(weightsFileName, IN);
    if(!weightsFile)
    {
      std::cerr << "Error opening external weights file " << weightsFileName << std::endl
                << "Exiting..." << std::endl;
      return -1;
    }
    while(!weightsFile.eof())
    {
      Real currentWeight;
      weightsFile >> currentWeight;
      if(weightsFile.eof()) break;
      if(!weightsFile)
      {
        std::cerr << "Error reading external weights file " << weightsFileName << std::endl
                  << "Exiting..." << std::endl;
        return -1;
      }
      if(currentWeight < 0.0)
      {
        std::cerr << "Error: External weights must be nonnegative" << std::endl
                  << "Exiting..." << std::endl;
        return -1;
      }
      weights.push_back(currentWeight);
    }
    if(weights.size() != numInputImgs)
    {
      std::cerr << "Error: Inconsistent number of images in external weights file " << weightsFileName << std::endl
                << "Exiting..."  << std::endl;
      return -1;
    }
   }
  else 
  {
    // All weights are one
    for(size_t i = 0; i < numInputImgs; ++i)
      weights.push_back(1.0);
  }


  // Reparameterize to constant (weighted) arc length or distance to initial point

  std::vector<Real> parametersAfterEuler;

  if(found_chords) // Use constant distance to initial point
  {
    parametersAfterEuler.push_back(0.0);
    for(size_t i = 1; i < numInputImgs; ++i)
    {
      Real distSq = 0.0;
      for(size_t j = 0; j < numParameters; ++j)
      {
        Real delta = weights[i]*imageOPsAfterEuler[i][j] - 
                     weights[0]*imageOPsAfterEuler[0][j];
        distSq += delta*delta;
      }
      parametersAfterEuler.push_back(sqrt(distSq));
    }
  }
  else // Use constant arc length
  {
    parametersAfterEuler.push_back(0.0);
    for(size_t i = 1; i < numInputImgs; ++i)
    {
      Real distSq = 0.0;
      for(size_t j = 0; j < numParameters; ++j)
      {
        Real delta = weights[i]*imageOPsAfterEuler[i][j] - 
                     weights[i - 1]*imageOPsAfterEuler[i - 1][j];
        distSq += delta*delta;
      }
      parametersAfterEuler.push_back(parametersAfterEuler[i - 1] + sqrt(distSq));
    }
  }

  Real totalLength = parametersAfterEuler[numInputImgs - 1];
//  if (found_climbing && numOutputImgs%5!=0) {std::cerr << "Number of output images should be divisible by 5 for climbing images to be used properly." << std::endl; return -1;}
  // Build array of output parameter values
  std::vector<Real> outputPars;
  if(found_extpars)
  {
    // Read them from the external file
    IOFile parsFile(parsFileName, IN);
    if(!parsFile)
    {
      std::cerr << "Error opening external parameters file " << parsFileName << std::endl
                << "Exiting..." << std::endl;
      return -1;
    }
    while(!parsFile.eof())
    {
      Real currentPar;
      parsFile >> currentPar;
      if(parsFile.eof()) break;
      if(!parsFile)
      {
        std::cerr << "Error reading external parameters file " << parsFileName << std::endl
                  << "Exiting..." << std::endl;
        return -1;
      }
      if(currentPar < 0.0 || currentPar > 1.0)
      {
        std::cerr << "Error: Relative parameters must be between zero and one" << std::endl
                  << "Exiting..." << std::endl;
        return -1;
      }
      outputPars.push_back(currentPar*totalLength);
    }
    if(outputPars.size() != numOutputImgs)
    {
      std::cerr << "Error: Inconsistent number of images in external parameters file " << parsFileName << std::endl
                << "Exiting..." << std::endl;
      return -1;
    }
  }
  else if(!found_climbing) // No external file given, make equispaced
    for(size_t i = 0; i < numOutputImgs; ++i)
      outputPars.push_back(totalLength*(Real)i/(Real)(numOutputImgs - 1));

  #define climbing_image_multiplier 0.0001
  //for now, I'm going to assume the climbing image is the very end of the string, this commented piece makes the 9th image the climbing image (8th index)
  else//JAKE create parameters for climbing image string 
    {//this works unless you want over 100 images, then you need to change the multiplier
//    for (size_t i=0; i<numOutputImgs; i++)
//      {
//      if (i < 2*(numOutputImgs/5)-1)         outputPars.push_back((parametersAfterEuler[climbing_image]/Real(2*numOutputImgs/5))*Real(i));
//      else if (i == 2*(numOutputImgs/5)-1)   outputPars.push_back(parametersAfterEuler[climbing_image] - totalLength*climbing_image_multiplier);
//      else if (i==2*(numOutputImgs/5))       outputPars.push_back(parametersAfterEuler[climbing_image]);
//      else if (i == 2*(numOutputImgs/5)+1)   outputPars.push_back(parametersAfterEuler[climbing_image] + totalLength*climbing_image_multiplier);
//      else {      outputPars.push_back(parametersAfterEuler[climbing_image] + 
//                                             Real(i-(2*numOutputImgs/5))*(totalLength-parametersAfterEuler[climbing_image])/Real(numOutputImgs-(2*numOutputImgs/5)-1));
//           }     
//                                        
//      }
//this is just to see whats going on
   for (size_t i=0; i<numOutputImgs; i++) outputPars.push_back(parametersAfterEuler[climbing_image]/Real(numOutputImgs-1)*Real(i));
   for (size_t i=0; i<numOutputImgs; i++) std::cerr << outputPars[i]/totalLength << std::endl;
   }


  COPString newImageOPs;
  for(size_t i = 0; i < numOutputImgs; ++i)
  {
    if(found_bspline)
      newImageOPs.push_back(evaluateBSpline<CrystalOrderParameters>(
                             parametersAfterEuler, 
                             imageOPsAfterEuler, 
                             outputPars[i]));
    else if(found_linear)
      newImageOPs.push_back(evaluateLinearSpline<CrystalOrderParameters>(
                             parametersAfterEuler, 
                             imageOPsAfterEuler, 
                             outputPars[i]));
    else // Should never get here! 
    {
      std::cerr << "Internal error - bad interpolation rules" << std::endl
                << "Exiting..." << std::endl;
      return -1;
    }
  }




////JAKE
//  std::vector<Real> parametersBeforeEuler;
//  if(found_chords) // Use constant distance to initial point
//  {
//    parametersBeforeEuler.push_back(0.0);
//    for(size_t i = 1; i < numInputImgs; ++i)
//    {
//      Real distSq = 0.0;
//      for(size_t j = 0; j < numParameters; ++j)
//      {
//        Real delta = weights[i]*imageOPs[i][j] -
//                     weights[0]*imageOPs[0][j];
//        distSq += delta*delta;
//      }
//      parametersBeforeEuler.push_back(sqrt(distSq));
//    }
//  }
//  else // Use constant arc length
//  {
//    parametersBeforeEuler.push_back(0.0);
//    for(size_t i = 1; i < numInputImgs; ++i)
//    {
//      Real distSq = 0.0;
//      for(size_t j = 0; j < numParameters; ++j)
//      {
//        Real delta = weights[i]*imageOPs[i][j] -
//                     weights[i - 1]*imageOPs[i - 1][j];
//        distSq += delta*delta;
//      }
//      parametersBeforeEuler.push_back(parametersBeforeEuler[i - 1] + sqrt(distSq));
//    }
//  }

//  Real totalLength = parametersAfterEuler[numInputImgs - 1];


//utter shite JAKE dec 2023
//  COPString newImageOPsDense;
//  COPString imageOPsDense;
//  std::vector<Real> bigPar;
//  Real fakeNumOutput = 100.0;
//  Real count=0.0;
//  for(size_t i = 0; i < 101; ++i)
//  {
//    bigPar.push_back(count/fakeNumOutput); 
//    if(found_bspline){
//     bigPar.push_back(count/fakeNumOutput); 
//     newImageOPsDense.push_back(evaluateBSpline<CrystalOrderParameters>(
//                             parametersAfterEuler, 
//                             imageOPsAfterEuler, 
//                             bigPar[i]));
//      imageOPsDense.push_back(evaluateBSpline<CrystalOrderParameters>(
//                             parametersBeforeEuler, 
//                             imageOPs, 
//                             bigPar[i]));
//     count++;
//     }
//    else if(found_linear){
//      newImageOPsDense.push_back(evaluateLinearSpline<CrystalOrderParameters>(
//                             parametersAfterEuler, 
//                             imageOPsAfterEuler, 
//                             bigPar[i]));
//
// imageOPsDense.push_back(evaluateBSpline<CrystalOrderParameters>(
//                             parametersBeforeEuler,
//                             imageOPs,
//                             bigPar[i]));
//
//     count++;   
// }
//    else // Should never get here! 
//    {
//      std::cerr << "Internal error - bad interpolation rules" << std::endl
//                << "Exiting..." << std::endl;
//      return -1;
//    }
//  }

//JAKE





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

//JAKE Nov 24, getting to the bottom of issues with reinterpolating, frechet distance not changing wit timestep size
//  for (size_t i=0; i<numOutputImgs; ++i) {
//    double maximum = 0;
//    for (size_t j=0; j<64 ; ++j){
//      if  ((newImageOPs[i][j+64]-imageOPs[i][j+64])*( newImageOPs[i][j+64]-imageOPs[i][j+64]) > maximum) maximum = sqrt((newImageOPs[i][j+64]-imageOPs[i][j+64])*( newImageOPs[i][j+64]-imageOPs[i][j+64]));
//    }
//  std::cout << i << " " <<maximum << std::endl;
//  }
for (size_t i=0; i<parametersAfterEuler.size(); ++i) std::cout << outputPars[i] << std::endl;

  std::cout << "Frechet distance between old and new strings: "
            << frechetDistance<COPString, CrystalOrderParameters, Real>(imageOPs, newImageOPs)
            << std::endl;

//  std::cout << "Fixed Frechet distance between old and new strings: "
//            << frechetDistance<COPString, CrystalOrderParameters, Real>(imageOPsDense, newImageOPsDense)
//            << std::endl;

  

  
  return 0;
}

