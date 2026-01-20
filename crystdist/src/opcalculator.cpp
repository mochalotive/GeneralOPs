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
** Crystal Order Parameters class
**
** Notes:
**
** - The implementation of the "per-molecule" case is not very elegant,
**   and it wastes some memory. It might be a good idea to change it for
**   a later version.
** - Note that p-norms are not applied to local densities
*/


#include <iostream>
#include <fstream>
#include <string>
#include "crystdist/include/opcalculator.h"
#include <cstdlib>
//#define TanhModifier 5.0
#define tanhtest true
// Constructors

Real tanh_rational_approximation(Real x)
  {
  if (x > 3) return (1.0);
  else if (x > -3.0 && x< 3.0) return (x*(105.0 + 10.0*x*x) / (105.0 + 45.0*x*x + x*x*x*x));
  else return( -1.0);
  }



OrderParameterCalculator::OrderParameterCalculator()
:
ops_()
{}

OrderParameterCalculator::OrderParameterCalculator(System<PointMolecule> const &system, 
                                                   CrystalDistributionParameters const &parameters, 
                                                   OrderParameterGrid const &grid,
                                                   Real const &switchWidth,
                                                   Real const &cutoff, 
                                                   bool const useLogs,
                                                   Real const &pValue,
                                                   Real const &TanhModifier,
                                                   bool const useNorm)
:
ops_()
{ 
  ops_.switchWidth = switchWidth;
  ops_.cutoffSq = cutoff*cutoff;
  ops_.grid = grid;
  calculate(system, parameters, useLogs, pValue, TanhModifier, useNorm); 
}

// Interface

void OrderParameterCalculator::calculate(System<PointMolecule> const &system, 
                                         CrystalDistributionParameters const &parameters,
                                         bool const useLogs,
                                         Real const &pValue,
                                         Real const &TanhModifier,
                                         bool const useNorm)
{
  // Note that this assumes (without verifying) that the system is homogeneous (e.g. internal DOFs are
  // the same for all molecules)
  if(system.size() < 1)
  {
    std::cerr << "Error in CrystalOrderParameters::calculate(): Empty system" << std::endl;
    return;
  }

  // Store number of groups
  size_t const numGroups = parameters.means.size();
  if(numGroups < 1)
  {
    std::cerr << "Error in CrystalOrderParameters::calculate: Empty parameter set" << std::endl;
    clearOPs();
    return;
  }

  // Sanity check for p-norms
  if(pValue > 0.0 && pValue < 1.0)
  {
    std::cerr << "Error in CrystalOrderParameters::calculate: p value must be at least 1" << std::endl;
    clearOPs();
    return;
  }

  // Initialize variables
  clearOPs();

  size_t const numCells = ops_.grid.numCells();
  bool const useGrid = (numCells > 0);
  std::vector<Real> totalWeights((numCells > 0)?numCells:system.size(), 0);
  size_t const numIntDOFs = system[0].numInternalDOFs();
  if(useGrid) ops_.resize(numCells, numIntDOFs);
  else ops_.resize(system.size(), numIntDOFs);  // Do OPs per molecule
  //JAKE 2025
  //for compatitbility with plumed
  std::vector<Real> coordination;
  for(size_t i = 0; i < system.size(); ++i)
    {
    //one can interpret this as
    //every molecule is  "coordinated with itself" 
    coordination.push_back(1.0);
    for(size_t j = 0; j < system.size(); ++j)
      { //if nopbc flag, within cutoff, and a non-zero switch width (it'll crash later anyway if switch width is zero)
      if (system.lattice().fakeNOPBC && norm2(system.lattice().difference(system[i].position, system[j].position)) < ops_.cutoffSq   && i!= j && ops_.switchWidth != 0.0)
        coordination[i] = coordination[i]+ 1;
                         // - tanh((-sqrt(ops_.cutoffSq) + sqrt(norm2( //in my experience, this is close enough without a switching function
                         //   system.lattice().difference(system[i].position, system[j].position))))/ops_.switchWidth );   
      }
    }


  for(size_t i = 0; i < system.size(); ++i)
  {
    // Find cell index and update number of molecules
    CellData iData;
    if(useGrid)
      { 
      iData = cellData(system[i].position, system.lattice());
      }
    else  // Per-molecule case
    {
      iData.clear();
      iData.cellIndices.push_back(i);
      iData.cellWeights.push_back(1.0);
    }
    size_t const iDataSize = iData.size();
    
    for(size_t iIndex = 0; iIndex < iDataSize; ++iIndex)
      totalWeights[iData.cellIndices[iIndex]] += iData.cellWeights[iIndex];
    for(size_t j= i+1; j < system.size(); ++j) // Loop over neighbors
    {
    // Check neighbor is within cutoff
      if(ops_.cutoffSq && 
         norm2(system.lattice().difference(system[i].position, system[j].position)) > ops_.cutoffSq)
         continue;
      // Find cell index and relative configurations
      CellData jData;
      if(useGrid) jData = cellData(system[j].position, system.lattice());
      else
      {
        jData.clear();
        jData.cellIndices.push_back(j);
        jData.cellWeights.push_back(1.0);
      }
      size_t const jDataSize = jData.size();
      RelativeConfiguration const confij(system[i], system[j], system.lattice());
      RelativeConfiguration const confji(system[j], system[i], system.lattice());
      for(size_t iGroup = 0; iGroup < numGroups; ++iGroup)
      {
        size_t numPeaks = parameters.means[iGroup].size();
        for(size_t iPeak = 0; iPeak < numPeaks; ++iPeak)
        {
          // Distance OPs
          Real deltaDistance = confij.distance - parameters.means[iGroup][iPeak].distance; 
          Real distanceOP = parameters.normalizationFactors[iGroup][iPeak].distance*
                            exp(-0.5*parameters.concentrationParameters[iGroup][iPeak].distance*
                                deltaDistance*deltaDistance);
          if(pValue > 0.0)
          {
            for(size_t iIndex = 0; iIndex < iDataSize; ++iIndex)
              ops_.distance[iData.cellIndices[iIndex]] += 
                pow(iData.cellWeights[iIndex]*distanceOP, pValue);
            for(size_t jIndex = 0; jIndex < jDataSize; ++jIndex)
              ops_.distance[jData.cellIndices[jIndex]] += 
                pow(jData.cellWeights[jIndex]*distanceOP, pValue);
          }
          else
          {
            for(size_t iIndex = 0; iIndex < iDataSize; ++iIndex)
           {
              ops_.distance[iData.cellIndices[iIndex]] += iData.cellWeights[iIndex]*distanceOP/coordination[i];
           }
            for(size_t jIndex = 0; jIndex < jDataSize; ++jIndex)
              ops_.distance[jData.cellIndices[jIndex]] += jData.cellWeights[jIndex]*distanceOP/coordination[j];
          }

          // Bond orientation and relative orientation OPs
          if(confij.types.first == GENERAL && confij.types.second == GENERAL)
          {
            // General case: Bond orientation = Vector3D, Relative orientation = Quaternion

            Real bondDotij = vector(confij.bondOrientation)*
                             vector(parameters.means[iGroup][iPeak].bondOrientation);
            Real bondDotji = vector(confji.bondOrientation)*
                             vector(parameters.means[iGroup][iPeak].bondOrientation);
            Real bondOPij = parameters.normalizationFactors[iGroup][iPeak].bondOrientation*
                            exp(parameters.concentrationParameters[iGroup][iPeak].bondOrientation*
                                bondDotij);
            Real bondOPji = parameters.normalizationFactors[iGroup][iPeak].bondOrientation*
                            exp(parameters.concentrationParameters[iGroup][iPeak].bondOrientation*
                                bondDotji);
            if(pValue > 0.0)
            {
              for(size_t iIndex = 0; iIndex < iDataSize; ++iIndex)
                ops_.bondOrientation[iData.cellIndices[iIndex]] += 
                  pow(iData.cellWeights[iIndex]*distanceOP*bondOPij, pValue);
              for(size_t jIndex = 0; jIndex < jDataSize; ++jIndex)
                ops_.bondOrientation[jData.cellIndices[jIndex]] += 
                  pow(jData.cellWeights[jIndex]*distanceOP*bondOPji, pValue);
            }
            else
            {
              for(size_t iIndex = 0; iIndex < iDataSize; ++iIndex)
                ops_.bondOrientation[iData.cellIndices[iIndex]] += 
                  iData.cellWeights[iIndex]*distanceOP*bondOPij/coordination[i];
              for(size_t jIndex = 0; jIndex < jDataSize; ++jIndex)
                ops_.bondOrientation[jData.cellIndices[jIndex]] += 
                  jData.cellWeights[jIndex]*distanceOP*bondOPji/coordination[j];
            }
            Real relDotij = dot(parameters.means[iGroup][iPeak].relativeOrientation,
                                confij.relativeOrientation);
            Real relDotji = dot(parameters.means[iGroup][iPeak].relativeOrientation,
                                confji.relativeOrientation);
            Real relOPij = parameters.normalizationFactors[iGroup][iPeak].relativeOrientation*
                           exp(parameters.concentrationParameters[iGroup][iPeak].relativeOrientation*
                               relDotij*relDotij);
            Real relOPji = parameters.normalizationFactors[iGroup][iPeak].relativeOrientation*
                           exp(parameters.concentrationParameters[iGroup][iPeak].relativeOrientation*
                               relDotji*relDotji);
            if(pValue > 0.0)
            {
              for(size_t iIndex = 0; iIndex < iDataSize; ++iIndex)
                ops_.relativeOrientation[iData.cellIndices[iIndex]] += 
                  pow(iData.cellWeights[iIndex]*distanceOP*relOPij, pValue);
              for(size_t jIndex = 0; jIndex < jDataSize; ++jIndex)
                ops_.relativeOrientation[jData.cellIndices[jIndex]] += 
                  pow(jData.cellWeights[jIndex]*distanceOP*relOPji, pValue);
            }
            else
            {
              for(size_t iIndex = 0; iIndex < iDataSize; ++iIndex)
                ops_.relativeOrientation[iData.cellIndices[iIndex]] += 
                   iData.cellWeights[iIndex]*distanceOP*relOPij/coordination[i];
              for(size_t jIndex = 0; jIndex < jDataSize; ++jIndex)
                ops_.relativeOrientation[jData.cellIndices[jIndex]] += 
                  jData.cellWeights[jIndex]*distanceOP*relOPji/coordination[i];
            }
            // Total OPs added 07-14-10
            if(pValue > 0.0)
            {
              for(size_t iIndex = 0; iIndex < iDataSize; ++iIndex)
                ops_.total[iData.cellIndices[iIndex]] +=
                  pow(iData.cellWeights[iIndex]*distanceOP*bondOPij*relOPij, pValue);
              for(size_t jIndex = 0; jIndex < jDataSize; ++jIndex)
                ops_.total[jData.cellIndices[jIndex]] +=
                  pow(jData.cellWeights[jIndex]*distanceOP*bondOPji*relOPji, pValue);
            }
            else
            {
                          for(size_t iIndex = 0; iIndex < iDataSize; ++iIndex){
                ops_.total[iData.cellIndices[iIndex]] +=
                  iData.cellWeights[iIndex]*distanceOP*bondOPij*relOPij/coordination[i];
}
              for(size_t jIndex = 0; jIndex < jDataSize; ++jIndex)
                ops_.total[jData.cellIndices[jIndex]] +=
                  jData.cellWeights[jIndex]*distanceOP*bondOPji*relOPji/coordination[j];
            }
          }
          else if(confij.types.first == LINEAR_ASYMMETRIC && confij.types.second == LINEAR_ASYMMETRIC)
          {
            // Axisymmetric case: Bond orientation and relative orientation are both angles
            Real deltaBondij = confij.bondOrientation.w - 
                               parameters.means[iGroup][iPeak].bondOrientation.w;
            Real deltaBondji = confji.bondOrientation.w - 
                               parameters.means[iGroup][iPeak].bondOrientation.w;
            Real bondOPij = parameters.normalizationFactors[iGroup][iPeak].bondOrientation*
                            exp(parameters.concentrationParameters[iGroup][iPeak].bondOrientation*
                                cos(deltaBondij));
            Real bondOPji = parameters.normalizationFactors[iGroup][iPeak].bondOrientation*
                            exp(parameters.concentrationParameters[iGroup][iPeak].bondOrientation*
                                cos(deltaBondji));
            if(pValue > 0.0)
            {
              for(size_t iIndex = 0; iIndex < iDataSize; ++iIndex)
                ops_.bondOrientation[iData.cellIndices[iIndex]] += 
                  pow(iData.cellWeights[iIndex]*distanceOP*bondOPij, pValue);
              for(size_t jIndex = 0; jIndex < jDataSize; ++jIndex)
                ops_.bondOrientation[jData.cellIndices[jIndex]] += 
                  pow(jData.cellWeights[jIndex]*distanceOP*bondOPji, pValue);
            }
            else
            {
              for(size_t iIndex = 0; iIndex < iDataSize; ++iIndex)
                ops_.bondOrientation[iData.cellIndices[iIndex]] += 
                  iData.cellWeights[iIndex]*distanceOP*bondOPij;
              for(size_t jIndex = 0; jIndex < jDataSize; ++jIndex)
                ops_.bondOrientation[jData.cellIndices[jIndex]] += 
                  jData.cellWeights[jIndex]*distanceOP*bondOPji;
            }
            Real deltaRelij = confij.relativeOrientation.w -
                              parameters.means[iGroup][iPeak].relativeOrientation.w;
            Real deltaRelji = confji.relativeOrientation.w -
                              parameters.means[iGroup][iPeak].relativeOrientation.w;
            Real relOPij = parameters.normalizationFactors[iGroup][iPeak].relativeOrientation*
                           exp(parameters.concentrationParameters[iGroup][iPeak].relativeOrientation*
                               cos(deltaRelij));
            Real relOPji = parameters.normalizationFactors[iGroup][iPeak].relativeOrientation*
                           exp(parameters.concentrationParameters[iGroup][iPeak].relativeOrientation*
                               cos(deltaRelji));
            if(pValue > 0.0)
            {
              for(size_t iIndex = 0; iIndex < iDataSize; ++iIndex)
                ops_.relativeOrientation[iData.cellIndices[iIndex]] += 
                  pow(iData.cellWeights[iIndex]*distanceOP*relOPij, pValue);
              for(size_t jIndex = 0; jIndex < jDataSize; ++jIndex)
                ops_.relativeOrientation[jData.cellIndices[jIndex]] += 
                  pow(jData.cellWeights[jIndex]*distanceOP*relOPji, pValue);
            }
            else
            {
              for(size_t iIndex = 0; iIndex < iDataSize; ++iIndex)
                ops_.relativeOrientation[iData.cellIndices[iIndex]] += 
                  iData.cellWeights[iIndex]*distanceOP*relOPij;
              for(size_t jIndex = 0; jIndex < jDataSize; ++jIndex)
                ops_.relativeOrientation[jData.cellIndices[jIndex]] += 
                  jData.cellWeights[jIndex]*distanceOP*relOPji;
            }



            // Total OPs added 07-14-10
            if(pValue > 0.0)
            {
              for(size_t iIndex = 0; iIndex < iDataSize; ++iIndex)
                ops_.total[iData.cellIndices[iIndex]] +=
                  pow(iData.cellWeights[iIndex]*distanceOP*bondOPij*relOPij, pValue);
              for(size_t jIndex = 0; jIndex < jDataSize; ++jIndex)
                ops_.total[jData.cellIndices[jIndex]] +=
                  pow(jData.cellWeights[jIndex]*distanceOP*bondOPji*relOPji, pValue);
            }
            else
            {
              for(size_t iIndex = 0; iIndex < iDataSize; ++iIndex)
                ops_.total[iData.cellIndices[iIndex]] +=
                  iData.cellWeights[iIndex]*distanceOP*bondOPij*relOPij;
              for(size_t jIndex = 0; jIndex < jDataSize; ++jIndex)
                ops_.total[jData.cellIndices[jIndex]] +=
                  jData.cellWeights[jIndex]*distanceOP*bondOPji*relOPji;
            }
          }
          else if((confij.types.first == LINEAR_SYMMETRIC && confij.types.second == LINEAR_SYMMETRIC) ||
                  (confij.types.first == PLANAR_SYMMETRIC && confij.types.second == PLANAR_SYMMETRIC))
          {
            // Axisymmetric case with symmetry plane: Bond orientation and relative orientation are both angles
            Real deltaBondij = confij.bondOrientation.w - 
                               parameters.means[iGroup][iPeak].bondOrientation.w;
            Real deltaBondji = confji.bondOrientation.w - 
                               parameters.means[iGroup][iPeak].bondOrientation.w;
            Real bondOPij = parameters.normalizationFactors[iGroup][iPeak].bondOrientation*
                            exp(parameters.concentrationParameters[iGroup][iPeak].bondOrientation*
                                cos(2.0*deltaBondij));  // 2 due to symmetry
            Real bondOPji = parameters.normalizationFactors[iGroup][iPeak].bondOrientation*
                            exp(parameters.concentrationParameters[iGroup][iPeak].bondOrientation*
                                cos(2.0*deltaBondji));  // 2 due to symmetry
            if(pValue > 0.0)
            {
              for(size_t iIndex = 0; iIndex < iDataSize; ++iIndex)
                ops_.bondOrientation[iData.cellIndices[iIndex]] += 
                  pow(iData.cellWeights[iIndex]*distanceOP*bondOPij, pValue);
              for(size_t jIndex = 0; jIndex < jDataSize; ++jIndex)
                ops_.bondOrientation[jData.cellIndices[jIndex]] += 
                  pow(jData.cellWeights[jIndex]*distanceOP*bondOPji, pValue);
            }
            else
            {
              for(size_t iIndex = 0; iIndex < iDataSize; ++iIndex)
                ops_.bondOrientation[iData.cellIndices[iIndex]] += 
                  iData.cellWeights[iIndex]*distanceOP*bondOPij;
              for(size_t jIndex = 0; jIndex < jDataSize; ++jIndex)
                ops_.bondOrientation[jData.cellIndices[jIndex]] += 
                  jData.cellWeights[jIndex]*distanceOP*bondOPji;
            }

            Real deltaRelij = confij.relativeOrientation.w -
                              parameters.means[iGroup][iPeak].relativeOrientation.w;
            Real deltaRelji = confji.relativeOrientation.w -
                              parameters.means[iGroup][iPeak].relativeOrientation.w;
            Real relOPij = parameters.normalizationFactors[iGroup][iPeak].relativeOrientation*
                           exp(parameters.concentrationParameters[iGroup][iPeak].relativeOrientation*
                               cos(2.0*deltaRelij));  // 2 due to symmetry
            Real relOPji = parameters.normalizationFactors[iGroup][iPeak].relativeOrientation*
                           exp(parameters.concentrationParameters[iGroup][iPeak].relativeOrientation*
                               cos(2.0*deltaRelji));  // 2 due to symmetry
            if(pValue > 0.0)
            {
              for(size_t iIndex = 0; iIndex < iDataSize; ++iIndex)
                ops_.relativeOrientation[iData.cellIndices[iIndex]] += 
                  pow(iData.cellWeights[iIndex]*distanceOP*relOPij, pValue);
              for(size_t jIndex = 0; jIndex < jDataSize; ++jIndex)
                ops_.relativeOrientation[jData.cellIndices[jIndex]] += 
                  pow(jData.cellWeights[jIndex]*distanceOP*relOPji, pValue);
            }
            else
            {
              for(size_t iIndex = 0; iIndex < iDataSize; ++iIndex)
                ops_.relativeOrientation[iData.cellIndices[iIndex]] += 
                  iData.cellWeights[iIndex]*distanceOP*relOPij;
              for(size_t jIndex = 0; jIndex < jDataSize; ++jIndex)
                ops_.relativeOrientation[jData.cellIndices[jIndex]] += 
                  jData.cellWeights[jIndex]*distanceOP*relOPji;
            }
            // Total OPs added 07-14-10
            if(pValue > 0.0)
            {
              for(size_t iIndex = 0; iIndex < iDataSize; ++iIndex)
                ops_.total[iData.cellIndices[iIndex]] +=
                  pow(iData.cellWeights[iIndex]*distanceOP*bondOPij*relOPij, pValue);
              for(size_t jIndex = 0; jIndex < jDataSize; ++jIndex)
                ops_.total[jData.cellIndices[jIndex]] +=
                  pow(jData.cellWeights[jIndex]*distanceOP*bondOPji*relOPji, pValue);
            }
            else {
              for(size_t iIndex = 0; iIndex < iDataSize; ++iIndex)
                ops_.total[iData.cellIndices[iIndex]] +=
                  iData.cellWeights[iIndex]*distanceOP*bondOPij*relOPij;
              for(size_t jIndex = 0; jIndex < jDataSize; ++jIndex)
                ops_.total[jData.cellIndices[jIndex]] +=
                  jData.cellWeights[jIndex]*distanceOP*bondOPji*relOPji;
            }
          }
          else
          {
            std::cerr << "Error in CrystalOrderParameters::calculate(): "
                      << "Combination of molecule types not implemented" << std::endl;
            return;
          }
          // Internal DOF order parameters
          // Note that this assumes that both molecules have the same types and number
          // of internal DOFs

          for(size_t iDOF = 0; iDOF < numIntDOFs; ++iDOF)
          {
            switch(confij.internalDOFs.first[iDOF].type)
            {
              case DISTANCE:
              {
                // Gaussian
                Real deltaIntDOFij1 = confij.internalDOFs.first[iDOF].value -
                                      parameters.means[iGroup][iPeak].internalDOFs.first[iDOF].value;
                Real deltaIntDOFij2 = confij.internalDOFs.second[iDOF].value -
                                      parameters.means[iGroup][iPeak].internalDOFs.second[iDOF].value;
                Real intOPij1 = parameters.normalizationFactors[iGroup][iPeak].internalDOFs.first[iDOF]*
                    exp(-0.5*parameters.concentrationParameters[iGroup][iPeak].internalDOFs.first[iDOF]*
                         deltaIntDOFij1*deltaIntDOFij1);
                Real intOPij2 = parameters.normalizationFactors[iGroup][iPeak].internalDOFs.second[iDOF]*
                    exp(-0.5*parameters.concentrationParameters[iGroup][iPeak].internalDOFs.second[iDOF]*
                         deltaIntDOFij2*deltaIntDOFij2);
                if(pValue > 0.0)
                {
                  for(size_t iIndex = 0; iIndex < iDataSize; ++iIndex)
                    ops_.internal[iData.cellIndices[iIndex]][iDOF] += 
                      pow(iData.cellWeights[iIndex]*distanceOP*intOPij1*intOPij2, pValue);
                }
                else
                {
                  for(size_t iIndex = 0; iIndex < iDataSize; ++iIndex)
                    ops_.internal[iData.cellIndices[iIndex]][iDOF] += 
                      iData.cellWeights[iIndex]*distanceOP*intOPij1*intOPij2;
                }
                
                Real deltaIntDOFji1 = confji.internalDOFs.first[iDOF].value -
                                      parameters.means[iGroup][iPeak].internalDOFs.first[iDOF].value;
                Real deltaIntDOFji2 = confji.internalDOFs.second[iDOF].value -
                                      parameters.means[iGroup][iPeak].internalDOFs.second[iDOF].value;
                Real intOPji1 = parameters.normalizationFactors[iGroup][iPeak].internalDOFs.first[iDOF]*
                    exp(-0.5*parameters.concentrationParameters[iGroup][iPeak].internalDOFs.first[iDOF]*
                         deltaIntDOFji1*deltaIntDOFji1);
                Real intOPji2 = parameters.normalizationFactors[iGroup][iPeak].internalDOFs.second[iDOF]*
                    exp(-0.5*parameters.concentrationParameters[iGroup][iPeak].internalDOFs.second[iDOF]*
                         deltaIntDOFji2*deltaIntDOFji2);
                if(pValue > 0.0)
                {
                  for(size_t jIndex = 0; jIndex < jDataSize; ++jIndex)
                    ops_.internal[jData.cellIndices[jIndex]][iDOF] += 
                      pow(jData.cellWeights[jIndex]*distanceOP*intOPji1*intOPji2, pValue);
                }
                else
                {
                  for(size_t jIndex = 0; jIndex < jDataSize; ++jIndex)
                    ops_.internal[jData.cellIndices[jIndex]][iDOF] += 
                      jData.cellWeights[jIndex]*distanceOP*intOPji1*intOPji2;
                }
              }
                break;
              case ANGLE:
              {
                // von Mises
                Real deltaIntDOFij1 = confij.internalDOFs.first[iDOF].value -
                                      parameters.means[iGroup][iPeak].internalDOFs.first[iDOF].value;
                Real deltaIntDOFij2 = confij.internalDOFs.second[iDOF].value -
                                      parameters.means[iGroup][iPeak].internalDOFs.second[iDOF].value;
                Real intOPij1 = parameters.normalizationFactors[iGroup][iPeak].internalDOFs.first[iDOF]*
                    exp(parameters.concentrationParameters[iGroup][iPeak].internalDOFs.first[iDOF]*
                        cos(deltaIntDOFij1));
                Real intOPij2 = parameters.normalizationFactors[iGroup][iPeak].internalDOFs.second[iDOF]*
                    exp(parameters.concentrationParameters[iGroup][iPeak].internalDOFs.second[iDOF]*
                        cos(deltaIntDOFij2));
                if(pValue > 0.0)
                {
                  for(size_t iIndex = 0; iIndex < iDataSize; ++iIndex)
                    ops_.internal[iData.cellIndices[iIndex]][iDOF] += 
                      pow(iData.cellWeights[iIndex]*distanceOP*intOPij1*intOPij2, pValue);
                }
                else
                {
                  for(size_t iIndex = 0; iIndex < iDataSize; ++iIndex)
                    ops_.internal[iData.cellIndices[iIndex]][iDOF] += 
                      iData.cellWeights[iIndex]*distanceOP*intOPij1*intOPij2;
                }

                Real deltaIntDOFji1 = confji.internalDOFs.first[iDOF].value -
                                      parameters.means[iGroup][iPeak].internalDOFs.first[iDOF].value;
                Real deltaIntDOFji2 = confji.internalDOFs.second[iDOF].value -
                                      parameters.means[iGroup][iPeak].internalDOFs.second[iDOF].value;
                Real intOPji1 = parameters.normalizationFactors[iGroup][iPeak].internalDOFs.first[iDOF]*
                    exp(parameters.concentrationParameters[iGroup][iPeak].internalDOFs.first[iDOF]*
                        cos(deltaIntDOFji1));
                Real intOPji2 = parameters.normalizationFactors[iGroup][iPeak].internalDOFs.second[iDOF]*
                    exp(parameters.concentrationParameters[iGroup][iPeak].internalDOFs.second[iDOF]*
                        cos(deltaIntDOFji2));
                if(pValue > 0.0)
                {
                  for(size_t jIndex = 0; jIndex < jDataSize; ++jIndex)
                    ops_.internal[jData.cellIndices[jIndex]][iDOF] += 
                      pow(jData.cellWeights[jIndex]*distanceOP*intOPji1*intOPji2, pValue);
                }
                else
                {
                  for(size_t jIndex = 0; jIndex < jDataSize; ++jIndex)
                    ops_.internal[jData.cellIndices[jIndex]][iDOF] += 
                      jData.cellWeights[jIndex]*distanceOP*intOPji1*intOPji2;
                }
              }
                break;
              case DIHEDRAL:
              {
                // von Mises with symmetry factor
                Real deltaIntDOFij1 = confij.internalDOFs.first[iDOF].value -
                                      parameters.means[iGroup][iPeak].internalDOFs.first[iDOF].value;
                deltaIntDOFij1 *= (Real)confij.internalDOFs.first[iDOF].symmetryNumber;
                Real deltaIntDOFij2 = confij.internalDOFs.second[iDOF].value -
                                      parameters.means[iGroup][iPeak].internalDOFs.second[iDOF].value;
                deltaIntDOFij2 *= (Real)confij.internalDOFs.second[iDOF].symmetryNumber;
                Real intOPij1 = parameters.normalizationFactors[iGroup][iPeak].internalDOFs.first[iDOF]*
                    exp(parameters.concentrationParameters[iGroup][iPeak].internalDOFs.first[iDOF]*
                        cos(deltaIntDOFij1));
                Real intOPij2 = parameters.normalizationFactors[iGroup][iPeak].internalDOFs.second[iDOF]*
                    exp(parameters.concentrationParameters[iGroup][iPeak].internalDOFs.second[iDOF]*
                        cos(deltaIntDOFij2));
                if(pValue > 0.0)
                {
                  for(size_t iIndex = 0; iIndex < iDataSize; ++iIndex)
                    ops_.internal[iData.cellIndices[iIndex]][iDOF] += 
                      pow(iData.cellWeights[iIndex]*distanceOP*intOPij1*intOPij2, pValue);
                }
                else
                {
                  for(size_t iIndex = 0; iIndex < iDataSize; ++iIndex)
                    ops_.internal[iData.cellIndices[iIndex]][iDOF] += 
                      iData.cellWeights[iIndex]*distanceOP*intOPij1*intOPij2;
                }

                Real deltaIntDOFji1 = confji.internalDOFs.first[iDOF].value -
                                      parameters.means[iGroup][iPeak].internalDOFs.first[iDOF].value;
                deltaIntDOFji1 *= (Real)confji.internalDOFs.first[iDOF].symmetryNumber;
                Real deltaIntDOFji2 = confji.internalDOFs.second[iDOF].value -
                                      parameters.means[iGroup][iPeak].internalDOFs.second[iDOF].value;
                deltaIntDOFji2 *= (Real)confji.internalDOFs.second[iDOF].symmetryNumber;
                Real intOPji1 = parameters.normalizationFactors[iGroup][iPeak].internalDOFs.first[iDOF]*
                    exp(parameters.concentrationParameters[iGroup][iPeak].internalDOFs.first[iDOF]*
                        cos(deltaIntDOFji1));
                Real intOPji2 = parameters.normalizationFactors[iGroup][iPeak].internalDOFs.second[iDOF]*
                    exp(parameters.concentrationParameters[iGroup][iPeak].internalDOFs.second[iDOF]*
                        cos(deltaIntDOFji2));
                if(pValue > 0.0)
                {
                  for(size_t jIndex = 0; jIndex < jDataSize; ++jIndex)
                    ops_.internal[jData.cellIndices[jIndex]][iDOF] += 
                      pow(jData.cellWeights[jIndex]*distanceOP*intOPji1*intOPji2, pValue);
                }
                else
                {
                  for(size_t jIndex = 0; jIndex < jDataSize; ++jIndex)
                    ops_.internal[jData.cellIndices[jIndex]][iDOF] += 
                      jData.cellWeights[jIndex]*distanceOP*intOPji1*intOPji2;
                }
              }
                break;
            }
          }
        } // End of loop over peaks
      } // End of loop over groups
    } // End of loop over neighbors
  } // End of main loop over molecules
  
  //defined via cellweight / cellVolume, averaged over all 64
  //#define rho_liq 0.000145363235
  //defined by total molecules over total volume, should be the same but its not. 
  //#define rho_liq 0.00058138
  Real cellVolume = useGrid?(system.lattice().volume()/ops_.grid.numCells()):1.0;
  for(size_t i = 0; i < numCells; i++) // Won't run if numCells = 0 (per molecule case)
  {
    Real totalWeight = totalWeights[i];
    if(totalWeight == 0) continue;
    if(pValue > 0.0)
    {
      ops_.distance[i] = pow(ops_.distance[i], 1.0/pValue);
      ops_.bondOrientation[i] = pow(ops_.bondOrientation[i], 1.0/pValue);
      ops_.relativeOrientation[i] = pow(ops_.relativeOrientation[i], 1.0/pValue);
      for(size_t j = 0; j < ops_.internal[i].size(); ++j)
        ops_.internal[i][j] = pow(ops_.internal[i][j], 1.0/pValue);
      ops_.total[i] = pow(ops_.total[i], 1.0/pValue);
    }
    else if (system.lattice().fakeNOPBC)
    {
    continue; //nothing to do!
    }
    else
    {
     ops_.distance[i] /= totalWeight;
     ops_.bondOrientation[i] /= totalWeight;
     ops_.relativeOrientation[i] /= totalWeight;
     ops_.total[i] /= totalWeight;
     for(size_t j = 0; j < ops_.internal[i].size(); ++j)
       ops_.internal[i][j] /= totalWeight;
    }
  }
  if(useLogs)
  {
  for(size_t i = 0; i < ops_.distance.size(); ++i)
        ops_.distance[i] = tanh(TanhModifier*ops_.distance[i]);
    for(size_t i = 0; i < ops_.bondOrientation.size(); ++i)
      ops_.bondOrientation[i] = tanh(TanhModifier*ops_.bondOrientation[i]);
    for(size_t i = 0; i < ops_.total.size(); ++i)
      ops_.total[i] = tanh(TanhModifier*ops_.total[i]);
    for(size_t i = 0; i < ops_.localDensity.size(); ++i) //not implemented
    for(size_t i = 0; i < ops_.internal.size(); ++i)
      for(size_t j = 0; j < ops_.internal[i].size(); ++j)
        ops_.internal[i][j] = tanh(TanhModifier*ops_.internal[i][j]);

//old log definition
//    for(size_t i = 0; i < ops_.distance.size(); ++i)
//      ops_.distance[i] = (ops_.distance[i] > 0.0)?log(ops_.distance[i]):-99999;
//    for(size_t i = 0; i < ops_.bondOrientation.size(); ++i)
//      ops_.bondOrientation[i] = (ops_.bondOrientation[i] > 0.0)?log(ops_.bondOrientation[i]):-99999;
//    for(size_t i = 0; i < ops_.relativeOrientation.size(); ++i)
//      ops_.relativeOrientation[i] = (ops_.relativeOrientation[i] > 0.0)?log(ops_.relativeOrientation[i]):-99999;
//    for(size_t i = 0; i < ops_.total.size(); ++i)
//      ops_.total[i] = (ops_.total[i] > 0.0)?log(ops_.total[i]):-99999;
//    for(size_t i = 0; i < ops_.localDensity.size(); ++i)
//      ops_.localDensity[i] = (ops_.localDensity[i] > 0.0)?log(ops_.localDensity[i]):-99999;
//    for(size_t i = 0; i < ops_.internal.size(); ++i)
//      for(size_t j = 0; j < ops_.internal[i].size(); ++j)
//        ops_.internal[i][j] = (ops_.internal[i][j] > 0.0)?log(ops_.internal[i][j]):-99999;
  }

//ridiculousness 
////JAKE 2025 NEW OP TESTING
////double shellBond[system.size()][system.size()];
////double shellBondPnorm = 0.0;
////double shellBondPnormSmall = 0.0;
////double pnormOP = 0.0;
////double pnormOPSmall = 0.0;
//////Vector scaledLatice = system.lattice();
////for (size_t i=0; i<system.size(); ++i)
////for (size_t j=0; j<system.size(); ++j)
////  shellBond[i][j] = 0.0;
////for (size_t i=0; i< system.size(); ++i)
////{
////  pnormOP += pow(ops_.bondOrientation[i], 10.0);
////  pnormOPSmall += pow(ops_.bondOrientation[i], 2.0);
////  for (size_t j = i+1; j< system.size(); ++j)
////    {
////         Vector3D temp = system.lattice().difference(system[i].position, system[j].position);
////         double LatticeXNorm = sqrt(norm2(system.lattice().latticeVector(0)));
////         double LatticeYNorm = sqrt(norm2(system.lattice().latticeVector(1)));
////         double LatticeZNorm = sqrt(norm2(system.lattice().latticeVector(2)));
////         temp[0] /= LatticeXNorm;
////         temp[1] /= LatticeYNorm;
////         temp[2] /= LatticeZNorm;
////         shellBond[i][j] = shellBond[j][i] = norm2(temp)*ops_.bondOrientation[i]*ops_.bondOrientation[j];
//////         std::cout << i << " " << j << " " << shellBond[i][j] << std::endl;
////         shellBondPnorm += 2*pow( shellBond[i][j] ,10.0 );
////         shellBondPnormSmall += 2*pow(shellBond[i][j], 2.0);
////    }
////}
//////for (int i =0; i<system.size(); ++i)
//////  for (int j=0; j<system.size(); ++j)
////         Vector3D temp = system.lattice().difference(system[0].position, system[51].position);
////         double LatticeXNorm = sqrt(norm2(system.lattice().latticeVector(0)));
////         double LatticeYNorm = sqrt(norm2(system.lattice().latticeVector(1)));
////         double LatticeZNorm = sqrt(norm2(system.lattice().latticeVector(2)));
////         temp[0] /= LatticeXNorm;
////         temp[1] /= LatticeYNorm;
////         temp[2] /= LatticeZNorm;
////
//        
////std::cout << "i " << ops_.bondOrientation[0]  << std::endl;
////std::cout << "i " << ops_.bondOrientation[51]  << std::endl;
////std::cout << "dij2 " << norm2(temp) << std::endl;
////std::cout << pow(shellBondPnorm, 1.0/10.0) << std::endl;
////std::cout << pow(shellBondPnormSmall, 1.0/2.0) << std::endl;
////std::cout << pow(pnormOP, 1.0/10.0) << std::endl;
////std::cout << pow(pnormOPSmall, 1.0/2.0) << std::endl;
}
//JAKE
void OrderParameterCalculator::printQuatProduct(System<PointMolecule> const &system, CrystalDistributionParameters const &parameters){

size_t const numGroups = parameters.means.size();

//std::ifstream xs; xs.open("/rs1/researchers/e/eesantis/jpmckibb/x");
//std::ifstream ys; ys.open("/rs1/researchers/e/eesantis/jpmckibb/y");
//std::ifstream zs; zs.open("/rs1/researchers/e/eesantis/jpmckibb/z");

for (size_t i=0;i<system.size();++i){
  Real finalBop=0;
for (size_t j=0;j<system.size();++j){
  Real bopSum=0;
  Vector3D fakeDistance = system[j].position - system[i].position;
  Quaternion fakeDistanceQuat(fakeDistance);
  Quaternion fakeOr = ~(system[i].orientation)*fakeDistanceQuat*system[i].orientation;
  if (i!=j) fakeOr/=norm(fakeDistance);
  for (size_t iGroup=0; iGroup<numGroups; ++iGroup){
    size_t numPeaks = parameters.means[iGroup].size(); 
  for (size_t iPeak=0; iPeak<numPeaks; iPeak++){
    Real deltaDistance = norm(fakeDistance) - parameters.means[iGroup][iPeak].distance;
    //Real bondDotij = vector(fakeOr)*vector(parameters.means[iGroup][iPeak].bondOrientation);
    Real bondDotij = dot(fakeOr,parameters.means[iGroup][iPeak].bondOrientation);
    Real bondOPij = parameters.normalizationFactors[iGroup][iPeak].bondOrientation*exp(parameters.concentrationParameters[iGroup][iPeak].bondOrientation*bondDotij);
//std::cout << bondOPij << std::endl;
Real distanceOP = parameters.normalizationFactors[iGroup][iPeak].distance*exp(-0.5*parameters.concentrationParameters[iGroup][iPeak].distance*deltaDistance*deltaDistance);
//    std::cout << bopSum << std::endl;
    bopSum+=(distanceOP*bondOPij);


      
}}
finalBop+=bopSum;
}
//std::cout << finalBop << std::endl;
}

//    std::string x_val;    
//    std::string y_val;
//    std::string z_val;
//    if ( xs.is_open() ) { // always check whether the file is open
//        xs >> x_val; // pipe file's content into stream
//}
//    if ( ys.is_open() ) { // always check whether the file is open
//        ys >> y_val; // pipe file's content into stream
//}
//    if ( zs.is_open() ) { // always check whether the file is open
//        zs >> z_val; // pipe file's content into stream
//}
//
//
//     Real x_dist = std::strtod(x_val.c_str(),0);  
//     Real y_dist = std::strtod(y_val.c_str(),0);  
//     Real z_dist = std::strtod(z_val.c_str(),0);  




////RelativeConfiguration const confij(system[i], system[j], system.lattice());
//Vector3D fakeDistance =  system[j].position - system[i].position;
////Vector3D fakeDistance(x_dist, y_dist, z_dist);
////Quaternion bond(0.0, x_dist, y_dist, z_dist);
////Quaternion rhat(0.0, fakeDistance.x, fakeDistance.y, fakeDistance.z);
//Quaternion rhat(fakeDistance);
////rotateVector(~system[i].orientation, fakeDistance);
//Quaternion bondProd = ~(system[i].orientation)*rhat*system[i].orientation / norm(fakeDistance); 
////Real normFakeDistance = norm(fakeDistance);
////Quaternion bondProd = fakeDistance / normFakeDistance;
////Quaternion bondProd = ~(system[i].orientation)*rhat;
//
//std::cout << bondProd.w << " " << bondProd.x << " " << bondProd.y << " " << bondProd.z <<std::endl;
//}
//else std::cout << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 <<std::endl;
//
//// Vector3D bond = lattice.difference(first.position, second.position); fake distance 
////    rotateVector(~first.orientation, bond); // Project onto first molecule frame
////    distance = norm(bond);
////    bondOrientation = bond/distance;
//
//
//
//
//
//}
//}

//for (size_t i=0; i<system.size(); ++i){
//std::cout << system[i].orientation.w << " " << system[i].orientation.x << " " << system[i].orientation.y << " " << system[i].orientation.z << std::endl; 
//
//}
//





//compare quats
//for (size_t i=0;i<system.size();++i){
//for (size_t j=0;j<system.size();++j){
//  if (i!=j){
//  Quaternion quatprod = ~(system[i].orientation)*system[j].orientation; 
////  RelativeConfiguration const confij(system[i], system[j], system.lattice());
//  std::cout << quatprod.w << " " << quatprod.x << " " << quatprod.y  << " " << quatprod.z << std::endl;
//}
//else std::cout << 1.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 <<std::endl;
//}
//}

//
//number of groups and peaks in the groups 
//size_t const numGroups = parameters.means.size();
//////dirty dirty dirty
//    std::ifstream myfile; myfile.open("/rs1/researchers/e/eesantis/jpmckibb/temp");
////
//     
//  for(size_t i=0;i<system.size();++i){
//    Real finalRelOP=0;
//    for(size_t j=0;j<system.size();j++){
//    std::string mystring;
//    if ( myfile.is_open() ) { // always check whether the file is open
//    myfile >> mystring; // pipe file's content into stream
//}
//        Real mydistance = std::strtod(mystring.c_str(),0)*10.0;       
//        Real fakeDistance =  norm(system[j].position - system[i].position);
//     if (i!=j){
//     Real relOP=0;
//    ////replace confij.distance with ones I'm reading in
//    RelativeConfiguration const confij(system[i], system[j], system.lattice());
//       for (int iGroup=0; iGroup<numGroups; ++iGroup){
//        size_t numPeaks = parameters.means[iGroup].size(); 
//        for (int iPeak=0; iPeak<numPeaks; iPeak++){
//            Real deltaDistance = fakeDistance - parameters.means[iGroup][iPeak].distance;
//
//            Real relDotij = dot(parameters.means[iGroup][iPeak].relativeOrientation,
//                                confij.relativeOrientation);
//            Real relOPij = parameters.normalizationFactors[iGroup][iPeak].relativeOrientation*
//                           exp(parameters.concentrationParameters[iGroup][iPeak].relativeOrientation*
//                               relDotij*relDotij); 
//            Real distanceOP = parameters.normalizationFactors[iGroup][iPeak].distance*exp(-0.5*parameters.concentrationParameters[iGroup][iPeak].distance*deltaDistance*deltaDistance);
//            relOP+=(relOPij*distanceOP);
//
//
//
//
////            distanceOP += parameters.normalizationFactors[iGroup][iPeak].distance*exp(-0.5*parameters.concentrationParameters[iGroup][iPeak].distance*deltaDistance*deltaDistance);
//            }
//            }
//finalRelOP+=relOP;
////      std::cout << relOP << std::endl;
//        //std::cout << mydistance << " "  << norm(system[j].position - system[i].position) << std::endl;
//}
//         else int i = 0;//std::cout << i << " " << 0.000 << std::endl;
////  std::cout << i << " " << distanceOP << std::endl;
////
//}
//std::cout << i << " " << finalRelOP << std::endl;
////
//}
}


// Private functions

CellData const OrderParameterCalculator::cellData(Vector3D const &position, Lattice const &lattice)
{
  CellData cellData;
  if(lattice.numDimensions() != 3)
  {
    std::cerr << "Error in OrderParameterCalculator::cellData - invalid number of dimensions " << std::endl;
    return cellData;
  }

//don't account for pbc, just pass every cell as an index and weight, we only run this script sparingly, doesn't need to be memory efficient at all
if (lattice.fakeNOPBC)
{
  if (ops_.switchWidth == 0.0)
  {
  std::cerr << "Error in OrderParameterCalculator::cellData - need finite switch width the pbc disabled." << std::endl;
  return cellData;

  }
  //won't work at all if not using a cubic system
  Vector3D xlent =lattice.latticeVector(0);
  Vector3D ylent =lattice.latticeVector(1);
  Vector3D zlent =lattice.latticeVector(2);  
  //extract box length
  Real xlen = xlent[0];
  Real ylen = ylent[1];
  Real zlen = zlent[2];
  //size of grid
  size_t const nx = ops_.grid.x;
  size_t const ny = ops_.grid.y;
  size_t const nz = ops_.grid.z;
  //size of each cube
  Real m=xlen/Real(nx);
  Real n=ylen/Real(ny);
  Real o=zlen/Real(nz);

  for(size_t inx = 0; inx < nx; ++inx)
  {
  for(size_t iny = 0; iny < ny; ++iny)
  {
  for(size_t inz = 0; inz < nz; ++inz)
    {
    cellData.cellIndices.push_back(inx + nx*(iny + ny*inz));
    cellData.cellWeights.push_back(
             (1-tanh((-xlen/2.0 + inx*m -position.x)/ops_.switchWidth))*(1-tanh((position.x+xlen/2.0 - (inx+1)*m )/ops_.switchWidth))*
             (1-tanh((-ylen/2.0 + iny*n -position.y)/ops_.switchWidth))*(1-tanh((position.y+ylen/2.0 - (iny+1)*n )/ops_.switchWidth))*
             (1-tanh((-zlen/2.0 + inz*o -position.z)/ops_.switchWidth))*(1-tanh((position.z+zlen/2.0 - (inz+1)*o )/ops_.switchWidth))/64);
    }
  }
  }   
return cellData;
}
  Vector3D const pos = position - ops_.grid.origin; // For compatibility with NAMD
  Vector3D scaled(lattice.reciprocalVector(0)*pos,
                  lattice.reciprocalVector(1)*pos,
                  lattice.reciprocalVector(2)*pos);
  // Numbers of cells

  size_t const nx = ops_.grid.x;
  size_t const ny = ops_.grid.y;
  size_t const nz = ops_.grid.z;

  // Indices that contribute to average
  std::vector<size_t> indx, indy, indz;
  // In case molecule has drifted out of simulation box
  scaled.x -= 1 + floor(scaled.x - 0.5);
  scaled.y -= 1 + floor(scaled.y - 0.5);
  scaled.z -= 1 + floor(scaled.z - 0.5);

//JAKE 2025
  //my explanation of this function
  //if u = scaled.w + 0.5
  //simplifies to  floor(nw*(u - floor(u)))
  //so essentially we're finding the matissa of u with u-floor(u).  We just need it to be posistion, hence the +0.5 business. 
  //So then when me multiply by length, we can figure out the cell number that way by flooring it. 
  // Indices of nearest cell center
  // 
  size_t ixc = (size_t)floor(nx*(scaled.x + 0.5 - floor(scaled.x + 0.5)));
  size_t iyc = (size_t)floor(ny*(scaled.y + 0.5 - floor(scaled.y + 0.5)));
  size_t izc = (size_t)floor(nz*(scaled.z + 0.5 - floor(scaled.z + 0.5)));
  if(ixc > nx - 1) ixc = nx - 1;
  if(iyc > ny - 1) iyc = ny - 1;
  if(izc > nz - 1) izc = nz - 1;
  indx.push_back(ixc);
  indy.push_back(iyc);
  indz.push_back(izc);
  // Nearest cell center in scaled coordinates and difference
  Vector3D const center(((Real)ixc+0.5)/(Real)nx - 0.5,
                        ((Real)iyc+0.5)/(Real)ny - 0.5, 
                        ((Real)izc+0.5)/(Real)nz - 0.5);

  Vector3D delta = scaled - center;

  // Partial contributions to weights
  std::vector<Real> fx, fy, fz;

  Real scaledDiff = (ops_.switchWidth > 0.0)?(fabs(nx*delta.x) - 0.5)/ops_.switchWidth:-2.0;
  if(scaledDiff > -1.0)
  {
    if(delta.x > 0.0)
    {
      if(ixc == nx - 1 ) indx.push_back(0); //this is where it wraps around
      else indx.push_back(ixc + 1);
    }
    else
    {
      if(ixc == 0 && !lattice.fakeNOPBC) indx.push_back(nx - 1);//also here
      else indx.push_back(ixc - 1);
    }
    fx.push_back(0.5 + 0.25*scaledDiff*(scaledDiff*scaledDiff - 3));
    fx.push_back(1.0 - fx[0]);
  }
  else{
    fx.push_back(1.0);
   }
  scaledDiff = (ops_.switchWidth > 0.0)?(fabs(ny*delta.y) - 0.5)/ops_.switchWidth:-2.0;
  if(scaledDiff > -1.0)
  {
    if(delta.y > 0.0)
    {
      if(iyc == ny - 1 ) indy.push_back(0);
      else indy.push_back(iyc + 1);
    }
    else
    {
      if(iyc == 0 ) indy.push_back(ny - 1);
      else indy.push_back(iyc - 1);
    }
    fy.push_back(0.5 + 0.25*scaledDiff*(scaledDiff*scaledDiff - 3));
    fy.push_back(1.0 - fy[0]);
  }
  else
    fy.push_back(1.0);

  scaledDiff = (ops_.switchWidth > 0.0)?(fabs(nz*delta.z) - 0.5)/ops_.switchWidth:-2.0;
  if(scaledDiff > -1.0)
  {
    if(delta.z > 0.0)
    {
      if(izc == nz - 1 ) indz.push_back(0);
      else indz.push_back(izc + 1);
    }
    else
    {
      if(izc == 0 ) indz.push_back(nz - 1);
      else indz.push_back(izc - 1);
    }
    fz.push_back(0.5 + 0.25*scaledDiff*(scaledDiff*scaledDiff - 3));
    fz.push_back(1.0 - fz[0]);
  }
  else
    fz.push_back(1.0);


  for(size_t inx = 0; inx < indx.size(); ++inx)
  for(size_t iny = 0; iny < indy.size(); ++iny)
  for(size_t inz = 0; inz < indz.size(); ++inz)
    {
    cellData.cellIndices.push_back(indx[inx] + nx*(indy[iny] + ny*indz[inz]));
    cellData.cellWeights.push_back(fx[inx]*fy[iny]*fz[inz]);
    }
  return cellData;

}

