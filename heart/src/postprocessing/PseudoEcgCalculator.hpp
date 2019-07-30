/*

Copyright (c) 2005-2019, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef _PSEUDOECGALCALCULATOR_HPP_
#define _PSEUDOECGALCALCULATOR_HPP_

#include "AbstractFunctionalCalculator.hpp"
#include "AbstractTetrahedralMesh.hpp"
#include "ChastePoint.hpp"
#include "UblasCustomFunctions.hpp"
#include "Hdf5DataReader.hpp"

/**
 * This class implements a pseudo-ECG calculator. It is a concrete class of AbstractFunctionalCalculator.
 * The pseudo-ECG is defined as the integral over the mesh of the following integrand:
 *
 * - D * grad (solution) dot grad (1/r)
 *
 *  where D is a diffusion coefficient and r is the distance between a recording electrode
 *  and a given point in the mesh.
 *
 *  References for the formula that defines the pseudo-ECG:
 *
 *  Gima K, Rudy Y Circ Res (2002) 90:889-896 (equation 1)
 *  Baher et al. Am J Physiol (2007) 292:H180-H189 (equation 5)
 */
/*
 * NB I don't think that this class needs to be templated over PROBLEM_DIM,
 *  since even when working with Bidomain we only integrate one thing at once,
 *  and hence still need to call AbstractFunctionalCalculator<ELEMENT_DIM, SPACE_DIM, 1>.
 *  See PostProcessingWriter.cpp line ~141. The fact that it is only explicitly
 *  instantiated (at the bottom of the cpp) for cases with 1, would agree with this.
 *
 *  But at the same time it probably isn't worth breaking users' code by changing it.
 *  Might also make it easier to copy and paste for any new class that does a
 *  fancier integral with more quantities.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
class PseudoEcgCalculator : public AbstractFunctionalCalculator<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>
{
private:

    friend class TestPseudoEcgCalculator;//for testing

    Hdf5DataReader* mpDataReader; /**< An HDF5 reader from which to get the solution*/
    unsigned mNumberOfNodes; /**< Number of nodes in the mesh (got from the data reader)*/
    unsigned mNumTimeSteps;/**< Number of time steps in the simulation (got from the data reader)*/
    AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>& mrMesh;/**< A mesh used by the calculator*/
    ChastePoint<SPACE_DIM> mProbeElectrode; /**<The point from where we want to calculate the pseudoECG*/
    double mDiffusionCoefficient;/**<The diffusion coefficient D*/
    std::string mVariableName;/**< the variable for which we want to calculate the pseudo ecg, defaults to "V"*/
    unsigned mTimestepStride; /**< The number of timesteps in a stride (so that we don't have to compute all the ECGs).  This defaults to 1.*/
    /**
     * @return the integrand.
     * The pseudo-ECG is defined as the integral over the mesh of the following integrand:
     *
     * - D * grad (solution) dot grad (1/r)
     *
     * @param rX The point in space
     * @param rU The unknown as a vector, u(i) = u_i
     * @param rGradU The gradient of the unknown as a matrix, rGradU(i,j) = d(u_i)/d(X_j)
     */
    double GetIntegrand(ChastePoint<SPACE_DIM>& rX,
                        c_vector<double,PROBLEM_DIM>& rU,
                        c_matrix<double,PROBLEM_DIM,SPACE_DIM>& rGradU);


    /**
     * Calculates the pseudo-ECG and returns its value at the given time step.
     *
     * @param timeStep the time step where we want to calculate the pseudo-ecg
     * @return the pseudo ECG at the given time step.
     *
     */
    double ComputePseudoEcgAtOneTimeStep(unsigned timeStep);

    /**
     * Whether we should not calculate the Pseudo ECG on this element.
     *
     * This method returns true if we are integrating voltage and the
     * tissue element is in the bath.
     *
     * @param rElement  the element of interest
     * @return  whether we should skip this element.
     */
    bool ShouldSkipThisElement(Element<ELEMENT_DIM,SPACE_DIM>& rElement);

public:

    /**
     * Constructor
     *
     * @param rMesh A reference to the mesh
     * @param rProbeElectrode The location of the recording electrode
     * @param rDirectory  The directory where the simulation results are stored
     * @param rHdf5FileName The file name  where the simulation results are stored
     * @param rVariableName  The name of the voltage variable (is V by default)
     * @param timestepStride The number of timesteps in a stride (so that we don't
     *        have to compute all the ECGs).  This defaults to 1.
     *
     */
    PseudoEcgCalculator(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>& rMesh,
                        const ChastePoint<SPACE_DIM>& rProbeElectrode,
                        const FileFinder& rDirectory,
                        const std::string& rHdf5FileName,
                        const std::string& rVariableName = "V",
                        unsigned timestepStride = 1);

    /**
     * Destructor.
     */
    ~PseudoEcgCalculator();


    /**
     * Sets the value of the diffusion coefficient (D)
     *
     * @param diffusionCoefficient The desired value of the diffusion coefficient
     *
     */
    void SetDiffusionCoefficient(double diffusionCoefficient);

    /**
     *
     * Calculates and writes the pseudo-ECG to file. the file will be named PseudoEcgFromElectrodeAt_x_y_z.dat,
     * where x,y,z are replaced by the location of the electrode.
     * It will contain one column of numbers, each being the pseudoECG at each time step.
     * It will be created  by the master processor into /output relaitive to where the output
     * directory is set (by the HeartConfig or by default)
     *
     */
    void WritePseudoEcg();
};

#endif //_PSEUDOECGALCALCULATOR_HPP_
