/*

Copyright (c) 2005-2025, University of Oxford.
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

#ifndef VOLTAGEINTERPOLATERONTOMECHANICSMESH_HPP_
#define VOLTAGEINTERPOLATERONTOMECHANICSMESH_HPP_

#include <vector>
#include <string>
#include "UblasIncludes.hpp"
#include "ReplicatableVector.hpp"
#include "TetrahedralMesh.hpp"
#include "QuadraticMesh.hpp"
#include "FineCoarseMeshPair.hpp"

/**
 *  A class which can be used to convert the voltage from an electrics
 *  (or electromechanics) simulation onto a coarser mechanics mesh, by interpolation. The
 *  class provides a method to interpolate values from the fine mesh 
 *  onto the coarse mesh and the main meythod to output a HDF5 file corresponding to nodes 
 *  on the mechanics mesh, and convert it to CMGUI output.
 */
template<unsigned DIM>
class VoltageInterpolaterOntoMechanicsMesh
{

private:

    /**< The fine mesh (values will be interpolated from this one)*/
    TetrahedralMesh<DIM,DIM>& mrElectricsMesh;
    /**< The coarse mesh (values will be interpolated onto this one)*/
    QuadraticMesh<DIM>& mrMechanicsMesh;
    /**< Mesh pair object. Initialized by the constructor using the two meshes passed in*/
    FineCoarseMeshPair<DIM>* mpMeshPair;

public:
    /**
     *  Constructor. Sets up the Mesh Pair object used for the interpolation later.
     *
     *  @param rElectricsMesh The electrics mesh
     *  @param rMechanicsMesh The mechanics mesh
     */
    VoltageInterpolaterOntoMechanicsMesh(TetrahedralMesh<DIM,DIM>& rElectricsMesh,
                                         QuadraticMesh<DIM>& rMechanicsMesh);

    /**
     * Destructor. 
     */
    ~VoltageInterpolaterOntoMechanicsMesh();
    
    /**
     * Does the actual interpolation. 
     * 
     * @param rValuesOnCoarseMesh It will be filled with the interpolated values. Memory must 
     *                            have been allocated before (this method checks): the vector must be of 
     *                            the same size as the number of nodes of the coarser mesh
     * @param rValuesOnFineMesh The values on the fine mesh to be interpolated onto the coarser mesh
     */
    void InterpolateOnCoarseMesh(std::vector<double>& rValuesOnCoarseMesh, ReplicatableVector& rValuesOnFineMesh);


    /**
     * The main method of the class. It reads the h5 files in directory (filename prefix inputFileNamePrefix)
     * It then considers all the variables whose name is specified in rVariableNames, does
     * the interplation using  InterpolateOnCoarseMesh, and creates output for all the variables over 
     * all the time steps included in the h5 file. It does so by first writing out to h5 
     * and then converting to CMGUI.
     *  
     *  @param rVariableNames vector of names of variables contained in the input h5 file and
     *                        that you want to be interpolated.
     *  @param directory Directory the voltage file is in
     *  @param inputFileNamePrefix Filename (without ".h5") of the electrics solution HDF5 file
     */
    void OutputToCmgui(std::vector<std::string>& rVariableNames, std::string directory, std::string inputFileNamePrefix);
};

#endif /*VOLTAGEINTERPOLATERONTOMECHANICSMESH_HPP_*/
