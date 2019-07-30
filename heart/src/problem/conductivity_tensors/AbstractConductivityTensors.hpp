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
#ifndef ABSTRACTCONDUCTIVITYTENSORS_HPP_
#define ABSTRACTCONDUCTIVITYTENSORS_HPP_

#include <vector>
#include <string>
#include <memory>
#include "UblasIncludes.hpp"
#include "AbstractTetrahedralMesh.hpp"
#include "FibreReader.hpp"
#include "Exception.hpp"
#include "FileFinder.hpp"

/**
 * Base class for different representations of conductivity tensors.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class AbstractConductivityTensors
{
protected:
    AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>* mpMesh; /**< Mesh on which to apply*/
    bool mUseNonConstantConductivities; /**< Whether conductivities can be non-constant*/
    bool mUseFibreOrientation; /**< Set by SetFibreOrientationFile so that fibre orientation can be read*/

    /** Single constant conductivities for all space (when mUseNonConstantConductivities==false)*/
    c_vector<double, SPACE_DIM> mConstantConductivities; // mS/cm

    /**
     * Non-constant conductivities for each element (when mUseNonConstantConductivities==true)
     * The size of this vector should match the number of local elements in the mesh
     */
    std::vector<c_vector<double, SPACE_DIM> >* mpNonConstantConductivities; // mS/cm

    /** Container for conductivity tensors
     * (single, size=1 [one for all space] or multiple, size=num local elements [one for each local element]) */
    std::vector< c_matrix<double,SPACE_DIM,SPACE_DIM> > mTensors;

    /** Set by Init() in the base classes*/
    bool mInitialised;

    /** Fibre orientation file (see SetFibreOrientationFile)*/
    FileFinder mFibreOrientationFile;

    /** Fibre file reader */
    std::shared_ptr<FibreReader<SPACE_DIM> > mFileReader;

public:

    AbstractConductivityTensors();

    virtual ~AbstractConductivityTensors();

    /**
     *  Sets a file for reading the fibre orientation from.
     *
     *  @param rFibreOrientationFile  the file defining the fibre orientation
     */
    void SetFibreOrientationFile(const FileFinder &rFibreOrientationFile);

    /**
     *  Sets constant conductivities for all the elements of the mesh.
     *  @param constantConductivities Longitudinal, Transverse (y axis) and Normal conductivity (z axis)
     *
     *  We need explicit instantiation of this method to make sure that c_vector length matches SPACE_DIM.
     *  Compiler won't detect mismatches.
     */
    void SetConstantConductivities(c_vector<double, 1> constantConductivities);

    /**
     *  Sets constant conductivities for all the elements of the mesh.
     *  @param constantConductivities Longitudinal, Transverse (y axis) and Normal conductivity (z axis)
     *
     *  We need explicit instantiation of this method to make sure that c_vector length matches SPACE_DIM.
     *  Compiler won't detect mismatches.
     */
    void SetConstantConductivities(c_vector<double, 2> constantConductivities);

    /**
     *  Sets constant conductivities for all the elements of the mesh.
     *  @param constantConductivities Longitudinal, Transverse (y axis) and Normal conductivity (z axis)
     *
     *  We need explicit instantiation of this method to make sure that c_vector length matches SPACE_DIM.
     *  Compiler won't detect mismatches.
     */
    virtual void SetConstantConductivities(c_vector<double, 3> constantConductivities);


    /**
     *
     *  @param pNonConstantConductivities pointer to vector of conductivities (one per local element)
     */
    void SetNonConstantConductivities(std::vector<c_vector<double, SPACE_DIM> >* pNonConstantConductivities);

    /**
     *  Computes the tensors based in all the info set
     * @param pMesh a pointer to the mesh on which these tensors are to be used
     */
    virtual void Init(AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM> *pMesh)  = 0;

    /**
     *  @return the diffussion tensor of the element number "index"
     *
     *  @param global_index Global index of the element of the mesh
     */
    c_matrix<double,SPACE_DIM,SPACE_DIM>& operator[](const unsigned global_index);
};

#endif /*ABSTRACTCONDUCTIVITYTENSORS_HPP_*/
