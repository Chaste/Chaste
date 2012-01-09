/*

Copyright (C) University of Oxford, 2005-2012

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

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
    std::auto_ptr<FibreReader<SPACE_DIM> > mFileReader;

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
    virtual void Init(AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM> *pMesh) throw (Exception) = 0;

    /**
     *  Returns the diffussion tensor of the element number "index"
     *
     *  @param global_index Global index of the element of the mesh
     */
    c_matrix<double,SPACE_DIM,SPACE_DIM>& operator[](const unsigned global_index);
};



#endif /*ABSTRACTCONDUCTIVITYTENSORS_HPP_*/
