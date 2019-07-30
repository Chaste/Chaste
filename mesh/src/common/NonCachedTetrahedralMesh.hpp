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

#ifndef NONCACHEDTETRAHEDRALMESH_HPP_
#define NONCACHEDTETRAHEDRALMESH_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "TetrahedralMesh.hpp"

/**
 * A drop-in replacement for TetrahedralMesh that doesn't cache any
 * jacobian-related data.
 *
 * It thus provides essentially a serial version of the memory-efficient
 * DistributedTetrahedralMesh, enabling the use of larger meshes on
 * single-processor machines.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class NonCachedTetrahedralMesh : public TetrahedralMesh< ELEMENT_DIM, SPACE_DIM>
{
private:
    /** Needed for serialization.*/
    friend class boost::serialization::access;
    /**
     * Serialize the mesh.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
       archive & boost::serialization::base_object<TetrahedralMesh<ELEMENT_DIM, SPACE_DIM> >(*this);
    }
public:

    /** Reimplemented to do no caching */
    void RefreshJacobianCachedData();

    /**
     * Get the Jacobian matrix and its determinant for a given element.
     *
     * @param elementIndex index of the element in the mesh
     * @param rJacobian the Jacobian matrix
     * @param rJacobianDeterminant the determinant of the Jacobian matrix
     */
    void GetJacobianForElement(unsigned elementIndex, c_matrix<double, SPACE_DIM, SPACE_DIM>& rJacobian, double& rJacobianDeterminant) const;

    /**
     * Get the Jacobian matrix, its inverse and its determinant for a given element.
     *
     * @param elementIndex index of the element in the mesh
     * @param rJacobian the Jacobian matrix
     * @param rJacobianDeterminant the determinant of the Jacobian matrix
     * @param rInverseJacobian the inverse Jacobian matrix
     */
    void GetInverseJacobianForElement(unsigned elementIndex, c_matrix<double, SPACE_DIM, ELEMENT_DIM>& rJacobian, double& rJacobianDeterminant, c_matrix<double, ELEMENT_DIM, SPACE_DIM>& rInverseJacobian) const;

    /**
     * Get the weighted direction and the determinant of the Jacobian for a given element.
     *
     * @param elementIndex index of the element in the mesh
     * @param rWeightedDirection the weighted direction
     * @param rJacobianDeterminant the determinant of the Jacobian matrix
     */
    void GetWeightedDirectionForElement(unsigned elementIndex, c_vector<double, SPACE_DIM>& rWeightedDirection, double& rJacobianDeterminant) const;

    /**
     * Get the weighted direction and the determinant of the Jacobian for a given boundary element.
     *
     * @param elementIndex index of the element in the mesh
     * @param rWeightedDirection the weighted direction
     * @param rJacobianDeterminant the determinant of the Jacobian matrix
     */
    void GetWeightedDirectionForBoundaryElement(unsigned elementIndex, c_vector<double, SPACE_DIM>& rWeightedDirection, double& rJacobianDeterminant) const;
};

EXPORT_TEMPLATE_CLASS2(NonCachedTetrahedralMesh, 1, 1)
EXPORT_TEMPLATE_CLASS2(NonCachedTetrahedralMesh, 2, 2)
EXPORT_TEMPLATE_CLASS2(NonCachedTetrahedralMesh, 3, 3)

#endif /*NONCACHEDTETRAHEDRALMESH_HPP_*/
