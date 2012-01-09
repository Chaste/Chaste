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
