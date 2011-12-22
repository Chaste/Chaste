/*

Copyright (C) University of Oxford, 2005-2011

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


#ifndef _ABSTRACTTETRAHEDRALELEMENT_HPP_
#define _ABSTRACTTETRAHEDRALELEMENT_HPP_

#include <vector>

#include "UblasMatrixInclude.hpp"
#include "UblasVectorInclude.hpp"

#include "AbstractElement.hpp"
#include "Node.hpp"

/**
 * This abstract class defines a tetrahedral element for use in the Finite Element Method.
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class AbstractTetrahedralElement : public AbstractElement<ELEMENT_DIM,SPACE_DIM>
{
protected:

    /**
     * Refresh the Jacobian for this element.
     *
     * @param rJacobian  the Jacobian matrix
     */
    void RefreshJacobian(c_matrix<double, SPACE_DIM, ELEMENT_DIM>& rJacobian);

public:

    /**
     * Constructor which takes in a vector of nodes.
     *
     * @param index  the index of the element in the mesh
     * @param rNodes  the nodes owned by the element
     */
    AbstractTetrahedralElement(unsigned index, const std::vector<Node<SPACE_DIM>*>& rNodes);

    /**
     * Default constructor, which doesn't fill in any nodes.
     * The nodes must be added later.
     *
     * @param index  the index of the element in the mesh (defaults to INDEX_IS_NOT_USED)
     */
    AbstractTetrahedralElement(unsigned index=INDEX_IS_NOT_USED);

    /**
     * Virtual destructor, since this class has virtual methods.
     */
    virtual ~AbstractTetrahedralElement()
    {}

    /**
     * Get the location of the centroid of the element.
     */
    c_vector<double, SPACE_DIM> CalculateCentroid() const;

    /**
     * Compute the Jacobian for this element.
     *
     * @param rJacobian  the Jacobian matrix
     * @param rJacobianDeterminant  the determinant of the Jacobian
     */
    void CalculateJacobian(c_matrix<double, SPACE_DIM, ELEMENT_DIM>& rJacobian,
                           double& rJacobianDeterminant);

    /**
     * Compute the weighted direction for this element.
     *
     * @param rWeightedDirection  the weighted direction vector
     * @param rJacobianDeterminant  the determinant of the Jacobian
     */
    void CalculateWeightedDirection(c_vector<double, SPACE_DIM>& rWeightedDirection, double& rJacobianDeterminant);

    /**
     * Compute a unit vector normal to this element, if possible.
     */
    c_vector<double, SPACE_DIM> CalculateNormal();

    /**
     * Compute the Jacobian and the inverse Jacobian for this element.
     *
     * @param rJacobian  the Jacobian matrix (output)
     * @param rJacobianDeterminant  the determinant of the Jacobian (output)
     * @param rInverseJacobian  the inverse Jacobian matrix (output)
     */
    void CalculateInverseJacobian(c_matrix<double, SPACE_DIM, ELEMENT_DIM>& rJacobian,
                                  double& rJacobianDeterminant,
                                  c_matrix<double, ELEMENT_DIM, SPACE_DIM>& rInverseJacobian);

    /** Get the volume of an element (or area in 2d, or length in 1d)
     * @param determinant a pre-calculated Jacobian determinant for this element.
     * @return volume (which is simply the determinant weighted by the SPACE_DIM)
     */
    double GetVolume(double determinant) const;

    /**
     * Place in the pIndices array, the global indices (within the stiffness matrix)
     * of the degrees of freedom associated with this element.
     *
     * @param problemDim the problem dimension e.g. 2 for Bidomain.
     * @param pIndices where to store results: an unsigned array with ELEMENT_DIM+1 entries.
     *
     */
    void GetStiffnessMatrixGlobalIndices(unsigned problemDim, unsigned* pIndices) const;
};

//////////////////////////////////////////////////////////////////////
//                  Specialization for 0d elements                  //
//////////////////////////////////////////////////////////////////////

/**
 * Specialization for 0d elements so we don't get errors from Boost on some
 * compilers.
 */
template<unsigned SPACE_DIM>
class AbstractTetrahedralElement<0, SPACE_DIM> : public AbstractElement<0,SPACE_DIM>
{
public:

    /**
     * Constructor which takes in a vector of nodes.
     *
     * @param index  the index of the element in the mesh
     * @param rNodes  the nodes owned by the element
     */
    AbstractTetrahedralElement(unsigned index, const std::vector<Node<SPACE_DIM>*>& rNodes);

    /**
     * Default constructor, which doesn't fill in any nodes.
     * The nodes must be added later.
     *
     * @param index  the index of the element in the mesh (defaults to INDEX_IS_NOT_USED)
     */
    AbstractTetrahedralElement(unsigned index=INDEX_IS_NOT_USED);

    /**
     * Virtual destructor, since this class has virtual methods.
     */
    virtual ~AbstractTetrahedralElement()
    {}

    /**
     * Get the location of the centroid of the element.
     */
    c_vector<double, SPACE_DIM> CalculateCentroid() const;

    /**
     * Compute the weighted direction for this element.
     *
     * @param rWeightedDirection  the weighted direction vector
     * @param rJacobianDeterminant  the determinant of the Jacobian
     */
    void CalculateWeightedDirection(c_vector<double, SPACE_DIM>& rWeightedDirection, double& rJacobianDeterminant);

    /**
     * Compute a unit vector normal to this element, if possible.
     */
    c_vector<double, SPACE_DIM> CalculateNormal();

    /**
     * Place in the pIndices array, the global indices (within the stiffness matrix)
     * of the degrees of freedom associated with this element.
     *
     * @param problemDim the problem dimension e.g. 2 for Bidomain.
     * @param pIndices where to store results: an unsigned array with ELEMENT_DIM+1 entries.
     *
     */
    void GetStiffnessMatrixGlobalIndices(unsigned problemDim, unsigned* pIndices) const;
};

#endif //_ABSTRACTTETRAHEDRALELEMENT_HPP_
