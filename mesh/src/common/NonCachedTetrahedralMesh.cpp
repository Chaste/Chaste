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

#include "NonCachedTetrahedralMesh.hpp"
#include "Exception.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NonCachedTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::RefreshJacobianCachedData()
{
    // Don't do any caching
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NonCachedTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetJacobianForElement(unsigned elementIndex, c_matrix<double, SPACE_DIM, SPACE_DIM>& rJacobian, double& rJacobianDeterminant) const
{
    EXCEPTION("Use GetInverseJacobianForElement to retrieve Jacobian data instead.");
}


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NonCachedTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetInverseJacobianForElement(
        unsigned elementIndex,
        c_matrix<double, SPACE_DIM, ELEMENT_DIM>& rJacobian,
        double& rJacobianDeterminant,
        c_matrix<double, ELEMENT_DIM, SPACE_DIM>& rInverseJacobian) const
{
    this->mElements[this->SolveElementMapping(elementIndex)]->CalculateInverseJacobian(rJacobian, rJacobianDeterminant, rInverseJacobian);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NonCachedTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetWeightedDirectionForElement(unsigned elementIndex, c_vector<double, SPACE_DIM>& rWeightedDirection, double& rJacobianDeterminant) const
{
    // See comment in AbstractTetrahedralMesh::GetWeightedDirectionForBoundaryElement()
    EXCEPTION("Probably redundant method.");
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void NonCachedTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetWeightedDirectionForBoundaryElement(
        unsigned elementIndex,
        c_vector<double, SPACE_DIM>& rWeightedDirection,
        double& rJacobianDeterminant) const
{
    this->mBoundaryElements[this->SolveBoundaryElementMapping(elementIndex)]->CalculateWeightedDirection(rWeightedDirection, rJacobianDeterminant );
}


/////////////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////////////

template class NonCachedTetrahedralMesh<3,3>;
template class NonCachedTetrahedralMesh<2,2>;
template class NonCachedTetrahedralMesh<1,1>;


// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS2(NonCachedTetrahedralMesh, 1, 1)
EXPORT_TEMPLATE_CLASS2(NonCachedTetrahedralMesh, 2, 2)
EXPORT_TEMPLATE_CLASS2(NonCachedTetrahedralMesh, 3, 3)
