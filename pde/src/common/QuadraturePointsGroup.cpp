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

#include "QuadraturePointsGroup.hpp"

template<unsigned DIM>
QuadraturePointsGroup<DIM>::QuadraturePointsGroup(TetrahedralMesh<DIM,DIM>& rMesh,
                          GaussianQuadratureRule<DIM>& rQuadRule)
{
    mNumElements = rMesh.GetNumElements();
    mNumQuadPointsPerElement = rQuadRule.GetNumQuadPoints();
    data.resize(mNumElements*mNumQuadPointsPerElement, zero_vector<double>(DIM));

    // Loop over elements
    for (unsigned elem_index=0; elem_index<rMesh.GetNumElements(); elem_index++)
    {
        Element<DIM,DIM>& r_elem = *(rMesh.GetElement(elem_index));

        c_vector<double, DIM+1> linear_phi;
        for (unsigned quad_index=0; quad_index<rQuadRule.GetNumQuadPoints(); quad_index++)
        {
            const ChastePoint<DIM>& quadrature_point = rQuadRule.rGetQuadPoint(quad_index);

            LinearBasisFunction<DIM>::ComputeBasisFunctions(quadrature_point, linear_phi);

            // Interpolate to calculate quadrature point
            c_vector<double,DIM> X = zero_vector<double>(DIM);
            for (unsigned node_index=0; node_index<DIM+1; node_index++)
            {
                X += linear_phi(node_index)*rMesh.GetNode( r_elem.GetNodeGlobalIndex(node_index) )->rGetLocation();
            }

            // Save the quadrature point
            assert(elem_index<mNumElements);
            assert(quad_index<mNumQuadPointsPerElement);
            data[ elem_index*mNumQuadPointsPerElement + quad_index ] = X;
        }
    }
}

template<unsigned DIM>
c_vector<double,DIM>& QuadraturePointsGroup<DIM>::Get(unsigned elementIndex, unsigned quadIndex)
{
    assert(elementIndex<mNumElements);
    assert(quadIndex<mNumQuadPointsPerElement);
    return data[ elementIndex*mNumQuadPointsPerElement + quadIndex ];
}

template<unsigned DIM>
c_vector<double,DIM>& QuadraturePointsGroup<DIM>::Get(unsigned i)
{
    assert(i < mNumElements*mNumQuadPointsPerElement);
    return data[i];
}

template<unsigned DIM>
unsigned QuadraturePointsGroup<DIM>::GetNumElements() const
{
    return mNumElements;
}

template<unsigned DIM>
unsigned QuadraturePointsGroup<DIM>::GetNumQuadPointsPerElement() const
{
    return mNumQuadPointsPerElement;
}

template<unsigned DIM>
unsigned QuadraturePointsGroup<DIM>::Size() const
{
    return mNumElements*mNumQuadPointsPerElement;
}

////////////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
////////////////////////////////////////////////////////////////////////////////////

template class QuadraturePointsGroup<1>;
template class QuadraturePointsGroup<2>;
template class QuadraturePointsGroup<3>;
