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

#include "QuadraturePointsGroup.hpp"
#include "Exception.hpp"

template<unsigned DIM>
QuadraturePointsGroup<DIM>::QuadraturePointsGroup(AbstractTetrahedralMesh<DIM,DIM>& rMesh,
                                                  GaussianQuadratureRule<DIM>& rQuadRule)
{
    c_vector<double, DIM> unset;
    for (unsigned i=0; i<DIM; i++)
    {
        unset(i)=DOUBLE_UNSET;
    }
    mNumElements = rMesh.GetNumElements();
    mNumQuadPointsPerElement = rQuadRule.GetNumQuadPoints();
    data.resize(mNumElements*mNumQuadPointsPerElement, unset);

    // Loop over elements
    for (typename AbstractTetrahedralMesh<DIM, DIM>::ElementIterator iter = rMesh.GetElementIteratorBegin();
                  iter != rMesh.GetElementIteratorEnd();
                  ++iter)
    {
//        Element<DIM,DIM>& r_elem = *(rMesh.GetElement(elem_index));
        unsigned elem_index = iter->GetIndex();

        c_vector<double, DIM+1> linear_phi;
        for (unsigned quad_index=0; quad_index<rQuadRule.GetNumQuadPoints(); quad_index++)
        {
            const ChastePoint<DIM>& quadrature_point = rQuadRule.rGetQuadPoint(quad_index);

            LinearBasisFunction<DIM>::ComputeBasisFunctions(quadrature_point, linear_phi);

            // Interpolate to calculate quadrature point
            c_vector<double,DIM> physical_quad_point = zero_vector<double>(DIM);
            for (unsigned node_index=0; node_index<DIM+1; node_index++)
            {
                physical_quad_point += linear_phi(node_index)*(iter->GetNode(node_index))->rGetLocation();
            }

            // Save the quadrature point
            assert(elem_index<mNumElements);
            assert(quad_index<mNumQuadPointsPerElement);
            data[ elem_index*mNumQuadPointsPerElement + quad_index ] = physical_quad_point;
        }
    }
}

template<unsigned DIM>
c_vector<double,DIM>& QuadraturePointsGroup<DIM>::rGet(unsigned elementIndex, unsigned quadIndex)
{
    assert(elementIndex<mNumElements);
    assert(quadIndex<mNumQuadPointsPerElement);
    return data[ elementIndex*mNumQuadPointsPerElement + quadIndex ];
}

template<unsigned DIM>
c_vector<double,DIM>& QuadraturePointsGroup<DIM>::rGet(unsigned i)
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

//////////////// Explicit instantiation//////////////

template class QuadraturePointsGroup<1>;
template class QuadraturePointsGroup<2>;
template class QuadraturePointsGroup<3>;
