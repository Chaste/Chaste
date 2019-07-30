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

#include "CuboidMeshConstructor.hpp"
#include "MathsCustomFunctions.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CuboidMeshConstructor<ELEMENT_DIM, SPACE_DIM>::Construct(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>& rMesh, unsigned meshRefinementNum, double meshWidth)
{
    // The class has only been tested for ELEMENT_DIM == SPACE_DIM or ELEMENT_DIM == 1 && SPACE_DIM == 3
    assert(ELEMENT_DIM == SPACE_DIM || (ELEMENT_DIM == 1 && SPACE_DIM == 3));

    mMeshWidth = meshWidth;
    assert(meshRefinementNum < 30); // Sanity

    // Create the mesh
    unsigned mesh_size = SmallPow(2u, meshRefinementNum+2); // number of elements in each dimension
    double scaling = mMeshWidth/(double) mesh_size;
    rMesh.ConstructRegularSlabMesh(scaling, mMeshWidth, mMeshWidth, mMeshWidth);
    mNumElements = rMesh.GetNumElements();
    mNumNodes = rMesh.GetNumNodes();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double CuboidMeshConstructor<ELEMENT_DIM, SPACE_DIM>::GetWidth()
{
    return mMeshWidth;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned CuboidMeshConstructor<ELEMENT_DIM, SPACE_DIM>::GetNumElements()
{
    return mNumElements;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned CuboidMeshConstructor<ELEMENT_DIM, SPACE_DIM>::GetNumNodes()
{
    return mNumNodes;
}

////// Explicit instantiation////

template class CuboidMeshConstructor<1>;
template class CuboidMeshConstructor<2>;
template class CuboidMeshConstructor<3>;
template class CuboidMeshConstructor<1,3>;
