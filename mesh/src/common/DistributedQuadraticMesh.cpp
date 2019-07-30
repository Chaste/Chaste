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


#include "DistributedQuadraticMesh.hpp"

template<unsigned DIM>
DistributedQuadraticMesh<DIM>::DistributedQuadraticMesh(DistributedTetrahedralMeshPartitionType::type partitioningMethod)
    : DistributedTetrahedralMesh<DIM, DIM>(partitioningMethod)
{

}

template<unsigned DIM>
DistributedQuadraticMesh<DIM>::~DistributedQuadraticMesh()
{
}

template<unsigned DIM>
void DistributedQuadraticMesh<DIM>::ConstructFromMeshReader(AbstractMeshReader<DIM, DIM>& rAbsMeshReader)
{
    TrianglesMeshReader<DIM, DIM>* p_mesh_reader=dynamic_cast<TrianglesMeshReader<DIM, DIM>*>(&rAbsMeshReader);
    assert(p_mesh_reader != nullptr);

    unsigned order_of_elements = 1;
    if (p_mesh_reader)
    {
        //A triangles mesh reader will let you read with non-linear elements
        order_of_elements = p_mesh_reader->GetOrderOfElements();
    }

    // If it is a linear TrianglesMeshReader or any other reader (which are all linear)
    if (order_of_elements == 1)
    {
        EXCEPTION("Cannot convert a (linear) tetrahedral mesh directly to a DistributedQuadraticMesh.  Please convert to QuadraticMesh and save in that format first.");
    }
    this->mMeshIsLinear = false;
    DistributedTetrahedralMesh<DIM,DIM>::ConstructFromMeshReader(*p_mesh_reader);
    assert(this->GetNumBoundaryElements() > 0u);
    QuadraticMeshHelper<DIM>::AddInternalNodesToElements(this, p_mesh_reader);
    QuadraticMeshHelper<DIM>::AddInternalNodesToBoundaryElements(this, p_mesh_reader);
    QuadraticMeshHelper<DIM>::CheckBoundaryElements(this);
}

// Explicit instantiation
template class DistributedQuadraticMesh<1>;
template class DistributedQuadraticMesh<2>;
template class DistributedQuadraticMesh<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(DistributedQuadraticMesh)
