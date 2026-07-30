
/*

Copyright (c) 2005-2026, University of Oxford.
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

#include "SemSingleElementMeshGenerator.hpp"

#include "Exception.hpp"
#include "SemEnumerations.hpp"
#include "SemLatticeGeometry.hpp"
#include "SemMesh.hpp"


template <unsigned DIM> SemSingleElementMeshGenerator<DIM>::SemSingleElementMeshGenerator(
    const std::array<unsigned, DIM>& numNodes,
    double scaleFactor,
    SemLatticeType nodeLattice)
    : mpMesh{std::make_shared<SemMesh<DIM>>()},
      mNumNodes{ numNodes },
      mScaleFactor{ scaleFactor },
      mNodeLattice{ nodeLattice }
{
    for (unsigned i = 0; i < DIM; ++i)
    {
        if (mNumNodes[i] == 0u)
        {
            EXCEPTION("SemSingleElementMeshGenerator: each entry of numNodes must be >= 1");
        }
    }

    if (scaleFactor <= 0.0)
    {
        EXCEPTION("SemSingleElementMeshGenerator: scaleFactor must be positive");
    }

    mNodeSpacing = mScaleFactor / static_cast<double>(mNumNodes[0]);

    std::vector<c_vector<double, DIM>> positions = GenerateNodePositions();
    GenerateMesh(positions);
}

template <unsigned DIM>
std::shared_ptr<SemMesh<DIM>> SemSingleElementMeshGenerator<DIM>::GetMesh()
{
    return mpMesh;
}

template <unsigned DIM>
std::vector<c_vector<double, DIM>> SemSingleElementMeshGenerator<DIM>::GenerateNodePositions() const
{
    std::array<double, DIM> spacing;
    spacing.fill(mNodeSpacing);

    return GenerateSemLatticePositions<DIM>(mNumNodes, spacing, mNodeLattice);
}

template <unsigned DIM>
bool SemSingleElementMeshGenerator<DIM>::IsBoundaryNode(unsigned flatIndex) const
{
    if constexpr (DIM == 1)
    {
        return (flatIndex == 0u) || (flatIndex == mNumNodes[0] - 1u);
    }
    else if constexpr (DIM == 2)
    {
        const unsigned i = flatIndex % mNumNodes[0];
        const unsigned j = flatIndex / mNumNodes[0];
        return (i == 0u) || (i == mNumNodes[0] - 1u) || (j == 0u) || (j == mNumNodes[1] - 1u);
    }
    else
    {
        const unsigned i = flatIndex % mNumNodes[0];
        const unsigned j = (flatIndex / mNumNodes[0]) % mNumNodes[1];
        const unsigned k = flatIndex / (mNumNodes[0] * mNumNodes[1]);
        return (i == 0u) || (i == mNumNodes[0] - 1u) ||
               (j == 0u) || (j == mNumNodes[1] - 1u) ||
               (k == 0u) || (k == mNumNodes[2] - 1u);
    }
}

template <unsigned DIM> void SemSingleElementMeshGenerator<DIM>::GenerateMesh(std::vector<c_vector<double, DIM>> positions)
{
    unsigned int new_element_id = mpMesh->GetNumElements();
    auto new_element = new SemElement<DIM>(new_element_id, {});

    for (unsigned i = 0; i < positions.size(); ++i)
    {
        unsigned new_node_index = mpMesh->GetNumNodes();
        Node<DIM>* new_node = new Node<DIM>(new_node_index, positions[i]);
        new_node->SetRegion(IsBoundaryNode(i) ? SEM_BOUNDARY_REGION : SEM_INTERIOR_REGION);
        new_node->SetRadius(0.05);
        new_node->AddElement(new_element_id);

        mpMesh->AddNode(new_node);
        new_element->AddNode(new_node);
    }

    mpMesh->AddElement(new_element);
    NodeMap map(mpMesh->GetNumNodes());
    mpMesh->ReMesh(map);
}

template class SemSingleElementMeshGenerator<1>;
template class SemSingleElementMeshGenerator<2>;
template class SemSingleElementMeshGenerator<3>;
