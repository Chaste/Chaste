
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

#include "SemMultiElementMeshGenerator.hpp"

#include "Exception.hpp"
#include "SemEnumerations.hpp"
#include "SemLatticeGeometry.hpp"
#include "SemMesh.hpp"


template <unsigned DIM> SemMultiElementMeshGenerator<DIM>::SemMultiElementMeshGenerator(
    const std::array<unsigned, DIM>& numNodesPerElem, const std::array<unsigned, DIM>& numElems, double scaleFactor,
    SemLatticeType::Value nodeLattice, SemLatticeType::Value elementLattice)
    : mpMesh{std::make_shared<SemMesh<DIM>>()},
      mNumNodesPerElem{ numNodesPerElem },
      mNumElems{ numElems },
      mScaleFactor{ scaleFactor },
      mNodeLattice{ nodeLattice },
      mElementLattice{ elementLattice }
{
    for (unsigned i = 0; i < DIM; ++i)
    {
        if (mNumNodesPerElem[i] == 0u)
        {
            EXCEPTION("SemMultiElementMeshGenerator: each entry of numNodesPerElem must be >= 1");
        }

        if (mNumElems[i] == 0u)
        {
            EXCEPTION("SemMultiElementMeshGenerator: each entry of numElems must be >= 1");
        }
    }

    if (scaleFactor <= 0.0)
    {
        EXCEPTION("SemMultiElementMeshGenerator: scaleFactor must be positive");
    }

    mNodeSpacing = mScaleFactor / static_cast<double>(mNumNodesPerElem[0]);

    std::vector<c_vector<double, DIM>> positions = GenerateNodePositions();

    // Space elements so that nodes in adjacent elements are no closer than nodes within an
    // element, whichever lattice either is placed on. On a cubic lattice each dimension can be
    // spaced by its own extent plus one node spacing, which for cubic nodes recovers the original
    // spacing mNodeSpacing * mNumNodesPerElem[i]; a close-packed lattice staggers its rows, so it
    // needs a single spacing accounting for the directions in which elements interlock.
    const c_vector<double, DIM> element_extent = GetSemLatticeExtent<DIM>(positions);
    if (mElementLattice == SemLatticeType::SEM_LATTICE_CLOSE_PACKED)
    {
        mElemSpacing.fill(GetSemClosePackedSpacing<DIM>(element_extent, mNodeSpacing));
    }
    else
    {
        for (unsigned i = 0; i < DIM; ++i)
        {
            mElemSpacing[i] = element_extent[i] + mNodeSpacing;
        }
    }

    std::vector<c_vector<double, DIM>> offsets = GenerateElementOffsets();
    GenerateMesh(positions, offsets);
}

template <unsigned DIM>
std::shared_ptr<SemMesh<DIM>> SemMultiElementMeshGenerator<DIM>::GetMesh()
{
    return mpMesh;
}

template <unsigned DIM>
std::vector<c_vector<double, DIM>> SemMultiElementMeshGenerator<DIM>::GenerateNodePositions() const
{
    std::array<double, DIM> spacing;
    spacing.fill(mNodeSpacing);

    return GenerateSemLatticePositions<DIM>(mNumNodesPerElem, spacing, mNodeLattice);
}

template <unsigned DIM>
std::vector<c_vector<double, DIM>> SemMultiElementMeshGenerator<DIM>::GenerateElementOffsets() const
{
    return GenerateSemLatticePositions<DIM>(mNumElems, mElemSpacing, mElementLattice);
}

template <unsigned DIM>
bool SemMultiElementMeshGenerator<DIM>::IsBoundaryNode(unsigned flatIndex) const
{
    if constexpr (DIM == 1)
    {
        return (flatIndex == 0u) || (flatIndex == mNumNodesPerElem[0] - 1u);
    }
    else if constexpr (DIM == 2)
    {
        const unsigned i = flatIndex % mNumNodesPerElem[0];
        const unsigned j = flatIndex / mNumNodesPerElem[0];
        return (i == 0u) || (i == mNumNodesPerElem[0] - 1u) || (j == 0u) || (j == mNumNodesPerElem[1] - 1u);
    }
    else
    {
        const unsigned i = flatIndex % mNumNodesPerElem[0];
        const unsigned j = (flatIndex / mNumNodesPerElem[0]) % mNumNodesPerElem[1];
        const unsigned k = flatIndex / (mNumNodesPerElem[0] * mNumNodesPerElem[1]);
        return (i == 0u) || (i == mNumNodesPerElem[0] - 1u) ||
               (j == 0u) || (j == mNumNodesPerElem[1] - 1u) ||
               (k == 0u) || (k == mNumNodesPerElem[2] - 1u);
    }
}

template <unsigned DIM> void SemMultiElementMeshGenerator<DIM>::GenerateMesh(
    std::vector<c_vector<double, DIM>> positions, std::vector<c_vector<double, DIM>> offsets)
{
    for (const auto& elem_offset : offsets)
    {
        unsigned int new_element_id = mpMesh->GetNumElements();
        auto new_element = new SemElement<DIM>(new_element_id, {});

        for (unsigned i = 0; i < positions.size(); ++i)
        {
            unsigned new_node_index = mpMesh->GetNumNodes();
            Node<DIM>* new_node = new Node<DIM>(new_node_index, positions[i] + elem_offset);
            new_node->SetRegion(IsBoundaryNode(i) ? SEM_BOUNDARY_REGION : SEM_INTERIOR_REGION);
            new_node->AddElement(new_element_id);

            mpMesh->AddNode(new_node);
            new_element->AddNode(new_node);
        }

        mpMesh->AddElement(new_element);
        NodeMap map(mpMesh->GetNumNodes());
        mpMesh->ReMesh(map);
    }
}

template class SemMultiElementMeshGenerator<1>;
template class SemMultiElementMeshGenerator<2>;
template class SemMultiElementMeshGenerator<3>;
