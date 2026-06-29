
/*

Copyright (c) 2005-2025, University of Oxford.
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
#include "SemMesh.hpp"


template <unsigned DIM> SemMultiElementMeshGenerator<DIM>::SemMultiElementMeshGenerator(
    const std::array<unsigned, DIM>& numNodesPerElem, const std::array<unsigned, DIM>& numElems, double scaleFactor)
    : mpMesh{std::make_shared<SemMesh<DIM>>()},
      mNumNodesPerElem{ numNodesPerElem },
      mNumElems{ numElems },
      mScaleFactor{ scaleFactor }
{
    unsigned num_all_nodes_per_elem = 1;
    unsigned num_all_elems = 1;
    for (unsigned i = 0; i < DIM; ++i)
    {
        if (mNumNodesPerElem[i] == 0u)
        {
            EXCEPTION("SemMultiElementMeshGenerator: each entry of numNodesPerElem must be >= 1");
        }
        num_all_nodes_per_elem *= mNumNodesPerElem[i];

        if (mNumElems[i] == 0u)
        {
            EXCEPTION("SemMultiElementMeshGenerator: each entry of numElems must be >= 1");
        }
        num_all_elems *= mNumElems[i];
    }

    mNumAllNodesPerElem = num_all_nodes_per_elem;
    mNumAllElems = num_all_elems;

    if (scaleFactor <= 0.0)
    {
        EXCEPTION("SemMultiElementMeshGenerator: scaleFactor must be positive");
    }

    mNodeSpacing = mScaleFactor / static_cast<double>(mNumNodesPerElem[0]);

    for (unsigned i = 0; i < DIM; ++i)
    {
        mElemSpacing[i] = mNodeSpacing * static_cast<double>(mNumNodesPerElem[i]);
    }

    std::vector<c_vector<double, DIM>> positions = GenerateNodePositions();
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
    std::vector<c_vector<double, DIM>> positions;
    positions.reserve(mNumAllNodesPerElem);

    if constexpr (DIM == 1)
    {
        for (unsigned i = 0; i < mNumNodesPerElem[0]; ++i)
        {
            c_vector<double,1> v;
            v[0] = static_cast<double>(i) * mNodeSpacing;
            positions.push_back(v);
        }
    }
    else if constexpr (DIM == 2)
    {
        for (unsigned j = 0; j < mNumNodesPerElem[1]; ++j)
        {
            for (unsigned i = 0; i < mNumNodesPerElem[0]; ++i)
            {
                c_vector<double,2> v;
                v[0] = static_cast<double>(i) * mNodeSpacing;
                v[1] = static_cast<double>(j) * mNodeSpacing;
                positions.push_back(v);
            }
        }
    }
    else if constexpr (DIM == 3)
    {
        for (unsigned k = 0; k < mNumNodesPerElem[2]; ++k)
        {
            for (unsigned j = 0; j < mNumNodesPerElem[1]; ++j)
            {
                for (unsigned i = 0; i < mNumNodesPerElem[0]; ++i)
                {
                    c_vector<double,3> v;
                    v[0] = static_cast<double>(i) * mNodeSpacing;
                    v[1] = static_cast<double>(j) * mNodeSpacing;
                    v[2] = static_cast<double>(k) * mNodeSpacing;
                    positions.push_back(v);
                }
            }
        }
    }
    else
    {
        NEVER_REACHED;
    }

    return positions;
}

template <unsigned DIM>
std::vector<c_vector<double, DIM>> SemMultiElementMeshGenerator<DIM>::GenerateElementOffsets() const
{
    std::vector<c_vector<double, DIM>> offsets;
    offsets.reserve(mNumAllElems);

    if constexpr (DIM == 1)
    {
        for (unsigned i = 0; i < mNumElems[0]; ++i)
        {
            c_vector<double,1> v;
            v[0] = static_cast<double>(i) * mElemSpacing[i];
            offsets.push_back(v);
        }
    }
    else if constexpr (DIM == 2)
    {
        for (unsigned j = 0; j < mNumElems[1]; ++j)
        {
            for (unsigned i = 0; i < mNumElems[0]; ++i)
            {
                c_vector<double,2> v;
                v[0] = static_cast<double>(i) * mElemSpacing[0];
                v[1] = static_cast<double>(j) * mElemSpacing[1];
                offsets.push_back(v);
            }
        }
    }
    else if constexpr (DIM == 3)
    {
        for (unsigned k = 0; k < mNumElems[2]; ++k)
        {
            for (unsigned j = 0; j < mNumElems[1]; ++j)
            {
                for (unsigned i = 0; i < mNumElems[0]; ++i)
                {
                    c_vector<double,3> v;
                    v[0] = static_cast<double>(i) * mElemSpacing[0];
                    v[1] = static_cast<double>(j) * mElemSpacing[1];
                    v[2] = static_cast<double>(k) * mElemSpacing[2];
                    offsets.push_back(v);
                }
            }
        }
    }
    else
    {
        NEVER_REACHED;
    }

    return offsets;
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
        new_element->SetCellId(new_element_id);

        for (unsigned i = 0; i < positions.size(); ++i)
        {
            unsigned new_node_index = mpMesh->GetNumNodes();
            Node<DIM>* new_node = new Node<DIM>(new_node_index, positions[i] + elem_offset);
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
}

template class SemMultiElementMeshGenerator<1>;
template class SemMultiElementMeshGenerator<2>;
template class SemMultiElementMeshGenerator<3>;
