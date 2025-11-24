
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

#include "SemSingleElementMeshGenerator.hpp"

#include "Exception.hpp"
#include "SemMesh.hpp"


template <unsigned DIM> SemSingleElementMeshGenerator<DIM>::SemSingleElementMeshGenerator(
    const std::array<unsigned, DIM>& numNodes,
    double scaleFactor)
    : mpMesh{std::make_shared<SemMesh<DIM>>()},
      mNumNodes{ numNodes },
      mScaleFactor{ scaleFactor }
{
    unsigned num_all_nodes = 1;
    for (unsigned i = 0; i < DIM; ++i)
    {
        if (mNumNodes[i] == 0u)
        {
            EXCEPTION("SemSingleElementMeshGenerator: each entry of numNodes must be >= 1");
        }
        num_all_nodes *= mNumNodes[i];
    }

    mNumAllNodes = num_all_nodes;

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
    std::vector<c_vector<double, DIM>> positions;
    positions.reserve(mNumAllNodes);

    if constexpr (DIM == 1)
    {
        for (unsigned i = 0; i < mNumNodes[0]; ++i)
        {
            c_vector<double,1> v;
            v[0] = static_cast<double>(i) * mNodeSpacing;
            positions.push_back(v);
        }
    }
    else if constexpr (DIM == 2)
    {
        for (unsigned j = 0; j < mNumNodes[1]; ++j)
        {
            for (unsigned i = 0; i < mNumNodes[0]; ++i)
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
        for (unsigned k = 0; k < mNumNodes[2]; ++k)
        {
            for (unsigned j = 0; j < mNumNodes[1]; ++j)
            {
                for (unsigned i = 0; i < mNumNodes[0]; ++i)
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

template <unsigned DIM> void SemSingleElementMeshGenerator<DIM>::GenerateMesh(std::vector<c_vector<double, DIM>> positions)
{
    unsigned int new_element_id = mpMesh->GetNumElements();
    auto new_element = new SemElement<DIM>(new_element_id, {});
    new_element->SetCellId(new_element_id);

    std::vector<unsigned> all_node_indices;
    all_node_indices.reserve(mNumAllNodes);

    for (unsigned i = 0; i < positions.size(); ++i)
    {
        unsigned new_node_index = mpMesh->GetNumNodes();
        Node<DIM>* new_node = new Node<DIM>(new_node_index, positions[i]);
        new_node->SetRegion(0u);
        new_node->SetRadius(0.05);
        new_node->AddElement(new_element_id);

        // Add the node to the mesh
        mpMesh->AddNode(new_node);
        new_element->AddNode(new_node);

        // Add the node indices to the relevant layer groups
        all_node_indices.push_back(new_node_index);
    }

    // Add the layer groups to the element
    new_element->AddInteractionLayer("all-nodes", all_node_indices);

    // Add the element to the mesh
    mpMesh->AddElement(new_element);
    NodeMap map(mpMesh->GetNumNodes());
    mpMesh->ReMesh(map);
}

template class SemSingleElementMeshGenerator<1>;
template class SemSingleElementMeshGenerator<2>;
template class SemSingleElementMeshGenerator<3>;
