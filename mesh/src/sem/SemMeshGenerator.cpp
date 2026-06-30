
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

#include "SemMeshGenerator.hpp"

#include "SemEnumerations.hpp"
#include <boost/foreach.hpp>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
#include "TrianglesMeshReader.hpp"
#include "OutputFileHandler.hpp"
#include "RandomNumberGenerator.hpp"
#include "MathsCustomFunctions.hpp"
#include "ChasteSyscalls.hpp"

SemMeshGenerator::SemMeshGenerator(unsigned numNodesAlongWidth, unsigned numNodesAlongLength, unsigned ghosts, double scaleFactor)
  : mpMesh(std::make_shared<SemMesh<2>>()),
    mMeshFilename("mesh"),
    mDomainWidth(numNodesAlongWidth * scaleFactor),
    mNumCellWidth(numNodesAlongWidth),
    mNumCellLength(numNodesAlongLength)
{
}

SemMeshGenerator::SemMeshGenerator()
    : mpMesh(std::make_shared<SemMesh<2>>())
{
}

std::shared_ptr<SemMesh<2> > SemMeshGenerator::GetMesh()
{
    return mpMesh;
}

std::vector<unsigned> SemMeshGenerator::GetCellLocationIndices()
{
    std::vector<unsigned> location_indices;

    for (unsigned i=0; i<mpMesh->GetNumNodes(); i++)
    {
        if (mGhostNodeIndices.find(i)==mGhostNodeIndices.end())
        {
            location_indices.push_back(i);
        }
    }
    return location_indices;
}

std::shared_ptr<SemMesh<2> > SemMeshGenerator::GetCircularMesh(double radius)
{
    if (!mGhostNodeIndices.empty())
    {
        EXCEPTION("Cannot call GetCircularMesh on a SemMesh with ghost nodes");
    }

    // Centre the mesh at (0,0)
    c_vector<double,2> centre = zero_vector<double>(2);
    for (unsigned i=0; i<mpMesh->GetNumNodes(); i++)
    {
        centre += mpMesh->GetNode(i)->rGetLocation();
    }
    centre /= (double)mpMesh->GetNumNodes();

    mpMesh->Translate(-centre[0], -centre[1]);

    // Iterate over nodes, deleting any that lie more than the specified radius from (0,0)
    for (unsigned i=0; i<mpMesh->GetNumAllNodes(); i++)
    {
        if (norm_2(mpMesh->GetNode(i)->rGetLocation()) >= radius)
        {
            mpMesh->DeleteNodePriorToReMesh(i);
        }
        else
        {
            // Jiggle the data
            c_vector<double,2>& r_location = mpMesh->GetNode(i)->rGetModifiableLocation();
            c_vector<double,2> shift;
            RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
            double max_jiggle = radius*5e-6;
            shift[0] = max_jiggle*(p_gen->ranf()-0.5);
            shift[1] = max_jiggle*(p_gen->ranf()-0.5);
            r_location += shift;
        }
    }

    // Remesh
    NodeMap map(mpMesh->GetNumNodes());
    mpMesh->ReMesh(map);

    return mpMesh;
}

double SemMeshGenerator::GetDomainDepth()
{
    return mDomainDepth;
}

double SemMeshGenerator::GetDomainWidth()
{
    return mDomainWidth;
}

void SemMeshGenerator::GenerateSingleCell(std::array<double, 2> center, std::array<double, 2> dimensions, std::array<double, 2> num_nodes_on_axis)
{
    // Single cell is represented by a single SemElement
    unsigned int new_element_id = mpMesh->GetNumElements();
    auto new_element = new SemElement<2>(new_element_id, {});
    new_element->SetCellId(new_element_id);

    // Generate the nodes
    for (unsigned x = 0; x < num_nodes_on_axis[0]; x++) {
        for (unsigned y = 0; y < num_nodes_on_axis[1]; y++) {
            double x_location_normed = double(x) / double(num_nodes_on_axis[0]);
            double y_location_normed = double(y) / double(num_nodes_on_axis[1]);

            // Compute the node location
            double x_location = center[0] + (x_location_normed - 0.5) * dimensions[0];
            double y_location = center[1] + (y_location_normed - 0.5) * dimensions[1];

            const unsigned max_x_index = num_nodes_on_axis[0] - 1;
            const unsigned max_y_index = num_nodes_on_axis[1] - 1;
            const bool is_boundary_node = x == 0 || y == 0 || x == max_x_index || y == max_y_index;

            unsigned int new_node_index = mpMesh->GetNumNodes();
            std::vector<double> node_location = {x_location, y_location};
            Node<2>* new_node = new Node<2>(new_node_index, node_location, is_boundary_node);

            new_node->SetRegion(is_boundary_node ? SEM_BOUNDARY_REGION : SEM_INTERIOR_REGION);
            new_node->SetRadius(0.05);

            mpMesh->AddNode(new_node);
            new_element->AddNode(new_node);
        }
    }

    new_element->RegisterWithNodes();

    // Add the element to the mesh
    mpMesh->AddElement(new_element);
    NodeMap map(mpMesh->GetNumNodes());
    mpMesh->ReMesh(map);
}
