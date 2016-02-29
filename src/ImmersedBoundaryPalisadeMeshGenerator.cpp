
/*

Copyright (c) 2005-2015, University of Oxford.
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

#include "ImmersedBoundaryPalisadeMeshGenerator.hpp"
#include "RandomNumberGenerator.hpp"

#include "Debug.hpp"

ImmersedBoundaryPalisadeMeshGenerator::ImmersedBoundaryPalisadeMeshGenerator(unsigned numCellsWide, unsigned numNodesPerCell, double ellipseExponent, double cellAspectRatio, double randomYMult, bool membrane)
  : mpMesh(NULL),
    mNumCellsWide(numCellsWide),
    mNumNodesPerCell(numNodesPerCell),
    mEllipseExponent(ellipseExponent),
    mCellAspectRatio(cellAspectRatio),
    mRandomYMult(randomYMult),
    mMembrane(membrane)
{
    // Check for sensible input
    assert(mNumCellsWide > 0);
    assert(mNumNodesPerCell > 3);
    assert(mEllipseExponent > 0.0);
    assert(mCellAspectRatio > 0.0); // aspect ratio is cell height / cell width
    assert(fabs(mRandomYMult < 2.0));

    // Helper vectors
    unit_vector<double> x_unit(2,0);
    unit_vector<double> y_unit(2,1);

    // Calculate cell sizes and rescale if necessary
    double cell_width = 1.0 / double(mNumCellsWide);
    double cell_height = mCellAspectRatio * cell_width;

    if (cell_height > 0.8)
    {
        cell_height = 0.8;
        cell_width = cell_height / mCellAspectRatio;
    }

    // Generate a reference superellipse
    SuperellipseGenerator* p_gen = new SuperellipseGenerator(mNumNodesPerCell, mEllipseExponent, cell_width, cell_height, 0.0, 0.0);
    std::vector<c_vector<double, 2> > locations = p_gen->GetPointsAsVectors();

    // Calculate which locations will be the element corners
    double top_height = p_gen->GetHeightOfTopSurface();
    double bot_height = cell_height - top_height;

    if (top_height < 0.5 * cell_height || top_height > cell_height)
    {
        EXCEPTION("Something went wrong calculating the height of top surface of the cell.");
    }

    // Corner 0 is left-apical, 1 is right-apical, 2 is right-basal, 3 is left-basal
    std::vector<unsigned> corner_indices(4, UINT_MAX);

    // First anticlockwise corner is corner 1
    for(unsigned location = 0 ; location < locations.size() ; location++)
    {
        if (locations[location][1] > top_height)
        {
            corner_indices[1] = location;
            break;
        }
    }

    // Second anticlockwise corner is corner 0
    for(unsigned location = unsigned(0.25 * mNumNodesPerCell) ; location < locations.size() ; location++)
    {
        if (locations[location][1] < top_height)
        {
            corner_indices[0] = location - 1;
            break;
        }
    }

    // Third anticlockwise corner is corner 3
    for(unsigned location = unsigned(0.5 * mNumNodesPerCell) ; location < locations.size() ; location++)
    {
        if (locations[location][1] < bot_height)
        {
            corner_indices[3] = location;
            break;
        }
    }

    // Fourth anticlockwise corner is corner 2
    for(unsigned location = unsigned(0.75 * mNumNodesPerCell) ; location < locations.size() ; location++)
    {
        if (locations[location][1] > bot_height)
        {
            corner_indices[2] = location - 1;
            break;
        }
    }

    if ( (corner_indices[0] == UINT_MAX) ||
         (corner_indices[1] == UINT_MAX) ||
         (corner_indices[2] == UINT_MAX) ||
         (corner_indices[3] == UINT_MAX) )
    {
        EXCEPTION("At least one corner not tagged properly");
    }

    if ( (corner_indices[0] < corner_indices[1]) ||
         (corner_indices[3] < corner_indices[0]) ||
         (corner_indices[2] < corner_indices[3]) )
    {
        EXCEPTION("Something went wrong when tagging corner locations");
    }

    if ( corner_indices[0] - corner_indices[1] != corner_indices[2] - corner_indices[3] )
    {
        EXCEPTION("Apical and basal surfaces are different sizes");
    }

    // Create a vector of immersed boundary elements and vector of nodes
    std::vector<ImmersedBoundaryElement<2,2>*> ib_elements;
    std::vector<Node<2>*> nodes;

    // Helper c_vector for offsets in x and y
    c_vector<double, 2> x_offset = x_unit * cell_width;
    c_vector<double, 2> y_offset = y_unit * (1.0 - cell_height) / 2.0;

    // Add the membrane element, if there is one
    if (mMembrane)
    {
        // Aim to give the membrane roughly the same node-spacing as the other cells
        double perimeter = 2.0 * (cell_height + cell_width);
        double node_spacing = perimeter / (double)mNumNodesPerCell;

        unsigned num_membrane_nodes = (unsigned)floor(1.0 / node_spacing);

        std::vector<Node<2>*> nodes_this_elem;

        for (unsigned mem_node_idx = 0 ; mem_node_idx < num_membrane_nodes ; mem_node_idx++)
        {
            // Calculate location of node
            c_vector<double, 2> location = 0.97 * y_offset + ( 0.5 / num_membrane_nodes + double(mem_node_idx) / num_membrane_nodes ) * x_unit;

            // Create the new node
            nodes.push_back(new Node<2>(mem_node_idx, location, true));
            nodes_this_elem.push_back(nodes.back());
        }

        // Create the membrane element
        ib_elements.push_back(new ImmersedBoundaryElement<2,2>(0, nodes_this_elem));
    }

    RandomNumberGenerator* p_rand_gen = RandomNumberGenerator::Instance();

    for (unsigned cell_idx = 0 ; cell_idx < mNumCellsWide ; cell_idx++)
    {
        std::vector<Node<2>*> nodes_this_elem;
        double temp_rand = p_rand_gen->ranf();

        for(unsigned location = 0 ; location < locations.size() ; location++)
        {
            unsigned node_index = nodes.size();
            c_vector<double, 2> scaled_location = 0.95 * locations[location] + x_offset * (cell_idx + 0.5);

            scaled_location[0] = fmod(scaled_location[0], 1.0);
            scaled_location[1] *= (1.0 + temp_rand * mRandomYMult);
            scaled_location[1] += y_offset[1];

            nodes.push_back(new Node<2>(node_index, scaled_location, true));
            nodes_this_elem.push_back(nodes.back());
        }

        ib_elements.push_back(new ImmersedBoundaryElement<2,2>(ib_elements.size(), nodes_this_elem));

        // Pass in the correct corners
        std::vector<Node<2>*>& r_elem_corners = ib_elements.back()->rGetCornerNodes();
        for (unsigned corner = 0 ; corner < 4 ; corner++)
        {
            r_elem_corners.push_back(ib_elements.back()->GetNode(corner_indices[corner]));
        }
    }

    // Create a mesh with the nodes and elements vectors
    if (mMembrane)
    {
        mpMesh = new ImmersedBoundaryMesh<2,2>(nodes, ib_elements, 256, 256, 0);
    }
    else
    {
        mpMesh = new ImmersedBoundaryMesh<2,2>(nodes, ib_elements, 256, 256);
    }
}

ImmersedBoundaryPalisadeMeshGenerator::~ImmersedBoundaryPalisadeMeshGenerator()
{
    delete mpMesh;
}

ImmersedBoundaryMesh<2,2>* ImmersedBoundaryPalisadeMeshGenerator::GetMesh()
{
    return mpMesh;
}

void ImmersedBoundaryPalisadeMeshGenerator::SetRandomYMult(double mult)
{
    assert(fabs(mult) < 1.0);
    mRandomYMult = mult;
}


