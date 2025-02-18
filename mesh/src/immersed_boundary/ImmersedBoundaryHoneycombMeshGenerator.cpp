
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

#include "ImmersedBoundaryHoneycombMeshGenerator.hpp"

ImmersedBoundaryHoneycombMeshGenerator::ImmersedBoundaryHoneycombMeshGenerator(unsigned numElementsX,
                                                                               unsigned numElementsY,
                                                                               unsigned numNodesPerEdge,
                                                                               double proportionalGap,
                                                                               double padding)
        : mpMesh(nullptr)
{
    // Check for sensible input
    assert(numElementsX > 0);
    assert(numElementsY > 0);
    assert(numNodesPerEdge > 2);
    assert(proportionalGap > 0.0);
    assert(padding > 0.0 && padding < 0.5);

    // Helper vectors
    unit_vector<double> x_unit(2,0);
    unit_vector<double> y_unit(2,1);

    std::vector<Node<2>*> nodes;
    std::vector<ImmersedBoundaryElement<2,2>*>  elements;

    double width = 0.5 + 1.5 * numElementsX;
    double height = sqrt(3.0) * numElementsY + 0.5 * sqrt(3.0) * (numElementsY > 1);

    double max = width > height ? width : height;

    double radius = (1.0 - 2.0 * padding) / max;

    // Reset width and height to their scaled values
    width *= radius;
    height *= radius;

    // We need to calculate the required centre of each hexagon
    std::vector<c_vector<double, 2> > offsets(numElementsY * numElementsX);

    // Specify the centre of the bottom-left most hexagon
    c_vector<double, 2> global_offset;
    if(width > height)
    {
        global_offset[0] = padding + radius;
        global_offset[1] = 0.5 - 0.5 * height + 0.5 * radius * sqrt(3.0);
    }
    else
    {
        global_offset[0] = 0.5 - 0.5 * width + radius;
        global_offset[1] = padding + 0.5 * sqrt(3.0) * radius;
    }

    // Calculate the centres
    for (unsigned x = 0; x < numElementsX; x++)
    {
        for (unsigned y = 0; y < numElementsY; y++)
        {
            unsigned idx = x * numElementsY + y;

            offsets[idx][0] = global_offset[0] + 1.5 * radius * x;
            offsets[idx][1] = global_offset[1] + 0.5 * sqrt(3.0) * radius * (2.0 * y + (double)(x % 2 == 1));
        }
    }

    // Get locations for the reference hexagon
    std::vector<c_vector<double, 2> > node_locations = GetUnitHexagon(numNodesPerEdge);

    double scale = 1.0 - proportionalGap;

    // For each calculated centre, create the nodes representing each location around that hexagon
    for (unsigned offset = 0; offset < offsets.size(); offset++)
    {
        std::vector<Node<2>*> nodes_this_elem;

        for (unsigned location = 0; location < node_locations.size(); location++)
        {
            unsigned index = offset * node_locations.size() + location;
            Node<2>* p_node = new Node<2>(index, offsets[offset] + scale * radius * node_locations[location], true);

            nodes_this_elem.push_back(p_node);
            nodes.push_back(p_node);
        }

        ImmersedBoundaryElement<2,2>* p_elem = new ImmersedBoundaryElement<2,2>(offset, nodes_this_elem);

        // Set whether it's on the boundary
        if (offset < numElementsY ||                        // left edge
            offset >= numElementsY * (numElementsX - 1) ||  // right edge
            offset % numElementsY == 0 ||                   // bottom edge
            (offset + 1) % numElementsY == 0)               // top edge
        {
            p_elem->SetIsBoundaryElement(true);
        }

        elements.push_back(p_elem);
    }

    // Create the mesh, cells, cell population, and simulation
    mpMesh = new ImmersedBoundaryMesh<2,2>(nodes, elements);
}

ImmersedBoundaryHoneycombMeshGenerator::~ImmersedBoundaryHoneycombMeshGenerator()
{
    delete mpMesh;
}

ImmersedBoundaryMesh<2,2>* ImmersedBoundaryHoneycombMeshGenerator::GetMesh()
{
    return mpMesh;
}

std::vector<c_vector<double, 2> > ImmersedBoundaryHoneycombMeshGenerator::GetUnitHexagon(unsigned numPtsPerSide)
{
    std::vector<c_vector<double, 2> > locations(numPtsPerSide * 6);

    // Find locations of the six vertices (and the seventh is the same as the first)
    std::vector<c_vector<double, 2> > vertices(7);
    double sixty_degrees = M_PI / 3.0;
    for (unsigned vertex = 0; vertex < vertices.size(); vertex++)
    {
        vertices[vertex][0] = cos((double)vertex * sixty_degrees);
        vertices[vertex][1] = sin((double)vertex * sixty_degrees);
    }

    // Between each pair of vertices, fill in the specified number of locations evenly spaced along the edge
    for (unsigned vertex = 0; vertex < 6; vertex++)
    {
        c_vector<double, 2> this_vertex = vertices[vertex];
        c_vector<double, 2> next_vertex = vertices[vertex + 1];
        c_vector<double, 2> vec_between = next_vertex - this_vertex;

        for (unsigned i = 0; i < numPtsPerSide; i++)
        {
            locations[vertex * numPtsPerSide + i] = this_vertex + i * vec_between / (double)numPtsPerSide;
        }
    }
    return locations;
}
