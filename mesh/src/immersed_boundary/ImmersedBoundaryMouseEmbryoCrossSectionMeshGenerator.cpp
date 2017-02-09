/*

Copyright (c) 2005-2017, University of Oxford.
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

// Created by bartosz on 29/01/17.

#include "ImmersedBoundaryMouseEmbryoCrossSectionMeshGenerator.hpp"
#include "ImmersedBoundaryEnumerations.hpp"

ImmersedBoundaryMouseEmbryoCrossSectionMeshGenerator::ImmersedBoundaryMouseEmbryoCrossSectionMeshGenerator(unsigned numCells,
                                                                                                           unsigned numNodesMembrane,
                                                                                                           double ellipseExponent,
                                                                                                           double membraneWidth,
                                                                                                           double membraneHeight,
                                                                                                           double cellHeight)
        : mpMesh(NULL)
{
    // Check for sensible input
    assert(numCells > 0);
    assert(numNodesMembrane > 3);
    assert(membraneWidth > 0.0);
    assert(membraneHeight > 0.0);
    assert(cellHeight > 0.0);

    if (membraneWidth > 1.0 || membraneHeight > 1.0)
    {
        EXCEPTION("Membrane Width exceeds domain bounds [0, 1]x[0, 1]");
    }

    if (membraneWidth + 2.0 * cellHeight > 1.0 || membraneHeight + 2.0 * cellHeight > 1.0)
    {
        EXCEPTION("The cells are too large! Decrease the cellHeight parameter.");
    }

    // Helper vectors
    unit_vector<double> x_unit(2, 0);
    unit_vector<double> y_unit(2, 1);

    // Test the CalculateProportionOfPointsInSuperellipse member function
    double test_arc_length = CalculateTotalArcLengthOfSuperellipse(1.0, 1.0, 1.0);
    if (fabs(M_PI - test_arc_length) > 1e-5)
    {
        EXCEPTION("The CalculateTotalArcLengthOfSuperellipse function does not work properly!");
    }

    // Calculate the number of nodes in part of the membrane
    double total_arc_length = CalculateTotalArcLengthOfSuperellipse(membraneWidth, 2.0 * membraneHeight, ellipseExponent);
    double proportion = total_arc_length / (total_arc_length + 2.0 * membraneWidth);
    unsigned num_nodes_in_curved_part = unsigned(numNodesMembrane * proportion);
    unsigned num_nodes_in_line_part = numNodesMembrane - num_nodes_in_curved_part;

    // Initialise the membrane
    SuperellipseGenerator* p_gen = new SuperellipseGenerator(2 * num_nodes_in_curved_part, ellipseExponent,
                                                             membraneWidth, 2.0 * membraneHeight,
                                                             0.5 * (1.0 - membraneWidth), 0.5 * (1.0 - membraneHeight));
    std::vector<c_vector<double, 2> > locations = p_gen->GetPointsAsVectors();
    double node_spacing = p_gen->GetTargetNodeSpacing();

    // Create vectors of immersed boundary elements, laminas, nodes
    std::vector<ImmersedBoundaryElement<2, 2>*> ib_elements;
    std::vector<ImmersedBoundaryElement<1, 2>*> ib_laminas;
    std::vector<Node<2>*> nodes;

    // Position nodes on a line that represents the top of the cross-section (supposed to be stationary)
    for (unsigned i = 0; i < num_nodes_in_line_part; i++)
    {
        unsigned node_index = nodes.size();
        c_vector<double, 2> position = locations[0] - i * node_spacing * x_unit;
        Node<2> *p_node = new Node<2>(node_index, position, true);
        nodes.push_back(p_node);
    }

    // Position nodes in a half of the superellipse
    for (unsigned location = 0.5 * locations.size() + 1; location < locations.size(); location++)
    {
        unsigned node_index = nodes.size();
        Node<2> *p_node = new Node<2>(node_index, locations[location], true);
        nodes.push_back(p_node);
    }

    ib_elements.push_back(new ImmersedBoundaryElement<2, 2>(0, nodes));

    // Calculate the number of nodes per cell
    unsigned num_mem_nodes_per_cell = (0.5 * locations.size()) / double(numCells);
    unsigned shift = 0.5*((0.5 * locations.size()) - num_mem_nodes_per_cell * numCells);

    // Place the cells around the membrane
    for (unsigned cell_index = 0; cell_index < numCells; cell_index++)
    {
        // Now add the cells around the membrane, first add all the base points
        std::vector<c_vector<double, 2> > base;
        for (unsigned i = 0; i < num_mem_nodes_per_cell; i++)
        {
            base.push_back(locations[shift + (0.5 * locations.size() + 1) + (cell_index * num_mem_nodes_per_cell) + i]);
        }

        // Calculate the normals to the membrane
        c_vector<double, 2> normal_1 = CalculateNormal(locations[shift + (0.5 * locations.size() + 1) + cell_index * num_mem_nodes_per_cell],
                                                       locations[shift + (0.5 * locations.size() + 2) + cell_index * num_mem_nodes_per_cell]);
        c_vector<double, 2> normal_2 = CalculateNormal(locations[shift + (0.5 * locations.size() + 1) + (cell_index + 1) * num_mem_nodes_per_cell],
                                                       locations[shift + (0.5 * locations.size() + 2) + (cell_index + 1) * num_mem_nodes_per_cell]);

        std::vector<c_vector<double, 2> > cell_locations = BuildACell(node_spacing, base, normal_1, normal_2, cellHeight);

        unsigned index_of_last_node = nodes.size();

        // Add these to the list of nodes and IB elements
        std::vector<Node<2>*> nodes_this_elem;
        for (unsigned node_index = 0; node_index < cell_locations.size(); node_index++)
        {
            Node<2> *p_node = new Node<2>(index_of_last_node + node_index, cell_locations[node_index], true);
            nodes_this_elem.push_back(p_node);
            nodes.push_back(nodes_this_elem.back());
        }

        ib_elements.push_back(new ImmersedBoundaryElement<2, 2>(ib_elements.size(), nodes_this_elem));

    }

    mpMesh = new ImmersedBoundaryMesh<2,2>(nodes, ib_elements, ib_laminas);
}

ImmersedBoundaryMouseEmbryoCrossSectionMeshGenerator::~ImmersedBoundaryMouseEmbryoCrossSectionMeshGenerator()
{
    delete mpMesh;
}

ImmersedBoundaryMesh<2, 2>* ImmersedBoundaryMouseEmbryoCrossSectionMeshGenerator::GetMesh()
{
    return mpMesh;
}

double ImmersedBoundaryMouseEmbryoCrossSectionMeshGenerator::CalculateTotalArcLengthOfSuperellipse(double membraneWidth,
                                                                                                    double membraneHeight,
                                                                                                    double ellipseExponent)
{

    unsigned dense_pts = 50000;

    // scale the value of membraneHeight, as we need half of the superellipse

    // Vectors to store information from the loop
    std::vector<double> cumulative_arc_length(dense_pts);
    c_vector<double, 2> dense_locations;

    // Helper variables for the loop
    double dense_spacing = (2.0 * M_PI) / (double)dense_pts;
    double cumulative_dist = 0.0;
    double temp_sin, temp_cos;

    // Fill in first location by hand
    cumulative_arc_length[0] = 0.0;
    dense_locations[0] = membraneWidth * 0.5;
    dense_locations[1] = 0.0;
    c_vector<double, 2> first_point = dense_locations;

    // Fill in all other locations
    for (unsigned point = 1; point < dense_pts; point++)
    {
        // Update cumulative distance
        cumulative_dist += dense_spacing;

        // Get the current sin and cos values
        temp_cos = cos(cumulative_dist);
        temp_sin = sin(cumulative_dist);

        c_vector<double, 2> previous_dense_locations = dense_locations;

        // x and y are powers of cos and sin, with the sign of sin and cos
        dense_locations[0] = copysign(pow(fabs(temp_cos), ellipseExponent), temp_cos) * membraneWidth * 0.5;
        dense_locations[1] = copysign(pow(fabs(temp_sin), ellipseExponent), temp_sin) * membraneHeight * 0.5;

        // Check that the new point created lies within [-width/2, width/2] x [-height/2, height/2]
        assert(fabs(dense_locations[0]) < membraneWidth * 0.5 + 1e-10);
        assert(fabs(dense_locations[1]) < membraneHeight * 0.5 + 1e-10);

        // Fill in arc length property
        cumulative_arc_length[point] = cumulative_arc_length[point - 1] + norm_2(dense_locations - previous_dense_locations);
    }

    double total_arc_length = cumulative_arc_length[dense_pts - 1] + norm_2(dense_locations - first_point);

    return total_arc_length;
}

c_vector<double, 2> ImmersedBoundaryMouseEmbryoCrossSectionMeshGenerator::CalculateNormal(c_vector<double, 2> point_1,
                                                                                          c_vector<double, 2> point_2)
{
    // Calculate the tangent
    c_vector<double, 2> tangent = point_2 - point_1;

    c_vector<double, 2> normal;

    // Calculate the vector perpendicular to the tangent
    normal[0] = tangent[1];
    normal[1] = -tangent[0];

    // Normalize this vector
    normal /= norm_2(normal);

    return normal;
}

const std::vector<c_vector<double, 2> > ImmersedBoundaryMouseEmbryoCrossSectionMeshGenerator::BuildACell(double targetNodeSpacing,
                                                                                                         std::vector<c_vector<double, 2> > base,
                                                                                                         c_vector<double, 2> normal_1,
                                                                                                         c_vector<double, 2> normal_2,
                                                                                                         double cell_height) const
{
    // Check that both vector are of unit length
    assert(fabs(norm_2(normal_1) - 1.0) < 1.0e-5);
    assert(fabs(norm_2(normal_2) - 1.0) < 1.0e-5);

    // Create a container to hold all the points
    std::vector<c_vector<double, 2> > cell;

    // First move the base in the direction normal to the mean of the two normals
    c_vector<double, 2> avg_normal = 0.5 * (normal_1 + normal_2);
    for (unsigned i = 0; i < base.size(); i++)
    {
        base[i] += (targetNodeSpacing * avg_normal);
        cell.push_back(base[i]);
    }

    // Calculate the number of nodes for the sides of a cell
    unsigned num_nodes_side = round(cell_height/targetNodeSpacing);

    std::vector<c_vector<double, 2> > side_1;
    std::vector<c_vector<double, 2> > side_2;

    // Create the sides of the cell
    for (unsigned i = 0; i < num_nodes_side; i++)
    {
        side_1.push_back(base[base.size()-1] + normal_2 * i * targetNodeSpacing);
        side_2.push_back(base[0] + normal_1 * i * targetNodeSpacing);
    }

    // Calculate the number of nodes for the top of a cell
    c_vector<double, 2> top_vec = side_2[num_nodes_side-1] - side_1[num_nodes_side-1];
    double top_length = norm_2(top_vec);
    top_vec /= norm_2(top_vec);
    unsigned num_nodes_top = round(top_length/targetNodeSpacing);

    // Add the points of the first side
    for (unsigned i = 1; i < num_nodes_side; i++)
    {
        cell.push_back(side_1[i]);
    }

    // Add the points of the top of the cell
    for (unsigned i = 1; i < num_nodes_top; i++)
    {
        c_vector<double, 2> vec = side_1[num_nodes_side-1] + top_vec * i * targetNodeSpacing;
        cell.push_back(vec);
    }

    // Add the points of the other side of the cell
    for (unsigned i = 1; i < num_nodes_side; i++)
    {
        cell.push_back(side_2[num_nodes_side - i]);
    }

    return cell;

}
