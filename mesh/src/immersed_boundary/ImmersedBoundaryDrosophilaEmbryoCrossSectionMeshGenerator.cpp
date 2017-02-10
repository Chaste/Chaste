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

// Created by bartosz on 09/02/17.

#include "ImmersedBoundaryDrosophilaEmbryoCrossSectionMeshGenerator.hpp"
#include "ImmersedBoundaryEnumerations.hpp"

ImmersedBoundaryDrosophilaEmbryoCrossSectionMeshGenerator::ImmersedBoundaryDrosophilaEmbryoCrossSectionMeshGenerator(unsigned numCells,
                                                                                                                     unsigned numNodesPerCell,
                                                                                                                     double cellHeight,
                                                                                                                     double radius,
                                                                                                                     double aspectRatio,
                                                                                                                     double ellipseExponent) : mpMesh(NULL)
{
    // Check for sensible input
    assert(numCells > 0);
    assert(numNodesPerCell > 3);
    assert(cellHeight > 0);
    assert(radius >0);
    assert(aspectRatio > 0);
    assert(ellipseExponent > 0);

    // Calculate the necessary parameters
    double membrane_width_base = 2.0 * aspectRatio * radius;
    double membrane_height_base = 2.0 * radius;
    double membrane_width_top = 2.0 * aspectRatio * radius + 2.0 * cellHeight;
    double membrane_height_top = 2.0 * radius + 2.0 * cellHeight;

    double membrane_base_arc_length = CalculateTotalArcLengthOfSuperellipse(membrane_width_base, membrane_height_base, ellipseExponent);
    double membrane_top_arc_length = CalculateTotalArcLengthOfSuperellipse(membrane_width_top, membrane_height_top, ellipseExponent);

    double total_cell_height_length = 2.0 * numCells * cellHeight;
    double proportion_base = membrane_base_arc_length / (membrane_base_arc_length + membrane_top_arc_length + total_cell_height_length);
    double proportion_top = membrane_top_arc_length / (membrane_base_arc_length + membrane_top_arc_length + total_cell_height_length);
    unsigned total_num_nodes = numCells * numNodesPerCell;

    unsigned membrane_base_num_nodes = proportion_base * total_num_nodes;
    unsigned membrane_top_num_nodes = proportion_top * total_num_nodes;

    // Create the circles around which the cells will be placed
    SuperellipseGenerator* p_gen_base = new SuperellipseGenerator(membrane_base_num_nodes, ellipseExponent, membrane_width_base, membrane_height_base, 0.5 * (1.0 - membrane_width_base), 0.5 * (1.0 - membrane_height_base));
    std::vector<c_vector<double, 2> > base_locations = p_gen_base->GetPointsAsVectors();
    delete p_gen_base;

    SuperellipseGenerator* p_gen_top = new SuperellipseGenerator(membrane_top_num_nodes, ellipseExponent, membrane_width_top, membrane_height_top, 0.5 * (1.0 - membrane_width_top), 0.5 * (1.0 - membrane_height_top));
    std::vector<c_vector<double, 2> > top_locations = p_gen_top->GetPointsAsVectors();
    delete p_gen_top;

    // Create the containers for the nodes, immersed boundary elements and laminas
    std::vector<Node<2>*> nodes;
    std::vector<ImmersedBoundaryElement<2,2>*> ib_elements;
    std::vector<ImmersedBoundaryElement<1,2>*> ib_laminas;

    for (unsigned i = 0; i < top_locations.size(); i++)
    {
        Node<2>* p_node = new Node<2>(i, top_locations[i], true);
        nodes.push_back(p_node);
    }

    ib_elements.push_back(new ImmersedBoundaryElement<2,2>(0, nodes));

    mpMesh = new ImmersedBoundaryMesh<2,2>(nodes, ib_elements, ib_laminas);

}

ImmersedBoundaryDrosophilaEmbryoCrossSectionMeshGenerator::~ImmersedBoundaryDrosophilaEmbryoCrossSectionMeshGenerator()
{
    delete mpMesh;
}

ImmersedBoundaryMesh<2,2>* ImmersedBoundaryDrosophilaEmbryoCrossSectionMeshGenerator::GetMesh()
{
    return mpMesh;
}

double ImmersedBoundaryDrosophilaEmbryoCrossSectionMeshGenerator::CalculateTotalArcLengthOfSuperellipse(double membraneWidth,
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

c_vector<double, 2> ImmersedBoundaryDrosophilaEmbryoCrossSectionMeshGenerator::CalculateNormal(c_vector<double, 2> point_1,
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

const std::vector<c_vector<double, 2> > ImmersedBoundaryDrosophilaEmbryoCrossSectionMeshGenerator::BuildACell(double targetNodeSpacing,
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