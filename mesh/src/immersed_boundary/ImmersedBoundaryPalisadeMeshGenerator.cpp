
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

#include "ImmersedBoundaryPalisadeMeshGenerator.hpp"
#include "ImmersedBoundaryEnumerations.hpp"
#include "RandomNumberGenerator.hpp"

ImmersedBoundaryPalisadeMeshGenerator::ImmersedBoundaryPalisadeMeshGenerator(unsigned numCellsWide,
                                                                             unsigned numNodesPerCell,
                                                                             double ellipseExponent,
                                                                             double cellAspectRatio,
                                                                             double randomYMult,
                                                                             bool basalLamina,
                                                                             bool apicalLamina,
                                                                             bool leakyLaminas,
                                                                             unsigned numFluidMeshPoints,
                                                                             double absoluteGap)
    : mpMesh(NULL)
{
    // Check for sensible input
    assert(numCellsWide > 0);
    assert(numNodesPerCell > 3);
    assert(ellipseExponent > 0.0);
    assert(cellAspectRatio > 0.0); // aspect ratio is cell height / cell width
    assert(fabs(randomYMult) < 2.0);
    assert( (leakyLaminas && numFluidMeshPoints < UINT_MAX) || (!leakyLaminas) );

    // If the number of fluid mesh points is specified, calculate an appropriate number of nodes per cell automatically
    bool override_nodes_per_cell = numFluidMeshPoints < UINT_MAX;

    // Helper vectors
    unit_vector<double> x_unit(2,0);
    unit_vector<double> y_unit(2,1);

    // Calculate cell sizes and rescale if necessary
    double cell_width = 1.0 / double(numCellsWide);
    double cell_height = cellAspectRatio * cell_width;

    if (cell_height > 0.8)
    {
        cell_height = 0.8;
        cell_width = cell_height / cellAspectRatio;
    }

    // If we are overriding the number of nodes per cell, alter it so that the spacing ratio is approx 0.5
    if (override_nodes_per_cell)
    {
        double cell_perimeter = 2.0 * (cell_height + cell_width);
        double fluid_mesh_spacing = 1.0 / static_cast<double>(numFluidMeshPoints);
        double target_node_spacing = 0.5 * fluid_mesh_spacing;
        numNodesPerCell = static_cast<unsigned>(cell_perimeter / target_node_spacing);
    }

    const double width_scale_fac = absoluteGap == DOUBLE_UNSET ? 0.95 : (cell_width - absoluteGap) / cell_width;

    // Generate a reference superellipse
    std::unique_ptr<SuperellipseGenerator> p_gen(new SuperellipseGenerator(numNodesPerCell, ellipseExponent, cell_width, cell_height, 0.0, 0.0));
    std::vector<c_vector<double, 2> > locations = p_gen->GetPointsAsVectors();

    /*
     * The top and bottom heights are the heights at which there is maximum curvature in the superellipse.
     * These define the apical and basal domains.
     *       _____
     *    __/_____\__ Apical cutoff
     *    __|_____|__ Periapical cutoff
     *      |     |
     *      |     |
     *    __|_____|__ Basal cutoff
     *      \_____/__ y = 0
     */
    double apical_cutoff = p_gen->GetHeightOfTopSurface();
    double basal_cutoff = cell_height - apical_cutoff;
    double periapical_cutoff = apical_cutoff - basal_cutoff;

    if (apical_cutoff < 0.5 * cell_height || apical_cutoff > cell_height || periapical_cutoff < basal_cutoff)
    {
        EXCEPTION("Something went wrong calculating the height of top surface of the cell."); //LCOV_EXCL_LINE
    }

    // Corner 0 is left-apical, 1 is right-apical, 2 is right-basal, 3 is left-basal
    std::vector<unsigned> corner_indices(4, UINT_MAX);

    // Node 0 is middle of right lateral side

    // First anticlockwise is RIGHT_APICAL_CORNER
    for (unsigned location = 0; location < locations.size(); location++)
    {
        if (locations[location][1] > apical_cutoff)
        {
            corner_indices[RIGHT_APICAL_CORNER] = location;
            break;
        }
    }

    // Second anticlockwise is LEFT_APICAL_CORNER
    for (unsigned location = unsigned(0.25 * numNodesPerCell); location < locations.size(); location++)
    {
        if (locations[location][1] < apical_cutoff)
        {
            corner_indices[LEFT_APICAL_CORNER] = location - 1;
            break;
        }
    }

    // Third anticlockwise is LEFT_BASAL_CORNER
    for (unsigned location = unsigned(0.5 * numNodesPerCell); location < locations.size(); location++)
    {
        if (locations[location][1] < basal_cutoff)
        {
            corner_indices[LEFT_BASAL_CORNER] = location;
            break;
        }
    }

    // Fourth anticlockwise is RIGHT_BASAL_CORNER
    for (unsigned location = unsigned(0.75 * numNodesPerCell); location < locations.size(); location++)
    {
        if (locations[location][1] > basal_cutoff)
        {
            corner_indices[RIGHT_BASAL_CORNER] = location - 1;
            break;
        }
    }

    if ( (corner_indices[0] == UINT_MAX) ||
         (corner_indices[1] == UINT_MAX) ||
         (corner_indices[2] == UINT_MAX) ||
         (corner_indices[3] == UINT_MAX) )
    {
        EXCEPTION("At least one corner not tagged"); //LCOV_EXCL_LINE
    }

    // Should have RIGHT_APICAL_CORNER < LEFT_APICAL_CORNER < LEFT_BASAL_CORNER < RIGHT_BASAL_CORNER
    if ( (corner_indices[RIGHT_APICAL_CORNER] > corner_indices[LEFT_APICAL_CORNER]) ||
         (corner_indices[LEFT_APICAL_CORNER] > corner_indices[LEFT_BASAL_CORNER]) ||
         (corner_indices[LEFT_BASAL_CORNER] > corner_indices[RIGHT_BASAL_CORNER]) )
    {
        EXCEPTION("Something went wrong when tagging corner locations"); //LCOV_EXCL_LINE
    }

    if ( corner_indices[LEFT_APICAL_CORNER] - corner_indices[RIGHT_APICAL_CORNER] !=
         corner_indices[RIGHT_BASAL_CORNER] - corner_indices[LEFT_BASAL_CORNER] )
    {
        EXCEPTION("Apical and basal surfaces are different sizes"); //LCOV_EXCL_LINE
    }

    // Create vectors of immersed boundary elements, laminas, nodes
    std::vector<ImmersedBoundaryElement<2,2>*> ib_elements;
    std::vector<ImmersedBoundaryElement<1,2>*> ib_laminas;
    std::vector<Node<2>*> nodes;

    // Helper c_vector for offsets in x and y
    c_vector<double, 2> x_offset = x_unit * cell_width;
    c_vector<double, 2> y_offset = y_unit * (1.0 - cell_height) / 2.0;

    // Aim to give the laminas roughly the same node-spacing as the other cells...
    double lamina_node_spacing = 2.0 * (cell_height + cell_width) / numNodesPerCell;
    //... unless the laminas are to be leaky, in which case the spacing should be bigger
    if (leakyLaminas && override_nodes_per_cell)
    {
        // If laminas are leaky, increase the spacing so fluid can flow through
        lamina_node_spacing = 2.0 / static_cast<double>(numFluidMeshPoints);
    }

    // Add the basal lamina, if there is one
    if (basalLamina)
    {
        // The height of the lamina is offset by a proportion of the cell height
        double lam_height = absoluteGap == DOUBLE_UNSET ? y_offset[1] - 0.05 * cell_height : y_offset[1] - absoluteGap;

        unsigned num_lamina_nodes = static_cast<unsigned>(floor(1.0 / lamina_node_spacing));

        std::vector<Node<2>*> nodes_this_elem;

        for (unsigned node_idx = 0; node_idx < num_lamina_nodes; node_idx++)
        {
            // Calculate location of node
            c_vector<double, 2> location = lam_height * y_unit + ( 0.5 / num_lamina_nodes + double(node_idx) / num_lamina_nodes ) * x_unit;

            // Create the new node
            Node<2>* p_node = new Node<2>(nodes.size(), location, true);
            p_node->SetRegion(LAMINA_REGION);

            nodes.push_back(p_node);
            nodes_this_elem.push_back(p_node);
        }

        // Create the lamina element
        ib_laminas.push_back(new ImmersedBoundaryElement<1,2>(0, nodes_this_elem));
    }

    // Add the apical lamina, if there is one
    if (apicalLamina)
    {
        if (randomYMult != 0.0)
        {
            for (auto& p_node : nodes)
            {
                delete p_node;
            }
            for (auto& p_lam : ib_laminas)
            {
                delete p_lam;
            }
            EXCEPTION("Currently no random y variation allowed with an apical lamina");
        }

        // The height of the lamina is offset by a proportion of the cell height
        double lam_height = y_offset[1] + cell_height;

        unsigned num_lamina_nodes = static_cast<unsigned>(floor(1.0 / lamina_node_spacing));

        std::vector<Node<2>*> nodes_this_elem;

        for (unsigned node_idx = 0; node_idx < num_lamina_nodes; node_idx++)
        {
            // Calculate location of node
            c_vector<double, 2> location = lam_height * y_unit + ( 0.5 / num_lamina_nodes + double(node_idx) / num_lamina_nodes ) * x_unit;

            // Create the new node
            Node<2>* p_node = new Node<2>(nodes.size(), location, true);
            p_node->SetRegion(LAMINA_REGION);

            nodes.push_back(p_node);
            nodes_this_elem.push_back(p_node);
        }

        // Create the lamina element
        ib_laminas.push_back(new ImmersedBoundaryElement<1,2>(1, nodes_this_elem));
    }

    RandomNumberGenerator* p_rand_gen = RandomNumberGenerator::Instance();

    for (unsigned elem_idx = 0; elem_idx < numCellsWide; elem_idx++)
    {
        std::vector<Node<2>*> nodes_this_elem;
        double random_variation = p_rand_gen->ranf() * randomYMult;
        c_vector<double, 2> scaled_location;

        for (unsigned location = 0; location < locations.size(); location++)
        {
            unsigned node_index = nodes.size();
            scaled_location = width_scale_fac * locations[location] + x_offset * (elem_idx + 0.5);
            scaled_location[0] = fmod(scaled_location[0], 1.0);
            scaled_location[1] *= (1.0 + random_variation);
            scaled_location[1] += y_offset[1];

            // Create the new node
            Node<2>* p_node = new Node<2>(node_index, scaled_location, true);

            // Tag the node region
            if (locations[location][1] <= basal_cutoff)
            {
                if (locations[location][0] < 0.5 * cell_width)
                {
                    p_node->SetRegion(LEFT_BASAL_REGION);
                }
                else
                {
                    p_node->SetRegion(RIGHT_BASAL_REGION);
                }
            }
            else if (locations[location][1] <= periapical_cutoff)
            {
                if (locations[location][0] < 0.5 * cell_width)
                {
                    p_node->SetRegion(LEFT_LATERAL_REGION);
                }
                else
                {
                    p_node->SetRegion(RIGHT_LATERAL_REGION);
                }
            }
            else if (locations[location][1] < apical_cutoff)
            {
                if (locations[location][0] < 0.5 * cell_width)
                {
                    p_node->SetRegion(LEFT_PERIAPICAL_REGION);
                }
                else
                {
                    p_node->SetRegion(RIGHT_PERIAPICAL_REGION);
                }
            }
            else
            {
                if (locations[location][0] < 0.5 * cell_width)
                {
                    p_node->SetRegion(LEFT_APICAL_REGION);
                }
                else
                {
                    p_node->SetRegion(RIGHT_APICAL_REGION);
                }
            }

            nodes.push_back(p_node);
            nodes_this_elem.push_back(p_node);
        }

        ib_elements.push_back(new ImmersedBoundaryElement<2,2>(elem_idx, nodes_this_elem));

        // Pass in the corners calculated above
        std::vector<Node<2>*>& r_elem_corners = ib_elements.back()->rGetCornerNodes();
        for (unsigned corner = 0; corner < 4; corner++)
        {
            r_elem_corners.push_back(ib_elements.back()->GetNode(corner_indices[corner]));
        }
    }

    if (override_nodes_per_cell)
    {
        mpMesh = new ImmersedBoundaryMesh<2, 2>(nodes, ib_elements, ib_laminas, numFluidMeshPoints, numFluidMeshPoints);
    }
    else
    {
        mpMesh = new ImmersedBoundaryMesh<2, 2>(nodes, ib_elements, ib_laminas);
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