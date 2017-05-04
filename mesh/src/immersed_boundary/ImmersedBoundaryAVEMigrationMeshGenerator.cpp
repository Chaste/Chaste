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

// Created by bartosz on 07/02/17.

#include "ImmersedBoundaryAVEMigrationMeshGenerator.hpp"
#include "ImmersedBoundaryEnumerations.hpp"

ImmersedBoundaryAVEMigrationMeshGenerator::ImmersedBoundaryAVEMigrationMeshGenerator(unsigned numCellsHorizontally,
                                                                                     unsigned numCellsVertically,
                                                                                     unsigned numNodesPerCell,
                                                                                     double ellipseExponent,
                                                                                     double spacing,
                                                                                     double shift)
{
    // Check for sensible input
    assert(numCellsHorizontally > 0);
    assert(numCellsVertically > 0);
    assert(numNodesPerCell > 3);
    assert(ellipseExponent > 0);
    assert(spacing > 0);

    // Helper vectors
    unit_vector<double> x_unit(2, 0);
    unit_vector<double> y_unit(2, 1);

    // Calculate the size of cells
    double cell_width =  (1.0 - (numCellsHorizontally) * spacing)/numCellsHorizontally;
    double cell_height =  (1.0 - (numCellsVertically) * spacing)/numCellsVertically;

    // Set up a reference cell
    SuperellipseGenerator* p_gen = new SuperellipseGenerator(numNodesPerCell, ellipseExponent, cell_width, cell_height, 0.5*spacing, 0.5*spacing);
    std::vector<c_vector<double, 2> > locations =  p_gen->GetPointsAsVectors();

    // Set up the containers for holding the nodes and ib_elements
    std::vector<Node<2>*> nodes;
    std::vector<ImmersedBoundaryElement<2, 2>*> ib_elements;
    std::vector<ImmersedBoundaryElement<1, 2>*> ib_laminas;

    for (unsigned vert_index = 0; vert_index < numCellsVertically; vert_index++)
    {
        for (unsigned horz_index = 0; horz_index < numCellsHorizontally; horz_index++)
        {
            std::vector<Node<2>*> nodes_this_elem;
            // Shift the reference cell to an appropriate place
            for ( unsigned location_index = 0; location_index < locations.size(); location_index++)
            {
//                c_vector<double, 2> position = locations[location_index] + horz_index * (cell_width + spacing) * x_unit + vert_index * (cell_height + spacing) * y_unit;

                c_vector<double, 2> position;
                if (vert_index % 2 != 0)
                {
                    position = locations[location_index] + horz_index * (cell_width + spacing) * x_unit + vert_index * (cell_height + spacing) * y_unit;
                }
                else
                {
                    position = locations[location_index] + horz_index * (cell_width + spacing) * x_unit + vert_index * (cell_height + spacing) * y_unit;
                    position[0] = fmod(position[0] + shift * (cell_width + spacing), 1.0);
                }

                Node<2>* p_node = new Node<2>(nodes.size(), position, true);
                nodes.push_back(p_node);
                nodes_this_elem.push_back(nodes.back());
            }
            ib_elements.push_back(new ImmersedBoundaryElement<2, 2>(ib_elements.size(), nodes_this_elem));
        }
    }

    mpMesh = new ImmersedBoundaryMesh<2,2>(nodes, ib_elements, ib_laminas);
}

ImmersedBoundaryAVEMigrationMeshGenerator::~ImmersedBoundaryAVEMigrationMeshGenerator()
{
    delete mpMesh;
}

ImmersedBoundaryMesh<2, 2>* ImmersedBoundaryAVEMigrationMeshGenerator::GetMesh()
{
    return mpMesh;
}