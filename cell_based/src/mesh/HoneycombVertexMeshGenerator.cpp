/*

Copyright (C) University of Oxford, 2005-2011

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/

#include "HoneycombVertexMeshGenerator.hpp"

HoneycombVertexMeshGenerator::HoneycombVertexMeshGenerator(unsigned numElementsAcross,
                                                           unsigned numElementsUp,
                                                           bool isFlatBottom,
                                                           double cellRearrangementThreshold,
                                                           double t2Threshold)
{
    assert(numElementsAcross > 0);
    assert(numElementsUp > 0);
    assert(cellRearrangementThreshold > 0.0);
    assert(t2Threshold > 0.0);

    std::vector<Node<2>*> nodes;
    std::vector<VertexElement<2,2>*>  elements;

    unsigned node_index = 0;
    unsigned node_indices[6];
    unsigned element_index;

    // Create the nodes, row by row, from the bottom up

    // On the first row we have numElementsAcross nodes, all of which are boundary nodes
    for (unsigned i=0; i<numElementsAcross; i++)
    {
        Node<2>* p_node = new Node<2>(node_index, true, i+0.5, 0);
        nodes.push_back(p_node);
        node_index++;
    }

    /*
     * On each interior row we have numElementsAcross+1 nodes. On the second and penultimate
     * row all nodes are boundary nodes. On other rows the first and last nodes only
     * are boundary nodes.
     */
    for (unsigned j=1; j<2*numElementsUp+1; j++)
    {
        for (unsigned i=0; i<=numElementsAcross; i++)
        {
            double x_coord = ((j%4 == 0)||(j%4 == 3)) ? i+0.5 : i;
            double y_coord = (1.5*j - 0.5*(j%2))*0.5/sqrt(3);
            bool is_boundary_node = (j==1 || j==2*numElementsUp || i==0 || i==numElementsAcross) ? true : false;

            Node<2>* p_node = new Node<2>(node_index, is_boundary_node, x_coord, y_coord);
            nodes.push_back(p_node);
            node_index++;
        }
    }

    /*
     * On the last row we have numElementsAcross nodes, all of which are boundary nodes.
     */
    double y_coord = (1.5*(2*numElementsUp+1) - 0.5*((2*numElementsUp+1)%2))*0.5/sqrt(3);
    if (((2*numElementsUp+1)%4 == 0)||((2*numElementsUp+1)%4 == 3))
    {
        Node<2>* p_node = new Node<2>(node_index, true, 0.5, y_coord);
        nodes.push_back(p_node);
        node_index++;
    }
    for (unsigned i=1; i<numElementsAcross; i++)
    {
        double x_coord = (((2*numElementsUp+1)%4 == 0)||((2*numElementsUp+1)%4 == 3)) ? i+0.5 : i;

        Node<2>* p_node = new Node<2>(node_index, true, x_coord, y_coord);
        nodes.push_back(p_node);
        node_index++;
    }
    if (((2*numElementsUp+1)%4 == 1)||((2*numElementsUp+1)%4 == 2))
    {
        Node<2>* p_node = new Node<2>(node_index, true, numElementsAcross, y_coord);
        nodes.push_back(p_node);
        node_index++;
    }

    /*
     * Create the elements. The array node_indices contains the
     * global node indices from bottom, going anticlockwise.
     */
    for (unsigned j=0; j<numElementsUp; j++)
    {
        for (unsigned i=0; i<numElementsAcross; i++)
        {
            if (j==0)
            {
                node_indices[0] = i;
            }
            else
            {
                node_indices[0] = 2*j*(numElementsAcross+1) - 1*(j%2==0) + i; // different for even/odd rows
            }
            node_indices[1] = node_indices[0] + numElementsAcross + 1 + 1*(j%2==0 && j>0);
            node_indices[2] = node_indices[1] + numElementsAcross + 1;
            node_indices[3] = node_indices[2] + numElementsAcross + 1*(j%2==1 && j<numElementsUp-1);
            node_indices[4] = node_indices[2] - 1;
            node_indices[5] = node_indices[1] - 1;

            std::vector<Node<2>*> element_nodes;
            for (unsigned k=0; k<6; k++)
            {
               element_nodes.push_back(nodes[node_indices[k]]);
            }

            element_index = j*numElementsAcross + i;
            VertexElement<2,2>* p_element = new VertexElement<2,2>(element_index, element_nodes);
            elements.push_back(p_element);
        }
    }

    mpMesh = new MutableVertexMesh<2,2>(nodes, elements, cellRearrangementThreshold, t2Threshold);
}

HoneycombVertexMeshGenerator::~HoneycombVertexMeshGenerator()
{
    delete mpMesh;
}

MutableVertexMesh<2,2>* HoneycombVertexMeshGenerator::GetMesh()
{
    return mpMesh;
}
