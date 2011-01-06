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

#include "PottsMeshGenerator.hpp"


PottsMeshGenerator::PottsMeshGenerator(unsigned numNodesAcross, unsigned numNodesUp,
                                       unsigned numElementsAcross, unsigned numElementsUp,
                                       unsigned elementWidth, unsigned elementHeight)
{
    assert(numElementsAcross > 0);
    assert(numElementsUp > 0);
    assert(elementWidth > 0);
    assert(elementHeight > 0);
    assert(numElementsAcross*elementWidth==numNodesAcross);
    assert(numElementsUp*elementHeight<=numNodesUp);

    std::vector<Node<2>*> nodes;
    std::vector<PottsElement*>  elements;

    unsigned node_index = 0;
    unsigned node_indices[elementWidth*elementHeight];
    unsigned element_index;

    /*
     * Create the nodes, row by row, from the bottom up
     * On the first and last row we have numNodesAcross nodes, all of which are boundary
     * nodes. On each interior row we have numNodesAcross nodes, the first and last nodes
     * are boundary nodes.
     */
    for (unsigned j=0; j<numNodesUp; j++)
    {
        for (unsigned i=0; i<numNodesAcross; i++)
        {
            bool is_boundary_node = (j==0 || j==numNodesUp-1 || i==0 || i==numNodesAcross-1) ? true : false;

            Node<2>* p_node = new Node<2>(node_index, is_boundary_node, i, j);
            nodes.push_back(p_node);
            node_index++;
        }
    }


    /*
     * Create the elements. The array node_indices contains the
     * global node indices, in increasing order.
     */
    for (unsigned j=0; j<numElementsUp; j++)
    {
        for (unsigned i=0; i<numElementsAcross; i++)
        {
            for (unsigned l=0; l<elementHeight; l++)
            {
                for (unsigned k=0; k<elementWidth; k++)
                {
                    node_indices[l*elementWidth + k] = j*elementHeight*numNodesAcross + i*elementWidth
                                                       + l*numNodesAcross + k;
                }
            }

            std::vector<Node<2>*> element_nodes;
            for (unsigned k=0; k<elementHeight*elementWidth; k++)
            {
               element_nodes.push_back(nodes[node_indices[k]]);
            }

            element_index = j*numElementsAcross + i;
            PottsElement* p_element = new PottsElement(element_index, element_nodes);
            elements.push_back(p_element);
        }
    }
    mpMesh = new PottsMesh(nodes, elements);
}


PottsMeshGenerator::~PottsMeshGenerator()
{
    delete mpMesh;
}

PottsMesh* PottsMeshGenerator::GetMesh()
{
    return mpMesh;
}



