/*

Copyright (C) University of Oxford, 2005-2012

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

#include "CylindricalHoneycombVertexMeshGenerator.hpp"

CylindricalHoneycombVertexMeshGenerator::CylindricalHoneycombVertexMeshGenerator(unsigned numElementsAcross,
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

    assert(numElementsAcross > 1);
    assert(numElementsAcross%2==0); // numElementsAcross must be even for cylindrical meshes

    // Create the nodes
    for (unsigned j=0; j<=2*numElementsUp+1; j++)
    {
        if (isFlatBottom && (j==1))
        {
            // Flat bottom to cylindrical mesh
            for (unsigned i=0; i<=numElementsAcross-1; i++)
            {
                Node<2>* p_node = new Node<2>(node_index, true, i, 0.0);
                nodes.push_back(p_node);
                node_index++;
            }
        }
        /*
         * On each interior row we have numElementsAcross+1 nodes. On the first second, penultimate,
         *  and last rows all nodes are boundary nodes. On other rows no nodes are boundary nodes.
         */
        else
        {
            for (unsigned i=0; i<=numElementsAcross-1; i++)
            {
                double x_coord = ((j%4 == 0)||(j%4 == 3)) ? i+0.5 : i;
                double y_coord = (1.5*j - 0.5*(j%2))*0.5/sqrt(3);
                bool is_boundary_node = (j==0 || j==1 || j==2*numElementsUp || j==2*numElementsUp+1) ? true : false;

                Node<2>* p_node = new Node<2>(node_index, is_boundary_node , x_coord, y_coord);
                nodes.push_back(p_node);
                node_index++;
            }
        }
    }

    /*
     * Create the elements. The array node_indices contains the
     * global node indices from bottom, going anticlockwise.
     */
    for (unsigned j=0; j<numElementsUp; j++)
    {
        for (unsigned i=0; i<numElementsAcross; i++)
        {
            element_index = j*numElementsAcross + i;

            node_indices[0] = 2*j*(numElementsAcross) + i + 1*(j%2==1);
            node_indices[1] = node_indices[0] + numElementsAcross + 1*(j%2==0);
            node_indices[2] = node_indices[0] + 2*numElementsAcross + 1*(j%2==0);
            node_indices[3] = node_indices[0] + 3*numElementsAcross;
            node_indices[4] = node_indices[0] + 2*numElementsAcross - 1*(j%2==1);
            node_indices[5] = node_indices[0] + numElementsAcross - 1*(j%2==1);

            if (i==numElementsAcross-1) // on far right
            {
                node_indices[0] -= numElementsAcross*(j%2==1);
                node_indices[1] -= numElementsAcross;
                node_indices[2] -= numElementsAcross;
                node_indices[3] -= numElementsAcross*(j%2==1);
            }

            std::vector<Node<2>*> element_nodes;
            for (unsigned k=0; k<6; k++)
            {
               element_nodes.push_back(nodes[node_indices[k]]);
            }
            VertexElement<2,2>* p_element = new VertexElement<2,2>(element_index, element_nodes);
            elements.push_back(p_element);
        }
    }

    // If the mesh has an imposed flat bottom delete unnessesary nodes
    if (isFlatBottom)
    {
        for (unsigned i=0; i<numElementsAcross; i++)
        {
            // Move node 0 to the same position as node 5 then merge the nodes
            nodes[i]->SetPoint(nodes[i+numElementsAcross]->GetPoint());
//                SetNode(i, nodes[i+numElementsAcross]->GetPoint());
//                PerformNodeMerge(mNodes[i], mNodes[i+numElementsAcross]);
        }
    }
    mpMesh = new Cylindrical2dVertexMesh(numElementsAcross, nodes, elements, cellRearrangementThreshold, t2Threshold);
}

MutableVertexMesh<2,2>* CylindricalHoneycombVertexMeshGenerator::GetMesh()
{
    EXCEPTION("A cylindrical mesh was created but a normal mesh is being requested.");
    return mpMesh; // Not really
}

Cylindrical2dVertexMesh* CylindricalHoneycombVertexMeshGenerator::GetCylindricalMesh()
{
    return (Cylindrical2dVertexMesh*) mpMesh;
}
