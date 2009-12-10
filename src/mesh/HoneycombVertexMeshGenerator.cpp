/*

Copyright (C) University of Oxford, 2005-2009

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
                                                           bool isCylindrical,
                                                           bool isFlatBottom,
                                                           double cellRearrangementThreshold,
                                                           double edgeDivisionThreshold,
                                                           double t2Threshold)
    : mIsCylindrical(isCylindrical)
{
    assert(numElementsAcross > 0);
    assert(numElementsUp > 0);
    assert(cellRearrangementThreshold > 0.0);
    assert(edgeDivisionThreshold > 0.0);
    assert(t2Threshold > 0.0);

    std::vector<Node<2>*> nodes;
    std::vector<VertexElement<2,2>*>  elements;

    unsigned node_index = 0;
    unsigned node_indices[6];
    unsigned element_index;

    if (mIsCylindrical)
    {
        assert(numElementsAcross > 1);
        assert(numElementsAcross%2==0); // numElementsAcross must be even for cylindrical meshes ///\todo why?

        // Create the nodes
        for (unsigned j=0; j<=2*numElementsUp+1; j++)
        {
            if (isFlatBottom && (j==1))
            {
                // Flat bottom to cylindrical mesh
                for (unsigned i=0; i<=numElementsAcross-1; i++)
                {
                    Node<2>* p_node = new Node<2>(node_index, false, i, 0.0);
                    nodes.push_back(p_node);
                    node_index++;
                }
            }
            else
            {
                for (unsigned i=0; i<=numElementsAcross-1; i++)
                {
                    double x_coord = ((j%4 == 0)||(j%4 == 3)) ? i+0.5 : i;
                    double y_coord = (1.5*j - 0.5*(j%2))*0.5/sqrt(3);
    
                    Node<2>* p_node = new Node<2>(node_index, false, x_coord, y_coord);
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
        mpMesh = new Cylindrical2dVertexMesh(numElementsAcross, nodes, elements, cellRearrangementThreshold, edgeDivisionThreshold, t2Threshold);
    }
    else
    {
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
    
        mpMesh = new VertexMesh<2,2>(nodes, elements, cellRearrangementThreshold, edgeDivisionThreshold, t2Threshold);
    }
}


HoneycombVertexMeshGenerator::~HoneycombVertexMeshGenerator()
{
    delete mpMesh;
}


VertexMesh<2,2>* HoneycombVertexMeshGenerator::GetMesh()
{
    if (mIsCylindrical)
    {
        EXCEPTION("A cylindrical mesh was created but a normal mesh is being requested.");
    }
    return mpMesh;
}


Cylindrical2dVertexMesh* HoneycombVertexMeshGenerator::GetCylindricalMesh()
{
    if (!mIsCylindrical)
    {
        EXCEPTION("A normal mesh was created but a cylindrical mesh is being requested.");
    }
    return (Cylindrical2dVertexMesh*) mpMesh;
}
