/*

Copyright (c) 2005-2019, University of Oxford.
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

#include "CylindricalHoneycombVertexMeshGenerator.hpp"

CylindricalHoneycombVertexMeshGenerator::CylindricalHoneycombVertexMeshGenerator(unsigned numElementsAcross,
                                                           unsigned numElementsUp,
                                                           bool isFlatBottom,
                                                           double cellRearrangementThreshold,
                                                           double t2Threshold)
{
    // numElementsAcross must be even for cylindrical meshes
    ///\todo This should be an exception (#2401)
    assert(numElementsAcross > 1);
    assert(numElementsAcross%2 == 0);

    assert(numElementsUp > 0);
    assert(cellRearrangementThreshold > 0.0);
    assert(t2Threshold > 0.0);

    std::vector<Node<2>*> nodes;
    std::vector<VertexElement<2,2>*>  elements;

    unsigned node_index = 0;
    unsigned node_indices[6];
    unsigned element_index;

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
         * On each interior row we have numElementsAcross+1 nodes. All nodes on the first, second,
         * penultimate and last rows are boundary nodes. On other rows no nodes are boundary nodes.
         */
        else
        {
            for (unsigned i=0; i<=numElementsAcross-1; i++)
            {
                double x_coord = ((j%4 == 0)||(j%4 == 3)) ? i+0.5 : i;
                double y_coord = (1.5*j - 0.5*(j%2))*0.5/sqrt(3.0);
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

            node_indices[0] = 2*j*numElementsAcross + i + 1*(j%2==1);
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

    // If imposing a flat bottom, delete unnecessary nodes from the mesh
    if (isFlatBottom)
    {
        for (unsigned i=0; i<numElementsAcross; i++)
        {
            nodes[i]->SetPoint(nodes[i+numElementsAcross]->GetPoint());
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
