/*

Copyright (c) 2005-2014, University of Oxford.
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

#include "HexagonalPrism3dVertexMeshGenerator.hpp"
#include "Debug.hpp"

HexagonalPrism3dVertexMeshGenerator::HexagonalPrism3dVertexMeshGenerator(unsigned numElementsInXDirection,
    unsigned numElementsInYDirection,
    double elementSideLength,
    double elementHeight)
{
    assert(numElementsInXDirection > 0);
    assert(numElementsInYDirection > 0);
    assert(elementSideLength > 0.0);
    assert(elementHeight > 0.0);

    // First, create the lower and upper nodes, row by row in increasing x then y
    unsigned lower_node_index = 0;
    unsigned num_lower_nodes = 2*(numElementsInXDirection + numElementsInXDirection*numElementsInYDirection + numElementsInYDirection);

    std::vector<Node<3>*> nodes(2*num_lower_nodes);

    // On each first row we have numElementsInXDirection nodes, all of which are boundary nodes
    for (unsigned i=0; i<numElementsInXDirection; i++)
    {
        double x = i+0.5;
        double y = 0.0;

        // Lower node
        Node<3>* p_lower_node = new Node<3>(lower_node_index, true, x, y, 0.0);
        nodes[lower_node_index] = p_lower_node;

        // Upper node
        unsigned upper_node_index = num_lower_nodes + lower_node_index;
        Node<3>* p_upper_node = new Node<3>(upper_node_index, true, x, y, elementHeight);
        p_upper_node->AddNodeAttribute(1.0);
        nodes[upper_node_index] = p_upper_node;

        lower_node_index++;
    }
    // On each interior row we have numElementsInXDirection+1 nodes
    for (unsigned j=1; j<2*numElementsInYDirection+1; j++)
    {
        for (unsigned i=0; i<=numElementsInXDirection; i++)
        {
            double x = ((j%4 == 0)||(j%4 == 3)) ? i+0.5 : i;
            double y = (1.5*j - 0.5*(j%2))*0.5/sqrt(3.0);

            // Lower node
            Node<3>* p_lower_node = new Node<3>(lower_node_index, true, x, y, 0.0);
            nodes[lower_node_index] = p_lower_node;

            // Upper node
            unsigned upper_node_index = num_lower_nodes + lower_node_index;
            Node<3>* p_upper_node = new Node<3>(upper_node_index, true, x, y, elementHeight);
            p_upper_node->AddNodeAttribute(1.0);
            nodes[upper_node_index] = p_upper_node;

            lower_node_index++;
        }
    }
    // On the last row we have numElementsInXDirection nodes
    if (((2*numElementsInYDirection+1)%4 == 0)||((2*numElementsInYDirection+1)%4 == 3))
    {
        double x = 0.5;
        double y = (1.5*(2*numElementsInYDirection+1) - 0.5*((2*numElementsInYDirection+1)%2))*0.5/sqrt(3.0);

        // Lower node
        Node<3>* p_lower_node = new Node<3>(lower_node_index, true, x, y, 0.0);
        nodes[lower_node_index] = p_lower_node;

        // Upper node
        unsigned upper_node_index = num_lower_nodes + lower_node_index;
        Node<3>* p_upper_node = new Node<3>(upper_node_index, true, x, y, elementHeight);
        p_upper_node->AddNodeAttribute(1.0);
        nodes[upper_node_index] = p_upper_node;

        lower_node_index++;
    }
    for (unsigned i=1; i<numElementsInXDirection; i++)
    {
        double x = (((2*numElementsInYDirection+1)%4 == 0)||((2*numElementsInYDirection+1)%4 == 3)) ? i+0.5 : i;
        double y = (1.5*(2*numElementsInYDirection+1) - 0.5*((2*numElementsInYDirection+1)%2))*0.5/sqrt(3.0);

        // Lower node
        Node<3>* p_lower_node = new Node<3>(lower_node_index, true, x, y, 0.0);
        nodes[lower_node_index] = p_lower_node;

        // Upper node
        unsigned upper_node_index = num_lower_nodes + lower_node_index;
        Node<3>* p_upper_node = new Node<3>(upper_node_index, true, x, y, elementHeight);
        p_upper_node->AddNodeAttribute(1.0);
        nodes[upper_node_index] = p_upper_node;

        lower_node_index++;
    }
    if (((2*numElementsInYDirection+1)%4 == 1)||((2*numElementsInYDirection+1)%4 == 2))
    {
        double x = numElementsInXDirection;
        double y = (1.5*(2*numElementsInYDirection+1) - 0.5*((2*numElementsInYDirection+1)%2))*0.5/sqrt(3.0);

        // Lower node
        Node<3>* p_lower_node = new Node<3>(lower_node_index, true, x, y, 0.0);
        nodes[lower_node_index] = p_lower_node;

        // Upper node
        unsigned upper_node_index = num_lower_nodes + lower_node_index;
        Node<3>* p_upper_node = new Node<3>(upper_node_index, true, x, y, elementHeight);
        p_upper_node->AddNodeAttribute(1.0);
        nodes[upper_node_index] = p_upper_node;

        lower_node_index++; ///\todo remove this line?
    }

    unsigned num_lower_faces = numElementsInXDirection*numElementsInYDirection;
    unsigned num_lateral_faces = 3*numElementsInXDirection*numElementsInYDirection + 2*(numElementsInXDirection + numElementsInYDirection) - 1;

    std::vector<VertexElement<2,3>*> faces(2*num_lower_faces + num_lateral_faces);

    // Create the lower and upper faces of elements (global node indices for each face are given anticlockwise)
    unsigned node_indices_for_face[6];
    for (unsigned j=0; j<numElementsInYDirection; j++)
    {
        for (unsigned i=0; i<numElementsInXDirection; i++)
        {
            if (j==0)
            {
                node_indices_for_face[0] = i;
            }
            else
            {
                node_indices_for_face[0] = 2*j*(numElementsInXDirection+1) - 1*(j%2==0) + i;
            }
            node_indices_for_face[1] = node_indices_for_face[0] + numElementsInXDirection + 1 + 1*(j%2==0 && j>0);
            node_indices_for_face[2] = node_indices_for_face[1] + numElementsInXDirection + 1;
            node_indices_for_face[3] = node_indices_for_face[2] + numElementsInXDirection + 1*(j%2==1 && j<numElementsInYDirection-1);
            node_indices_for_face[4] = node_indices_for_face[2] - 1;
            node_indices_for_face[5] = node_indices_for_face[1] - 1;

            // Create lower face
            std::vector<Node<3>*> nodes_for_lower_face;
            for (unsigned k=0; k<6; k++)
            {
                nodes_for_lower_face.push_back(nodes[node_indices_for_face[k]]);
            }
            unsigned lower_face_index = j*numElementsInXDirection + i;
            VertexElement<2,3>* p_lower_face = new VertexElement<2,3>(lower_face_index, nodes_for_lower_face);
            faces[lower_face_index] = p_lower_face;

            // Create upper face
            std::vector<Node<3>*> nodes_for_upper_face;
            for (unsigned k=0; k<6; k++)
            {
                nodes_for_upper_face.push_back(nodes[node_indices_for_face[k] + num_lower_nodes]);
            }
            unsigned upper_face_index = j*numElementsInXDirection + i + num_lower_faces;
            VertexElement<2,3>* p_upper_face = new VertexElement<2,3>(upper_face_index, nodes_for_upper_face);
            faces[upper_face_index] = p_upper_face;
        }
    }

    // Create the lateral faces of elements (global node indices for each face are given anticlockwise)

    // Create the lateral faces of the 'bottom left' element (in the xy plane)
    unsigned lateral_node_indices_for_element[6];
    lateral_node_indices_for_element[0] = 0;
    lateral_node_indices_for_element[1] = numElementsInXDirection + 1;
    lateral_node_indices_for_element[2] = 2*numElementsInXDirection + 2;
    lateral_node_indices_for_element[3] = 3*numElementsInXDirection + 2;
    lateral_node_indices_for_element[4] = 2*numElementsInXDirection + 1;
    lateral_node_indices_for_element[5] = numElementsInXDirection;

    for (unsigned local_face_index=0; local_face_index<6; local_face_index++)
    {
        unsigned node_indices_for_face[4];
        node_indices_for_face[0] = lateral_node_indices_for_element[local_face_index];
        node_indices_for_face[1] = lateral_node_indices_for_element[(local_face_index+1)%6];
        node_indices_for_face[2] = node_indices_for_face[1] + num_lower_nodes;
        node_indices_for_face[3] = node_indices_for_face[0] + num_lower_nodes;

        std::vector<Node<3>*> nodes_for_lateral_face;
        for (unsigned k=0; k<4; k++)
        {
            nodes_for_lateral_face.push_back(nodes[node_indices_for_face[k]]);
        }

        unsigned lateral_face_index = 2*num_lower_faces + local_face_index;
        VertexElement<2,3>* p_lateral_face = new VertexElement<2,3>(lateral_face_index, nodes_for_lateral_face);
        faces[lateral_face_index] = p_lateral_face;
    }

    // Create the lateral faces of the other elements in this row (in the x direction)
    for (unsigned local_elem_index=1; local_elem_index<numElementsInXDirection; local_elem_index++)
    {
        unsigned lateral_node_indices_for_element[6];
        lateral_node_indices_for_element[0] = numElementsInXDirection + local_elem_index;
        lateral_node_indices_for_element[1] = local_elem_index;
        lateral_node_indices_for_element[2] = numElementsInXDirection + 1 + local_elem_index;
        lateral_node_indices_for_element[3] = 2*numElementsInXDirection + 2 + local_elem_index;
        lateral_node_indices_for_element[4] = 3*numElementsInXDirection + 2 + local_elem_index;
        lateral_node_indices_for_element[5] = 2*numElementsInXDirection + 1 + local_elem_index;

        for (unsigned local_face_index=0; local_face_index<5; local_face_index++)
        {
            unsigned node_indices_for_face[4];
            node_indices_for_face[0] = lateral_node_indices_for_element[local_face_index];
            node_indices_for_face[1] = lateral_node_indices_for_element[local_face_index+1];
            node_indices_for_face[2] = node_indices_for_face[1] + num_lower_nodes;
            node_indices_for_face[3] = node_indices_for_face[0] + num_lower_nodes;

            std::vector<Node<3>*> nodes_for_lateral_face;
            for (unsigned k=0; k<4; k++)
            {
                nodes_for_lateral_face.push_back(nodes[node_indices_for_face[k]]);
            }

            unsigned lateral_face_index = 2*num_lower_faces + 6 + 5*(local_elem_index-1) + local_face_index;
            VertexElement<2,3>* p_lateral_face = new VertexElement<2,3>(lateral_face_index, nodes_for_lateral_face);
            faces[lateral_face_index] = p_lateral_face;
        }
    }

    // Create all other lateral faces for elements
    for (unsigned j=1; j<numElementsInYDirection; j++)
    {
        unsigned smallest_lateral_face_index_this_row = 2*num_lower_faces + 2*numElementsInXDirection - 1 + j*(3*numElementsInXDirection + 2);

        if (j%2 == 1)
        {
            // Create the lateral faces of the left-most element in this row (in the x direction)
            unsigned lateral_node_indices_for_leftmost_element[5];
            lateral_node_indices_for_leftmost_element[0] = (numElementsInXDirection + 1)*(1 + 2*j);
            lateral_node_indices_for_leftmost_element[1] = lateral_node_indices_for_leftmost_element[0] + numElementsInXDirection + 1;
            lateral_node_indices_for_leftmost_element[2] = lateral_node_indices_for_leftmost_element[1] + numElementsInXDirection + 1*(j < numElementsInYDirection-1);
            lateral_node_indices_for_leftmost_element[3] = lateral_node_indices_for_leftmost_element[1] - 1;
            lateral_node_indices_for_leftmost_element[4] = lateral_node_indices_for_leftmost_element[0] - 1;

            for (unsigned local_face_index=0; local_face_index<4; local_face_index++)
            {
                unsigned node_indices_for_face[4];
                node_indices_for_face[0] = lateral_node_indices_for_leftmost_element[local_face_index];
                node_indices_for_face[1] = lateral_node_indices_for_leftmost_element[local_face_index+1];
                node_indices_for_face[2] = node_indices_for_face[1] + num_lower_nodes;
                node_indices_for_face[3] = node_indices_for_face[0] + num_lower_nodes;

                std::vector<Node<3>*> nodes_for_lateral_face;
                for (unsigned k=0; k<4; k++)
                {
                    nodes_for_lateral_face.push_back(nodes[node_indices_for_face[k]]);
                }

                unsigned lateral_face_index = smallest_lateral_face_index_this_row + local_face_index;
                VertexElement<2,3>* p_lateral_face = new VertexElement<2,3>(lateral_face_index, nodes_for_lateral_face);
                faces[lateral_face_index] = p_lateral_face;
            }

            // Create the lateral faces of the intermediate elements in this row (in the x direction)
            for (unsigned i=1; i<numElementsInXDirection-1; i++)
            {
                unsigned lateral_node_indices_for_element[4];
                lateral_node_indices_for_element[0] = (numElementsInXDirection + 1)*(1 + 2*j) + i;
                lateral_node_indices_for_element[1] = lateral_node_indices_for_element[0] + numElementsInXDirection + 1;
                lateral_node_indices_for_element[2] = lateral_node_indices_for_element[1] + numElementsInXDirection + 1*(j < numElementsInYDirection-1);
                lateral_node_indices_for_element[3] = lateral_node_indices_for_element[1] - 1;

                for (unsigned local_face_index=0; local_face_index<3; local_face_index++)
                {
                    unsigned node_indices_for_face[4];
                    node_indices_for_face[0] = lateral_node_indices_for_element[local_face_index];
                    node_indices_for_face[1] = lateral_node_indices_for_element[local_face_index+1];
                    node_indices_for_face[2] = node_indices_for_face[1] + num_lower_nodes;
                    node_indices_for_face[3] = node_indices_for_face[0] + num_lower_nodes;

                    std::vector<Node<3>*> nodes_for_lateral_face;
                    for (unsigned k=0; k<4; k++)
                    {
                        nodes_for_lateral_face.push_back(nodes[node_indices_for_face[k]]);
                    }

                    unsigned lateral_face_index = smallest_lateral_face_index_this_row + 1 + 3*i + local_face_index;
                    VertexElement<2,3>* p_lateral_face = new VertexElement<2,3>(lateral_face_index, nodes_for_lateral_face);
                    faces[lateral_face_index] = p_lateral_face;
                }
            }

            // Create the lateral faces of the right-most element in this row (in the x direction)
            unsigned lateral_node_indices_for_rightmost_element[5];
            lateral_node_indices_for_rightmost_element[0] = (numElementsInXDirection + 1)*(1 + 2*j) - 2;
            lateral_node_indices_for_rightmost_element[1] = lateral_node_indices_for_rightmost_element[0] + numElementsInXDirection + 1;
            lateral_node_indices_for_rightmost_element[2] = lateral_node_indices_for_rightmost_element[1] + numElementsInXDirection + 1;
            lateral_node_indices_for_rightmost_element[3] = lateral_node_indices_for_rightmost_element[2] + numElementsInXDirection + 1*(j < numElementsInYDirection-1);
            lateral_node_indices_for_rightmost_element[4] = lateral_node_indices_for_rightmost_element[2] - 1;

            for (unsigned local_face_index=0; local_face_index<4; local_face_index++)
            {
                unsigned node_indices_for_face[4];
                node_indices_for_face[0] = lateral_node_indices_for_rightmost_element[local_face_index];
                node_indices_for_face[1] = lateral_node_indices_for_rightmost_element[local_face_index+1];
                node_indices_for_face[2] = node_indices_for_face[1] + num_lower_nodes;
                node_indices_for_face[3] = node_indices_for_face[0] + num_lower_nodes;

                std::vector<Node<3>*> nodes_for_lateral_face;
                for (unsigned k=0; k<4; k++)
                {
                    nodes_for_lateral_face.push_back(nodes[node_indices_for_face[k]]);
                }

                unsigned lateral_face_index = smallest_lateral_face_index_this_row + 1 + 3*(numElementsInXDirection-1) + local_face_index;
                VertexElement<2,3>* p_lateral_face = new VertexElement<2,3>(lateral_face_index, nodes_for_lateral_face);
                faces[lateral_face_index] = p_lateral_face;
            }
        }
        else // j%2 == 0
        {
            // Create the lateral faces of the left-most element in this row (in the x direction)
            unsigned lateral_node_indices_for_leftmost_element[6];
            lateral_node_indices_for_leftmost_element[0] = (numElementsInXDirection + 1)*(2*j + 1);
            lateral_node_indices_for_leftmost_element[1] = lateral_node_indices_for_leftmost_element[0] + numElementsInXDirection + 1;
            lateral_node_indices_for_leftmost_element[2] = lateral_node_indices_for_leftmost_element[1] + numElementsInXDirection;
            lateral_node_indices_for_leftmost_element[3] = lateral_node_indices_for_leftmost_element[1] - 1;
            lateral_node_indices_for_leftmost_element[4] = lateral_node_indices_for_leftmost_element[0] - 1;
            lateral_node_indices_for_leftmost_element[5] = lateral_node_indices_for_leftmost_element[4] - numElementsInXDirection - 1;

            for (unsigned local_face_index=0; local_face_index<5; local_face_index++)
            {
                unsigned node_indices_for_face[4];
                node_indices_for_face[0] = lateral_node_indices_for_leftmost_element[local_face_index];
                node_indices_for_face[1] = lateral_node_indices_for_leftmost_element[local_face_index+1];
                node_indices_for_face[2] = node_indices_for_face[1] + num_lower_nodes;
                node_indices_for_face[3] = node_indices_for_face[0] + num_lower_nodes;

                std::vector<Node<3>*> nodes_for_lateral_face;
                for (unsigned k=0; k<4; k++)
                {
                    nodes_for_lateral_face.push_back(nodes[node_indices_for_face[k]]);
                }

                unsigned lateral_face_index = smallest_lateral_face_index_this_row + local_face_index;
                VertexElement<2,3>* p_lateral_face = new VertexElement<2,3>(lateral_face_index, nodes_for_lateral_face);
                faces[lateral_face_index] = p_lateral_face;
            }

            // Create the lateral faces of the other elements in this row (in the x direction)
            for (unsigned i=1; i<numElementsInXDirection; i++)
            {
                unsigned lateral_node_indices_for_element[4];
                lateral_node_indices_for_element[0] = (numElementsInXDirection + 1)*(1 + 2*j) + i;
                lateral_node_indices_for_element[1] = lateral_node_indices_for_element[0] + numElementsInXDirection + 1;
                lateral_node_indices_for_element[2] = lateral_node_indices_for_element[1] + numElementsInXDirection;
                lateral_node_indices_for_element[3] = lateral_node_indices_for_element[1] - 1;

                for (unsigned local_face_index=0; local_face_index<3; local_face_index++)
                {
                    unsigned node_indices_for_face[4];
                    node_indices_for_face[0] = lateral_node_indices_for_element[local_face_index];
                    node_indices_for_face[1] = lateral_node_indices_for_element[local_face_index+1];
                    node_indices_for_face[2] = node_indices_for_face[1] + num_lower_nodes;
                    node_indices_for_face[3] = node_indices_for_face[0] + num_lower_nodes;

                    std::vector<Node<3>*> nodes_for_lateral_face;
                    for (unsigned k=0; k<4; k++)
                    {
                        nodes_for_lateral_face.push_back(nodes[node_indices_for_face[k]]);
                    }

                    unsigned lateral_face_index = smallest_lateral_face_index_this_row + 2 + 3*i + local_face_index;
                    VertexElement<2,3>* p_lateral_face = new VertexElement<2,3>(lateral_face_index, nodes_for_lateral_face);
                    faces[lateral_face_index] = p_lateral_face;
                }
            }
        }
    }

    // Create vector of face orientations
    ///\todo think carefully about whether all faces are oriented anticlockwise
    std::vector<bool> face_orientations;
    for (unsigned i=0; i<faces.size(); i++)
    {
        face_orientations.push_back(true);
    }

    std::vector<VertexElement<3,3>*> elements;

    // Store the indices of the faces that will form the 'bottom left' element (in the xy plane)
    unsigned face_indices[8];
    face_indices[0] = 0;
    face_indices[1] = num_lower_faces;
    face_indices[2] = 2*num_lower_faces;
    face_indices[3] = 2*num_lower_faces + 1;
    face_indices[4] = 2*num_lower_faces + 2;
    face_indices[5] = 2*num_lower_faces + 3;
    face_indices[6] = 2*num_lower_faces + 4;
    face_indices[7] = 2*num_lower_faces + 5;

    // Create and populate vectors of the faces and face orientations for this element
    std::vector<VertexElement<2,3>*> faces_for_element;
    std::vector<bool> orientations_for_element;
    for (unsigned k=0; k<8; k++)
    {
        faces_for_element.push_back(faces[face_indices[k]]);
        orientations_for_element.push_back(face_orientations[face_indices[k]]);
    }

    elements.push_back(new VertexElement<3,3>(0, faces_for_element, orientations_for_element));

    // Create the other elements in this row (in the x direction)
    for (unsigned i=1; i<numElementsInXDirection; i++)
    {
        // Store the indices of the faces that will form this element
        unsigned face_indices[8];
        face_indices[0] = i;
        face_indices[1] = i + num_lower_faces;
        face_indices[2] = 2*num_lower_faces + 5*i + 1;
        face_indices[3] = 2*num_lower_faces + 5*i + 2;
        face_indices[4] = 2*num_lower_faces + 5*i + 3;
        face_indices[5] = 2*num_lower_faces + 5*i + 4;
        face_indices[6] = 2*num_lower_faces + 5*i + 5;
        face_indices[7] = 2*num_lower_faces + 5*i - 2 - 2*(i==1);

        // Create and populate vectors of the faces and face orientations for this element
        std::vector<VertexElement<2,3>*> faces_for_element;
        std::vector<bool> orientations_for_element;
        for (unsigned k=0; k<8; k++)
        {
            faces_for_element.push_back(faces[face_indices[k]]);
            orientations_for_element.push_back(face_orientations[face_indices[k]]);
        }

        elements.push_back(new VertexElement<3,3>(i, faces_for_element, orientations_for_element));
    }

    // Create all other elements
    if (numElementsInYDirection != 1)
    {
        for (unsigned j=1; j<numElementsInYDirection; j++)
        {
            for (unsigned i=0; i<numElementsInXDirection; i++)
            {
                // Store the indices of the faces that will form this element
                unsigned face_indices[8];
                face_indices[0] = i + j*numElementsInXDirection;
                face_indices[1] = i + j*numElementsInXDirection + num_lower_faces;

                if (j%2 == 0)
                {
                    face_indices[2] = 2*numElementsInXDirection*numElementsInYDirection + 5*numElementsInXDirection + 1 + (j-2)*(3*numElementsInXDirection + 2) + (3*i + 1)*(i > 0) + 1*(i == numElementsInXDirection-1) + 2;
                }
                else
                {
                    if (j == 1)
                    {
                        if (i == numElementsInXDirection - 1)
                        {
                            face_indices[2] = 2*numElementsInXDirection*numElementsInYDirection + 5*numElementsInXDirection + 1 + (j-1)*(3*numElementsInXDirection + 2) + (3*i - 1) + 2;
                        }
                        else
                        {
                            face_indices[2] = 2*numElementsInXDirection*numElementsInYDirection + 5*(i+1) + 4 + 1;
                        }
                    }
                    else
                    {
                        if (i == numElementsInXDirection - 1)
                        {
                            face_indices[2] = 2*numElementsInXDirection*numElementsInYDirection + 5*numElementsInXDirection + 1 + (j-1)*(3*numElementsInXDirection + 2) + (3*(i-1) + 1)*(i > 0) + 1*(i == numElementsInXDirection) + 2 + 1;
                        }
                        else
                        {
                            face_indices[2] = 2*numElementsInXDirection*numElementsInYDirection + 5*numElementsInXDirection + 1 + (j-2)*(3*numElementsInXDirection + 2) + (3*i + 5) + 2;
                        }
                    }
                }

                /**
                 * If j=0, then f_3(i,j) = 2*n_x*n_y + 1 + (5*i + 2)*(i > 0).
                 * If j>0 is even, then f_3(i,j) = 2*n_x*n_y + 5*n_x + 1 + (j-1)*(3*n_x + 2) + (3*i + 2)*(i > 0).
                 * If j>0 is odd, then f_3(i,j) = 2*n_x*n_y + 5*n_x + 1 + (j-1)*(3*n_x + 2) + (3*i + 1)*(i > 0) + 1*(i == n_x - 1).
                 */
                face_indices[3] = 2*numElementsInXDirection*numElementsInYDirection + 5*numElementsInXDirection + 1 + (j-1)*(3*numElementsInXDirection + 2);
                if (j%2 == 0)
                {
                    face_indices[3] += (3*i + 2)*(i > 0);
                }
                else
                {
                    face_indices[3] += (3*i + 1)*(i > 0) + 1*(i == numElementsInXDirection-1);
                }

                face_indices[4] = face_indices[3] + 1;
                face_indices[5] = face_indices[4] + 1;

                /**
                 * If i=0, then f_6(i,j) = f_5(i,j) + 1.
                 * If i>0, then f_6(i,j) = f_3(i-1,j).
                 */
                if (i == 0)
                {
                    face_indices[6] = face_indices[5] + 1;
                }
                else
                {
                    face_indices[6] = 2*numElementsInXDirection*numElementsInYDirection + 5*numElementsInXDirection + 1 + (j-1)*(3*numElementsInXDirection + 2);
                    if (j%2 == 0)
                    {
                        face_indices[6] += (3*i -1)*(i > 1);
                    }
                    else
                    {
                        face_indices[6] += (3*i -2)*(i > 1);
                    }
                }

                /**
                 * Here we have j>0.
                 * If i=0, then f_7(i,j) = f_6(i,j) + 1.
                 * If i>0, then f_7(i,j) = f_4(i,j-1).
                 */
                if (j == 1)
                {
                    face_indices[7] = 2*numElementsInXDirection*numElementsInYDirection + 2 + 7*(i > 0) + 5*(i-1)*(i > 1);
                }
                else if (j%2 == 0)
                {
                    if (i == 0)
                    {
                        face_indices[7] = face_indices[6] + 1;
                    }
                    else
                    {
                        face_indices[7] = 1 + 2*numElementsInXDirection*numElementsInYDirection + 5*numElementsInXDirection + 1 + (j-2)*(3*numElementsInXDirection + 2) + (3*i -2)*(i > 1);
                    }
                }
                else
                {
                    face_indices[7] = 1 + 2*numElementsInXDirection*numElementsInYDirection + 5*numElementsInXDirection + 1 + (j-2)*(3*numElementsInXDirection + 2) + (3*i + 2)*(i > 0);
                }

                // Create and populate vectors of the faces and face orientations for this element
                std::vector<VertexElement<2,3>*> faces_for_element;
                std::vector<bool> orientations_for_element;
                for (unsigned k=0; k<8; k++)
                {
                    faces_for_element.push_back(faces[face_indices[k]]);
                    orientations_for_element.push_back(face_orientations[face_indices[k]]);
                }

                unsigned elem_index = i + j*numElementsInXDirection;
                elements.push_back(new VertexElement<3,3>(elem_index, faces_for_element, orientations_for_element));
            }
        }
    }

    /*
     * Finally, reorder the local indices of the nodes within each element,
     * so that the local indices 0 to 5 correspond to the lower nodes in
     * the element ordered anticlockwise and the local indices 6 to 11
     * correspond to the upper nodes in the element ordered anticlockwise.
     */
//    for (unsigned elem_index=0; elem_index<elements.size(); elem_index++)
//    { //commented to try the order before reordering, because I suspect that the set sorting might have problem
//        VertexElement<3,3>* p_element = elements[elem_index];
//
//        std::vector<Node<3>*> temp_nodes;
//        for (unsigned temp_index=0; temp_index<p_element->GetNumNodes(); temp_index++)
//        {
//            temp_nodes.push_back(p_element->GetNode(temp_index));
//        }
//
//        // Node 0 already has the correct local index
//        p_element->UpdateNode(1, temp_nodes[4]);
//        p_element->UpdateNode(2, temp_nodes[8]);
//        p_element->UpdateNode(3, temp_nodes[10]);
//        p_element->UpdateNode(4, temp_nodes[6]);
//        p_element->UpdateNode(5, temp_nodes[2]);
//        p_element->UpdateNode(6, temp_nodes[1]);
//        p_element->UpdateNode(7, temp_nodes[5]);
//        p_element->UpdateNode(8, temp_nodes[9]);
//        p_element->UpdateNode(9, temp_nodes[11]);
//        p_element->UpdateNode(10, temp_nodes[7]);
//        p_element->UpdateNode(11, temp_nodes[3]);
//    }

//for (unsigned l=0; l<elements.size(); l++)
//{
//    PRINT_VARIABLE(elements[l]->GetIndex());
//    PRINT_VARIABLE(elements[l]->GetNumFaces());
//    for (unsigned m=0; m<elements[l]->GetNumFaces(); m++)
//    {
//        PRINT_VARIABLE(elements[l]->GetFace(m)->GetIndex());
//    }
//}

    mpMesh = new MutableVertexMesh<3,3>(nodes, elements);//cellRearrangementThreshold, t2Threshold); ///\todo

    // Scale the mesh so that each element's area takes the value elementArea
    mpMesh->Scale(sqrt(3.0)*elementSideLength, sqrt(3.0)*elementSideLength);
}

HexagonalPrism3dVertexMeshGenerator::~HexagonalPrism3dVertexMeshGenerator()
{
    delete mpMesh;
}

MutableVertexMesh<3,3>* HexagonalPrism3dVertexMeshGenerator::GetMesh()
{
    return mpMesh;
}
