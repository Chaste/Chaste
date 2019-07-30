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

#ifndef TESTVERTEXELEMENT_HPP_
#define TESTVERTEXELEMENT_HPP_

#include <cxxtest/TestSuite.h>

#include <cmath> //for M_PI

#include "VertexElement.hpp"
#include "Element.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestVertexElement : public CxxTest::TestSuite
{
public:

    void Test1dVertexElementIn2d()
    {
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, false, 1.0, 1.0));

        VertexElement<1,2> element(0, nodes);

        // Test RegisterWithNodes()
        element.RegisterWithNodes();

        for (unsigned node_index=0; node_index<element.GetNumNodes(); node_index++)
        {
            TS_ASSERT_EQUALS(element.GetNode(node_index)->GetNumContainingElements(), 1u);
        }

        TS_ASSERT_EQUALS(element.GetNumNodes(), 2u);
        TS_ASSERT_EQUALS(element.GetNode(0)->GetIndex(), 0u);
        TS_ASSERT_EQUALS(element.GetNode(1)->GetIndex(), 1u);

        // Test GetNodeLocalIndex()
        TS_ASSERT_EQUALS(element.GetNodeLocalIndex(0), 0u);
        TS_ASSERT_EQUALS(element.GetNodeLocalIndex(1), 1u);

        // Test Face methods
        TS_ASSERT_EQUALS(element.GetNumFaces(), 0u);
        VertexElement<0,2>* p_face = element.GetFace(0);
        TS_ASSERT(!p_face);
        TS_ASSERT_EQUALS(element.FaceIsOrientatedClockwise(0), false);

        // Test UpdateNode()
        Node<2>* p_node_2 = new Node<2>(2, false, 1.2, 1.3);
        element.UpdateNode(0, p_node_2);

        TS_ASSERT_DELTA(element.GetNode(0)->rGetLocation()[0], 1.2, 1e-12);
        TS_ASSERT_DELTA(element.GetNode(0)->rGetLocation()[1], 1.3, 1e-12);

        // Test ResetIndex()
        TS_ASSERT_EQUALS(element.GetIndex(), 0u);
        element.ResetIndex(5);
        TS_ASSERT_EQUALS(element.GetIndex(), 5u);

        // Test DeleteNode() and AddNode()
        TS_ASSERT_EQUALS(element.GetNumNodes(), 2u);
        element.DeleteNode(1);
        TS_ASSERT_EQUALS(element.GetNumNodes(), 1u);

        Node<2>* p_node_3 = new Node<2>(3, false, 0.1, 0.4);
        element.AddNode( p_node_3, 0);
        TS_ASSERT_EQUALS(element.GetNumNodes(), 2u);

        TS_ASSERT_DELTA(element.GetNode(0)->rGetLocation()[0], 1.2, 1e-12);
        TS_ASSERT_DELTA(element.GetNode(0)->rGetLocation()[1], 1.3, 1e-12);
        TS_ASSERT_DELTA(element.GetNode(1)->rGetLocation()[0], 0.1, 1e-12);
        TS_ASSERT_DELTA(element.GetNode(1)->rGetLocation()[1], 0.4, 1e-12);

        // Test GetNodeLocalIndex()
        TS_ASSERT_EQUALS(element.GetNodeLocalIndex(0), UINT_MAX);
        TS_ASSERT_EQUALS(element.GetNodeLocalIndex(1), UINT_MAX);

        // Test MarkAsDeleted()
        element.MarkAsDeleted();

        for (unsigned node_index=0; node_index<element.GetNumNodes(); node_index++)
        {
            TS_ASSERT_EQUALS(element.GetNode(node_index)->GetNumContainingElements(), 0u);
        }

        // Tidy up
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
        delete p_node_2;
        delete p_node_3;
    }

    void TestCreateVertexElement()
    {
        // Make 8 nodes to assign to a cube element
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, false, 0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, false, 1.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(2, false, 0.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(3, false, 0.0, 0.0, 1.0));
        nodes.push_back(new Node<3>(4, false, 1.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(5, false, 0.0, 1.0, 1.0));
        nodes.push_back(new Node<3>(6, false, 1.0, 0.0, 1.0));
        nodes.push_back(new Node<3>(7, false, 1.0, 1.0, 1.0));

        std::vector<Node<3>*> nodes_face_0, nodes_face_1, nodes_face_2, nodes_face_3, nodes_face_4, nodes_face_5;

        // Make 6 square faces out of these nodes
        nodes_face_0.push_back(nodes[0]);
        nodes_face_0.push_back(nodes[2]);
        nodes_face_0.push_back(nodes[4]);
        nodes_face_0.push_back(nodes[1]);

        nodes_face_1.push_back(nodes[4]);
        nodes_face_1.push_back(nodes[7]);
        nodes_face_1.push_back(nodes[5]);
        nodes_face_1.push_back(nodes[2]);

        nodes_face_2.push_back(nodes[7]);
        nodes_face_2.push_back(nodes[6]);
        nodes_face_2.push_back(nodes[1]);
        nodes_face_2.push_back(nodes[4]);

        nodes_face_3.push_back(nodes[0]);
        nodes_face_3.push_back(nodes[3]);
        nodes_face_3.push_back(nodes[5]);
        nodes_face_3.push_back(nodes[2]);

        nodes_face_4.push_back(nodes[1]);
        nodes_face_4.push_back(nodes[6]);
        nodes_face_4.push_back(nodes[3]);
        nodes_face_4.push_back(nodes[0]);

        nodes_face_5.push_back(nodes[7]);
        nodes_face_5.push_back(nodes[6]);
        nodes_face_5.push_back(nodes[3]);
        nodes_face_5.push_back(nodes[5]);

        std::vector<VertexElement<2,3>*> faces;
        faces.push_back(new VertexElement<2,3>(0, nodes_face_0));
        faces.push_back(new VertexElement<2,3>(1, nodes_face_1));
        faces.push_back(new VertexElement<2,3>(2, nodes_face_2));
        faces.push_back(new VertexElement<2,3>(3, nodes_face_3));
        faces.push_back(new VertexElement<2,3>(4, nodes_face_4));
        faces.push_back(new VertexElement<2,3>(5, nodes_face_5));

        std::vector<bool> orientations(faces.size());
        for (unsigned i=0; i<faces.size(); i++)
        {
            orientations[i] = true;
        }

        // Make a cube element out of these faces
        VertexElement<3,3> element(0, faces, orientations);

        TS_ASSERT_EQUALS(element.GetNumNodes(), 8u);
        TS_ASSERT_EQUALS(element.GetNumFaces(), 6u);

        TS_ASSERT_EQUALS(element.GetIndex(), 0u);

        // Test the position of some random nodes
        TS_ASSERT_DELTA(element.GetFace(0)->GetNode(0)->rGetLocation()[0], 0.0, 1e-6);
        TS_ASSERT_DELTA(element.GetFace(0)->GetNode(0)->rGetLocation()[1], 0.0, 1e-6);
        TS_ASSERT_DELTA(element.GetFace(0)->GetNode(0)->rGetLocation()[2], 0.0, 1e-6);

        TS_ASSERT_DELTA(element.GetFace(5)->GetNode(2)->rGetLocation()[0], 0.0, 1e-6);
        TS_ASSERT_DELTA(element.GetFace(5)->GetNode(2)->rGetLocation()[1], 0.0, 1e-6);
        TS_ASSERT_DELTA(element.GetFace(5)->GetNode(2)->rGetLocation()[2], 1.0, 1e-6);

        // Test orientations
        for (unsigned face_index=0; face_index<element.GetNumFaces(); face_index++)
        {
            TS_ASSERT_EQUALS(element.FaceIsOrientatedClockwise(face_index), true);
        }

        // Tidy up
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
        for (unsigned i=0; i<faces.size(); i++)
        {
            delete faces[i];
        }
    }

    void TestVertexElementFaceConstructor()
    {
        // Create a regular hexagon
        std::vector<Node<2>*> nodes;
        std::vector<VertexElement<1,2>*> faces;
        std::vector<Node<2>*> face_nodes;
        std::vector<bool> orientations;
        unsigned num_nodes = 6;
        for (unsigned i=0; i<num_nodes; i++)
        {
            double theta = 2.0*M_PI*(double)(i)/(double)(num_nodes);
            nodes.push_back(new Node<2>(i, false, cos(theta), sin(theta)));
        }

        for (unsigned i=0; i<num_nodes-1; i++)
        {
            face_nodes.clear();
            face_nodes.push_back(nodes[i]);
            face_nodes.push_back(nodes[i+1]);

            faces.push_back(new VertexElement<1,2>(i, face_nodes));
            orientations.push_back(true);
        }

        // Create a face with negative orientation
        face_nodes.clear();
        face_nodes.push_back(nodes[0]);
        face_nodes.push_back(nodes[num_nodes-1]);
        faces.push_back(new VertexElement<1,2>(num_nodes-1, face_nodes));
        orientations.push_back(false);

        // Create element
        VertexElement<2,2> vertex_element(0, faces, orientations);

        TS_ASSERT_EQUALS(vertex_element.GetNumNodes(), 6u);
        TS_ASSERT_EQUALS(vertex_element.GetNumFaces(), 6u);

        // Test that each face has the correct orientation and number of nodes
        for (unsigned face_index=0; face_index<6; face_index++)
        {
            TS_ASSERT_EQUALS(vertex_element.GetFace(face_index)->GetNumNodes(), 2u);

            bool is_clockwise = (face_index==5) ? false : true;
            TS_ASSERT_EQUALS(vertex_element.FaceIsOrientatedClockwise(face_index), is_clockwise);
        }

        // Tidy up
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
            delete faces[i];
        }
    }

    void TestAltenativeConstructor()
    {
        // Create element
        VertexElement<2,2> vertex_element(5);

        // Test member variables
        TS_ASSERT_EQUALS(vertex_element.GetIndex(), 5u);
        TS_ASSERT_EQUALS(vertex_element.GetNumNodes(), 0u);
        TS_ASSERT_EQUALS(vertex_element.GetNumFaces(), 0u);
    }

    void TestVertexElementDeleteAndAddNode()
    {
        // Create nodes and Faces
        std::vector<Node<2>*> nodes;
        unsigned num_nodes = 6;
        for (unsigned i=0; i<num_nodes; i++)
        {
            double theta = 2.0*M_PI*(double)(i)/(double)(num_nodes);
            nodes.push_back(new Node<2>(i, false, cos(theta), sin(theta)));
        }

        // Create element
        VertexElement<2,2> vertex_element(0, nodes);

        TS_ASSERT_EQUALS(vertex_element.GetNumNodes(), 6u);

        vertex_element.DeleteNode(3); // Removes (-1,0) node
        vertex_element.DeleteNode(0); // Removes (1,0) node

        // Test node is removed
        TS_ASSERT_EQUALS(vertex_element.GetNumNodes(), 4u);

        // Test other nodes are updated
        TS_ASSERT_DELTA(vertex_element.GetNode(0)->GetPoint()[0], 0.5, 1e-9);
        TS_ASSERT_DELTA(vertex_element.GetNode(0)->GetPoint()[1], 0.5*sqrt(3.0), 1e-9);

        TS_ASSERT_DELTA(vertex_element.GetNode(1)->GetPoint()[0], -0.5, 1e-9);
        TS_ASSERT_DELTA(vertex_element.GetNode(1)->GetPoint()[1], 0.5*sqrt(3.0), 1e-9);

        TS_ASSERT_DELTA(vertex_element.GetNode(2)->GetPoint()[0], -0.5, 1e-9);
        TS_ASSERT_DELTA(vertex_element.GetNode(2)->GetPoint()[1], -0.5*sqrt(3.0), 1e-9);

        TS_ASSERT_DELTA(vertex_element.GetNode(3)->GetPoint()[0], 0.5, 1e-9);
        TS_ASSERT_DELTA(vertex_element.GetNode(3)->GetPoint()[1], -0.5*sqrt(3.0), 1e-9);

        // Add new node
        Node<2>* p_new_node = new Node<2>(4, false, 0.0, 0.0);
        vertex_element.AddNode(p_new_node, 3); // Add node at (0,0) between nodes 3 and 0

        // Test node is added
        TS_ASSERT_EQUALS(vertex_element.GetNumNodes(), 5u);

        // Test other nodes are updated
        TS_ASSERT_DELTA(vertex_element.GetNode(0)->GetPoint()[0], 0.5, 1e-9);
        TS_ASSERT_DELTA(vertex_element.GetNode(0)->GetPoint()[1], 0.5*sqrt(3.0), 1e-9);

        TS_ASSERT_DELTA(vertex_element.GetNode(1)->GetPoint()[0], -0.5, 1e-9);
        TS_ASSERT_DELTA(vertex_element.GetNode(1)->GetPoint()[1], 0.5*sqrt(3.0), 1e-9);

        TS_ASSERT_DELTA(vertex_element.GetNode(2)->GetPoint()[0], -0.5, 1e-9);
        TS_ASSERT_DELTA(vertex_element.GetNode(2)->GetPoint()[1], -0.5*sqrt(3.0), 1e-9);

        TS_ASSERT_DELTA(vertex_element.GetNode(3)->GetPoint()[0], 0.5, 1e-9);
        TS_ASSERT_DELTA(vertex_element.GetNode(3)->GetPoint()[1], -0.5*sqrt(3.0), 1e-9);

        TS_ASSERT_DELTA(vertex_element.GetNode(4)->GetPoint()[0], 0.0, 1e-9);
        TS_ASSERT_DELTA(vertex_element.GetNode(4)->GetPoint()[1], 0.0, 1e-9);

        // Tidy up
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
        delete p_new_node;
    }

    void TestMarkAsDeleted()
    {
        // Create nodes
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
        nodes.push_back(new Node<2>(2, false, 1.0, 1.0));
        nodes.push_back(new Node<2>(3, false, 0.0, 1.0));

        // Create element
        VertexElement<2,2> vertex_element(0, nodes);
        vertex_element.RegisterWithNodes();

        for (unsigned i=0; i<nodes.size(); i++)
        {
            TS_ASSERT_EQUALS(vertex_element.GetNode(i)->GetNumContainingElements(), 1u);
        }

        vertex_element.MarkAsDeleted();

        for (unsigned i=0; i<nodes.size(); i++)
        {
            TS_ASSERT_EQUALS(vertex_element.GetNode(i)->GetNumContainingElements(), 0u);
        }

        // Tidy up
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

    void TestUpdateNode()
    {
        // Create nodes
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
        nodes.push_back(new Node<2>(2, false, 1.0, 1.0));
        nodes.push_back(new Node<2>(3, false, 0.0, 1.0));

        // Create element
        VertexElement<2,2> vertex_element(0, nodes);
        vertex_element.RegisterWithNodes();

        TS_ASSERT_DELTA(vertex_element.GetNode(2)->rGetLocation()[0], 1.0, 1e-12);
        TS_ASSERT_DELTA(vertex_element.GetNode(2)->rGetLocation()[0], 1.0, 1e-12);

        // Update location of node 2
        Node<2>* p_node = new Node<2>(4, false, 1.2, 1.3);
        vertex_element.UpdateNode(2, p_node);

        TS_ASSERT_DELTA(vertex_element.GetNode(2)->rGetLocation()[0], 1.2, 1e-12);
        TS_ASSERT_DELTA(vertex_element.GetNode(2)->rGetLocation()[1], 1.3, 1e-12);

        // Tidy up
        for (unsigned i=0; i<nodes.size(); ++i)
        {
            delete nodes[i];
        }
        delete p_node;
    }

    void TestAddFace()
    {
        // Create nodes
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, true, 0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, true, 2.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(2, true, 2.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(3, true, 0.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(4, true, 0.0, 0.0, 3.0));
        nodes.push_back(new Node<3>(5, true, 2.0, 0.0, 3.0));

        // Create faces
        std::vector<VertexElement<2,3>*> faces;

        std::vector<Node<3>*> nodes_face_0;
        nodes_face_0.push_back(nodes[0]);
        nodes_face_0.push_back(nodes[1]);
        nodes_face_0.push_back(nodes[2]);
        nodes_face_0.push_back(nodes[3]);
        faces.push_back(new VertexElement<2,3>(0, nodes_face_0));

        std::vector<Node<3>*> nodes_face_1;
        nodes_face_1.push_back(nodes[0]);
        nodes_face_1.push_back(nodes[3]);
        nodes_face_1.push_back(nodes[4]);
        faces.push_back(new VertexElement<2,3>(1, nodes_face_1));

        std::vector<Node<3>*> nodes_face_2;
        nodes_face_2.push_back(nodes[1]);
        nodes_face_2.push_back(nodes[2]);
        nodes_face_2.push_back(nodes[5]);
        faces.push_back(new VertexElement<2,3>(2, nodes_face_2));

        std::vector<Node<3>*> nodes_face_3;
        nodes_face_3.push_back(nodes[4]);
        nodes_face_3.push_back(nodes[5]);
        nodes_face_3.push_back(nodes[2]);
        nodes_face_3.push_back(nodes[3]);
        faces.push_back(new VertexElement<2,3>(3, nodes_face_3));

        // Create element
        VertexElement<3,3> element(0);
        TS_ASSERT_EQUALS(element.GetNumNodes(), 0u);
        TS_ASSERT_EQUALS(element.GetNumFaces(), 0u);

        // Add first face to element
        element.AddFace(faces[0]);
        TS_ASSERT_EQUALS(element.GetNumNodes(), 4u);
        TS_ASSERT_EQUALS(element.GetNumFaces(), 1u);

        // Add second face to element
        element.AddFace(faces[1]);
        TS_ASSERT_EQUALS(element.GetNumNodes(), 5u);
        TS_ASSERT_EQUALS(element.GetNumFaces(), 2u);

        // Add third face to element
        element.AddFace(faces[2]);
        TS_ASSERT_EQUALS(element.GetNumNodes(), 6u);
        TS_ASSERT_EQUALS(element.GetNumFaces(), 3u);

        // Add fourth face to element
        element.AddFace(faces[3]);
        TS_ASSERT_EQUALS(element.GetNumNodes(), 6u);
        TS_ASSERT_EQUALS(element.GetNumFaces(), 4u);

        // Tidy up
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
        for (unsigned i=0; i<faces.size(); i++)
        {
            delete faces[i];
        }
    }

    void TestGetNodeLocalIndex()
    {
        // Create nodes
        std::vector<Node<2>*> nodes;

        // This is a square
        nodes.push_back(new Node<2>(3, false, 0.0, 0.0));
        nodes.push_back(new Node<2>(2, false, 1.0, 0.0));
        nodes.push_back(new Node<2>(1, false, 1.0, 1.0));
        nodes.push_back(new Node<2>(0, false, 0.0, 1.0));

        // Create element
        VertexElement<2,2> vertex_element(INDEX_IS_NOT_USED, nodes);

        TS_ASSERT_EQUALS(vertex_element.GetNodeLocalIndex(0), 3u);
        TS_ASSERT_EQUALS(vertex_element.GetNodeLocalIndex(1), 2u);
        TS_ASSERT_EQUALS(vertex_element.GetNodeLocalIndex(2), 1u);
        TS_ASSERT_EQUALS(vertex_element.GetNodeLocalIndex(3), 0u);

        vertex_element.DeleteNode(3); // Removes (1,1) node

        TS_ASSERT_EQUALS(vertex_element.GetNodeLocalIndex(0), UINT_MAX);

        // Tidy up
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

    void TestIsElementOnBoundary()
    {
        // Test with a square element in 2D

        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(3, false, 0.0, 0.0));
        nodes.push_back(new Node<2>(2, false, 1.0, 0.0));
        nodes.push_back(new Node<2>(1, false, 1.0, 1.0));
        nodes.push_back(new Node<2>(0, false, 0.0, 1.0));

        VertexElement<2,2> square_vertex_element(INDEX_IS_NOT_USED, nodes);

        TS_ASSERT_EQUALS(square_vertex_element.IsElementOnBoundary(), false);

        nodes[0]->SetAsBoundaryNode();

        TS_ASSERT_EQUALS(square_vertex_element.IsElementOnBoundary(), true);

        // Tidy up
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
        nodes.clear();

        // Now test with a line element in 2D

        nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, false, 1.0, 0.0));

        VertexElement<1,2> line_vertex_element(INDEX_IS_NOT_USED, nodes);

        TS_ASSERT_EQUALS(line_vertex_element.IsElementOnBoundary(), false);

        nodes[0]->SetAsBoundaryNode();

        TS_ASSERT_EQUALS(line_vertex_element.IsElementOnBoundary(), true);

        // Tidy up
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

};
#endif /*TESTVERTEXELEMENT_HPP_*/
