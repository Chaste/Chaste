/*

Copyright (c) 2005-2016, University of Oxford.
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

#ifndef TESTMONOLAYERVERTEXMESHGENERATOR_HPP_
#define TESTMONOLAYERVERTEXMESHGENERATOR_HPP_

#include <cxxtest/TestSuite.h>

#include "MonolayerVertexMeshGenerator.hpp"
#include "MonolayerVertexMeshCustomFunctions.hpp"

#include "VoronoiVertexMeshGenerator.hpp"       // Make Mesh from 2D Mesh
#include "HoneycombVertexMeshGenerator.hpp"     // Make Mesh from 2D Mesh

#include <algorithm>
#include "Debug.hpp"

#include "FakePetscSetup.hpp"


#include "GeodesicSphere23Generator.hpp"


class TestMonolayerVertexMeshGenerator : public CxxTest::TestSuite
{
public:

    void TestGeodesicShpere()
    {
        GeodesicSphere23Generator builder;
        {
            MutableVertexMesh<2, 3>* p_mesh = new MutableVertexMesh<2, 3>(builder.mNodes, builder.mFaces);
            VertexMeshWriter<2, 3> Writer("SphericalMesh", "Geodesic_0Division", false);
            Writer.WriteVtkUsingMeshWithCellId(*p_mesh);

//            MonolayerVertexMeshGenerator sBuilder("trii");
//            MutableVertexMesh<3,3>* pp_mesh33 = sBuilder.MakeSphericalMesh33(p_mesh, 5, 0.5);
//            sBuilder.WriteVtk("SphericalMesh", "0");

//            delete p_mesh;
        }
        PRINT_3_VARIABLES(builder.mNodes.size(), builder.mEdges.size(), builder.mFaces.size())

        for (unsigned i=0; i<builder.mNodes.size(); ++i)
        {
            TS_ASSERT_EQUALS(builder.mNodes[i]->GetNumContainingFaces(), 5);
        }
        for (unsigned i=0; i<builder.mEdges.size(); ++i)
        {
            TS_ASSERT_EQUALS(builder.mEdges[i]->FaceGetNumContainingElements(), 2);
            TS_ASSERT_EQUALS(builder.mEdges[i]->GetNumNodes(), 2);
        }
        for (unsigned i=0; i<builder.mFaces.size(); ++i)
        {
            TS_ASSERT_EQUALS(builder.mFaces[i]->GetNumNodes(), 3);
            TS_ASSERT_EQUALS(builder.mFaces[i]->GetNumFaces(), 3);
        }

        {
            MutableVertexMesh<2,3>* pp_mesh = builder.GetDual();
            VertexMeshWriter<2, 3> Writer("SphericalMesh", "Geodesic_0Dual", false);
            Writer.WriteVtkUsingMeshWithCellId(*pp_mesh);

            MonolayerVertexMeshGenerator sBuilder("ssss");
            MutableVertexMesh<3,3>* pp_mesh33 = sBuilder.MakeSphericalMesh33(pp_mesh, 5, 0.5);
            sBuilder.WriteVtk("SphericalMesh", "0");

//            delete pp_mesh;
//            delete pp_mesh33;
        }




        builder.SubDivide();
        PRINT_3_VARIABLES(builder.mNodes.size(), builder.mEdges.size(), builder.mFaces.size())
        {
            MutableVertexMesh<2, 3>* p_mesh = new MutableVertexMesh<2, 3>(builder.mNodes, builder.mFaces);
            VertexMeshWriter<2, 3> Writer("SphericalMesh", "Geodesic_1Divison", false);
            Writer.WriteVtkUsingMeshWithCellId(*p_mesh);

//            MonolayerVertexMeshGenerator sBuilder("trii", false);
//            MutableVertexMesh<3,3>* pp_mesh33 = sBuilder.MakeSphericalMesh33(p_mesh, 5, 0.5);
//            sBuilder.WriteVtk("SphericalMesh", "1");

//            delete p_mesh;
        }

        unsigned iii = 0;
        for (unsigned i=0; i<builder.mNodes.size(); ++i)
        {
            TS_ASSERT_DELTA(builder.mNodes[i]->GetNumContainingElements(), 5.5, 0.51);
            TS_ASSERT_DELTA(builder.mNodes[i]->GetNumContainingFaces(), 5.5, 0.51);
            TS_ASSERT_EQUALS(builder.mNodes[i]->GetNumContainingFaces(), builder.mNodes[i]->GetNumContainingElements());
            if (builder.mNodes[i]->GetNumContainingElements() == 5)
            {
                ++iii;
            }
        }
        TS_ASSERT_EQUALS(iii, 12);
        for (unsigned i=30; i<builder.mEdges.size(); ++i)
        {
            if (builder.mEdges[i]==NULL)
                continue;
            TS_ASSERT_EQUALS(builder.mEdges[i]->FaceGetNumContainingElements(), 2);
            TS_ASSERT_EQUALS(builder.mEdges[i]->GetNumNodes(), 2);
        }
        for (unsigned i=0; i<builder.mFaces.size(); ++i)
        {
            TS_ASSERT_EQUALS(builder.mFaces[i]->GetNumNodes(), 3);
            TS_ASSERT_EQUALS(builder.mFaces[i]->GetNumFaces(), 3);
        }
        {
            MutableVertexMesh<2,3>* pp_mesh = builder.GetDual();
            VertexMeshWriter<2, 3> Writer("SphericalMesh", "Geodesic_1Dual", false);
            Writer.WriteVtkUsingMeshWithCellId(*pp_mesh);

            MonolayerVertexMeshGenerator sBuilder("ssss");
            MutableVertexMesh<3,3>* pp_mesh33 = sBuilder.MakeSphericalMesh33(pp_mesh, 5, 0.5);
            sBuilder.WriteVtk("SphericalMesh", "1");

//            delete pp_mesh;
//            delete pp_mesh33;
        }




        builder.SubDivide();
        PRINT_3_VARIABLES(builder.mNodes.size(), builder.mEdges.size(), builder.mFaces.size())
        {
            MutableVertexMesh<2, 3>* p_mesh = new MutableVertexMesh<2, 3>(builder.mNodes, builder.mFaces);
            VertexMeshWriter<2, 3> Writer("SphericalMesh", "Geodesic_2Division", false);
            Writer.WriteVtkUsingMeshWithCellId(*p_mesh);

//            delete p_mesh;
        }

        iii = 0;
        for (unsigned i=0; i<builder.mNodes.size(); ++i)
        {
            TS_ASSERT_DELTA(builder.mNodes[i]->GetNumContainingElements(), 5.5, 0.51);
            TS_ASSERT_DELTA(builder.mNodes[i]->GetNumContainingFaces(), 5.5, 0.51);
            TS_ASSERT_EQUALS(builder.mNodes[i]->GetNumContainingFaces(), builder.mNodes[i]->GetNumContainingElements());
            if (builder.mNodes[i]->GetNumContainingElements() == 5)
            {
                ++iii;
            }
        }
        TS_ASSERT_EQUALS(iii, 12);
        for (unsigned i=30; i<builder.mEdges.size(); ++i)
        {
            if (builder.mEdges[i]==NULL)
                continue;
            TS_ASSERT_EQUALS(builder.mEdges[i]->FaceGetNumContainingElements(), 2);
            TS_ASSERT_EQUALS(builder.mEdges[i]->GetNumNodes(), 2);
        }
        for (unsigned i=0; i<builder.mFaces.size(); ++i)
        {
            TS_ASSERT_EQUALS(builder.mFaces[i]->GetNumNodes(), 3);
            TS_ASSERT_EQUALS(builder.mFaces[i]->GetNumFaces(), 3);
        }
        {
            MutableVertexMesh<2,3>* pp_mesh = builder.GetDual();
            VertexMeshWriter<2, 3> Writer("SphericalMesh", "Geodesic_2Dual", false);
            Writer.WriteVtkUsingMeshWithCellId(*pp_mesh);

            MonolayerVertexMeshGenerator sBuilder("ssss");
            MutableVertexMesh<3,3>* pp_mesh33 = sBuilder.MakeSphericalMesh33(pp_mesh, 5, 0.5);
            sBuilder.WriteVtk("SphericalMesh", "2");

//            delete pp_mesh;
//            delete pp_mesh33;
        }

MARK


        builder.SubDivide();
        PRINT_3_VARIABLES(builder.mNodes.size(), builder.mEdges.size(), builder.mFaces.size())
        {
            MutableVertexMesh<2, 3>* p_mesh = new MutableVertexMesh<2, 3>(builder.mNodes, builder.mFaces);
            VertexMeshWriter<2, 3> Writer("SphericalMesh", "Geodesic_3Divison", false);
            Writer.WriteVtkUsingMeshWithCellId(*p_mesh);

//            delete p_mesh;
        }

        iii = 0;
        for (unsigned i=0; i<builder.mNodes.size(); ++i)
        {
            TS_ASSERT_DELTA(builder.mNodes[i]->GetNumContainingElements(), 5.5, 0.51);
            TS_ASSERT_DELTA(builder.mNodes[i]->GetNumContainingFaces(), 5.5, 0.51);
            TS_ASSERT_EQUALS(builder.mNodes[i]->GetNumContainingFaces(), builder.mNodes[i]->GetNumContainingElements());
            if (builder.mNodes[i]->GetNumContainingElements() == 5)
            {
                ++iii;
            }
        }
        TS_ASSERT_EQUALS(iii, 12);
        for (unsigned i=30; i<builder.mEdges.size(); ++i)
        {
            if (builder.mEdges[i]==NULL)
                continue;
            TS_ASSERT_EQUALS(builder.mEdges[i]->FaceGetNumContainingElements(), 2);
            TS_ASSERT_EQUALS(builder.mEdges[i]->GetNumNodes(), 2);
        }
        for (unsigned i=0; i<builder.mFaces.size(); ++i)
        {
            TS_ASSERT_EQUALS(builder.mFaces[i]->GetNumNodes(), 3);
            TS_ASSERT_EQUALS(builder.mFaces[i]->GetNumFaces(), 3);
        }
        {
            MutableVertexMesh<2,3>* pp_mesh = builder.GetDual();
            VertexMeshWriter<2, 3> Writer("SphericalMesh", "Geodesic_3Dual", false);
            Writer.WriteVtkUsingMeshWithCellId(*pp_mesh);

            MonolayerVertexMeshGenerator sBuilder("ssss");
            MutableVertexMesh<3,3>* pp_mesh33 = sBuilder.MakeSphericalMesh33(pp_mesh, 5, 0.5);
            sBuilder.WriteVtk("SphericalMesh", "3");

//            delete pp_mesh;
//            delete pp_mesh33;
        }
    }

    void TestMakeNormalMesh()
    {
        /** Create a mesh comprising six nodes contained in two triangle and two rhomboid elements, as shown below.
        * We will test that that a T1 swap of the two central nodes is correctly implemented.
        *  _____
        * |\   /|
        * | \ / |
        * |  |  |
        * | / \ |
        * |/___\|
        */
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, true,  0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, true,  1.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(2, true,  1.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(3, true,  0.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(4, false, 0.5, 0.4, 0.0));
        nodes.push_back(new Node<3>(5, false, 0.5, 0.6, 0.0));

        unsigned node_indices_elem_0[3] = {2, 3, 5};
        unsigned node_indices_elem_1[4] = {4, 1, 2, 5};
        unsigned node_indices_elem_2[3] = {0, 1, 4};
        unsigned node_indices_elem_3[4] = {4, 5, 3, 0};

        MonolayerVertexMeshGenerator builder(nodes, "T1SwapWith4Elements");
        builder.BuildElementWith(3, node_indices_elem_0);
        builder.BuildElementWith(4, node_indices_elem_1);
        builder.BuildElementWith(3, node_indices_elem_2);
        builder.BuildElementWith(4, node_indices_elem_3);
        // A reference variable as mesh is noncopyable
        MutableVertexMesh<3, 3>* p_mesh = builder.GenerateMesh();
        builder.WriteVtkWithSubfolder("TestMonolayerGenerator");
        builder.WriteVtk("TestMonolayerGenerator");
        builder.PrintMesh();

        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 12u);
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 4u);
        TS_ASSERT_EQUALS(p_mesh->GetNumFaces(), 17u);

        const unsigned num_lateral_faces = p_mesh->GetNumFaces() - 2*p_mesh->GetNumElements();

        TS_ASSERT_EQUALS(p_mesh->GetNumNodes()/2+p_mesh->GetNumElements()+1-num_lateral_faces, 2u);

        for (unsigned i=0; i<p_mesh->GetNumElements(); ++i)
        {
            VertexElement<3, 3>* p_elem = p_mesh->GetElement(i);
            std::vector<unsigned> node_indices;
            for (unsigned j=0; j<p_elem->GetNumNodes()/2; ++j)
            {
                node_indices.push_back(p_elem->GetNodeGlobalIndex(j));
            }
            PRINT_CONTAINER(node_indices)
//            TS_ASSERT(std::equal(node_indices.begin(), node_indices.begin(), node_indices_elem_1));
        }
    }

    void TestCreateMeshFrom2dMesh()
    {
        const double z_height(3.14);
        HoneycombVertexMeshGenerator generator2 (5, 5, false, 1, 1, 7);
        MutableVertexMesh<2, 2>* p_mesh2 = generator2.GetMesh();
        MonolayerVertexMeshGenerator generator;
        MutableVertexMesh<3, 3>* p_mesh = generator.MakeMeshUsing2dMesh(*p_mesh2, z_height);

        TS_ASSERT_EQUALS(p_mesh2->GetNumNodes()*2, p_mesh->GetNumNodes());
        TS_ASSERT_EQUALS(p_mesh2->GetNumElements(), p_mesh->GetNumElements());
        const unsigned num_lateral_faces = p_mesh->GetNumFaces() - 2*p_mesh2->GetNumElements();

        TS_ASSERT_EQUALS(p_mesh2->GetNumNodes()+p_mesh2->GetNumElements()+1-num_lateral_faces, 2u);

        // Check if nodes are created properly (their location)
        for (unsigned i=0; i<p_mesh2->GetNumNodes(); ++i)
        {
            c_vector<double, 2> loc2 = p_mesh2->GetNode(i)->rGetLocation();
            c_vector<double, 3> loc_basal = p_mesh->GetNode(i)->rGetLocation();
            c_vector<double, 3> loc_apical = p_mesh->GetNode(i+p_mesh2->GetNumNodes())->rGetLocation();

            TS_ASSERT(IsBasalNode(p_mesh->GetNode(i)));
            TS_ASSERT(IsApicalNode(p_mesh->GetNode(i+p_mesh2->GetNumNodes())));
            TS_ASSERT_DELTA(loc_basal[0], loc2[0], 1e-6);
            TS_ASSERT_DELTA(loc_basal[1], loc2[1], 1e-6);
            TS_ASSERT_DELTA(loc_basal[2], 0, 1e-6);

            TS_ASSERT_DELTA(loc_apical[0], loc2[0], 1e-6);
            TS_ASSERT_DELTA(loc_apical[1], loc2[1], 1e-6);
            TS_ASSERT_DELTA(loc_apical[2], z_height, 1e-6);
        }

        for (unsigned i=0; i<p_mesh2->GetNumElements(); ++i)
        {
            VertexElement<2, 2>* p_elem2 = p_mesh2->GetElement(i);
            VertexElement<3, 3>* p_elem = p_mesh->GetElement(i);
            std::set<unsigned> node_indices2;
            std::set<unsigned> node_indices;
            for (unsigned j=0; j<p_elem2->GetNumNodes(); ++j)
            {
                node_indices2.insert(p_elem2->GetNodeGlobalIndex(j));
                node_indices.insert(p_elem->GetNodeGlobalIndex(j));
            }
            TS_ASSERT_EQUALS(node_indices, node_indices2);
        }
    }


//    void TestMonolayerRearrangement()
//    {
//MARK
//        // Make 8 nodes to assign to a cube element
//        std::vector<Node<3>*> nodes;
//        nodes.push_back(new Node<3>(0, false, 0.0, 0.0, 0.0));
//        nodes.push_back(new Node<3>(1, false, 1.0, 0.0, 0.0));
//        nodes.push_back(new Node<3>(2, false, 0.0, 1.0, 0.0));
//        nodes.push_back(new Node<3>(3, false, 1.0, 1.0, 0.0));
//        nodes.push_back(new Node<3>(4, false, 0.0, 0.0, 1.0));
//        nodes.push_back(new Node<3>(5, false, 1.0, 0.0, 1.0));
//        nodes.push_back(new Node<3>(6, false, 0.0, 1.0, 1.0));
//        nodes.push_back(new Node<3>(7, false, 1.0, 1.0, 1.0));
//        for (unsigned i=0; i<4; ++i)
//        {
//            SetNodeAsBasal(nodes[i]);
//            SetNodeAsApical(nodes[i+4]);
//        }
//
//        std::vector<Node<3>*> nodes_face_0, nodes_face_1, nodes_face_2, nodes_face_3, nodes_face_4, nodes_face_5;
//
//        // Make 6 square faces out of these nodes
//        nodes_face_0.push_back(nodes[0]);
//        nodes_face_0.push_back(nodes[1]);
//        nodes_face_0.push_back(nodes[3]);
//        nodes_face_0.push_back(nodes[2]);
//
//        nodes_face_1.push_back(nodes[3]);
//        nodes_face_1.push_back(nodes[7]);
//        nodes_face_1.push_back(nodes[6]);
//        nodes_face_1.push_back(nodes[2]);
//
//        nodes_face_2.push_back(nodes[7]);
//        nodes_face_2.push_back(nodes[5]);
//        nodes_face_2.push_back(nodes[1]);
//        nodes_face_2.push_back(nodes[3]);
//
//        nodes_face_3.push_back(nodes[0]);
//        nodes_face_3.push_back(nodes[4]);
//        nodes_face_3.push_back(nodes[6]);
//        nodes_face_3.push_back(nodes[2]);
//
//        nodes_face_4.push_back(nodes[1]);
//        nodes_face_4.push_back(nodes[5]);
//        nodes_face_4.push_back(nodes[4]);
//        nodes_face_4.push_back(nodes[0]);
//
//        nodes_face_5.push_back(nodes[4]);
//        nodes_face_5.push_back(nodes[5]);
//        nodes_face_5.push_back(nodes[7]);
//        nodes_face_5.push_back(nodes[6]);
//
//        std::vector<VertexElement<2,3>*> faces;
//        faces.push_back(new VertexElement<2,3>(0, nodes_face_0));
//        faces.push_back(new VertexElement<2,3>(1, nodes_face_1));
//        faces.push_back(new VertexElement<2,3>(2, nodes_face_2));
//        faces.push_back(new VertexElement<2,3>(3, nodes_face_3));
//        faces.push_back(new VertexElement<2,3>(4, nodes_face_4));
//        faces.push_back(new VertexElement<2,3>(5, nodes_face_5));
//
//        std::vector<bool> orientations(faces.size());
//        for (unsigned i=0; i<faces.size(); i++)
//        {
//            orientations[i] = true;
//        }
//
//        // Make a cube element out of these faces
//        VertexElement<3,3> element(0, faces, orientations);
//        SetElementAsMonolayer(&element);
//
//        TS_ASSERT_EQUALS(element.GetNumNodes(), 8u);
//        TS_ASSERT_EQUALS(element.GetNumFaces(), 6u);
//
//        TS_ASSERT_EQUALS(element.GetIndex(), 0u);
//
//        // Test the position of some random nodes
//        TS_ASSERT_DELTA(element.GetFace(0)->GetNode(0)->rGetLocation()[0], 0.0, 1e-6);
//        TS_ASSERT_DELTA(element.GetFace(0)->GetNode(0)->rGetLocation()[1], 0.0, 1e-6);
//        TS_ASSERT_DELTA(element.GetFace(0)->GetNode(0)->rGetLocation()[2], 0.0, 1e-6);
//
//        TS_ASSERT_DELTA(element.GetFace(5)->GetNode(2)->rGetLocation()[0], 0.0, 1e-6);
//        TS_ASSERT_DELTA(element.GetFace(5)->GetNode(2)->rGetLocation()[1], 0.0, 1e-6);
//        TS_ASSERT_DELTA(element.GetFace(5)->GetNode(2)->rGetLocation()[2], 1.0, 1e-6);
//
//        // Test orientations
//        for (unsigned face_index=0; face_index<element.GetNumFaces(); face_index++)
//        {
//            TS_ASSERT_EQUALS(element.FaceIsOrientatedAntiClockwise(face_index), true);
//        }
//
//        // Tidy up
//        for (unsigned i=0; i<nodes.size(); i++)
//        {
//            delete nodes[i];
//        }
//        for (unsigned i=0; i<faces.size(); i++)
//        {
//            delete faces[i];
//        }
//
//    }

};

#endif /* TESTMONOLAYERVERTEXMESHGENERATOR_HPP_ */
