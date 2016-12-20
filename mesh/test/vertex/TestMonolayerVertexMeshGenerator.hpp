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




/**
 * temporary for GeodesicSphere
 */
#include "Debug.hpp"
#include <vector>
#include <set>
#include <utility>              // for pair
#include <algorithm>            // sort
#include "UblasCustomFunctions.hpp"
#include "VertexElement.hpp"
#include "Node.hpp"
#include "MonolayerVertexMeshCustomFunctions.hpp"
#include "VertexMeshWriter.hpp"
#include "MutableVertexMesh.hpp"
#include <cmath>                                    // for atan2
#include "MonolayerVertexMeshGenerator.hpp"


class GeodesicSphereBuilder
{
private:
    unsigned mDepth;
    static const double X = 0.525731112119133606;
    static const double Z = 0.850650808352039932;

    class SortWithIndex
    {
    public:
        bool operator()(const std::pair<double, unsigned>& a, const std::pair<double, unsigned>& b)
        {
            return a.first < b.first;
        }
    };

public:
    std::vector<Node<3>*> mNodes;
    std::vector<VertexElement<1,3>*> mEdges;
    std::vector<VertexElement<2,3>*> mFaces;

    c_vector<double,3> normalise(const c_vector<double,3>& v)
                {
        return v/norm_2(v);
                }

    GeodesicSphereBuilder(const unsigned numDivision=0u)
    : mDepth(0)
    {
        const double vdata[12][3] = {
                {-X, 0.0, Z}, { X, 0.0, Z }, { -X, 0.0, -Z }, { X, 0.0, -Z },
                { 0.0, Z, X }, { 0.0, Z, -X }, { 0.0, -Z, X }, { 0.0, -Z, -X },
                { Z, X, 0.0 }, { -Z, X, 0.0 }, { Z, -X, 0.0 }, { -Z, -X, 0.0 }};
        for (unsigned i=0; i<12; ++i)
        {
            mNodes.push_back(new Node<3>(i, false, vdata[i][0], vdata[i][1], vdata[i][2]));
        }

        int tindices[20][3] = {
                {0, 4, 1}, { 0, 9, 4 }, { 9, 5, 4 }, { 4, 5, 8 }, { 4, 8, 1 },
                { 8, 10, 1 }, { 8, 3, 10 }, { 5, 3, 8 }, { 5, 2, 3 }, { 2, 7, 3 },
                { 7, 10, 3 }, { 7, 6, 10 }, { 7, 11, 6 }, { 11, 0, 6 }, { 0, 1, 6 },
                { 6, 1, 10 }, { 9, 0, 11 }, { 9, 11, 2 }, { 9, 2, 5 }, { 7, 2, 11 }};

        // Creating each face (edges at the same time.)
        for (unsigned i = 0; i<20; ++i)
        {
            std::vector<Node<3>*> this_face_nodes;
            std::vector<VertexElement<1,3>*> this_face_edges;
            // well, orientation doesn't really matter.
            std::vector<bool> this_face_orientations(3);

            for (unsigned ja=0; ja<3; ++ja)
            {
                this_face_nodes.push_back(mNodes[tindices[i][ja]]);

                const unsigned jb = (ja+1)%3;
                Node<3>* node_a = mNodes[tindices[i][ja]];
                Node<3>* node_b = mNodes[tindices[i][jb]];
//PRINT_3_VARIABLES(i, tindices[i][ja], tindices[i][jb])
//PRINT_CONTAINER(node_a_edges)
//PRINT_CONTAINER(node_b_edges)
                std::set<unsigned> shared_edge = GetSharedFaceIndices(node_a, node_b);

                unsigned current_edge_index = UINT_MAX;
                if (shared_edge.size() == 0)
                {
                    std::vector<Node<3>*> this_edge_nodes;
                    this_edge_nodes.push_back(node_a);
                    this_edge_nodes.push_back(node_b);
                    current_edge_index = mEdges.size();
                    mEdges.push_back(new VertexElement<1,3>(current_edge_index, this_edge_nodes));
                }
                else
                {
//TRACE("FOUND shared edge")
//PRINT_3_VARIABLES(i, tindices[i][ja], tindices[i][jb])
//PRINT_CONTAINER(node_a_edges)
//PRINT_CONTAINER(node_b_edges)
                    assert(shared_edge.size() == 1);
                    current_edge_index = *(shared_edge.begin());
                }

                this_face_edges.push_back(mEdges[current_edge_index]);
            }
            assert(mFaces.size() == i);
            mFaces.push_back(new VertexElement<2,3>(mFaces.size(), this_face_edges, this_face_orientations, this_face_nodes));
        }

PRINT_3_VARIABLES(mNodes.size(), mEdges.size(), mFaces.size());

        for (unsigned i=0; i<numDivision; ++i)
        {
            SubDivide();
PRINT_3_VARIABLES(mNodes.size(), mEdges.size(), mFaces.size());
        }
    }

    void SubDivide()
    {
        const unsigned num_old_nodes = mNodes.size();
        const unsigned num_old_edges = mEdges.size();
        const unsigned num_old_faces = mFaces.size();
        if (num_old_nodes+num_old_faces-num_old_edges != 2)
        {
            TRACE("EXCEPTION EULER IS THE BOSS")
                            EXCEPTION("Breaking Euler's Polyhedron Formula!");
        }

        // For sanity check later
        const unsigned num_new_total_nodes = num_old_nodes + num_old_edges;
        const unsigned num_new_edges = 2*num_old_edges + 3*num_old_faces;
        const unsigned num_new_faces = 4*num_old_faces;

        // Each edge will be "bisected", and each face will become 4 new faces.
        std::vector<Node<3>*> map_edge_to_new_nodes(mEdges.size());
        for (unsigned i=0; i<num_old_edges; ++i)
        {
            Node<3>* p_tmp_node = new Node<3>(mNodes.size(), normalise(mEdges[i]->GetCentroid()), false);
            map_edge_to_new_nodes[i] = p_tmp_node;
            mNodes.push_back(p_tmp_node);
        }
        assert(mNodes.size() == num_new_total_nodes);

        const unsigned cyc[3][3] = {{0,1,2},{1,2,0},{2,0,1}};
        // Creating new faces and edges
        for (unsigned old_face_index=0; old_face_index<num_old_faces; ++old_face_index)
        {
                /*
                 *           1
                 *           /\
                 *          /  \
                 *         /  b \
                 *      0'/______\2'
                 *       /\      /\
                 *      /  \  d /  \
                 *     /  c \  /  a \
                 *    /______\/______\
                 *   2        1'      0
                 */
            VertexElement<2, 3>* p_face = mFaces[old_face_index];
            assert(old_face_index == p_face->GetIndex());

            Node<3>* old_nodes[3] = {p_face->GetNode(0), p_face->GetNode(1), p_face->GetNode(2)};
            std::vector<Node<3>*> mid_nodes(3);
            for (unsigned i=0; i<3; ++i)
            {
                const std::set<unsigned> shared_edge =
                        GetSharedFaceIndices(old_nodes[cyc[i][0]], old_nodes[cyc[i][1]]);
                assert(shared_edge.size() == 1);
                assert(*(shared_edge.begin()) == p_face->GetFace(i)->GetIndex());
                mid_nodes[cyc[i][2]] = map_edge_to_new_nodes[*(shared_edge.begin())];
            }

            // create face a,b,c
            for (unsigned i=0; i<4; ++i)
            {
                std::vector<Node<3>*> this_face_nodes(3);
                if (i<3)
                {
                    this_face_nodes[0] = old_nodes[cyc[i][0]];
                    this_face_nodes[1] = mid_nodes[cyc[i][1]];
                    this_face_nodes[2] = mid_nodes[cyc[i][2]];
                }
                else
                {
                    this_face_nodes[0] = mid_nodes[0];
                    this_face_nodes[1] = mid_nodes[1];
                    this_face_nodes[2] = mid_nodes[2];
                }


                std::vector<VertexElement<1,3>*> this_face_edges(3);
                std::vector<bool> this_orientation(3);
                for (unsigned j=0; j<3; ++j)
                {
                    const std::set<unsigned> shared_edge = GetSharedFaceIndices(this_face_nodes[cyc[j][0]], this_face_nodes[cyc[j][1]]);
                    VertexElement<1,3>* p_tmp_edge;
                    if (shared_edge.size() == 0)
                    {
                        std::vector<Node<3>*> tmp_edge_nodes;
                        tmp_edge_nodes.push_back(this_face_nodes[cyc[j][0]]);
                        tmp_edge_nodes.push_back(this_face_nodes[cyc[j][1]]);

                        p_tmp_edge = new VertexElement<1,3>(mEdges.size(), tmp_edge_nodes);
                        p_tmp_edge->RegisterFaceWithNodes();
                        mEdges.push_back(p_tmp_edge);
                    }
                    else
                    {
                        assert(shared_edge.size() == 1);
                        p_tmp_edge = mEdges[*(shared_edge.begin())];
                    }
                    this_face_edges[j] = p_tmp_edge;
                }

                if (i<3)
                {
                    mFaces.push_back(new VertexElement<2,3>(mFaces.size(), this_face_edges,
                            this_orientation, this_face_nodes));
                }
                else
                {
                    p_face->MarkAsDeleted();
                    delete p_face;
                    mFaces[old_face_index] = new VertexElement<2,3>(old_face_index, this_face_edges,
                            this_orientation, this_face_nodes);
                }
            }
        }


        assert(mEdges.size() == num_old_edges+num_new_edges);
        for (unsigned i=0; i<num_old_edges; ++i)
        {
            mEdges[i]->MarkFaceAsDeleted();
            delete mEdges[i];

            mEdges[i] = mEdges[num_new_edges+i];
            mEdges[i]->FaceResetIndex(i);
        }

        for (unsigned i=num_new_edges, j=0; i<mEdges.size(); ++i, ++j)
        {
            assert(mEdges[i]->GetIndex() == j);
        }
        mEdges.resize(num_new_edges);

        assert(mFaces.size() == num_new_faces);
MARK;
    }

    MutableVertexMesh<2, 3>* GetDual()
    {
        std::vector<Node<3>*> dual_nodes(mFaces.size());
        for (unsigned dual_node_index=0; dual_node_index<mFaces.size(); ++dual_node_index)
        {
            dual_nodes[dual_node_index] = new Node<3>(dual_node_index, mFaces[dual_node_index]->GetCentroid(), false);
        }

        std::vector<VertexElement<2, 3>*> dual_faces(mNodes.size());
        for (unsigned dual_face_index=0; dual_face_index<mNodes.size(); ++dual_face_index)
        {
            const c_vector<double, 3> normal_v = normalise(mNodes[dual_face_index]->rGetLocation());

            const std::set<unsigned> s_tmp = mNodes[dual_face_index]->rGetContainingElementIndices();
            // Switch to vector for better access
            const std::vector<unsigned> this_dual_node_global_indices(s_tmp.begin(), s_tmp.end());
            const unsigned this_num_dual_nodes = this_dual_node_global_indices.size();


            std::vector<Node<3>*> unsorted_this_dual_face_nodes(this_num_dual_nodes);
            std::vector<std::pair<double, unsigned> > angles_with_index;
            unsorted_this_dual_face_nodes[0] = dual_nodes[this_dual_node_global_indices[0]];
            //            angles_with_index.push_back(std::make_pair(0.0, 0));

            // Create orthonormal basis with the face normal, first dual node location of the face.
            // e1 is obtained by Gram-Schmidt process.
            const c_vector<double, 3> e1 = normalise(unsorted_this_dual_face_nodes[0]->rGetLocation()
                    - inner_prod(unsorted_this_dual_face_nodes[0]->rGetLocation(), normal_v)*normal_v);
            // e2 simply by cross prod of normal_v and e1
            const c_vector<double, 3> e2 = VectorProduct(normal_v, e1);
            assert(abs(norm_2(e2)-1) < 1e-5);

            for (unsigned dual_node_local_index=0; dual_node_local_index<this_num_dual_nodes; ++dual_node_local_index)
            {
                Node<3>* p_node = dual_nodes[this_dual_node_global_indices[dual_node_local_index]];
                unsorted_this_dual_face_nodes[dual_node_local_index] = p_node;

                const c_vector<double, 3> r_node = p_node->rGetLocation();
                const double tmp_angle = atan2(inner_prod(r_node, e1), inner_prod(r_node, e2));
                angles_with_index.push_back(std::make_pair(tmp_angle, dual_node_local_index));
            }

            std::sort(angles_with_index.begin(), angles_with_index.end(), SortWithIndex());
            std::vector<Node<3>*> this_dual_face_nodes(this_num_dual_nodes);

            std::vector<double> Debug_anglesss;
            for (unsigned iii=0; iii<this_num_dual_nodes; ++iii)
            {
                this_dual_face_nodes[iii] = unsorted_this_dual_face_nodes[angles_with_index[iii].second];

                Debug_anglesss.push_back(angles_with_index[iii].first);
            }
//PRINT_VECTOR(Debug_anglesss);

            dual_faces[dual_face_index] = new VertexElement<2, 3>(dual_face_index, this_dual_face_nodes);
        }


        return new MutableVertexMesh<2,3>(dual_nodes, dual_faces);
    }
};





class TestMonolayerVertexMeshGenerator : public CxxTest::TestSuite
{
public:

    void TestGeodesicShpere()
    {
        GeodesicSphereBuilder builder;
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

    void TestCreateMeshFrom2dMesh() throw (Exception)
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
};

#endif /* TESTMONOLAYERVERTEXMESHGENERATOR_HPP_ */
