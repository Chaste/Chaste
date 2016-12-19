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

#ifndef TESTHEXAGONALPRISM3DVERTEXMESHGENERATOR_HPP_
#define TESTHEXAGONALPRISM3DVERTEXMESHGENERATOR_HPP_

#include <cxxtest/TestSuite.h>
#include "HexagonalPrism3dVertexMeshGenerator.hpp"
#include "FakePetscSetup.hpp"
#include <set>

class TestHexagonalPrism3dVertexMeshGenerator : public CxxTest::TestSuite
{
public:
    void TestSingle() throw (Exception)
    {
        HexagonalPrism3dVertexMeshGenerator generator(1, 1, 1.0, 1.0);
        MutableVertexMesh<3,3>* p_mesh = generator.GetMesh();
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 12u);
        TS_ASSERT_EQUALS(p_mesh->GetNumFaces(), 8u);
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 1u);

        VertexElement<3,3>* p_element_0 (p_mesh->GetElement(0));
        TS_ASSERT_EQUALS(p_element_0->GetNumNodes(), 12u);
        TS_ASSERT_EQUALS(p_element_0->GetNode(0)->GetIndex(), 0u);
        TS_ASSERT_EQUALS(p_element_0->GetNode(1)->GetIndex(), 2u);
        TS_ASSERT_EQUALS(p_element_0->GetNode(2)->GetIndex(), 4u);
        TS_ASSERT_EQUALS(p_element_0->GetNode(3)->GetIndex(), 5u);
        TS_ASSERT_EQUALS(p_element_0->GetNode(4)->GetIndex(), 3u);
        TS_ASSERT_EQUALS(p_element_0->GetNode(5)->GetIndex(), 1u);
        TS_ASSERT_EQUALS(p_element_0->GetNode(6)->GetIndex(), 6u);
        TS_ASSERT_EQUALS(p_element_0->GetNode(7)->GetIndex(), 8u);
        TS_ASSERT_EQUALS(p_element_0->GetNode(8)->GetIndex(), 10u);
        TS_ASSERT_EQUALS(p_element_0->GetNode(9)->GetIndex(), 11u);
        TS_ASSERT_EQUALS(p_element_0->GetNode(10)->GetIndex(), 9u);
        TS_ASSERT_EQUALS(p_element_0->GetNode(11)->GetIndex(), 7u);

        TS_ASSERT_DELTA(p_mesh->GetVolumeOfElement(0), 1, 1e-6);
    }

    void TestRowOfThreeElements() throw (Exception)
    {
        // Create a mesh comprising a row of three hexagonal prism elements in the x direction
        HexagonalPrism3dVertexMeshGenerator generator(3, 1, 5.0, 1.0);
        MutableVertexMesh<3,3>* p_mesh = generator.GetMesh();

        // Test that the mesh has the correct number of nodes, faces and elements
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 28u);
        TS_ASSERT_EQUALS(p_mesh->GetNumFaces(), 22u);
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 3u);

        // Test that each element has the correct number of nodes and faces
        for (unsigned elem_index=0; elem_index<p_mesh->GetNumElements(); elem_index++)
        {
            TS_ASSERT_EQUALS(p_mesh->GetElement(elem_index)->GetNumNodes(), 12u);
            TS_ASSERT_EQUALS(p_mesh->GetElement(elem_index)->GetNumFaces(), 8u);
            TS_ASSERT_DELTA(p_mesh->GetVolumeOfElement(elem_index), 5, 1e-6);
        }

        const unsigned a_tmp[6] = {0,1,8,9,15,16};
        const std::set<unsigned> s(a_tmp, a_tmp+18);
        // Test that each face has the correct number of nodes
        for (unsigned index=0; index<p_mesh->GetNumFaces(); index++)
        {
            unsigned face_index = p_mesh->GetFace(index)->GetIndex();
            if (s.find(face_index) != s.end())
            {
                TS_ASSERT_EQUALS(p_mesh->GetFace(index)->GetNumNodes(), 6u);
            }
            else
            {
                TS_ASSERT_EQUALS(p_mesh->GetFace(index)->GetNumNodes(), 4u);
            }
        }

        // Test all elements have the correct nodes and faces
        const unsigned node_indices[1][6] = {{0,4,8,11,7,3}};
        const unsigned face_indices[3][8] = {{0,1,2,3,4,5,6,7}, {8,9,10,11,12,13,3,14}, {15,16,17,18,19,20,11,21}};
        const unsigned half_num_nodes = p_mesh->GetNumNodes()/2;
        const unsigned max_y = 1;
        const unsigned max_x = 3;
        for (unsigned y_count=0, elem_index=0; y_count<max_y; ++y_count)
        {
            for (unsigned x_count=0; x_count<max_x; ++x_count, ++elem_index)
            {
                const VertexElement<3,3>* p_elem = p_mesh->GetElement(elem_index);
                for (unsigned i=0; i<6u; ++i)
                {
                    TS_ASSERT_EQUALS(p_elem->GetNode(i)->GetIndex(), node_indices[y_count][i] + x_count);
                    TS_ASSERT_EQUALS(p_elem->GetNode(i+6u)->GetIndex(), node_indices[y_count][i] + x_count + half_num_nodes);
                }
                for (unsigned i=0; i<8u; ++i)
                {
                    TS_ASSERT_EQUALS(p_elem->GetFace(i)->GetIndex(), face_indices[elem_index][i]);
                }
            }
        }

        // Create a vertex mesh writer
        VertexMeshWriter<3,3> vertex_mesh_writer("TestHexagonalPrism3dVertexMesh/RowOfThree", "vertex_mesh_3d_row_of_three");
        vertex_mesh_writer.WriteVtkUsingMeshWithCellId(*p_mesh);
    }

    void TestThreeByThreeElements() throw (Exception)
    {
        // Create a mesh comprising a row of nine hexagonal prism elements (three in the x direction, three in the y direction)
        HexagonalPrism3dVertexMeshGenerator generator(3, 3, 2/sqrt(3), 2.0);
        MutableVertexMesh<3,3>* p_mesh = generator.GetMesh();

        // Test that the mesh has the correct number of nodes, faces and elements
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 60u);
        TS_ASSERT_EQUALS(p_mesh->GetNumFaces(), 56u);
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 9u);

        // Test that each element has the correct number of nodes and faces
        for (unsigned elem_index=0; elem_index<p_mesh->GetNumElements(); elem_index++)
        {
            TS_ASSERT_EQUALS(p_mesh->GetElement(elem_index)->GetNumNodes(), 12u);
            TS_ASSERT_EQUALS(p_mesh->GetElement(elem_index)->GetNumFaces(), 8u);
        }

        const unsigned a_tmp[18] = {0,1,8,9,15,16,22,23,28,29,33,34,39,40,46,47,51,52};
        const std::set<unsigned> s(a_tmp, a_tmp+18);
        // Test that each face has the correct number of nodes
        for (unsigned index=0; index<p_mesh->GetNumFaces(); index++)
        {
            unsigned face_index = p_mesh->GetFace(index)->GetIndex();
            if (s.find(face_index) != s.end())
            {
                TS_ASSERT_EQUALS(p_mesh->GetFace(index)->GetNumNodes(), 6u);
            }
            else
            {
                TS_ASSERT_EQUALS(p_mesh->GetFace(index)->GetNumNodes(), 4u);
            }
        }

        // Test all elements have the correct nodes and faces
        const unsigned node_indices[3][6] = {{0,4,8,11,7,3},
                                             {8,12,16,20,15,11},
                                             {15,20,24,27,23,19}};
        const unsigned face_indices[9][8] = {{0,1,2,3,4,5,6,7}, {8,9,10,11,12,13,3,14}, {15,16,17,18,19,20,11,21},
                                             {22,23,13,24,25,26,27,4}, {28,29,20,30,31,32,24,12}, {33,34,35,36,37,38,30,19},
                                             {39,40,26,41,42,43,44,45}, {46,47,32,48,49,50,41,25}, {51,52,38,53,54,55,48,31}};
        const unsigned half_num_nodes = p_mesh->GetNumNodes()/2;
        const unsigned max_y = 3;
        const unsigned max_x = 3;
        for (unsigned y_count=0, elem_index=0; y_count<max_y; ++y_count)
        {
            for (unsigned x_count=0; x_count<max_x; ++x_count, ++elem_index)
            {
                const VertexElement<3,3>* p_elem = p_mesh->GetElement(elem_index);
                for (unsigned i=0; i<6u; ++i)
                {
                    TS_ASSERT_EQUALS(p_elem->GetNode(i)->GetIndex(), node_indices[y_count][i] + x_count);
                    TS_ASSERT_EQUALS(p_elem->GetNode(i+6u)->GetIndex(), node_indices[y_count][i] + x_count + half_num_nodes);
                }
                for (unsigned i=0; i<8u; ++i)
                {
                    TS_ASSERT_EQUALS(p_elem->GetFace(i)->GetIndex(), face_indices[elem_index][i]);
                }
            }
        }

        // Create a vertex mesh writer
        VertexMeshWriter<3,3> vertex_mesh_writer("TestHexagonalPrism3dVertexMesh/ThreeByThree", "vertex_mesh_3d_three_by_three");
        vertex_mesh_writer.WriteVtkUsingMeshWithCellId(*p_mesh);
    }

    void TestFourByFourElements() throw (Exception)
    {
        // Create a mesh
        HexagonalPrism3dVertexMeshGenerator generator(4, 4, 2.0, 2.0);
        MutableVertexMesh<3,3>* p_mesh = generator.GetMesh();

        // Test that the mesh has the correct number of nodes, faces and elements
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 96u);
        TS_ASSERT_EQUALS(p_mesh->GetNumFaces(), 95u);
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 16u);

        // Test that each face has the correct number of nodes
        const unsigned a_tmp[32] = {0,1, 8,9, 15,16, 22,23,
                                    29,30, 35,36, 40,41, 45,46,
                                    51,52, 58,59, 63,64, 68,69,
                                    73,74, 79,80, 84,85, 89,90};
        const std::set<unsigned> s(a_tmp, a_tmp+32);
        for (unsigned index=0; index<p_mesh->GetNumFaces(); ++index)
        {
            unsigned face_index = p_mesh->GetFace(index)->GetIndex();
            if (s.find(face_index) != s.end())
            {
                TS_ASSERT_EQUALS(p_mesh->GetFace(index)->GetNumNodes(), 6u);
            }
            else
            {
                TS_ASSERT_EQUALS(p_mesh->GetFace(index)->GetNumNodes(), 4u);
            }
        }

        // Test all elements have the correct nodes and faces
        const unsigned node_indices[4][6] =
        {{0,5,10,14,9,4}, {10,15,20,25,19,14}, {19,25,30,34,29,24}, {30,35,40,44,39,34}};

        const unsigned face_indices[16][8] = {
        {0,1,2,3,4,5,6,7}, {8,9,10,11,12,13,3,14},
        {15,16,17,18,19,20,11,21}, {22,23,24,25,26,27,18,28},
        {29,30,13,31,32,33,34,4}, {35,36,20,37,38,39,31,12},
        {40,41,27,42,43,44,37,19}, {45,46,47,48,49,50,42,26},
        {51,52,33,53,54,55,56,57}, {58,59,39,60,61,62,53,32},
        {63,64,44,65,66,67,60,38}, {68,69,50,70,71,72,65,43},
        {73,74,62,75,76,77,78,54}, {79,80,67,81,82,83,75,61},
        {84,85,72,86,87,88,81,66}, {89,90,91,92,93,94,86,71}};
        const unsigned half_num_nodes = p_mesh->GetNumNodes()/2;
        const unsigned max_y = 4;
        const unsigned max_x = 4;
        for (unsigned y_count=0, elem_index=0; y_count<max_y; ++y_count)
        {
            for (unsigned x_count=0; x_count<max_x; ++x_count, ++elem_index)
            {
                const VertexElement<3,3>* p_elem = p_mesh->GetElement(elem_index);
                for (unsigned i=0; i<6u; ++i)
                {
                    TS_ASSERT_EQUALS(p_elem->GetNode(i)->GetIndex(), node_indices[y_count][i] + x_count);
                    TS_ASSERT_EQUALS(p_elem->GetNode(i+6u)->GetIndex(), node_indices[y_count][i] + x_count + half_num_nodes);
                }
                for (unsigned i=0; i<8u; ++i)
                {
                    TS_ASSERT_EQUALS(p_elem->GetFace(i)->GetIndex(), face_indices[elem_index][i]);
                }
            }
        }

        // Create a vertex mesh writer
        VertexMeshWriter<3,3> vertex_mesh_writer3("TestHexagonalPrism3dVertexMesh/FourByFour", "vertex_mesh_3d_four_by_four");
        vertex_mesh_writer3.WriteVtkUsingMeshWithCellId(*p_mesh);
    }

    void TestLargeMesh() throw (Exception)
    {
        // Create a mesh
        unsigned num_rows = 5;
        unsigned num_columns = 5;
        HexagonalPrism3dVertexMeshGenerator generator(num_rows, num_columns, 1.0, 2.0);
        MutableVertexMesh<3,3>* p_mesh = generator.GetMesh();

        // Test that the mesh has the correct number of nodes, faces and elements
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 140u);
        TS_ASSERT_EQUALS(p_mesh->GetNumFaces(), 144u);
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), num_rows*num_columns);

        // Test that each face has the correct number of nodes
        const unsigned a_tmp[50] = {0,1,8,9,15,16,22,23,29,30,36,37,42,43,47,48,52,
                                    53,57,58,63,64,70,71,75,76,80,81,85,86,90,91,96,
                                    97,101,102,106,107,111,112,117,118,124,125,129,
                                    130,134,135,139,140};
        const std::set<unsigned> s (a_tmp, a_tmp+50);
        for (unsigned index=0; index<p_mesh->GetNumFaces(); ++index)
        {
            unsigned face_index = p_mesh->GetFace(index)->GetIndex();
            if (s.find(face_index) != s.end())
            {
                TS_ASSERT_EQUALS(p_mesh->GetFace(index)->GetNumNodes(), 6u);
            }
            else
            {
                TS_ASSERT_EQUALS(p_mesh->GetFace(index)->GetNumNodes(), 4u);
            }
        }
    }
};

#endif /*TESTHEXAGONALPRISM3DVERTEXMESHGENERATOR_HPP_*/
