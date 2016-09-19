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
#ifndef TESTHEXAGONALPRISM3DVERTEXMESHGENERATOR_HPP_
#define TESTHEXAGONALPRISM3DVERTEXMESHGENERATOR_HPP_

#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "FakePetscSetup.hpp"
#include "HexagonalPrism3dVertexMeshGenerator.hpp"
#include "Debug.hpp"

class TestHexagonalPrism3dVertexMeshGenerator : public AbstractCellBasedTestSuite
{
public:
    void TestRowOfThreeElements() throw (Exception)
    {
        // Create a mesh comprising a row of three hexagonal prism elements in the x direction
        HexagonalPrism3dVertexMeshGenerator generator(3, 1, 1.0, 1.0);
        MutableVertexMesh<3,3>* p_mesh = generator.GetMesh();

        // Test that the mesh has the correct number of nodes, faces and elements
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 28u);
        TS_ASSERT_EQUALS(p_mesh->GetNumFaces(), 22u);
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 3u);

        // Test that each face has the correct number of nodes
        for (unsigned index=0; index<p_mesh->GetNumFaces(); index++)
        {
            unsigned face_index = p_mesh->GetFace(index)->GetIndex();
            if (face_index < 6)
            {
                TS_ASSERT_EQUALS(p_mesh->GetFace(index)->GetNumNodes(), 6u);
            }
            else
            {
                TS_ASSERT_EQUALS(p_mesh->GetFace(index)->GetNumNodes(), 4u);
            }
        }

        // Test that each element has the correct number of nodes and faces
        for (unsigned elem_index=0; elem_index<p_mesh->GetNumElements(); elem_index++)
        {
            TS_ASSERT_EQUALS(p_mesh->GetElement(elem_index)->GetNumNodes(), 12u);
            TS_ASSERT_EQUALS(p_mesh->GetElement(elem_index)->GetNumFaces(), 8u);
        }

        // Test the element 0 has the correct nodes and faces
        VertexElement<3,3>* p_element_0 = p_mesh->GetElement(0);

        TS_ASSERT_EQUALS(p_element_0->GetNumNodes(), 12u);
        TS_ASSERT_EQUALS(p_element_0->GetNode(0)->GetIndex(), 0u);
        TS_ASSERT_EQUALS(p_element_0->GetNode(1)->GetIndex(), 14u);
        TS_ASSERT_EQUALS(p_element_0->GetNode(2)->GetIndex(), 3u);
        TS_ASSERT_EQUALS(p_element_0->GetNode(3)->GetIndex(), 17u);
        TS_ASSERT_EQUALS(p_element_0->GetNode(4)->GetIndex(), 4u);
        TS_ASSERT_EQUALS(p_element_0->GetNode(5)->GetIndex(), 18u);
        TS_ASSERT_EQUALS(p_element_0->GetNode(6)->GetIndex(), 7u);
        TS_ASSERT_EQUALS(p_element_0->GetNode(7)->GetIndex(), 21u);
        TS_ASSERT_EQUALS(p_element_0->GetNode(8)->GetIndex(), 8u);
        TS_ASSERT_EQUALS(p_element_0->GetNode(9)->GetIndex(), 22u);
        TS_ASSERT_EQUALS(p_element_0->GetNode(10)->GetIndex(), 11u);
        TS_ASSERT_EQUALS(p_element_0->GetNode(11)->GetIndex(), 25u);

        TS_ASSERT_EQUALS(p_element_0->GetNumFaces(), 8u);
        TS_ASSERT_EQUALS(p_element_0->GetFace(0)->GetIndex(), 0u);
        TS_ASSERT_EQUALS(p_element_0->GetFace(1)->GetIndex(), 3u);
        TS_ASSERT_EQUALS(p_element_0->GetFace(2)->GetIndex(), 6u);
        TS_ASSERT_EQUALS(p_element_0->GetFace(3)->GetIndex(), 7u);
        TS_ASSERT_EQUALS(p_element_0->GetFace(4)->GetIndex(), 8u);
        TS_ASSERT_EQUALS(p_element_0->GetFace(5)->GetIndex(), 9u);
        TS_ASSERT_EQUALS(p_element_0->GetFace(6)->GetIndex(), 10u);
        TS_ASSERT_EQUALS(p_element_0->GetFace(7)->GetIndex(), 11u);

        // Test the element 1 has the correct nodes and faces
        VertexElement<3,3>* p_element_1 = p_mesh->GetElement(1);

        TS_ASSERT_EQUALS(p_element_1->GetNumNodes(), 12u);
        TS_ASSERT_EQUALS(p_element_1->GetNode(0)->GetIndex(), 1u);
        TS_ASSERT_EQUALS(p_element_1->GetNode(1)->GetIndex(), 15u);
        TS_ASSERT_EQUALS(p_element_1->GetNode(2)->GetIndex(), 4u);
        TS_ASSERT_EQUALS(p_element_1->GetNode(3)->GetIndex(), 18u);
        TS_ASSERT_EQUALS(p_element_1->GetNode(4)->GetIndex(), 5u);
        TS_ASSERT_EQUALS(p_element_1->GetNode(5)->GetIndex(), 19u);
        TS_ASSERT_EQUALS(p_element_1->GetNode(6)->GetIndex(), 8u);
        TS_ASSERT_EQUALS(p_element_1->GetNode(7)->GetIndex(), 22u);
        TS_ASSERT_EQUALS(p_element_1->GetNode(8)->GetIndex(), 9u);
        TS_ASSERT_EQUALS(p_element_1->GetNode(9)->GetIndex(), 23u);
        TS_ASSERT_EQUALS(p_element_1->GetNode(10)->GetIndex(), 12u);
        TS_ASSERT_EQUALS(p_element_1->GetNode(11)->GetIndex(), 26u);

        TS_ASSERT_EQUALS(p_element_1->GetNumFaces(), 8u);
        TS_ASSERT_EQUALS(p_element_1->GetFace(0)->GetIndex(), 1u);
        TS_ASSERT_EQUALS(p_element_1->GetFace(1)->GetIndex(), 4u);
        TS_ASSERT_EQUALS(p_element_1->GetFace(2)->GetIndex(), 12u);
        TS_ASSERT_EQUALS(p_element_1->GetFace(3)->GetIndex(), 13u);
        TS_ASSERT_EQUALS(p_element_1->GetFace(4)->GetIndex(), 14u);
        TS_ASSERT_EQUALS(p_element_1->GetFace(5)->GetIndex(), 15u);
        TS_ASSERT_EQUALS(p_element_1->GetFace(6)->GetIndex(), 16u);
        TS_ASSERT_EQUALS(p_element_1->GetFace(7)->GetIndex(), 7u);

        // Test the element 2 has the correct nodes and faces
        VertexElement<3,3>* p_element_2 = p_mesh->GetElement(2);

        TS_ASSERT_EQUALS(p_element_2->GetNumNodes(), 12u);
        TS_ASSERT_EQUALS(p_element_2->GetNode(0)->GetIndex(), 2u);
        TS_ASSERT_EQUALS(p_element_2->GetNode(1)->GetIndex(), 16u);
        TS_ASSERT_EQUALS(p_element_2->GetNode(2)->GetIndex(), 5u);
        TS_ASSERT_EQUALS(p_element_2->GetNode(3)->GetIndex(), 19u);
        TS_ASSERT_EQUALS(p_element_2->GetNode(4)->GetIndex(), 6u);
        TS_ASSERT_EQUALS(p_element_2->GetNode(5)->GetIndex(), 20u);
        TS_ASSERT_EQUALS(p_element_2->GetNode(6)->GetIndex(), 9u);
        TS_ASSERT_EQUALS(p_element_2->GetNode(7)->GetIndex(), 23u);
        TS_ASSERT_EQUALS(p_element_2->GetNode(8)->GetIndex(), 10u);
        TS_ASSERT_EQUALS(p_element_2->GetNode(9)->GetIndex(), 24u);
        TS_ASSERT_EQUALS(p_element_2->GetNode(10)->GetIndex(), 13u);
        TS_ASSERT_EQUALS(p_element_2->GetNode(11)->GetIndex(), 27u);

        TS_ASSERT_EQUALS(p_element_2->GetNumFaces(), 8u);
        TS_ASSERT_EQUALS(p_element_2->GetFace(0)->GetIndex(), 2u);
        TS_ASSERT_EQUALS(p_element_2->GetFace(1)->GetIndex(), 5u);
        TS_ASSERT_EQUALS(p_element_2->GetFace(2)->GetIndex(), 17u);
        TS_ASSERT_EQUALS(p_element_2->GetFace(3)->GetIndex(), 18u);
        TS_ASSERT_EQUALS(p_element_2->GetFace(4)->GetIndex(), 19u);
        TS_ASSERT_EQUALS(p_element_2->GetFace(5)->GetIndex(), 20u);
        TS_ASSERT_EQUALS(p_element_2->GetFace(6)->GetIndex(), 21u);
        TS_ASSERT_EQUALS(p_element_2->GetFace(7)->GetIndex(), 14u);
    }

    void TestThreeByThreeElements() throw (Exception)
    {
        // Create a mesh comprising a row of nine hexagonal prism elements (three in the x direction, three in the y direction)
        HexagonalPrism3dVertexMeshGenerator generator(3, 3, 1.0, 2.0);
        MutableVertexMesh<3,3>* p_mesh = generator.GetMesh();

        // Test that the mesh has the correct number of nodes, faces and elements
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 60u);
        TS_ASSERT_EQUALS(p_mesh->GetNumFaces(), 56u);
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 9u);

        // Test that each face has the correct number of nodes
        for (unsigned index=0; index<p_mesh->GetNumFaces(); index++)
        {
            unsigned face_index = p_mesh->GetFace(index)->GetIndex();
            if (face_index < 18)
            {
                TS_ASSERT_EQUALS(p_mesh->GetFace(index)->GetNumNodes(), 6u);
            }
            else
            {
                TS_ASSERT_EQUALS(p_mesh->GetFace(index)->GetNumNodes(), 4u);
            }
        }

        // Test that each element has the correct number of nodes and faces
        for (unsigned elem_index=0; elem_index<p_mesh->GetNumElements(); elem_index++)
        {
            TS_ASSERT_EQUALS(p_mesh->GetElement(elem_index)->GetNumNodes(), 12u);
            TS_ASSERT_EQUALS(p_mesh->GetElement(elem_index)->GetNumFaces(), 8u);
        }

        // Test the element 0 has the correct nodes and faces
        VertexElement<3,3>* p_element_0 = p_mesh->GetElement(0);

        TS_ASSERT_EQUALS(p_element_0->GetNumNodes(), 12u);
        TS_ASSERT_EQUALS(p_element_0->GetNode(0)->GetIndex(), 0u);
        TS_ASSERT_EQUALS(p_element_0->GetNode(1)->GetIndex(), 30u);
        TS_ASSERT_EQUALS(p_element_0->GetNode(2)->GetIndex(), 3u);
        TS_ASSERT_EQUALS(p_element_0->GetNode(3)->GetIndex(), 33u);
        TS_ASSERT_EQUALS(p_element_0->GetNode(4)->GetIndex(), 4u);
        TS_ASSERT_EQUALS(p_element_0->GetNode(5)->GetIndex(), 34u);
        TS_ASSERT_EQUALS(p_element_0->GetNode(6)->GetIndex(), 7u);
        TS_ASSERT_EQUALS(p_element_0->GetNode(7)->GetIndex(), 37u);
        TS_ASSERT_EQUALS(p_element_0->GetNode(8)->GetIndex(), 8u);
        TS_ASSERT_EQUALS(p_element_0->GetNode(9)->GetIndex(), 38u);
        TS_ASSERT_EQUALS(p_element_0->GetNode(10)->GetIndex(), 11u);
        TS_ASSERT_EQUALS(p_element_0->GetNode(11)->GetIndex(), 41u);

        TS_ASSERT_EQUALS(p_element_0->GetNumFaces(), 8u);
        TS_ASSERT_EQUALS(p_element_0->GetFace(0)->GetIndex(), 0u);
        TS_ASSERT_EQUALS(p_element_0->GetFace(1)->GetIndex(), 9u);
        TS_ASSERT_EQUALS(p_element_0->GetFace(2)->GetIndex(), 18u);
        TS_ASSERT_EQUALS(p_element_0->GetFace(3)->GetIndex(), 19u);
        TS_ASSERT_EQUALS(p_element_0->GetFace(4)->GetIndex(), 20u);
        TS_ASSERT_EQUALS(p_element_0->GetFace(5)->GetIndex(), 21u);
        TS_ASSERT_EQUALS(p_element_0->GetFace(6)->GetIndex(), 22u);
        TS_ASSERT_EQUALS(p_element_0->GetFace(7)->GetIndex(), 23u);

        // Test the element 3 has the correct nodes and faces
        VertexElement<3,3>* p_element_3 = p_mesh->GetElement(3);

        TS_ASSERT_EQUALS(p_element_3->GetNumNodes(), 12u);
        TS_ASSERT_EQUALS(p_element_3->GetNode(0)->GetIndex(), 8u);
        TS_ASSERT_EQUALS(p_element_3->GetNode(1)->GetIndex(), 38u);
        TS_ASSERT_EQUALS(p_element_3->GetNode(2)->GetIndex(), 11u);
        TS_ASSERT_EQUALS(p_element_3->GetNode(3)->GetIndex(), 41u);
        TS_ASSERT_EQUALS(p_element_3->GetNode(4)->GetIndex(), 12u);
        TS_ASSERT_EQUALS(p_element_3->GetNode(5)->GetIndex(), 42u);
        TS_ASSERT_EQUALS(p_element_3->GetNode(6)->GetIndex(), 15u);
        TS_ASSERT_EQUALS(p_element_3->GetNode(7)->GetIndex(), 45u);
        TS_ASSERT_EQUALS(p_element_3->GetNode(8)->GetIndex(), 16u);
        TS_ASSERT_EQUALS(p_element_3->GetNode(9)->GetIndex(), 46u);
        TS_ASSERT_EQUALS(p_element_3->GetNode(10)->GetIndex(), 20u);
        TS_ASSERT_EQUALS(p_element_3->GetNode(11)->GetIndex(), 50u);

        TS_ASSERT_EQUALS(p_element_3->GetNumFaces(), 8u);
        TS_ASSERT_EQUALS(p_element_3->GetFace(0)->GetIndex(), 3u);
        TS_ASSERT_EQUALS(p_element_3->GetFace(1)->GetIndex(), 12u);
        TS_ASSERT_EQUALS(p_element_3->GetFace(2)->GetIndex(), 28u);
        TS_ASSERT_EQUALS(p_element_3->GetFace(3)->GetIndex(), 34u);
        TS_ASSERT_EQUALS(p_element_3->GetFace(4)->GetIndex(), 35u);
        TS_ASSERT_EQUALS(p_element_3->GetFace(5)->GetIndex(), 36u);
        TS_ASSERT_EQUALS(p_element_3->GetFace(6)->GetIndex(), 37u);
        TS_ASSERT_EQUALS(p_element_3->GetFace(7)->GetIndex(), 20u);

        // Test the element 4 has the correct nodes and faces
        VertexElement<3,3>* p_element_4 = p_mesh->GetElement(4);

        TS_ASSERT_EQUALS(p_element_4->GetNumNodes(), 12u);
        TS_ASSERT_EQUALS(p_element_4->GetNode(0)->GetIndex(), 9u);
        TS_ASSERT_EQUALS(p_element_4->GetNode(1)->GetIndex(), 39u);
        TS_ASSERT_EQUALS(p_element_4->GetNode(2)->GetIndex(), 12u);
        TS_ASSERT_EQUALS(p_element_4->GetNode(3)->GetIndex(), 42u);
        TS_ASSERT_EQUALS(p_element_4->GetNode(4)->GetIndex(), 13u);
        TS_ASSERT_EQUALS(p_element_4->GetNode(5)->GetIndex(), 43u);
        TS_ASSERT_EQUALS(p_element_4->GetNode(6)->GetIndex(), 16u);
        TS_ASSERT_EQUALS(p_element_4->GetNode(7)->GetIndex(), 46u);
        TS_ASSERT_EQUALS(p_element_4->GetNode(8)->GetIndex(), 17u);
        TS_ASSERT_EQUALS(p_element_4->GetNode(9)->GetIndex(), 47u);
        TS_ASSERT_EQUALS(p_element_4->GetNode(10)->GetIndex(), 21u);
        TS_ASSERT_EQUALS(p_element_4->GetNode(11)->GetIndex(), 51u);

        TS_ASSERT_EQUALS(p_element_4->GetNumFaces(), 8u);
        TS_ASSERT_EQUALS(p_element_4->GetFace(0)->GetIndex(), 4u);
        TS_ASSERT_EQUALS(p_element_4->GetFace(1)->GetIndex(), 13u);
        TS_ASSERT_EQUALS(p_element_4->GetFace(2)->GetIndex(), 33u);
        TS_ASSERT_EQUALS(p_element_4->GetFace(3)->GetIndex(), 38u);
        TS_ASSERT_EQUALS(p_element_4->GetFace(4)->GetIndex(), 39u);
        TS_ASSERT_EQUALS(p_element_4->GetFace(5)->GetIndex(), 40u);
        TS_ASSERT_EQUALS(p_element_4->GetFace(6)->GetIndex(), 34u);
        TS_ASSERT_EQUALS(p_element_4->GetFace(7)->GetIndex(), 27u);

        // Test the element 5 has the correct nodes and faces
        VertexElement<3,3>* p_element_5 = p_mesh->GetElement(5);

        TS_ASSERT_EQUALS(p_element_5->GetNumNodes(), 12u);
        TS_ASSERT_EQUALS(p_element_5->GetNode(0)->GetIndex(), 10u);
        TS_ASSERT_EQUALS(p_element_5->GetNode(1)->GetIndex(), 40u);
        TS_ASSERT_EQUALS(p_element_5->GetNode(2)->GetIndex(), 13u);
        TS_ASSERT_EQUALS(p_element_5->GetNode(3)->GetIndex(), 43u);
        TS_ASSERT_EQUALS(p_element_5->GetNode(4)->GetIndex(), 14u);
        TS_ASSERT_EQUALS(p_element_5->GetNode(5)->GetIndex(), 44u);
        TS_ASSERT_EQUALS(p_element_5->GetNode(6)->GetIndex(), 17u);
        TS_ASSERT_EQUALS(p_element_5->GetNode(7)->GetIndex(), 47u);
        TS_ASSERT_EQUALS(p_element_5->GetNode(8)->GetIndex(), 18u);
        TS_ASSERT_EQUALS(p_element_5->GetNode(9)->GetIndex(), 48u);
        TS_ASSERT_EQUALS(p_element_5->GetNode(10)->GetIndex(), 22u);
        TS_ASSERT_EQUALS(p_element_5->GetNode(11)->GetIndex(), 52u);

        TS_ASSERT_EQUALS(p_element_5->GetNumFaces(), 8u);
        TS_ASSERT_EQUALS(p_element_5->GetFace(0)->GetIndex(), 5u);
        TS_ASSERT_EQUALS(p_element_5->GetFace(1)->GetIndex(), 14u);
        TS_ASSERT_EQUALS(p_element_5->GetFace(2)->GetIndex(), 41u);
        TS_ASSERT_EQUALS(p_element_5->GetFace(3)->GetIndex(), 42u);
        TS_ASSERT_EQUALS(p_element_5->GetFace(4)->GetIndex(), 43u);
        TS_ASSERT_EQUALS(p_element_5->GetFace(5)->GetIndex(), 44u);
        TS_ASSERT_EQUALS(p_element_5->GetFace(6)->GetIndex(), 38u);
        TS_ASSERT_EQUALS(p_element_5->GetFace(7)->GetIndex(), 32u);

        // Test the element 6 has the correct nodes and faces
        VertexElement<3,3>* p_element_6 = p_mesh->GetElement(6);

        TS_ASSERT_EQUALS(p_element_6->GetNumNodes(), 12u);
        TS_ASSERT_EQUALS(p_element_6->GetNode(0)->GetIndex(), 15u);
        TS_ASSERT_EQUALS(p_element_6->GetNode(1)->GetIndex(), 45u);
        TS_ASSERT_EQUALS(p_element_6->GetNode(2)->GetIndex(), 19u);
        TS_ASSERT_EQUALS(p_element_6->GetNode(3)->GetIndex(), 49u);
        TS_ASSERT_EQUALS(p_element_6->GetNode(4)->GetIndex(), 20u);
        TS_ASSERT_EQUALS(p_element_6->GetNode(5)->GetIndex(), 50u);
        TS_ASSERT_EQUALS(p_element_6->GetNode(6)->GetIndex(), 23u);
        TS_ASSERT_EQUALS(p_element_6->GetNode(7)->GetIndex(), 53u);
        TS_ASSERT_EQUALS(p_element_6->GetNode(8)->GetIndex(), 24u);
        TS_ASSERT_EQUALS(p_element_6->GetNode(9)->GetIndex(), 54u);
        TS_ASSERT_EQUALS(p_element_6->GetNode(10)->GetIndex(), 27u);
        TS_ASSERT_EQUALS(p_element_6->GetNode(11)->GetIndex(), 57u);

        TS_ASSERT_EQUALS(p_element_6->GetNumFaces(), 8u);
        TS_ASSERT_EQUALS(p_element_6->GetFace(0)->GetIndex(), 6u);
        TS_ASSERT_EQUALS(p_element_6->GetFace(1)->GetIndex(), 15u);
        TS_ASSERT_EQUALS(p_element_6->GetFace(2)->GetIndex(), 36u);
        TS_ASSERT_EQUALS(p_element_6->GetFace(3)->GetIndex(), 45u);
        TS_ASSERT_EQUALS(p_element_6->GetFace(4)->GetIndex(), 46u);
        TS_ASSERT_EQUALS(p_element_6->GetFace(5)->GetIndex(), 47u);
        TS_ASSERT_EQUALS(p_element_6->GetFace(6)->GetIndex(), 48u);
        TS_ASSERT_EQUALS(p_element_6->GetFace(7)->GetIndex(), 49u);

        // Test the element 7 has the correct nodes and faces
        VertexElement<3,3>* p_element_7 = p_mesh->GetElement(7);

        TS_ASSERT_EQUALS(p_element_7->GetNumNodes(), 12u);
        TS_ASSERT_EQUALS(p_element_7->GetNode(0)->GetIndex(), 16u);
        TS_ASSERT_EQUALS(p_element_7->GetNode(1)->GetIndex(), 46u);
        TS_ASSERT_EQUALS(p_element_7->GetNode(2)->GetIndex(), 20u);
        TS_ASSERT_EQUALS(p_element_7->GetNode(3)->GetIndex(), 50u);
        TS_ASSERT_EQUALS(p_element_7->GetNode(4)->GetIndex(), 21u);
        TS_ASSERT_EQUALS(p_element_7->GetNode(5)->GetIndex(), 51u);
        TS_ASSERT_EQUALS(p_element_7->GetNode(6)->GetIndex(), 24u);
        TS_ASSERT_EQUALS(p_element_7->GetNode(7)->GetIndex(), 54u);
        TS_ASSERT_EQUALS(p_element_7->GetNode(8)->GetIndex(), 25u);
        TS_ASSERT_EQUALS(p_element_7->GetNode(9)->GetIndex(), 55u);
        TS_ASSERT_EQUALS(p_element_7->GetNode(10)->GetIndex(), 28u);
        TS_ASSERT_EQUALS(p_element_7->GetNode(11)->GetIndex(), 58u);

        TS_ASSERT_EQUALS(p_element_7->GetNumFaces(), 8u);
        TS_ASSERT_EQUALS(p_element_7->GetFace(0)->GetIndex(), 7u);
        TS_ASSERT_EQUALS(p_element_7->GetFace(1)->GetIndex(), 16u);
        TS_ASSERT_EQUALS(p_element_7->GetFace(2)->GetIndex(), 40u);
        TS_ASSERT_EQUALS(p_element_7->GetFace(3)->GetIndex(), 50u);
        TS_ASSERT_EQUALS(p_element_7->GetFace(4)->GetIndex(), 51u);
        TS_ASSERT_EQUALS(p_element_7->GetFace(5)->GetIndex(), 52u);
        TS_ASSERT_EQUALS(p_element_7->GetFace(6)->GetIndex(), 45u);
        TS_ASSERT_EQUALS(p_element_7->GetFace(7)->GetIndex(), 35u);

        // Test the element 8 has the correct nodes and faces
        VertexElement<3,3>* p_element_8 = p_mesh->GetElement(8);

        TS_ASSERT_EQUALS(p_element_8->GetNumNodes(), 12u);
        TS_ASSERT_EQUALS(p_element_8->GetNode(0)->GetIndex(), 17u);
        TS_ASSERT_EQUALS(p_element_8->GetNode(1)->GetIndex(), 47u);
        TS_ASSERT_EQUALS(p_element_8->GetNode(2)->GetIndex(), 21u);
        TS_ASSERT_EQUALS(p_element_8->GetNode(3)->GetIndex(), 51u);
        TS_ASSERT_EQUALS(p_element_8->GetNode(4)->GetIndex(), 22u);
        TS_ASSERT_EQUALS(p_element_8->GetNode(5)->GetIndex(), 52u);
        TS_ASSERT_EQUALS(p_element_8->GetNode(6)->GetIndex(), 25u);
        TS_ASSERT_EQUALS(p_element_8->GetNode(7)->GetIndex(), 55u);
        TS_ASSERT_EQUALS(p_element_8->GetNode(8)->GetIndex(), 26u);
        TS_ASSERT_EQUALS(p_element_8->GetNode(9)->GetIndex(), 56u);
        TS_ASSERT_EQUALS(p_element_8->GetNode(10)->GetIndex(), 29u);
        TS_ASSERT_EQUALS(p_element_8->GetNode(11)->GetIndex(), 59u);

        TS_ASSERT_EQUALS(p_element_8->GetNumFaces(), 8u);
        TS_ASSERT_EQUALS(p_element_8->GetFace(0)->GetIndex(), 8u);
        TS_ASSERT_EQUALS(p_element_8->GetFace(1)->GetIndex(), 17u);
        TS_ASSERT_EQUALS(p_element_8->GetFace(2)->GetIndex(), 44u);
        TS_ASSERT_EQUALS(p_element_8->GetFace(3)->GetIndex(), 53u);
        TS_ASSERT_EQUALS(p_element_8->GetFace(4)->GetIndex(), 54u);
        TS_ASSERT_EQUALS(p_element_8->GetFace(5)->GetIndex(), 55u);
        TS_ASSERT_EQUALS(p_element_8->GetFace(6)->GetIndex(), 50u);
        TS_ASSERT_EQUALS(p_element_8->GetFace(7)->GetIndex(), 39u);
    }

    void TestFourByFourElements() throw (Exception)
    {
        // Create a mesh
        HexagonalPrism3dVertexMeshGenerator generator(4, 4, 1.0, 2.0);
        MutableVertexMesh<3,3>* p_mesh = generator.GetMesh();

        // Test that the mesh has the correct number of nodes, faces and elements
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 96u);
        TS_ASSERT_EQUALS(p_mesh->GetNumFaces(), 95u);
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 16u);

        // Test that each face has the correct number of nodes
        for (unsigned index=0; index<p_mesh->GetNumFaces(); index++)
        {
            unsigned face_index = p_mesh->GetFace(index)->GetIndex();
            if (face_index < 32)
            {
                TS_ASSERT_EQUALS(p_mesh->GetFace(index)->GetNumNodes(), 6u);
            }
            else
            {
                TS_ASSERT_EQUALS(p_mesh->GetFace(index)->GetNumNodes(), 4u);
            }
        }

        // Test that each element has the correct number of nodes and faces
//        for (unsigned elem_index=0; elem_index<p_mesh->GetNumElements(); elem_index++)
//        {
//            if (p_mesh->GetElement(elem_index)->GetNumNodes() != 12)
//            {
//                PRINT_VARIABLE(elem_index);
//            }
//            TS_ASSERT_EQUALS(p_mesh->GetElement(elem_index)->GetNumNodes(), 12u);
//            TS_ASSERT_EQUALS(p_mesh->GetElement(elem_index)->GetNumFaces(), 8u);
//        }

        // Test the element 10 has the correct nodes and faces
        VertexElement<3,3>* p_element_10 = p_mesh->GetElement(10);

        TS_ASSERT_EQUALS(p_element_10->GetNumNodes(), 12u);
        TS_ASSERT_EQUALS(p_element_10->GetNode(0)->GetIndex(), 21u);
        TS_ASSERT_EQUALS(p_element_10->GetNode(1)->GetIndex(), 69u);
        TS_ASSERT_EQUALS(p_element_10->GetNode(2)->GetIndex(), 26u);
        TS_ASSERT_EQUALS(p_element_10->GetNode(3)->GetIndex(), 74u);
        TS_ASSERT_EQUALS(p_element_10->GetNode(4)->GetIndex(), 27u);
        TS_ASSERT_EQUALS(p_element_10->GetNode(5)->GetIndex(), 75u);
        TS_ASSERT_EQUALS(p_element_10->GetNode(6)->GetIndex(), 31u);
        TS_ASSERT_EQUALS(p_element_10->GetNode(7)->GetIndex(), 79u);
        TS_ASSERT_EQUALS(p_element_10->GetNode(8)->GetIndex(), 32u);
        TS_ASSERT_EQUALS(p_element_10->GetNode(9)->GetIndex(), 80u);
        TS_ASSERT_EQUALS(p_element_10->GetNode(10)->GetIndex(), 36u);
        TS_ASSERT_EQUALS(p_element_10->GetNode(11)->GetIndex(), 84u);

        TS_ASSERT_EQUALS(p_element_10->GetNumFaces(), 8u);
        TS_ASSERT_EQUALS(p_element_10->GetFace(0)->GetIndex(), 10u);
        TS_ASSERT_EQUALS(p_element_10->GetFace(1)->GetIndex(), 26u);
        TS_ASSERT_EQUALS(p_element_10->GetFace(2)->GetIndex(), 62u);
        TS_ASSERT_EQUALS(p_element_10->GetFace(3)->GetIndex(), 75u);
        TS_ASSERT_EQUALS(p_element_10->GetFace(4)->GetIndex(), 76u);
        TS_ASSERT_EQUALS(p_element_10->GetFace(5)->GetIndex(), 77u);
        TS_ASSERT_EQUALS(p_element_10->GetFace(6)->GetIndex(), 72u);
        TS_ASSERT_EQUALS(p_element_10->GetFace(7)->GetIndex(), 58u);

        // Test the element 14 has the correct nodes and faces
        VertexElement<3,3>* p_element_14 = p_mesh->GetElement(14);

        TS_ASSERT_EQUALS(p_element_14->GetNumNodes(), 12u);
        TS_ASSERT_EQUALS(p_element_14->GetNode(0)->GetIndex(), 32u);
        TS_ASSERT_EQUALS(p_element_14->GetNode(1)->GetIndex(), 80u);
        TS_ASSERT_EQUALS(p_element_14->GetNode(2)->GetIndex(), 36u);
        TS_ASSERT_EQUALS(p_element_14->GetNode(3)->GetIndex(), 84u);
        TS_ASSERT_EQUALS(p_element_14->GetNode(4)->GetIndex(), 37u);
        TS_ASSERT_EQUALS(p_element_14->GetNode(5)->GetIndex(), 85u);
        TS_ASSERT_EQUALS(p_element_14->GetNode(6)->GetIndex(), 41u);
        TS_ASSERT_EQUALS(p_element_14->GetNode(7)->GetIndex(), 89u);
        TS_ASSERT_EQUALS(p_element_14->GetNode(8)->GetIndex(), 42u);
        TS_ASSERT_EQUALS(p_element_14->GetNode(9)->GetIndex(), 90u);
        TS_ASSERT_EQUALS(p_element_14->GetNode(10)->GetIndex(), 46u);
        TS_ASSERT_EQUALS(p_element_14->GetNode(11)->GetIndex(), 94u);

        TS_ASSERT_EQUALS(p_element_14->GetNumFaces(), 8u);
        TS_ASSERT_EQUALS(p_element_14->GetFace(0)->GetIndex(), 14u);
        TS_ASSERT_EQUALS(p_element_14->GetFace(1)->GetIndex(), 30u);
        TS_ASSERT_EQUALS(p_element_14->GetFace(2)->GetIndex(), 80u);
        TS_ASSERT_EQUALS(p_element_14->GetFace(3)->GetIndex(), 88u);
        TS_ASSERT_EQUALS(p_element_14->GetFace(4)->GetIndex(), 89u);
        TS_ASSERT_EQUALS(p_element_14->GetFace(5)->GetIndex(), 90u);
        TS_ASSERT_EQUALS(p_element_14->GetFace(6)->GetIndex(), 85u);
        TS_ASSERT_EQUALS(p_element_14->GetFace(7)->GetIndex(), 76u);
    }

    void TestLargeMesh() throw (Exception)
    {
        // Create a mesh
        unsigned num_rows = 5;
        unsigned num_columns = 5;
        HexagonalPrism3dVertexMeshGenerator generator(5, 5, 1.0, 2.0);
        MutableVertexMesh<3,3>* p_mesh = generator.GetMesh();

        // Test that the mesh has the correct number of nodes, faces and elements
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 140u);
        TS_ASSERT_EQUALS(p_mesh->GetNumFaces(), 144u);
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 25u);

        // Test that each face has the correct number of nodes
        for (unsigned index=0; index<p_mesh->GetNumFaces(); index++)
        {
            unsigned face_index = p_mesh->GetFace(index)->GetIndex();
            if (face_index < 2*num_rows*num_columns)
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
