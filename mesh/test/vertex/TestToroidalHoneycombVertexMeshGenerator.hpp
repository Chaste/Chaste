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

#ifndef TESTTOROIDALHONEYCOMBVERTEXMESHGENERATOR_HPP_
#define TESTTOROIDALHONEYCOMBVERTEXMESHGENERATOR_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "ToroidalHoneycombVertexMeshGenerator.hpp"
#include "FileComparison.hpp"

//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

class TestToroidalHoneycombVertexMeshGenerator : public CxxTest::TestSuite
{
public:

    void TestToroidal2dVertexMeshGenerator()
    {
        ToroidalHoneycombVertexMeshGenerator generator(4, 4);

        // Coverage
        TS_ASSERT_THROWS_THIS(generator.GetMesh(),
                              "A toroidal mesh was created but a normal mesh is being requested.");

        // Create periodic mesh
        Toroidal2dVertexMesh* p_mesh = generator.GetToroidalMesh();

        // Test that the mesh has the correct numbers of nodes and elements
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 16u);
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 32u);

        // Test that some elements own the correct nodes
        VertexElement<2,2>* p_element0 = p_mesh->GetElement(0);
        TS_ASSERT_EQUALS(p_element0->GetNumNodes(), 6u);
        TS_ASSERT_EQUALS(p_element0->GetNodeGlobalIndex(0), 0u);
        TS_ASSERT_EQUALS(p_element0->GetNodeGlobalIndex(1), 5u);
        TS_ASSERT_EQUALS(p_element0->GetNodeGlobalIndex(2), 9u);
        TS_ASSERT_EQUALS(p_element0->GetNodeGlobalIndex(3), 12u);
        TS_ASSERT_EQUALS(p_element0->GetNodeGlobalIndex(4), 8u);
        TS_ASSERT_EQUALS(p_element0->GetNodeGlobalIndex(5), 4u);

        VertexElement<2,2>* p_element5 = p_mesh->GetElement(5);
        TS_ASSERT_EQUALS(p_element5->GetNumNodes(), 6u);
        TS_ASSERT_EQUALS(p_element5->GetNodeGlobalIndex(0), 10u);
        TS_ASSERT_EQUALS(p_element5->GetNodeGlobalIndex(1), 14u);
        TS_ASSERT_EQUALS(p_element5->GetNodeGlobalIndex(2), 18u);
        TS_ASSERT_EQUALS(p_element5->GetNodeGlobalIndex(3), 22u);
        TS_ASSERT_EQUALS(p_element5->GetNodeGlobalIndex(4), 17u);
        TS_ASSERT_EQUALS(p_element5->GetNodeGlobalIndex(5), 13u);

        VertexElement<2,2>* p_element7 = p_mesh->GetElement(7);
        TS_ASSERT_EQUALS(p_element7->GetNumNodes(), 6u);
        TS_ASSERT_EQUALS(p_element7->GetNodeGlobalIndex(0), 8u);
        TS_ASSERT_EQUALS(p_element7->GetNodeGlobalIndex(1), 12u);
        TS_ASSERT_EQUALS(p_element7->GetNodeGlobalIndex(2), 16u);
        TS_ASSERT_EQUALS(p_element7->GetNodeGlobalIndex(3), 20u);
        TS_ASSERT_EQUALS(p_element7->GetNodeGlobalIndex(4), 19u);
        TS_ASSERT_EQUALS(p_element7->GetNodeGlobalIndex(5), 15u);

        VertexElement<2,2>* p_element12 = p_mesh->GetElement(12);
        TS_ASSERT_EQUALS(p_element12->GetNumNodes(), 6u);
        TS_ASSERT_EQUALS(p_element12->GetNodeGlobalIndex(0), 25u);
        TS_ASSERT_EQUALS(p_element12->GetNodeGlobalIndex(1), 29u);
        TS_ASSERT_EQUALS(p_element12->GetNodeGlobalIndex(2), 1u);
        TS_ASSERT_EQUALS(p_element12->GetNodeGlobalIndex(3), 5u);
        TS_ASSERT_EQUALS(p_element12->GetNodeGlobalIndex(4), 0u);
        TS_ASSERT_EQUALS(p_element12->GetNodeGlobalIndex(5), 28u);

        VertexElement<2,2>* p_element15 = p_mesh->GetElement(15);
        TS_ASSERT_EQUALS(p_element15->GetNumNodes(), 6u);
        TS_ASSERT_EQUALS(p_element15->GetNodeGlobalIndex(0), 24u);
        TS_ASSERT_EQUALS(p_element15->GetNodeGlobalIndex(1), 28u);
        TS_ASSERT_EQUALS(p_element15->GetNodeGlobalIndex(2), 0u);
        TS_ASSERT_EQUALS(p_element15->GetNodeGlobalIndex(3), 4u);
        TS_ASSERT_EQUALS(p_element15->GetNodeGlobalIndex(4), 3u);
        TS_ASSERT_EQUALS(p_element15->GetNodeGlobalIndex(5), 31u);

        // Create a vertex mesh writer with toroidal mesh
        VertexMeshWriter<2,2> vertex_mesh_writer_1("TestToroidal2dVertexMesh", "toroidal_vertex_mesh");
        vertex_mesh_writer_1.WriteFilesUsingMesh(*p_mesh);

        OutputFileHandler handler_1("TestToroidal2dVertexMesh", false);
        std::string results_file1 = handler_1.GetOutputDirectoryFullPath() + "toroidal_vertex_mesh.node";
        std::string results_file2 = handler_1.GetOutputDirectoryFullPath() + "toroidal_vertex_mesh.cell";

        {
            FileComparison comparer(results_file1, "mesh/test/data/TestToroidal2dVertexMesh/toroidal_vertex_mesh.node");
            TS_ASSERT(comparer.CompareFiles());
        }
        {
            FileComparison comparer(results_file2, "mesh/test/data/TestToroidal2dVertexMesh/toroidal_vertex_mesh.cell");
            TS_ASSERT(comparer.CompareFiles());
        }
    }
};

#endif /*TESTTOROIDALHONEYCOMBVERTEXMESHGENERATOR_HPP_*/
