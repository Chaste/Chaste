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

#ifndef TESTCYLINDRICALHONEYCOMBVERTEXMESHGENERATOR_HPP_
#define TESTCYLINDRICALHONEYCOMBVERTEXMESHGENERATOR_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "CylindricalHoneycombVertexMeshGenerator.hpp"
#include "FileComparison.hpp"

//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

class TestCylindricalHoneycombVertexMeshGenerator : public CxxTest::TestSuite
{
public:

    void TestCylindrical2dVertexMeshGenerator()
    {
        CylindricalHoneycombVertexMeshGenerator generator(4, 4);

        // Coverage
        TS_ASSERT_THROWS_THIS(generator.GetMesh(),
                              "A cylindrical mesh was created but a normal mesh is being requested.");

        // Create periodic mesh
        Cylindrical2dVertexMesh* p_cylindrical_mesh = generator.GetCylindricalMesh();

        // The periodic mesh should have the same number of elements but fewer nodes
        TS_ASSERT_EQUALS(p_cylindrical_mesh->GetNumElements(), 16u);
        TS_ASSERT_EQUALS(p_cylindrical_mesh->GetNumNodes(), 40u);

        // Create a vertex mesh writer with cylindrical mesh
        VertexMeshWriter<2,2> vertex_mesh_writer_1("TestCylindrical2dVertexMesh", "cylindrical_vertex_mesh");
        vertex_mesh_writer_1.WriteFilesUsingMesh(*p_cylindrical_mesh);

        OutputFileHandler handler_1("TestCylindrical2dVertexMesh", false);
        std::string results_file1 = handler_1.GetOutputDirectoryFullPath() + "cylindrical_vertex_mesh.node";
        std::string results_file2 = handler_1.GetOutputDirectoryFullPath() + "cylindrical_vertex_mesh.cell";

        {
            FileComparison comparer(results_file1, "mesh/test/data/TestCylindrical2dVertexMesh/cylindrical_vertex_mesh.node");
            TS_ASSERT(comparer.CompareFiles());
        }
        {
            FileComparison comparer(results_file2, "mesh/test/data/TestCylindrical2dVertexMesh/cylindrical_vertex_mesh.cell");
            TS_ASSERT(comparer.CompareFiles());
        }

        // Create periodic mesh with flat bottom
        CylindricalHoneycombVertexMeshGenerator generator3(4, 4, true);
        Cylindrical2dVertexMesh* p_flat_cylindrical_mesh = generator3.GetCylindricalMesh();

        // The flat bottomed periodic mesh should have the same number of elements and nodes
        TS_ASSERT_EQUALS(p_flat_cylindrical_mesh->GetNumElements(), 16u);
        TS_ASSERT_EQUALS(p_flat_cylindrical_mesh->GetNumNodes(), 36u);

        // The bottom row of nodes should be boundary nodes
        TS_ASSERT_EQUALS(p_flat_cylindrical_mesh->GetNode(0)->IsBoundaryNode(), true);
        TS_ASSERT_EQUALS(p_flat_cylindrical_mesh->GetNode(1)->IsBoundaryNode(), true);
        TS_ASSERT_EQUALS(p_flat_cylindrical_mesh->GetNode(2)->IsBoundaryNode(), true);
        TS_ASSERT_EQUALS(p_flat_cylindrical_mesh->GetNode(3)->IsBoundaryNode(), true);
        TS_ASSERT_EQUALS(p_flat_cylindrical_mesh->GetNode(4)->IsBoundaryNode(), false);

        // Create a vertex mesh writer with cylindrical mesh
        VertexMeshWriter<2,2> vertex_mesh_writer_2("TestFlatCylindrical2dVertexMesh", "flat_cylindrical_vertex_mesh");
        vertex_mesh_writer_2.WriteFilesUsingMesh(*p_flat_cylindrical_mesh);

        OutputFileHandler handler_2("TestFlatCylindrical2dVertexMesh", false);
        results_file1 = handler_2.GetOutputDirectoryFullPath() + "flat_cylindrical_vertex_mesh.node";
        results_file2 = handler_2.GetOutputDirectoryFullPath() + "flat_cylindrical_vertex_mesh.cell";
        {
            FileComparison comparer(results_file1, "mesh/test/data/TestFlatCylindrical2dVertexMesh/flat_cylindrical_vertex_mesh.node");
            TS_ASSERT(comparer.CompareFiles());
        }
        {
            FileComparison comparer(results_file2, "mesh/test/data/TestFlatCylindrical2dVertexMesh/flat_cylindrical_vertex_mesh.cell");
            TS_ASSERT(comparer.CompareFiles());
        }
    }
};

#endif /*TESTCYLINDRICALHONEYCOMBVERTEXMESHGENERATOR_HPP_*/
