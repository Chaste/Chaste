/*

Copyright (C) University of Oxford, 2005-2010

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

#ifndef TESTCYLINDRICALHONEYCOMBVERTEXMESHGENERATOR_HPP_
#define TESTCYLINDRICALHONEYCOMBVERTEXMESHGENERATOR_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "CylindricalHoneycombVertexMeshGenerator.hpp"

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

        TS_ASSERT_EQUALS(system(("diff -I \"Created by Chaste\" " + results_file1 + " notforrelease_cell_based/test/data/TestCylindrical2dVertexMesh/cylindrical_vertex_mesh.node").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff -I \"Created by Chaste\" " + results_file2 + " notforrelease_cell_based/test/data/TestCylindrical2dVertexMesh/cylindrical_vertex_mesh.cell").c_str()), 0);

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

        TS_ASSERT_EQUALS(system(("diff -I \"Created by Chaste\" " + results_file1 + " notforrelease_cell_based/test/data/TestFlatCylindrical2dVertexMesh/flat_cylindrical_vertex_mesh.node").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff -I \"Created by Chaste\" " + results_file2 + " notforrelease_cell_based/test/data/TestFlatCylindrical2dVertexMesh/flat_cylindrical_vertex_mesh.cell").c_str()), 0);
    }
};

#endif /*TESTCYLINDRICALHONEYCOMBVERTEXMESHGENERATOR_HPP_*/
