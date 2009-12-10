/*

Copyright (C) University of Oxford, 2005-2009

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

#ifndef TESTHONEYCOMBVERTEXMESHGENERATOR_HPP_
#define TESTHONEYCOMBVERTEXMESHGENERATOR_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "HoneycombVertexMeshGenerator.hpp"

class TestHoneycombVertexMeshGenerator : public CxxTest::TestSuite
{
public:

    void TestVertexMeshGenerator() throw(Exception)
    {
        HoneycombVertexMeshGenerator generator(5, 3);

        // Coverage
        TS_ASSERT_THROWS_THIS(generator.GetCylindricalMesh(),
                              "A normal mesh was created but a cylindrical mesh is being requested.");

        // Create mesh
        VertexMesh<2,2>* p_mesh = generator.GetMesh();

        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 15u);
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 46u);

        // Test some random nodes are in the correct place
        TS_ASSERT_DELTA(p_mesh->GetNode(5)->rGetLocation()[0], 0.00000, 1e-3);
        TS_ASSERT_DELTA(p_mesh->GetNode(5)->rGetLocation()[1], 0.28867, 1e-3);

        TS_ASSERT_DELTA(p_mesh->GetNode(43)->rGetLocation()[0], 2.5, 1e-3);
        TS_ASSERT_DELTA(p_mesh->GetNode(43)->rGetLocation()[1], 2.8867, 1e-3);

        // Check that each node is contained in at least one element
        for (unsigned node_index=0; node_index<p_mesh->GetNumNodes(); node_index++)
        {
            std::set<unsigned> containing_elements = p_mesh->GetNode(node_index)->rGetContainingElementIndices();
            unsigned num_containing_elements = containing_elements.size();

            TS_ASSERT_LESS_THAN(0u, num_containing_elements);
        }

        // Check that some elements contain the correct nodes

        // Element 2 contains nodes 2, 8, 14, 19, 13 and 7
        TS_ASSERT_EQUALS(p_mesh->GetElement(2)->GetNumNodes(), 6u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(2)->GetNodeGlobalIndex(0), 2u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(2)->GetNodeGlobalIndex(1), 8u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(2)->GetNodeGlobalIndex(2), 14u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(2)->GetNodeGlobalIndex(3), 19u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(2)->GetNodeGlobalIndex(4), 13u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(2)->GetNodeGlobalIndex(5), 7u);

        // Element 9 contains nodes 16, 22, 28, 34, 27 and 21
        TS_ASSERT_EQUALS(p_mesh->GetElement(9)->GetNumNodes(), 6u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(9)->GetNodeGlobalIndex(0), 16u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(9)->GetNodeGlobalIndex(1), 22u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(9)->GetNodeGlobalIndex(2), 28u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(9)->GetNodeGlobalIndex(3), 34u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(9)->GetNodeGlobalIndex(4), 27u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(9)->GetNodeGlobalIndex(5), 21u);

        // Element 10 contains nodes 23, 30, 36, 41, 35 and 29
        TS_ASSERT_EQUALS(p_mesh->GetElement(10)->GetNumNodes(), 6u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(10)->GetNodeGlobalIndex(0), 23u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(10)->GetNodeGlobalIndex(1), 30u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(10)->GetNodeGlobalIndex(2), 36u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(10)->GetNodeGlobalIndex(3), 41u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(10)->GetNodeGlobalIndex(4), 35u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(10)->GetNodeGlobalIndex(5), 29u);
  
        // Element 14 contains nodes 27, 34, 40, 45, 39 and 33
        TS_ASSERT_EQUALS(p_mesh->GetElement(14)->GetNumNodes(), 6u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(14)->GetNodeGlobalIndex(0), 27u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(14)->GetNodeGlobalIndex(1), 34u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(14)->GetNodeGlobalIndex(2), 40u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(14)->GetNodeGlobalIndex(3), 45u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(14)->GetNodeGlobalIndex(4), 39u);
        TS_ASSERT_EQUALS(p_mesh->GetElement(14)->GetNodeGlobalIndex(5), 33u);

        // Check that some nodes are contained in the correct elements

        // Node 0 is only in element 3
        std::set<unsigned> temp_list;
        temp_list.insert(0u);
        TS_ASSERT_EQUALS(p_mesh->GetNode(0)->rGetContainingElementIndices(), temp_list);

        // Node 6 is in elements 0 and 1
        temp_list.insert(1u);
        TS_ASSERT_EQUALS(p_mesh->GetNode(6)->rGetContainingElementIndices(), temp_list);

        // Node 12 is in elements 0 and 1 and 5
        temp_list.insert(5u);
        TS_ASSERT_EQUALS(p_mesh->GetNode(12)->rGetContainingElementIndices(), temp_list);
    }

    void TestCylindrical2dVertexMeshGenerator()
    {
        HoneycombVertexMeshGenerator generator(4, 4, true);

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

        TS_ASSERT_EQUALS(system(("diff " + results_file1 + " notforrelease_cell_based/test/data/TestCylindrical2dVertexMesh/cylindrical_vertex_mesh.node").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_file2 + " notforrelease_cell_based/test/data/TestCylindrical2dVertexMesh/cylindrical_vertex_mesh.cell").c_str()), 0);

        // Create periodic mesh with flat bottom
        HoneycombVertexMeshGenerator generator3(4, 4, true, true);
        Cylindrical2dVertexMesh* p_flat_cylindrical_mesh = generator3.GetCylindricalMesh();

        // The flat bottomed periodic mesh should have the same number of elements and nodes
        TS_ASSERT_EQUALS(p_flat_cylindrical_mesh->GetNumElements(), 16u);
        TS_ASSERT_EQUALS(p_flat_cylindrical_mesh->GetNumNodes(), 36u);

        // Create a vertex mesh writer with cylindrical mesh
        VertexMeshWriter<2,2> vertex_mesh_writer_2("TestFlatCylindrical2dVertexMesh", "flat_cylindrical_vertex_mesh");
        vertex_mesh_writer_2.WriteFilesUsingMesh(*p_flat_cylindrical_mesh);

        OutputFileHandler handler_2("TestFlatCylindrical2dVertexMesh", false);
        results_file1 = handler_2.GetOutputDirectoryFullPath() + "flat_cylindrical_vertex_mesh.node";
        results_file2 = handler_2.GetOutputDirectoryFullPath() + "flat_cylindrical_vertex_mesh.cell";

        TS_ASSERT_EQUALS(system(("diff " + results_file1 + " notforrelease_cell_based/test/data/TestFlatCylindrical2dVertexMesh/flat_cylindrical_vertex_mesh.node").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_file2 + " notforrelease_cell_based/test/data/TestFlatCylindrical2dVertexMesh/flat_cylindrical_vertex_mesh.cell").c_str()), 0);
    }
};
#endif /*TESTHONEYCOMBVERTEXMESHGENERATOR_HPP_*/
