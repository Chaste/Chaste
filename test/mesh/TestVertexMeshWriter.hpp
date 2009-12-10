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
#ifndef TESTVERTEXMESHWRITER_HPP_
#define TESTVERTEXMESHWRITER_HPP_

#include <cxxtest/TestSuite.h>

#include <string>
#include <fstream>

#include "HoneycombVertexMeshGenerator.hpp"
#include "VertexMeshWriter.hpp"
#include "VtkWriter.hpp"
#include "VertexMeshReader.hpp"

class TestVertexMeshWriter : public CxxTest::TestSuite
{
public:

    void TestMeshWriter() throw(Exception)
    {
        std::vector<Node<2>*> basic_nodes;
        basic_nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        basic_nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
        basic_nodes.push_back(new Node<2>(2, false, 1.5, 1.0));
        basic_nodes.push_back(new Node<2>(3, false, 1.0, 2.0));
        basic_nodes.push_back(new Node<2>(4, false, 0.0, 1.0));
        basic_nodes.push_back(new Node<2>(5, false, 2.0, 0.0));
        basic_nodes.push_back(new Node<2>(6, false, 2.0, 3.0));

        std::vector<Node<2>*> nodes_elem_0, nodes_elem_1;

        // Make two elements out of these nodes
        nodes_elem_0.push_back(basic_nodes[0]);
        nodes_elem_0.push_back(basic_nodes[1]);
        nodes_elem_0.push_back(basic_nodes[2]);
        nodes_elem_0.push_back(basic_nodes[3]);
        nodes_elem_0.push_back(basic_nodes[4]);

        nodes_elem_1.push_back(basic_nodes[2]);
        nodes_elem_1.push_back(basic_nodes[5]);
        nodes_elem_1.push_back(basic_nodes[6]);

        std::vector<VertexElement<2,2>*> basic_vertex_elements;
        basic_vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem_0));
        basic_vertex_elements.push_back(new VertexElement<2,2>(1, nodes_elem_1));

        // Make a vertex mesh
        VertexMesh<2,2> basic_vertex_mesh(basic_nodes, basic_vertex_elements);

        // Create a vertex mesh writer
        VertexMeshWriter<2,2> vertex_mesh_writer("TestVertexMeshWriter", "vertex_mesh");
        vertex_mesh_writer.WriteFilesUsingMesh(basic_vertex_mesh);

        OutputFileHandler handler("TestVertexMeshWriter", false);
        std::string results_file1 = handler.GetOutputDirectoryFullPath() + "vertex_mesh.node";
        std::string results_file2 = handler.GetOutputDirectoryFullPath() + "vertex_mesh.cell";

        TS_ASSERT_EQUALS(system(("diff " + results_file1 + " notforrelease_cell_based/test/data/TestVertexMesh/vertex_mesh.node").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_file2 + " notforrelease_cell_based/test/data/TestVertexMesh/vertex_mesh.cell").c_str()), 0);

#ifdef CHASTE_VTK
        std::vector<double> cell_ids;
        cell_ids.push_back(0.0);
        cell_ids.push_back(1.0);

        vertex_mesh_writer.AddCellData("Cell IDs", cell_ids);
        // Add distance from origin into the node "point" data
        std::vector<double> distance;
        for (unsigned i=0; i<basic_vertex_mesh.GetNumNodes(); i++)
        {
            distance.push_back(norm_2(basic_vertex_mesh.GetNode(i)->rGetLocation()));
        }
        vertex_mesh_writer.AddPointData("Distance from origin", distance);

        vertex_mesh_writer.WriteVtkUsingMesh(basic_vertex_mesh);

        //1.5K uncompressed, 1.5K compressed
        std::string results_file3 = handler.GetOutputDirectoryFullPath() + "vertex_mesh.vtu";
        TS_ASSERT_EQUALS(system(("cmp  " + results_file3 + " notforrelease_cell_based/test/data/TestVertexMesh/vertex_mesh.vtu").c_str()), 0);

#endif //CHASTE_VTK
    }

    void TestMeshWriterWithDeletedNode() throw (Exception)
    {
        // Create mesh
        HoneycombVertexMeshGenerator generator(3, 3, false, false, 0.1, 2.0);
        VertexMesh<2,2>* p_mesh = generator.GetMesh();

        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 30u);
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 9u);

        /*
         * Delete element 0. This element contains 3 nodes that are
         * not contained in any other element and so will be marked
         * as deleted.
         */
        p_mesh->DeleteElementPriorToReMesh(0);

        // Write mesh to file
        VertexMeshWriter<2,2> mesh_writer("TestMeshWriterWithDeletedNode", "vertex_mesh");
        TS_ASSERT_THROWS_NOTHING(mesh_writer.WriteFilesUsingMesh(*p_mesh));

        // Read mesh back in from file
        std::string output_dir = mesh_writer.GetOutputDirectory();
        VertexMeshReader<2,2> mesh_reader(output_dir + "vertex_mesh");

        // We should have one less element and three less nodes
        TS_ASSERT_EQUALS(mesh_reader.GetNumNodes(), 27u);
        TS_ASSERT_EQUALS(mesh_reader.GetNumElements(), 8u);
    }

    void TestMeshVtkWriter3D() throw(Exception)
    {
        // Create 3D mesh
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, false, 0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, false, 1.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(2, false, 1.0, 2.0, 0.0));
        nodes.push_back(new Node<3>(3, false, 0.0, 2.0, 0.0));
        nodes.push_back(new Node<3>(4, false, 0.0, 2.0, 3.0));
        nodes.push_back(new Node<3>(5, false, 1.0, 0.0, 3.0));
        nodes.push_back(new Node<3>(6, false, 1.0, 2.0, 3.0));
        nodes.push_back(new Node<3>(7, false, 0.0, 2.0, 3.0));

        std::vector<VertexElement<3,3>*> elements;
        elements.push_back(new VertexElement<3,3>(0, nodes));

        VertexMesh<3,3> mesh3d(nodes, elements);
        TS_ASSERT_DELTA(mesh3d.GetWidth(0), 1.0, 1e-4);
        TS_ASSERT_DELTA(mesh3d.GetWidth(1), 2.0, 1e-4);
        TS_ASSERT_DELTA(mesh3d.GetWidth(2), 3.0, 1e-4);
        VertexMeshWriter<3,3> vertex_mesh_writer("TestVertexMeshWriter", "vertex_mesh_3d", false);

#ifdef CHASTE_VTK
        std::vector<double> cell_ids;
        cell_ids.push_back(0.0);
        vertex_mesh_writer.AddCellData("Cell IDs", cell_ids);

         // Add distance from origin into the node "point" data
        std::vector<double> distance;
        for (unsigned i=0; i<mesh3d.GetNumNodes(); i++)
        {
            distance.push_back(norm_2(mesh3d.GetNode(i)->rGetLocation()));
        }
        vertex_mesh_writer.AddPointData("Distance from origin", distance);

        vertex_mesh_writer.WriteVtkUsingMesh(mesh3d, "42");

        OutputFileHandler handler("TestVertexMeshWriter", false);

        //1.5K uncompressed, 1.5K compressed
        std::string results_file3 = handler.GetOutputDirectoryFullPath() + "vertex_mesh_3d_42.vtu";
        TS_ASSERT_EQUALS(system(("cmp  " + results_file3 + " notforrelease_cell_based/test/data/TestVertexMesh/vertex_mesh_3d.vtu").c_str()), 0);
#endif //CHASTE_VTK
    }

};


#endif /*TESTVERTEXMESHWRITER_HPP_*/
