/*

Copyright (C) University of Oxford, 2005-2011

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

#ifndef TESTPOTTSWRITER_HPP_
#define TESTPOTTSWRITER_HPP_

#include <cxxtest/TestSuite.h>

#include <string>
#include <fstream>

#include "MutableMesh.hpp"
#include "PottsMesh.hpp"
#include "PottsElement.hpp"
#include "PottsMeshWriter.hpp"
#include "VtkMeshWriter.hpp"

class TestPottsWriter : public CxxTest::TestSuite
{
public:
    void TestPottsMeshWriterIn2d() throw (Exception)
    {
        // Make 6 nodes to assign to two elements
        std::vector<Node<2>*> basic_nodes;
        basic_nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        basic_nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
        basic_nodes.push_back(new Node<2>(2, false, 2.0, 0.0));
        basic_nodes.push_back(new Node<2>(3, false, 0.0, 1.0));
        basic_nodes.push_back(new Node<2>(4, false, 1.0, 1.0));
        basic_nodes.push_back(new Node<2>(5, false, 2.0, 1.0));

        // Make two triangular elements out of these nodes
        std::vector<std::vector<Node<2>*> > nodes_elements(2);
        nodes_elements[0].push_back(basic_nodes[0]);
        nodes_elements[0].push_back(basic_nodes[1]);
        nodes_elements[0].push_back(basic_nodes[3]);

        nodes_elements[1].push_back(basic_nodes[2]);
        nodes_elements[1].push_back(basic_nodes[4]);
        nodes_elements[1].push_back(basic_nodes[5]);

        std::vector<PottsElement<2>*> basic_potts_elements;
        basic_potts_elements.push_back(new PottsElement<2>(0, nodes_elements[0]));
        basic_potts_elements.push_back(new PottsElement<2>(1, nodes_elements[1]));

        // Make a PottsMesh
        PottsMesh<2> basic_potts_mesh(basic_nodes, basic_potts_elements);

        // Create a potts mesh writer
        PottsMeshWriter<2> potts_mesh_writer("TestPottsMeshWriterIn2d", "potts_mesh_2d");

        // Write and check it's correct
        potts_mesh_writer.WriteFilesUsingMesh(basic_potts_mesh);

        OutputFileHandler handler("TestPottsMeshWriterIn2d", false);
        std::string results_file1 = handler.GetOutputDirectoryFullPath() + "potts_mesh_2d.node";
        std::string results_file2 = handler.GetOutputDirectoryFullPath() + "potts_mesh_2d.cell";

        // To ignore the provenance data we only go as far as
        TS_ASSERT_EQUALS(system(("diff -I \"Created by Chaste\" " + results_file1 + " notforrelease_cell_based/test/data/TestPottsMeshWriter/potts_mesh_2d.node").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff -I \"Created by Chaste\" " + results_file2 + " notforrelease_cell_based/test/data/TestPottsMeshWriter/potts_mesh_2d.cell").c_str()), 0);

    }
};

#endif // TESTPOTTSWRITER_HPP_
