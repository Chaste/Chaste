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
#ifndef TESTVERTEXMESHWRITERFORMUTABLE_HPP_
#define TESTVERTEXMESHWRITERFORMUTABLE_HPP_

#include <cxxtest/TestSuite.h>

#include <string>
#include <fstream>

#include "HoneycombMutableVertexMeshGenerator.hpp"
#include "VertexMeshWriter.hpp"
#include "VtkMeshWriter.hpp"
#include "VertexMeshReader.hpp"

class TestVertexMeshWriterForMutable : public CxxTest::TestSuite
{
public:
    void TestMeshWriterWithDeletedNode() throw (Exception)
    {
        // Create mesh
        HoneycombMutableVertexMeshGenerator generator(3, 3, false, 0.1, 2.0);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMutableMesh();

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
};


#endif /*TESTVERTEXMESHWRITERFORMUTABLE_HPP_*/
