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

#ifndef TESTPOTTSMESHWRITER_HPP_
#define TESTPOTTSMESHWRITER_HPP_

#include <cxxtest/TestSuite.h>

#include <string>
#include <fstream>

#include "MutableMesh.hpp"
#include "PottsMesh.hpp"
#include "PottsMeshGenerator.hpp"
#include "PottsElement.hpp"
#include "PottsMeshReader.hpp"
#include "PottsMeshWriter.hpp"
#include "VtkMeshWriter.hpp"
#include "FileComparison.hpp"
//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

class TestPottsMeshWriter : public CxxTest::TestSuite
{
public:

    void TestPottsMeshWriter2d()
    {
        // Create 2D mesh with 2 square elements
        PottsMeshGenerator<2> generator(4, 2, 2, 2, 1, 2);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create a Potts mesh writer
        PottsMeshWriter<2> potts_mesh_writer("TestPottsMeshWriter2d", "potts_mesh_2d");

        // Write and check it's correct
        potts_mesh_writer.WriteFilesUsingMesh(*p_mesh);

        OutputFileHandler handler("TestPottsMeshWriter2d", false);
        std::string results_file1 = handler.GetOutputDirectoryFullPath() + "potts_mesh_2d.node";
        std::string results_file2 = handler.GetOutputDirectoryFullPath() + "potts_mesh_2d.cell";

        {
            FileComparison comparer(results_file1,
                    "cell_based/test/data/TestPottsMeshWriter/potts_mesh_2d.node");
            TS_ASSERT(comparer.CompareFiles());
        }
        {
            FileComparison comparer(results_file2,
                    "cell_based/test/data/TestPottsMeshWriter/potts_mesh_2d.cell");
            TS_ASSERT(comparer.CompareFiles());
        }
    }

    void TestPottsMeshWriter3d()
    {
        // Create 3D mesh with 2 square elements
        PottsMeshGenerator<3> generator(2, 2, 1, 2, 1, 2, 2, 1, 2);
        PottsMesh<3>* p_mesh = generator.GetMesh();

        // Create a Potts mesh writer
        PottsMeshWriter<3> potts_mesh_writer("TestPottsMeshWriter3d", "potts_mesh_3d");

        // Write and check it's correct
        potts_mesh_writer.WriteFilesUsingMesh(*p_mesh);

        OutputFileHandler handler("TestPottsMeshWriter3d", false);
        std::string results_file1 = handler.GetOutputDirectoryFullPath() + "potts_mesh_3d.node";
        std::string results_file2 = handler.GetOutputDirectoryFullPath() + "potts_mesh_3d.cell";

        {
            FileComparison comparer(results_file1,
                    "cell_based/test/data/TestPottsMeshWriter/potts_mesh_3d.node");
            TS_ASSERT(comparer.CompareFiles());
        }
        {
            FileComparison comparer(results_file2,
                    "cell_based/test/data/TestPottsMeshWriter/potts_mesh_3d.cell");
            TS_ASSERT(comparer.CompareFiles());
        }
    }

    void TestReadingAndWritingElementAttributes()
    {
        // Read in a mesh with element attributes
        PottsMeshReader<2> mesh_reader("cell_based/test/data/TestPottsMeshReader2d/potts_mesh_with_element_attributes");
        TS_ASSERT_EQUALS(mesh_reader.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(mesh_reader.GetNumElementAttributes(), 1u);

        // Construct the mesh
        PottsMesh<2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetUnsignedAttribute(), 97u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetUnsignedAttribute(), 152u);

        // Write the mesh to file
        PottsMeshWriter<2> mesh_writer("TestReadingAndWritingElementAttributes", "potts_mesh_with_element_attributes");
        mesh_writer.WriteFilesUsingMesh(mesh);

        // Now read in the mesh that was written
        OutputFileHandler handler("TestReadingAndWritingElementAttributes", false);
        PottsMeshReader<2> mesh_reader2(handler.GetOutputDirectoryFullPath() + "potts_mesh_with_element_attributes");
        TS_ASSERT_EQUALS(mesh_reader2.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(mesh_reader2.GetNumElementAttributes(), 1u);

        // Construct the mesh again
        PottsMesh<2> mesh2;
        mesh2.ConstructFromMeshReader(mesh_reader);
        TS_ASSERT_EQUALS(mesh2.GetElement(0)->GetUnsignedAttribute(), 97u);
        TS_ASSERT_EQUALS(mesh2.GetElement(1)->GetUnsignedAttribute(), 152u);
    }

    void TestWriteFilesUsingMeshReader()
    {
        // Create a PottsMeshReader and use it to write mesh files
        PottsMeshReader<2> mesh_reader("cell_based/test/data/TestPottsMeshReader2d/potts_mesh_with_element_attributes");
        TS_ASSERT_EQUALS(mesh_reader.GetNumNodes(), 6u);
        TS_ASSERT_EQUALS(mesh_reader.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(mesh_reader.GetNumElementAttributes(), 1u);

        PottsMeshWriter<2> mesh_writer("TestWriteFilesUsingMeshReader", "potts_mesh");
        mesh_writer.WriteFilesUsingMeshReader(mesh_reader);

        // Now read in the mesh that was written
        OutputFileHandler handler("TestWriteFilesUsingMeshReader", false);
        PottsMeshReader<2> mesh_reader2(handler.GetOutputDirectoryFullPath() + "potts_mesh");
        TS_ASSERT_EQUALS(mesh_reader2.GetNumNodes(), 6u);
        TS_ASSERT_EQUALS(mesh_reader2.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(mesh_reader2.GetNumElementAttributes(), 1u);
    }
};

#endif // TESTPOTTSMESHWRITER_HPP_
