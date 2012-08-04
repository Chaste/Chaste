/*

Copyright (c) 2005-2012, University of Oxford.
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
#include "SemMesh.hpp"
#include "SemMeshGenerator.hpp"
#include "PottsElement.hpp"
#include "SemMeshReader.hpp"
#include "SemMeshWriter.hpp"
#include "VtkMeshWriter.hpp"
#include "FileComparison.hpp"

class TestSemMeshWriter : public CxxTest::TestSuite
{
public:

    void TestSemMeshWriter2d() throw (Exception)
    {
        // Create 2D mesh with 2 square elements
        SemMeshGenerator<2> generator(5, 5);
        SemMesh<2>* p_mesh = generator.GetMesh();

        // Create a Sem mesh writer
        SemMeshWriter<2> sem_mesh_writer("TestSemMeshWriter2d", "sem_mesh_2d");

        // Write and check it's correct
        sem_mesh_writer.WriteFilesUsingMesh(*p_mesh);

        OutputFileHandler handler("TestSemMeshWriter2d", false);
        std::string results_file1 = handler.GetOutputDirectoryFullPath() + "sem_mesh_2d.node";
        std::string results_file2 = handler.GetOutputDirectoryFullPath() + "sem_mesh_2d.cell";

        {
            FileComparison comparer(results_file1,
                    "cell_based/test/data/TestSemMeshWriter/sem_mesh_2d.node");
            TS_ASSERT(comparer.CompareFiles());
        }
        {
            FileComparison comparer(results_file2,
                    "cell_based/test/data/TestSemMeshWriter/sem_mesh_2d.cell");
            TS_ASSERT(comparer.CompareFiles());
        }
    }

    void TestSemMeshWriter3d() throw(Exception)
    {
        // Create 3D mesh with 2 square elements
        SemMeshGenerator<3> generator(5, 5, 5, 5, 5, 5);
        SemMesh<3>* p_mesh = generator.GetMesh();

        // Create a Sem mesh writer
        SemMeshWriter<3> sem_mesh_writer("TestSemMeshWriter3d", "sem_mesh_3d");

        // Write and check it's correct
        sem_mesh_writer.WriteFilesUsingMesh(*p_mesh);

        OutputFileHandler handler("TestSemMeshWriter3d", false);
        std::string results_file1 = handler.GetOutputDirectoryFullPath() + "sem_mesh_3d.node";
        std::string results_file2 = handler.GetOutputDirectoryFullPath() + "sem_mesh_3d.cell";

        {
            FileComparison comparer(results_file1,
                    "cell_based/test/data/TestSemMeshWriter/sem_mesh_3d.node");
            TS_ASSERT(comparer.CompareFiles());
        }
        {
            FileComparison comparer(results_file2,
                    "cell_based/test/data/TestSemMeshWriter/sem_mesh_3d.cell");
            TS_ASSERT(comparer.CompareFiles());
        }
    }

    void TestReadingAndWritingElementAttributes() throw(Exception)
    {
        // Read in a mesh with element attributes
        SemMeshReader<2> mesh_reader("cell_based/test/data/TestSemMeshReader2d/sem_mesh_2d");
        TS_ASSERT_EQUALS(mesh_reader.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(mesh_reader.GetNumElementAttributes(), 1u);

        // Construct the mesh
        SemMesh<2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        TS_ASSERT_EQUALS(mesh.GetElement(0)->GetUnsignedAttribute(), 97u);
        TS_ASSERT_EQUALS(mesh.GetElement(1)->GetUnsignedAttribute(), 98u);

        // Write the mesh to file
        SemMeshWriter<2> mesh_writer("TestReadingAndWritingElementAttributes", "sem_mesh_with_element_attributes");
        mesh_writer.WriteFilesUsingMesh(mesh);

        // Now read in the mesh that was written
        OutputFileHandler handler("TestReadingAndWritingElementAttributes", false);
        SemMeshReader<2> mesh_reader2(handler.GetOutputDirectoryFullPath() + "sem_mesh_with_element_attributes");
        TS_ASSERT_EQUALS(mesh_reader2.GetNumElements(), 2u);
        TS_ASSERT_EQUALS(mesh_reader2.GetNumElementAttributes(), 1u);

        // Construct the mesh again
        SemMesh<2> mesh2;
        mesh2.ConstructFromMeshReader(mesh_reader);
        TS_ASSERT_EQUALS(mesh2.GetElement(0)->GetUnsignedAttribute(), 97u);
        TS_ASSERT_EQUALS(mesh2.GetElement(1)->GetUnsignedAttribute(), 98u);
    }
};

#endif // TESTPOTTSMESHWRITER_HPP_
