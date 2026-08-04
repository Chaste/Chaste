/*

Copyright (c) 2005-2026, University of Oxford.
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
#ifndef TESTSEMMESHREADER_HPP_
#define TESTSEMMESHREADER_HPP_

#include <cxxtest/TestSuite.h>

#include <climits>
#include <fstream>
#include <string>

#include "SemMeshReader.hpp"
#include "SemMeshWriter.hpp"
#include "SemMesh.hpp"
#include "SemSingleElementMeshGenerator.hpp"
#include "SemMultiElementMeshGenerator.hpp"
#include "OutputFileHandler.hpp"
#include "PetscSetupAndFinalize.hpp"


class TestSemMeshReader : public CxxTest::TestSuite
{
private:

    /**
     * Write a valid 2D SEM mesh and return its base name (no extension). The reader's fixtures are
     * generated rather than checked in, so they cannot drift out of step with SemMeshWriter's
     * format -- which is still settling.
     *
     * @param rDirectory the output directory to write into
     *
     * @return the mesh base name, suitable for passing to a SemMeshReader
     */
    std::string WriteValidMesh(const std::string& rDirectory)
    {
        SemSingleElementMeshGenerator<2> generator({3, 3}, 0.5);
        auto p_mesh = generator.GetMesh();
        SemMeshWriter<2> writer(rDirectory, "sem_mesh", true);
        writer.WriteFilesUsingMesh(*p_mesh);

        return OutputFileHandler::GetChasteTestOutputDirectory() + rDirectory + "/sem_mesh";
    }

    /**
     * Copy a text file, optionally dropping one line and/or stopping early. Used to turn a valid
     * mesh file into a corrupt one.
     *
     * @param rSource path of the file to read
     * @param rDestination path of the file to write
     * @param lineToSkip zero-based index of a line to omit, or UINT_MAX to keep every line
     * @param maxLines stop after writing this many lines, or UINT_MAX for no limit
     */
    void CopyFileWithCorruption(const std::string& rSource,
                                const std::string& rDestination,
                                unsigned lineToSkip,
                                unsigned maxLines)
    {
        std::ifstream in_file(rSource.c_str());
        TS_ASSERT(in_file.is_open());
        std::ofstream out_file(rDestination.c_str());
        TS_ASSERT(out_file.is_open());

        std::string line;
        unsigned line_index = 0u;
        unsigned lines_written = 0u;
        while (std::getline(in_file, line) && lines_written < maxLines)
        {
            if (line_index != lineToSkip)
            {
                out_file << line << "\n";
                lines_written++;
            }
            line_index++;
        }
    }

public:

    void TestDimensionMismatchThrows()
    {
        EXIT_IF_PARALLEL;  // Writes and reads mesh files

        // Write a 2D SEM mesh to file
        SemSingleElementMeshGenerator<2> generator({3, 3}, 0.5);
        auto p_mesh = generator.GetMesh();
        SemMeshWriter<2> writer("TestSemMeshReaderDimension", "sem_mesh", true);
        writer.WriteFilesUsingMesh(*p_mesh);

        std::string mesh_base = OutputFileHandler::GetChasteTestOutputDirectory() + "TestSemMeshReaderDimension/sem_mesh";

        // Reading a 2D mesh file with a reader of a different dimension must throw
        TS_ASSERT_THROWS_CONTAINS((SemMeshReader<3>(mesh_base)), "does not match the reader's DIM");
        TS_ASSERT_THROWS_CONTAINS((SemMeshReader<1>(mesh_base)), "does not match the reader's DIM");

        // Reading with the matching dimension must succeed
        TS_ASSERT_THROWS_NOTHING((SemMeshReader<2>(mesh_base)));
    }

    void TestMissingFilesThrow()
    {
        EXIT_IF_PARALLEL;

        std::string missing_base = OutputFileHandler::GetChasteTestOutputDirectory()
                                   + "TestSemMeshReaderMissing/no_such_mesh";
        TS_ASSERT_THROWS_CONTAINS((SemMeshReader<2>(missing_base)), "Could not open data file");

        // A present .node file with no matching .cell file must still throw
        std::string valid_base = WriteValidMesh("TestSemMeshReaderMissingElements");
        OutputFileHandler handler("TestSemMeshReaderMissingElements", false);
        std::string partial_base = handler.GetOutputDirectoryFullPath() + "partial_mesh";
        CopyFileWithCorruption(valid_base + ".node", partial_base + ".node", UINT_MAX, UINT_MAX);

        TS_ASSERT_THROWS_CONTAINS((SemMeshReader<2>(partial_base)), "Could not open data file");
    }

    /**
     * SemMeshWriter never emits element attributes, so nothing produced by Chaste exercises the
     * reader's attribute column. The format supports one all the same, and a mesh file written by
     * hand or by another tool may carry it, so read one back.
     */
    void TestElementAttributeIsRead()
    {
        EXIT_IF_PARALLEL;

        // Borrow the node file from a valid mesh: 9 nodes in one element
        std::string valid_base = WriteValidMesh("TestSemMeshReaderElementAttribute");
        OutputFileHandler handler("TestSemMeshReaderElementAttribute", false);
        std::string attributed_base = handler.GetOutputDirectoryFullPath() + "attributed_mesh";
        CopyFileWithCorruption(valid_base + ".node", attributed_base + ".node", UINT_MAX, UINT_MAX);

        // Hand-write the element file, declaring one element attribute in the header
        {
            std::ofstream element_file((attributed_base + ".cell").c_str());
            TS_ASSERT(element_file.is_open());
            element_file << "1 1\n";
            element_file << "0 9 0 1 2 3 4 5 6 7 8 3\n";
        }

        SemMeshReader<2> reader(attributed_base);
        TS_ASSERT_EQUALS(reader.GetNumElementAttributes(), 1u);
        TS_ASSERT_EQUALS(reader.GetNumElements(), 1u);

        ElementData element_data = reader.GetNextElementData();
        TS_ASSERT_EQUALS(element_data.NodeIndices.size(), 9u);
        TS_ASSERT_EQUALS(element_data.NodeIndices[0], 0u);
        TS_ASSERT_EQUALS(element_data.NodeIndices[8], 8u);
        TS_ASSERT_DELTA(element_data.AttributeValue, 3.0, 1e-12);
    }

    void TestOutOfSequenceNodeDataThrows()
    {
        EXIT_IF_PARALLEL;

        std::string valid_base = WriteValidMesh("TestSemMeshReaderBadNodes");
        OutputFileHandler handler("TestSemMeshReaderBadNodes", false);
        std::string corrupt_base = handler.GetOutputDirectoryFullPath() + "corrupt_mesh";

        // Line 0 is the header and line 1 is node 0, so dropping line 3 removes node 2
        CopyFileWithCorruption(valid_base + ".node", corrupt_base + ".node", 3u, UINT_MAX);
        CopyFileWithCorruption(valid_base + ".cell", corrupt_base + ".cell", UINT_MAX, UINT_MAX);

        SemMeshReader<2> reader(corrupt_base);

        // The first two nodes are intact; the third is where the sequence breaks
        TS_ASSERT_THROWS_NOTHING(reader.GetNextNode());
        TS_ASSERT_THROWS_NOTHING(reader.GetNextNode());
        TS_ASSERT_THROWS_CONTAINS(reader.GetNextNode(), "Data for node 2 missing");
    }

    void TestOutOfSequenceElementDataThrows()
    {
        EXIT_IF_PARALLEL;

        // Two elements, so that dropping the first element line leaves a second to misread
        SemMultiElementMeshGenerator<2> generator({3, 3}, {2, 1}, 0.5);
        auto p_mesh = generator.GetMesh();
        SemMeshWriter<2> writer("TestSemMeshReaderBadElements", "sem_mesh", true);
        writer.WriteFilesUsingMesh(*p_mesh);

        OutputFileHandler handler("TestSemMeshReaderBadElements", false);
        std::string valid_base = handler.GetOutputDirectoryFullPath() + "sem_mesh";
        std::string corrupt_base = handler.GetOutputDirectoryFullPath() + "corrupt_mesh";

        CopyFileWithCorruption(valid_base + ".node", corrupt_base + ".node", UINT_MAX, UINT_MAX);

        // Line 0 is the header and line 1 is element 0, so dropping line 1 removes element 0
        CopyFileWithCorruption(valid_base + ".cell", corrupt_base + ".cell", 1u, UINT_MAX);

        SemMeshReader<2> reader(corrupt_base);
        TS_ASSERT_THROWS_CONTAINS(reader.GetNextElementData(), "Data for element 0 missing");
    }

    void TestTruncatedFileThrows()
    {
        EXIT_IF_PARALLEL;

        std::string valid_base = WriteValidMesh("TestSemMeshReaderTruncated");
        OutputFileHandler handler("TestSemMeshReaderTruncated", false);
        std::string truncated_base = handler.GetOutputDirectoryFullPath() + "truncated_mesh";

        // Keep the header and the first two nodes only, so the file runs out mid-read. This is a
        // different failure from an out-of-sequence index: the stream hits EOF rather than
        // reaching a line with the wrong index on it.
        CopyFileWithCorruption(valid_base + ".node", truncated_base + ".node", UINT_MAX, 3u);
        CopyFileWithCorruption(valid_base + ".cell", truncated_base + ".cell", UINT_MAX, UINT_MAX);

        SemMeshReader<2> reader(truncated_base);

        TS_ASSERT_THROWS_NOTHING(reader.GetNextNode());
        TS_ASSERT_THROWS_NOTHING(reader.GetNextNode());
        TS_ASSERT_THROWS_CONTAINS(reader.GetNextNode(), "incomplete data");
    }

    void TestResetRereadsFromTheStart()
    {
        EXIT_IF_PARALLEL;

        std::string mesh_base = WriteValidMesh("TestSemMeshReaderReset");
        SemMeshReader<2> reader(mesh_base);

        TS_ASSERT_EQUALS(reader.GetNumNodes(), 9u);
        TS_ASSERT_EQUALS(reader.GetNumElements(), 1u);

        std::vector<double> first_node = reader.GetNextNode();
        std::vector<double> second_node = reader.GetNextNode();
        ElementData first_element = reader.GetNextElementData();

        reader.Reset();

        // Both streams restart, and the header counts survive the reopen
        TS_ASSERT_EQUALS(reader.GetNumNodes(), 9u);
        TS_ASSERT_EQUALS(reader.GetNumElements(), 1u);

        std::vector<double> first_node_again = reader.GetNextNode();
        TS_ASSERT_EQUALS(first_node_again.size(), first_node.size());
        for (unsigned i = 0; i < first_node.size(); ++i)
        {
            TS_ASSERT_DELTA(first_node_again[i], first_node[i], 1e-9);
        }

        std::vector<double> second_node_again = reader.GetNextNode();
        for (unsigned i = 0; i < second_node.size(); ++i)
        {
            TS_ASSERT_DELTA(second_node_again[i], second_node[i], 1e-9);
        }

        ElementData first_element_again = reader.GetNextElementData();
        TS_ASSERT_EQUALS(first_element_again.NodeIndices, first_element.NodeIndices);
    }

    void TestFaceMethodsReportNoFaces()
    {
        EXIT_IF_PARALLEL;

        // SEM meshes have no faces. Both methods exist only to satisfy AbstractMeshReader.
        std::string mesh_base = WriteValidMesh("TestSemMeshReaderFaces");
        SemMeshReader<2> reader(mesh_base);

        TS_ASSERT_EQUALS(reader.GetNumFaces(), 0u);
        TS_ASSERT_EQUALS(reader.GetNumElementAttributes(), 0u);

        ElementData face_data = reader.GetNextFaceData();
        TS_ASSERT(face_data.NodeIndices.empty());
    }
};

#endif /*TESTSEMMESHREADER_HPP_*/
