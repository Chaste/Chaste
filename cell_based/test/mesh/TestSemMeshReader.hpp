/*

Copyright (c) 2005-2013, University of Oxford.
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

#ifndef TESTPOTTSMESHREADER_HPP_
#define TESTPOTTSMESHREADER_HPP_

#include <cxxtest/TestSuite.h>

#include <string>
#include <fstream>

#include "SemMeshReader.hpp"

class TestSemMeshReader : public CxxTest::TestSuite
{
public:

    /**
     * Check that input files are opened correctly.
     */
    void TestFilesOpen() throw(Exception)
    {
        SemMeshReader<2> mesh_reader("cell_based/test/data/TestSemMeshWriter/sem_mesh_2d");
    }

   /**
    * Check that the nodes are read correctly. Checks that the output vector
    * for a given input file is the correct length and that if the input file
    * is corrupted (missing nodes) then an exception is thrown.
    */
   void TestNodesDataRead() throw(Exception)
   {
       SemMeshReader<2> mesh_reader("cell_based/test/data/TestSemMeshWriter/sem_mesh_2d");

       TS_ASSERT_EQUALS(mesh_reader.GetNumNodes(), 2500u);

       SemMeshReader<2> mesh_reader2("cell_based/test/data/baddata/potts_mesh_bad_nodes");

       // Reads node 0 from file
       TS_ASSERT_THROWS_NOTHING(mesh_reader2.GetNextNode());

       // Reads node 3 from file when expecting number 1
       TS_ASSERT_THROWS_THIS(mesh_reader2.GetNextNode(), "Data for node 1 missing");
   }

   /**
    * Check that the elements are read correctly. Checks that the output vector
    * for a given input file is the correct length and that if the input file
    * is corrupted (missing elements) then an exception is thrown.
    */
   void TestElementsDataRead() throw(Exception)
   {
       SemMeshReader<2> mesh_reader("cell_based/test/data/TestSemMeshWriter/sem_mesh_2d");

       TS_ASSERT_EQUALS(mesh_reader.GetNumElements(), 25u);

       // Read element 0 from file
       ElementData data = mesh_reader.GetNextElementData();

       TS_ASSERT_EQUALS(data.NodeIndices.size(), 100u);
       TS_ASSERT_EQUALS(data.NodeIndices[0], 0u);
       TS_ASSERT_EQUALS(data.NodeIndices[1], 1u);
       TS_ASSERT_EQUALS(data.NodeIndices[2], 2u);
       TS_ASSERT_EQUALS(data.NodeIndices[3], 3u);

       // Read element 1 from file
       ElementData data2 = mesh_reader.GetNextElementData();

       TS_ASSERT_EQUALS(data2.NodeIndices.size(), 100u);
       TS_ASSERT_EQUALS(data2.NodeIndices[0], 100u);
       TS_ASSERT_EQUALS(data2.NodeIndices[1], 101u);
       TS_ASSERT_EQUALS(data2.NodeIndices[2], 102u);
       TS_ASSERT_EQUALS(data2.NodeIndices[3], 103u);

       TS_ASSERT_EQUALS(mesh_reader.GetNumElementAttributes(), 0u);

       mesh_reader.Reset();
       for (unsigned i=1; i<mesh_reader.GetNumElements(); i++)
       {
           ElementData data = mesh_reader.GetNextElementData();
           TS_ASSERT_EQUALS(data.AttributeValue, 0u);
       }

       SemMeshReader<2> mesh_reader2("cell_based/test/data/baddata/potts_mesh_bad_elements");

       // Reads element 0 from file
       TS_ASSERT_THROWS_NOTHING(mesh_reader2.GetNextElementData());

       // Reads element 2 from file when expecting number 1
       TS_ASSERT_THROWS_THIS(mesh_reader2.GetNextElementData(), "Data for element 1 missing");
   }

    /**
     * Check that GetNextNode() returns the coordinates of the correct node.
     * Compares the coordinates of the first two nodes with their known
     * values, checks that no errors are thrown for the remaining nodes and
     * that an error is thrown if we try to call the function too many times.
     */
    void TestGetNextNode() throw(Exception)
    {
        SemMeshReader<2> mesh_reader("cell_based/test/data/TestSemMeshWriter/sem_mesh_2d");

        std::vector<double> first_node;
        first_node = mesh_reader.GetNextNode();

        TS_ASSERT_DELTA(first_node[0], 0.0, 1e-6);
        TS_ASSERT_DELTA(first_node[1], 0.0, 1e-6)

        std::vector<double> next_node;
        next_node = mesh_reader.GetNextNode();

        TS_ASSERT_DELTA(next_node[0], 0.1904, 1e-4);
        TS_ASSERT_DELTA(next_node[1], 0.0, 1e-6);

        for (unsigned i=2; i<2500; i++)
        {
            TS_ASSERT_THROWS_NOTHING(next_node = mesh_reader.GetNextNode());
        }

        TS_ASSERT_THROWS_THIS(next_node = mesh_reader.GetNextNode(),
                "Cannot get the next line from node or element file due to incomplete data");
    }

    void TestReadingElementAttributes() throw(Exception)
    {
        SemMeshReader<2> mesh_reader("cell_based/test/data/TestSemMeshReader2d/sem_mesh_2d");

        TS_ASSERT_EQUALS(mesh_reader.GetNumElements(), 2u);

        TS_ASSERT_EQUALS(mesh_reader.GetNumElementAttributes(), 1u);

        ElementData next_element_info = mesh_reader.GetNextElementData();
        std::vector<unsigned> nodes = next_element_info.NodeIndices;
        TS_ASSERT_EQUALS(nodes.size(), 100u);
        TS_ASSERT_EQUALS(next_element_info.AttributeValue, 97u);

        next_element_info = mesh_reader.GetNextElementData();
        nodes = next_element_info.NodeIndices;
        TS_ASSERT_EQUALS(nodes.size(), 100u);
        TS_ASSERT_EQUALS(next_element_info.AttributeValue, 98u)

        /*
         * Coverage
         *
         * \todo The methods GetNextFaceData() and GetNumFaces() are not
         * fully implemented for SemMeshReader, but must be overridden
         * as they are pure virtual in the base class. When they are
         * implemented, these lines need to be replaced by proper tests.
         *
         * See also #1663.
         */
        ElementData face_data = mesh_reader.GetNextFaceData();
        TS_ASSERT_EQUALS(face_data.NodeIndices.empty(), true);
        TS_ASSERT_EQUALS(face_data.AttributeValue, 0u);

        TS_ASSERT_EQUALS(mesh_reader.GetNumFaces(), 0u);
        TS_ASSERT_EQUALS(mesh_reader.GetNumEdges(), 0u);
    }

    void TestSemMeshReader3d() throw(Exception)
    {
        SemMeshReader<3> mesh_reader("cell_based/test/data/TestSemMeshWriter/sem_mesh_3d");

        TS_ASSERT_EQUALS(mesh_reader.GetNumNodes(), 15625u);
        TS_ASSERT_EQUALS(mesh_reader.GetNumElements(), 125u);
        TS_ASSERT_EQUALS(mesh_reader.GetNumElementAttributes(), 0u);

        std::vector<double> next_node;
        next_node = mesh_reader.GetNextNode();

        TS_ASSERT_DELTA(next_node[0], 0.0, 1e-6);
        TS_ASSERT_DELTA(next_node[1], 0.0, 1e-6)
        TS_ASSERT_DELTA(next_node[2], 0.0, 1e-6)

        next_node = mesh_reader.GetNextNode();

        TS_ASSERT_DELTA(next_node[0], 0.3618, 1e-4);
        TS_ASSERT_DELTA(next_node[1], 0.0, 1e-6)
        TS_ASSERT_DELTA(next_node[2], 0.0, 1e-6)

        next_node = mesh_reader.GetNextNode();

        TS_ASSERT_DELTA(next_node[0], 0.7237, 1e-4);
        TS_ASSERT_DELTA(next_node[1], 0.0, 1e-6)
        TS_ASSERT_DELTA(next_node[2], 0.0, 1e-6)
    }

    void TestOtherExceptions() throw(Exception)
    {
        TS_ASSERT_THROWS_THIS(SemMeshReader<2> mesh_reader("cell_based/test/data/nonexistent_file"),
                "Could not open data file: cell_based/test/data/nonexistent_file.node");
        TS_ASSERT_THROWS_THIS(SemMeshReader<2> mesh_reader("cell_based/test/data/baddata/potts_mesh_without_element_file"),
                "Could not open data file: cell_based/test/data/baddata/potts_mesh_without_element_file.cell");
    }
};

#endif // TESTPOTTSMESHREADER_HPP_
