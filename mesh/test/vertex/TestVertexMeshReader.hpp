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
#ifndef TESTVERTEXMESHREADER2D_HPP_
#define TESTVERTEXMESHREADER2D_HPP_

#include <cxxtest/TestSuite.h>

#include "VertexMeshReader.hpp"
#include "PetscSetupAndFinalize.hpp"

/*
 * This typedef is just because we can't have lines such as
 * TS_ASSERT_THROWS_NOTHING(p_mesh_reader=new VertexMeshReader<2,2>(name));
 * because the macro thinks the comma separates two arguments
 */
typedef VertexMeshReader<2,2> READER_2D;

class TestVertexMeshReader : public CxxTest::TestSuite
{
public:

    /**
     * Check that input files are opened correctly.
     */
    void TestFilesOpen()
    {
        VertexMeshReader<2,2> mesh_reader("mesh/test/data/TestVertexMeshWriter/vertex_mesh_2d");
    }

    /**
     * Check that the nodes are read correctly. Checks that the output vector
     * for a given input file is the correct length and that if the input file
     * is corrupted (missing nodes) then an exception is thrown.
     */
    void TestNodesDataRead()
    {
        VertexMeshReader<2,2> mesh_reader("mesh/test/data/TestVertexMeshWriter/vertex_mesh_2d");

        TS_ASSERT_EQUALS(mesh_reader.GetNumNodes(), 7u);

        VertexMeshReader<2,2> mesh_reader2("mesh/test/data/baddata/vertex_mesh_bad_nodes");

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
    void TestElementsDataRead()
    {
        VertexMeshReader<2,2> mesh_reader("mesh/test/data/TestVertexMeshWriter/vertex_mesh_2d");

        TS_ASSERT_EQUALS(mesh_reader.GetNumElements(), 2u);

        // Read element 0 from file
        ElementData data1 = mesh_reader.GetNextElementData();

        TS_ASSERT_EQUALS(data1.NodeIndices.size(), 5u);
        TS_ASSERT_EQUALS(data1.NodeIndices[0], 0u);
        TS_ASSERT_EQUALS(data1.NodeIndices[1], 1u);
        TS_ASSERT_EQUALS(data1.NodeIndices[2], 2u);
        TS_ASSERT_EQUALS(data1.NodeIndices[3], 3u);
        TS_ASSERT_EQUALS(data1.NodeIndices[4], 4u);

        // Read element 1 from file
        ElementData data2 = mesh_reader.GetNextElementData();

        TS_ASSERT_EQUALS(data2.NodeIndices.size(), 3u);
        TS_ASSERT_EQUALS(data2.NodeIndices[0], 2u);
        TS_ASSERT_EQUALS(data2.NodeIndices[1], 5u);
        TS_ASSERT_EQUALS(data2.NodeIndices[2], 6u);

        TS_ASSERT_EQUALS(mesh_reader.GetNumElementAttributes(), 1u);

        mesh_reader.Reset();
        for (unsigned i=1; i<mesh_reader.GetNumElements(); i++)
        {
            ElementData data = mesh_reader.GetNextElementData();
            TS_ASSERT_EQUALS(data.AttributeValue, 0u);
        }

        VertexMeshReader<2,2> mesh_reader2("mesh/test/data/baddata/vertex_mesh_bad_elements");

        // Reads element 0 from file
        TS_ASSERT_THROWS_NOTHING(mesh_reader2.GetNextElementData());

        // Reads element 2 from file when expecting number 1
        TS_ASSERT_THROWS_THIS(mesh_reader2.GetNextElementData(), "Data for element 1 missing");
    }

    /**
     * Checks that nodes in the input data file are numbered sequentially.
     * (In the input file nodes must appear in increasing order since the node
     * number is only stored as the index of the vector in which the coordinates
     * are stored.)
     */
    void TestPermutedNodesFail()
    {
        VertexMeshReader<2,2> mesh_reader("mesh/test/data/baddata/vertex_mesh_permuted_nodes");
        TS_ASSERT_THROWS_THIS(for(unsigned i=0;i<mesh_reader.GetNumNodes();i++){mesh_reader.GetNextNode();}, "Data for node 3 missing")
    }

    /**
     * Check that GetNextNode() returns the coordinates of the correct node.
     * Compares the coordinates of the first two nodes with their known
     * values, checks that no errors are thrown for the remaining nodes and
     * that an error is thrown if we try to call the function too many times.
     */
    void TestGetNextNode()
    {
        VertexMeshReader<2,2> mesh_reader("mesh/test/data/TestVertexMeshWriter/vertex_mesh_2d");

        std::vector<double> first_node;
        first_node = mesh_reader.GetNextNode();

        TS_ASSERT_DELTA(first_node[0], 0.0, 1e-6);
        TS_ASSERT_DELTA(first_node[1], 0.0, 1e-6)

        std::vector<double> next_node;
        next_node = mesh_reader.GetNextNode();

        TS_ASSERT_DELTA(next_node[0], 1.0, 1e-6);
        TS_ASSERT_DELTA(next_node[1], 0.0, 1e-6);

        for (unsigned i=0; i<5; i++)
        {
            TS_ASSERT_THROWS_NOTHING(next_node = mesh_reader.GetNextNode());
        }

        TS_ASSERT_THROWS_THIS(next_node = mesh_reader.GetNextNode(),
                "Cannot get the next line from node or element file due to incomplete data");
    }

    /**
     * Check that GetNextElementData() works. Checks that no errors are thrown for
     * all of the elements and that an error is thrown if we try to call the
     * function too many times.
     */
    void TestGetNextElementData()
    {
        VertexMeshReader<2,2> mesh_reader("mesh/test/data/TestVertexMeshWriter/vertex_mesh_2d");

        std::vector<unsigned> next_element;
        for (unsigned i=0; i<mesh_reader.GetNumElements(); i++)
        {
            TS_ASSERT_THROWS_NOTHING(next_element = mesh_reader.GetNextElementData().NodeIndices);
        }

        TS_ASSERT_THROWS_THIS(next_element = mesh_reader.GetNextElementData().NodeIndices,
                "Cannot get the next line from node or element file due to incomplete data");
    }

    /**
     * Check that GetNextElementDataWithFaces() works. Checks that no errors are thrown for
     * all of the elements and that an error is thrown if we try to call the
     * function too many times.
     */
    void TestGetNextElementDataWithFaces()
    {
        // First test the case where there aren't actually any faces
        VertexMeshReader<3,3> mesh_reader1("mesh/test/data/TestVertexMeshWriter/vertex_mesh_3d");

        // This mesh should consist of a single cubic element with eight nodes and no faces
        TS_ASSERT_EQUALS(mesh_reader1.GetNumElements(), 1u);

        VertexElementData element_0_data = mesh_reader1.GetNextElementDataWithFaces();

        // Test node indices
        std::vector<unsigned> node_indices = element_0_data.NodeIndices;
        TS_ASSERT_EQUALS(node_indices.size(), 8u);
        for (unsigned i=0; i<8; i++)
        {
            TS_ASSERT_EQUALS(node_indices[i], i);
        }

        // Test there aren't any faces
        std::vector<ElementData> faces = element_0_data.Faces;
        TS_ASSERT_EQUALS(faces.size(), 0u);

        // Test an exception is thrown if we try to access the next element
        TS_ASSERT_THROWS_THIS(node_indices = mesh_reader1.GetNextElementDataWithFaces().NodeIndices,
                "Cannot get the next line from node or element file due to incomplete data");

        // Now test the case where there are faces
        VertexMeshReader<3,3> mesh_reader2("mesh/test/data/TestVertexMeshWriter/vertex_mesh_3d_with_faces");

        // This mesh should consist of a single tetrahedral element with four nodes and four faces
        TS_ASSERT_EQUALS(mesh_reader2.GetNumElements(), 1u);

        element_0_data = mesh_reader2.GetNextElementDataWithFaces();

        // Test there are four nodes owned by this element
        node_indices = element_0_data.NodeIndices;
        TS_ASSERT_EQUALS(node_indices.size(), 4u);

        // Test the node indices are correct (use a set comparison in case of funny business relating to tetgen versions)
        std::set<unsigned> node_indices_expected;
        std::set<unsigned> node_indices_returned;
        for (unsigned i=0; i<4; i++)
        {
            node_indices_expected.insert(i);
            node_indices_returned.insert(node_indices[i]);
        }
        TS_ASSERT_EQUALS(node_indices_expected, node_indices_returned);

        // Test there are four faces owned by this element
        faces = element_0_data.Faces;
        TS_ASSERT_EQUALS(faces.size(), 4u);

        // Test the first face has the correct index and owns the correct nodes
        ElementData face_0 = faces[0];
        TS_ASSERT_EQUALS(face_0.NodeIndices.size(), 3u);
        TS_ASSERT_EQUALS(face_0.NodeIndices[0], 3u);
        TS_ASSERT_EQUALS(face_0.NodeIndices[1], 0u);
        TS_ASSERT_EQUALS(face_0.NodeIndices[2], 2u);

        // Test an exception is thrown if we try to access the next element
        TS_ASSERT_THROWS_THIS(node_indices = mesh_reader2.GetNextElementDataWithFaces().NodeIndices,
                "Cannot get the next line from node or element file due to incomplete data");

        VertexMeshReader<3,3> mesh_reader3("mesh/test/data/baddata/vertex_mesh_3d_with_faces");

        // Reads element 1 from file when expecting number 0
        TS_ASSERT_THROWS_THIS(mesh_reader3.GetNextElementDataWithFaces(), "Data for element 0 missing");
    }

    void TestReadingWithNoElementAttributes()
    {
        VertexMeshReader<2,2> mesh_reader("mesh/test/data/TestVertexMeshReader2d/vertex_mesh_with_no_element_attributes");

        TS_ASSERT_EQUALS(mesh_reader.GetNumElements(), 2u);

        TS_ASSERT_EQUALS(mesh_reader.GetNumElementAttributes(), 0u);
        ElementData next_element_info = mesh_reader.GetNextElementData();
        TS_ASSERT_EQUALS(next_element_info.AttributeValue, 0u);
    }

    void TestReadingElementAttributes()
    {
        VertexMeshReader<2,2> mesh_reader("mesh/test/data/TestVertexMeshReader2d/vertex_mesh_with_element_attributes");

        TS_ASSERT_EQUALS(mesh_reader.GetNumElements(), 2u);

        TS_ASSERT_EQUALS(mesh_reader.GetNumElementAttributes(), 1u);

        ElementData next_element_info = mesh_reader.GetNextElementData();
        std::vector<unsigned> nodes = next_element_info.NodeIndices;
        TS_ASSERT_EQUALS(nodes.size(), 5u);
        TS_ASSERT_EQUALS(next_element_info.AttributeValue, 97u);

        next_element_info = mesh_reader.GetNextElementData();
        nodes = next_element_info.NodeIndices;
        TS_ASSERT_EQUALS(nodes.size(), 3u);
        TS_ASSERT_EQUALS(next_element_info.AttributeValue, 152u);

        /*
         * Coverage
         *
         * \todo The methods GetNextFaceData() and GetNumFaces() are not
         * fully implemented for VertexMeshReader, but must be overridden
         * as they are pure virtual in the base class. When they are
         * implemented, these lines need to be replaced by proper tests.
         *
         * See also #1001.
         */
        ElementData face_data = mesh_reader.GetNextFaceData();
        TS_ASSERT_EQUALS(face_data.NodeIndices.empty(), true);
        TS_ASSERT_EQUALS(face_data.AttributeValue, 0u);

        TS_ASSERT_EQUALS(mesh_reader.GetNumFaces(), 0u);
        TS_ASSERT_EQUALS(mesh_reader.GetNumEdges(), 0u);
    }

    void TestOtherExceptions()
    {
        TS_ASSERT_THROWS_THIS(READER_2D mesh_reader("mesh/test/data/nonexistent_file"),
                "Could not open data file: mesh/test/data/nonexistent_file.node");
        TS_ASSERT_THROWS_THIS(READER_2D mesh_reader("mesh/test/data/baddata/vertex_mesh_without_element_file"),
                "Could not open data file: mesh/test/data/baddata/vertex_mesh_without_element_file.cell");
    }
};

#endif /*TESTVERTEXMESHREADER2D_HPP_*/
