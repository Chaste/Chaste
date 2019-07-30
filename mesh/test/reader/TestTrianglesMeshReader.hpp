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


/**
 * Test suite for the TrianglesMeshReader class.
 *
 */

#ifndef _TESTTRIANGLESMESHREADER_HPP_
#define _TESTTRIANGLESMESHREADER_HPP_

#include <cxxtest/TestSuite.h>
#include "TrianglesMeshReader.hpp"
#include "GenericMeshReader.hpp"
#include "TrianglesMeshWriter.hpp"

// This is needed on Windows to ensure that meshes are looked up in the source folder.
#include "FakePetscSetup.hpp"

// these typedefs are just because can't have lines such as
//  TS_ASSERT_THROWS_NOTHING(p_mesh_reader=new TrianglesMeshReader<2,2>(name));
// because the macro thinks the comma separates two arguments
typedef TrianglesMeshReader<1,1> READER_1D;
typedef TrianglesMeshReader<2,2> READER_2D;
typedef TrianglesMeshReader<3,3> READER_3D;
typedef TrianglesMeshReader<0,1> READER_0D_IN_1D;


class TestTrianglesMeshReader : public CxxTest::TestSuite
{
public:

    /**
     * Check that input files are opened correctly.
     */
    void TestFilesOpen()
    {

        TS_ASSERT_THROWS_THIS(READER_2D mesh_reader1("no_mesh_files"),"Could not open data file: no_mesh_files.node");
        TS_ASSERT_THROWS_THIS(READER_2D mesh_reader1("mesh/test/data/square_no_ele_file"),
                           "Could not open data file: mesh/test/data/square_no_ele_file.ele");
        TS_ASSERT_THROWS_THIS(READER_2D mesh_reader1("mesh/test/data/square_no_edge_file"),
                           "Could not open data file: mesh/test/data/square_no_edge_file.edge");

        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_522_elements");
        // For coverage purposes, not sure how to test this functionality...
        mesh_reader.SetReadBufferSize(2*1024*1024); //2MB
    }

    /**
     * Check that the nodes are read correctly. Checks that the output vector
     * for a given input file is the correct length and that if the input file
     * is corrupted (missing nodes) then an exception is thrown.
     */
    void TestNodesDataRead()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_984_elements_indexed_from_1");

        TS_ASSERT_EQUALS( mesh_reader.GetNumNodes(), 543u);

        TrianglesMeshReader<2,2> mesh_reader2("mesh/test/data/baddata/bad_nodes_disk_522_elements");

        // Reads node 0 from file
        TS_ASSERT_THROWS_NOTHING(mesh_reader2.GetNextNode());
        // Reads node 3 from file when expecting number 1
        TS_ASSERT_THROWS_THIS(mesh_reader2.GetNextNode(),"Data for item 1 missing");
    }

    /**
     * Check that the elements are read correctly. Checks that the output vector
     * for a given input file is the correct length and that if the input file
     * is corrupted (missing elements) then an exception is thrown.
     */
    void TestElementsDataRead()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_984_elements_indexed_from_1");

        TS_ASSERT_EQUALS( mesh_reader.GetNumElements(), 984u);

        ElementData data1 = mesh_reader.GetNextElementData();
        TS_ASSERT_EQUALS(data1.NodeIndices.size(), 3u);
        TS_ASSERT_EQUALS(data1.NodeIndices[0], 309u);
        TS_ASSERT_EQUALS(data1.NodeIndices[1], 144u);
        TS_ASSERT_EQUALS(data1.NodeIndices[2], 310u);
        TS_ASSERT_EQUALS( mesh_reader.GetNumElementAttributes(), 0u);

        for (unsigned i=1; i<mesh_reader.GetNumElements(); i++)
        {
            ElementData data = mesh_reader.GetNextElementData();
            TS_ASSERT_EQUALS(data.AttributeValue, 0u);
        }

        TrianglesMeshReader<2,2> mesh_reader2("mesh/test/data/baddata/bad_elements_disk_522_elements");

        // Reads element 0 from file
        TS_ASSERT_THROWS_NOTHING(mesh_reader2.GetNextElementData());
        // Reads element 2 from file when expecting number 1
        TS_ASSERT_THROWS_THIS(mesh_reader2.GetNextElementData(),"Data for item 1 missing");
    }

    /**
     * Check that the faces are read correctly. Checks that the output vector
     * for a given input file is the correct length and that if the input file
     * is corrupted (missing faces) then an exception is thrown.
     */
    void TestFacesDataRead()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_984_elements_indexed_from_1");

        // TS_ASSERT_EQUALS( mesh_reader.GetNumFaces(), 1526u); // when all faces were read
        TS_ASSERT_EQUALS( mesh_reader.GetNumFaces(), 100u); // just boundary faces are read
        TS_ASSERT_EQUALS( mesh_reader.GetNumFaceAttributes(), 1u);

        for (unsigned i=1; i<mesh_reader.GetNumFaces(); i++)
        {
            ElementData data = mesh_reader.GetNextFaceData();
            TS_ASSERT_EQUALS(data.AttributeValue, 1u);

        }

        // First boundary face is #20, on its way through the file there's a gap between face 1 and face 10
        TS_ASSERT_THROWS_THIS(READER_2D mesh_reader2("mesh/test/data/baddata/bad_faces_disk_522_elements"),"Data for item 2 missing");
    }

    /**
     * Check that the faces are read correctly. Checks that the output vector
     * for a given input file is the correct length and that if the input file
     * is corrupted (missing faces) then an exception is thrown.
     */
    void TestFacesDataReadWithAttributes()
    {
        TrianglesMeshReader<3,3> mesh_reader("heart/test/data/box_shaped_heart/box_heart_nonnegative_flags");

        TS_ASSERT_EQUALS( mesh_reader.GetNumFaces(), 92u); // just boundary faces are read
        TS_ASSERT_EQUALS( mesh_reader.GetNumFaceAttributes(), 1u);

        bool read_zero_attribute = false;
        for (unsigned i=0; i<mesh_reader.GetNumFaces(); i++)
        {
            ElementData data = mesh_reader.GetNextFaceData();
            // Attributes are 0, 1, 2, or 3.
            TS_ASSERT_LESS_THAN(data.AttributeValue, 4u);
            if (data.AttributeValue == 0u)
            {
                read_zero_attribute = true;
            }
            TS_ASSERT(read_zero_attribute);
        }
    }

    /**
     * Checks that the reader can deal with (3-d) TetGen input files as well
     * as the previously considered (2-d) Triangles files. Checks that the
     * element output vector for a given input file is the correct length.
     */
    void Test3dDataRead()
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/slab_138_elements");
        TS_ASSERT_EQUALS(mesh_reader.GetNodeAttributes().size(), 0u);//no nodal attributes in this mesh
        TS_ASSERT_EQUALS(mesh_reader.GetNumElements(), 138u);
    }

    /**
     * Checks that nodes in the input data file are numbered sequentially.
     * (In the input file nodes must appear in increasing order since the node
     * number is only stored as the index of the vector in which the coordinates
     * are stored.)
     */
    void TestPermutedNodesFail()
    {
        TS_ASSERT_THROWS_THIS(READER_2D reader("mesh/test/data/baddata/permuted_nodes_disk_522_elements"),"Data for item 0 missing")
    }

    /**
     * Checks that elements have the correct number of nodes (i.e. one more
     * node than the dimension of the mesh). If quadratic basis functions are
     * required this should be dealt with elsewhere.
     */
    void TestOrder2ElementsFail()
    {
        TS_ASSERT_THROWS_THIS(new READER_2D("mesh/test/data/baddata/disk_522_order_2_elements"),
                "Number of nodes per elem, 6, does not match expected number, 3 "
                "(which is calculated given the order of elements chosen, 1 (1=linear, 2=quadratics)");
    }

    /**
     * Check that GetNextNode() returns the coordinates of the correct node and the correct node attributes.
     * Compares the coordinates of the first two nodes with their known
     * values, checks that no errors are thrown for the remaining nodes and
     * that an error is thrown if we try to call the function too many times.
     */
    void TestGetNextNode()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_984_elements");

        std::vector<double> first_node;
        first_node = mesh_reader.GetNextNode();

        TS_ASSERT_DELTA(first_node[0],  0.9980267283, 1e-6);
        TS_ASSERT_DELTA(first_node[1], -0.0627905195, 1e-6);

        // This mesh has zero attributes and one marker in the node file (as the header specifies).
        // we have to ensure that in such situation the last number in each node line is not mistakenly
        // read and interpreted as a node attribute (it is a node marker instead).
        TS_ASSERT_EQUALS(mesh_reader.GetNodeAttributes().size(), 0u);

        std::vector<double> next_node;
        next_node = mesh_reader.GetNextNode();

        TS_ASSERT_DELTA(next_node[0], 1.0, 1e-6);
        TS_ASSERT_DELTA(next_node[1], 0.0, 1e-6);

        for (int i=0; i<541; i++)
        {
            TS_ASSERT_THROWS_NOTHING(next_node = mesh_reader.GetNextNode());
        }

        TS_ASSERT_EQUALS(mesh_reader.GetNodeAttributes().size(), 0u);

        TS_ASSERT_THROWS_THIS(next_node = mesh_reader.GetNextNode(),
                              "File contains incomplete data: unexpected end of file.");
    }

    void TestGetNextNodeWithNodeAttributes()
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_2mm_12_elements_with_node_attributes");

        TS_ASSERT_EQUALS(mesh_reader.GetNodeAttributes().size(), 0u); // Check vector is empty at the beginning

        std::vector<double> first_node;
        first_node = mesh_reader.GetNextNode();//read a node

        // Check coordinates
        TS_ASSERT_DELTA(first_node[0], 0.0, 1e-6);
        TS_ASSERT_DELTA(first_node[1], 0.0, 1e-6);
        TS_ASSERT_DELTA(first_node[2], 0.0, 1e-6);

        // Check attributes
        TS_ASSERT_EQUALS(mesh_reader.GetNodeAttributes().size(), 2u);
        TS_ASSERT_DELTA(mesh_reader.GetNodeAttributes()[0], 25.2, 1e-6);
        TS_ASSERT_DELTA(mesh_reader.GetNodeAttributes()[1],   16.3, 1e-6);

        std::vector<double> second_node;
        second_node = mesh_reader.GetNextNode();//read another node

        // Check coordinates
        TS_ASSERT_DELTA(second_node[0], 0.2, 1e-6);
        TS_ASSERT_DELTA(second_node[1], 0.0, 1e-6);
        TS_ASSERT_DELTA(second_node[2], 0.0, 1e-6);

        // Check attributes
        TS_ASSERT_EQUALS(mesh_reader.GetNodeAttributes().size(), 2u); // Check that we remembered to clear it ( otherwise size would be 4 by now)
        TS_ASSERT_DELTA(mesh_reader.GetNodeAttributes()[0], 25.6, 1e-6);
        TS_ASSERT_DELTA(mesh_reader.GetNodeAttributes()[1],   15.0, 1e-6);

        // Read a few other nodes in succession
        std::vector<double> another_node;
        another_node = mesh_reader.GetNextNode();//3rd
        another_node = mesh_reader.GetNextNode();//4th
        another_node = mesh_reader.GetNextNode();//5th

        // Check the fifth node
        TS_ASSERT_DELTA(another_node[0], 0.0, 1e-6);
        TS_ASSERT_DELTA(another_node[1], 0.0, 1e-6);
        TS_ASSERT_DELTA(another_node[2], 0.2, 1e-6);

        TS_ASSERT_EQUALS(mesh_reader.GetNodeAttributes().size(), 2u);
        TS_ASSERT_DELTA(mesh_reader.GetNodeAttributes()[0], 2.3, 1e-6);
        TS_ASSERT_DELTA(mesh_reader.GetNodeAttributes()[1],   16.0, 1e-6);
    }

    /**
     * Check that GetNextElementData() works. Checks that no errors are thrown for
     * all of the elements and that an error is thrown if we try to call the
     * function too many times.
     */
    void TestGetNextElementData()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_984_elements");

        std::vector<unsigned> next_element;

        for (unsigned i=0; i<mesh_reader.GetNumElements(); i++)
        {
            TS_ASSERT_THROWS_NOTHING(next_element = mesh_reader.GetNextElementData().NodeIndices);
        }

        TS_ASSERT_THROWS_THIS(next_element = mesh_reader.GetNextElementData().NodeIndices,
                              "File contains incomplete data: unexpected end of file.");
    }

    /**
     * Check that GetNextEdgeData() works. Checks that no errors are thrown for
     * all of the edges and that an error is thrown if we try to call the
     * function too many times.
     */
    void TestGetNextEdgeData()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_984_elements");

        std::vector<unsigned> next_edge;

        TS_ASSERT_THROWS_NOTHING(next_edge = mesh_reader.GetNextFaceData().NodeIndices);
        TS_ASSERT_THROWS_NOTHING(next_edge = mesh_reader.GetNextFaceData().NodeIndices);

        for (unsigned i=2; i<mesh_reader.GetNumEdges(); i++)
        {
            TS_ASSERT_THROWS_NOTHING(next_edge = mesh_reader.GetNextEdgeData().NodeIndices);
        }

        TS_ASSERT_THROWS_THIS(next_edge = mesh_reader.GetNextEdgeData().NodeIndices,
                              "File contains incomplete data: unexpected end of file.");
    }

    /**
     * Check that the 1D data are read correctly. Check that the output vector
     * for a given input file is the correct length.
     */
    void Test1DMeshRead()
    {
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/trivial_1d_mesh");

        TS_ASSERT_EQUALS(mesh_reader.GetNumNodes(), 11u);
        TS_ASSERT_EQUALS(mesh_reader.GetNumElements(), 10u);
        TS_ASSERT_EQUALS(mesh_reader.GetNumFaces(), 2u);
    }

    void Test0DMeshIn1DSpaceFails()
    {
        TS_ASSERT_THROWS_THIS(new READER_0D_IN_1D("mesh/test/data/trivial_1d_mesh"),
                "Can\'t have a zero-dimensional mesh in a one-dimensional space");
    }

    void Test1DMeshIn2DSpace()
    {
        TrianglesMeshReader<1,2> mesh_reader("mesh/test/data/circle_outline");
        TS_ASSERT_EQUALS(mesh_reader.GetNumNodes(), 100u);
        TS_ASSERT_EQUALS(mesh_reader.GetNumElements(), 100u);

        // Note: don't test faces (end nodes), since they are culled later
        TrianglesMeshReader<1,2> mesh_reader2("mesh/test/data/semicircle_outline");
        TS_ASSERT_EQUALS(mesh_reader2.GetNumNodes(), 51u);
        TS_ASSERT_EQUALS(mesh_reader2.GetNumElements(), 50u);
    }

    void Test1DMeshIn3DSpace()
    {
        TrianglesMeshReader<1,3> mesh_reader("mesh/test/data/trivial_1d_in_3d_mesh");
        TS_ASSERT_EQUALS(mesh_reader.GetNumNodes(), 11u);
        TS_ASSERT_EQUALS(mesh_reader.GetNumElements(), 10u);
    }

    void Test2DMeshIn3DSpace()
    {
        TrianglesMeshReader<2,3> mesh_reader("mesh/test/data/slab_395_elements");
        TS_ASSERT_EQUALS(mesh_reader.GetNumNodes(), 132u);
        TS_ASSERT_EQUALS(mesh_reader.GetNumElements(), 224u);
        TS_ASSERT_EQUALS(mesh_reader.GetNumFaces(), 0u);

        TrianglesMeshReader<2,3> mesh_reader2("mesh/test/data/disk_in_3d");
        TS_ASSERT_EQUALS(mesh_reader2.GetNumNodes(), 312u);
        TS_ASSERT_EQUALS(mesh_reader2.GetNumElements(), 522u);
        // TS_ASSERT_EQUALS(mesh_reader2.GetNumFaces(), 833u); // when all faces were read
        TS_ASSERT_EQUALS(mesh_reader2.GetNumFaces(), 100u); // just boundary faces are read
    }

    void TestOtherExceptions()
    {
        // This should fail because SPACE_DIM doesn't match the dimension in the file
        TS_ASSERT_THROWS_THIS( READER_1D mesh_reader("mesh/test/data/disk_984_elements"),
                "SPACE_DIM  != dimension read from file ");
        //Indexed quadratic faces
        TS_ASSERT_THROWS_THIS( READER_1D mesh_reader2("mesh/test/data/baddata/bad_1D_0_to_1_10_elements_quadratic", 2, 2, true),
                "Boundary element file should not have containing element info if it is quadratic");

        //Exceptions due to unimplemented code
        TrianglesMeshReader<3,3> mesh_reader3("mesh/test/data/simple_cube");
        TS_ASSERT_THROWS_THIS(mesh_reader3.GetNode(0), "Random access is only implemented in mesh readers for binary mesh files.");
        TS_ASSERT_THROWS_THIS(mesh_reader3.GetElementData(0), "Random access is only implemented in mesh readers for binary mesh files.");
        TS_ASSERT_THROWS_THIS(mesh_reader3.GetFaceData(0), "Random access is only implemented in mesh readers for binary mesh files.");
        TS_ASSERT_THROWS_THIS(mesh_reader3.GetEdgeData(0), "Random access is only implemented in mesh readers for binary mesh files.");

        TS_ASSERT_THROWS_THIS(mesh_reader3.GetContainingElementIndices(0), "NCL file functionality is only implemented in mesh readers for binary mesh files.");
    }

    ////////////////////////////////////////////////////////
    // Quadratic tests
    ////////////////////////////////////////////////////////
    void TestReadingQuadraticMesh1d()
    {
        TS_ASSERT_THROWS_THIS( READER_1D wrong_reader("mesh/test/data/1D_0_to_1_10_elements_quadratics", 1),
                "Could not open data file: mesh/test/data/1D_0_to_1_10_elements_quadratics.node");

        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements_quadratic", 2);
        TS_ASSERT_EQUALS(mesh_reader.GetNumNodes(), 21u);
        TS_ASSERT_EQUALS(mesh_reader.GetNumElements(), 10u);

        std::vector<unsigned> next_element = mesh_reader.GetNextElementData().NodeIndices;

        TS_ASSERT_EQUALS(next_element.size(), 3u);

        TS_ASSERT_EQUALS(next_element[0], 0u);   // left node
        TS_ASSERT_EQUALS(next_element[1], 1u);   // right node
        TS_ASSERT_EQUALS(next_element[2], 11u);  // middle node

        for (unsigned i=1; i<10; i++)
        {
            next_element = mesh_reader.GetNextElementData().NodeIndices;
            TS_ASSERT_EQUALS(next_element.size(), 3u);
        }
    }

    void TestReadingQuadraticMesh2d()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements_quadratic", 2);
        TS_ASSERT_EQUALS(mesh_reader.GetNumNodes(), 17u*17u);
        TS_ASSERT_EQUALS(mesh_reader.GetNumElements(), 128u);

        std::vector<unsigned> next_element = mesh_reader.GetNextElementData().NodeIndices;

        TS_ASSERT_EQUALS(next_element.size(), 6u);

        TS_ASSERT_EQUALS(next_element[0], 53u);
        TS_ASSERT_EQUALS(next_element[1], 0u);
        TS_ASSERT_EQUALS(next_element[2], 54u);
        TS_ASSERT_EQUALS(next_element[3], 82u); // opposite to 53
        TS_ASSERT_EQUALS(next_element[4], 83u); // opposite to 0
        TS_ASSERT_EQUALS(next_element[5], 81u); // opposite to 54

        for (unsigned i=1; i<128; i++)
        {
            next_element = mesh_reader.GetNextElementData().NodeIndices;
            TS_ASSERT_EQUALS(next_element.size(), 6u);
        }
    }

    void TestReadingQuadraticMesh3d()
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/3D_Single_tetrahedron_element_quadratic", 2);
        TS_ASSERT_EQUALS(mesh_reader.GetNumNodes(), 10u);
        TS_ASSERT_EQUALS(mesh_reader.GetNumElements(), 1u);

        std::vector<unsigned> next_element = mesh_reader.GetNextElementData().NodeIndices;

        TS_ASSERT_EQUALS(next_element.size(), 10u);

        for (unsigned i=0; i<10; i++)
        {
            TS_ASSERT_EQUALS(next_element[i], i);
        }
    }

    void TestReadingElementAttributes()
    {
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements_with_attributes");

        TS_ASSERT_EQUALS(mesh_reader.GetNumElements(), 10u);
        TS_ASSERT_EQUALS(mesh_reader.GetNumElementAttributes(), 1u);

        for (unsigned i=0; i<10; i++)
        {
            ElementData next_element_info = mesh_reader.GetNextElementData();
            std::vector<unsigned> nodes = next_element_info.NodeIndices;
            TS_ASSERT_EQUALS(nodes.size(), 2u);
            if (i<5)
            {
                //unsigned
                TS_ASSERT_EQUALS(next_element_info.AttributeValue, i+1);
            }
            else
            {
                //floating point
                TS_ASSERT_EQUALS(next_element_info.AttributeValue, i%5+1.1);
            }
        }
    }

    void TestReadingContainingElementsOfBoundaryElements()
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_2mm_152_elements_v2", 1, 1, true);

        TS_ASSERT_EQUALS(mesh_reader.GetNumFaces(), 116u);
        TS_ASSERT_EQUALS(mesh_reader.GetNextFaceData().NodeIndices[0], 3u);     //face 0
        TS_ASSERT_EQUALS(mesh_reader.GetNextFaceData().ContainingElement, 36u); //face 1
        TS_ASSERT_EQUALS(mesh_reader.GetNextFaceData().NodeIndices[1], 36u);    //face 2
        TS_ASSERT_EQUALS(mesh_reader.GetNextFaceData().ContainingElement, 74u); //face 3
        TS_ASSERT_EQUALS(mesh_reader.GetNextFaceData().NodeIndices[2], 16u);    //face 4
        TS_ASSERT_EQUALS(mesh_reader.GetNextFaceData().ContainingElement, 4u);  //face 5
    }

    void TestReadingBinary1D()
    {
        TrianglesMeshReader<1,1> mesh_reader_ascii("mesh/test/data/trivial_1d_mesh");
        //TrianglesMeshWriter<1,1> writer("","trivial_1d_mesh_binary");
        //writer.SetWriteFilesAsBinary();
        //writer.WriteFilesUsingMeshReader(mesh_reader_ascii);
        TS_ASSERT_EQUALS(mesh_reader_ascii.GetNumFaces(), 2u);
        TS_ASSERT_EQUALS(mesh_reader_ascii.GetNumElements(), 10u);
        TS_ASSERT_EQUALS(mesh_reader_ascii.GetNumNodes(), 11u);

        TrianglesMeshReader<1,1> mesh_reader_bin("mesh/test/data/trivial_1d_mesh_binary");
        TS_ASSERT_EQUALS(mesh_reader_bin.GetNumFaces(), 2u);
        TS_ASSERT_EQUALS(mesh_reader_bin.GetNumElements(), 10u);
        TS_ASSERT_EQUALS(mesh_reader_bin.GetNumNodes(), 11u);
    }

    void TestReadingBinary()
    {
        TrianglesMeshReader<3,3> mesh_reader_ascii("mesh/test/data/squashed_cube");
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/squashed_cube_binary");

        TS_ASSERT_EQUALS(mesh_reader.GetNumFaces(), 12u);
        TS_ASSERT_EQUALS(mesh_reader.GetNumElements(), 12u);
        TS_ASSERT_EQUALS(mesh_reader.GetNumNodes(), 9u);
        TS_ASSERT_EQUALS(mesh_reader_ascii.GetNumFaces(), 12u);
        TS_ASSERT_EQUALS(mesh_reader_ascii.GetNumElements(), 12u);
        TS_ASSERT_EQUALS(mesh_reader_ascii.GetNumNodes(), 9u);

        TS_ASSERT( mesh_reader.IsFileFormatBinary() );
        TS_ASSERT( ! mesh_reader_ascii.IsFileFormatBinary() );

        /*
         * Check node locations
         */
        {
            std::vector<double> ascii_location(3u);
            std::vector<double> binary_location(3u);
            for (unsigned i=0; i<mesh_reader.GetNumNodes(); i++)
            {
                // Sequential reading in
                ascii_location = mesh_reader_ascii.GetNextNode();
                binary_location = mesh_reader.GetNextNode();
                TS_ASSERT_DELTA(ascii_location[0],binary_location[0],1e-12);
                TS_ASSERT_DELTA(ascii_location[1],binary_location[1],1e-12);
                TS_ASSERT_DELTA(ascii_location[2],binary_location[2],1e-12);
            }
            mesh_reader_ascii.Reset(); // You wouldn't believe how important this line is.
            for (unsigned i=0; i<mesh_reader.GetNumNodes(); i++)
            {
                // Random access
                ascii_location = mesh_reader_ascii.GetNextNode();
                binary_location = mesh_reader.GetNode(i);
                TS_ASSERT_DELTA(ascii_location[0],binary_location[0],1e-12);
                TS_ASSERT_DELTA(ascii_location[1],binary_location[1],1e-12);
                TS_ASSERT_DELTA(ascii_location[2],binary_location[2],1e-12);
            }

            // Now with the ASCII node iterator
            unsigned count = 0u;
            mesh_reader_ascii.Reset();
            mesh_reader.Reset();
            for (AbstractMeshReader<3,3>::NodeIterator node_it = mesh_reader_ascii.GetNodeIteratorBegin();
                    node_it != mesh_reader_ascii.GetNodeIteratorEnd();
                    ++node_it)
            {
                TS_ASSERT_EQUALS(count, node_it.GetIndex());
                binary_location = mesh_reader.GetNode(count);
                ascii_location = *node_it;
                TS_ASSERT_DELTA(ascii_location[0],binary_location[0],1e-12);
                TS_ASSERT_DELTA(ascii_location[1],binary_location[1],1e-12);
                TS_ASSERT_DELTA(ascii_location[2],binary_location[2],1e-12);
                count++;
            }
            TS_ASSERT_EQUALS(count, mesh_reader.GetNumNodes());
            // Now, again with the binary node iterator
            count = 0u;
            mesh_reader_ascii.Reset();
            mesh_reader.Reset();
            for (AbstractMeshReader<3,3>::NodeIterator node_it = mesh_reader.GetNodeIteratorBegin();
                    node_it != mesh_reader.GetNodeIteratorEnd();
                    ++node_it)
            {
                TS_ASSERT_EQUALS(count, node_it.GetIndex());
                ascii_location = mesh_reader_ascii.GetNextNode();
                binary_location = *node_it;
                TS_ASSERT_DELTA(ascii_location[0],binary_location[0],1e-12);
                TS_ASSERT_DELTA(ascii_location[1],binary_location[1],1e-12);
                TS_ASSERT_DELTA(ascii_location[2],binary_location[2],1e-12);
                count++;
            }
            TS_ASSERT_EQUALS(count, mesh_reader.GetNumNodes());

            // Now we test iterating the binary mesh reader over a subset of nodes only
            std::set<unsigned> even_indices;
            for (unsigned i=0; i<mesh_reader.GetNumNodes(); i += 2)
            {
                even_indices.insert(i);
            }
            mesh_reader_ascii.Reset();
            count = 0u;
            for (AbstractMeshReader<3,3>::NodeIterator node_it = mesh_reader.GetNodeIteratorBegin(even_indices);
                    node_it != mesh_reader.GetNodeIteratorEnd();
                    ++node_it)
            {
                TS_ASSERT_EQUALS(count*2u, node_it.GetIndex());
                count++;
                binary_location = *node_it;
                ascii_location = mesh_reader_ascii.GetNextNode();
                TS_ASSERT_DELTA(ascii_location[0],binary_location[0],1e-12);
                TS_ASSERT_DELTA(ascii_location[1],binary_location[1],1e-12);
                TS_ASSERT_DELTA(ascii_location[2],binary_location[2],1e-12);
                // Skip odd elements
                if (count <= mesh_reader_ascii.GetNumNodes() / 2u)
                {
                    mesh_reader_ascii.GetNextNode();
                }
            }
            TS_ASSERT_EQUALS(count, mesh_reader_ascii.GetNumNodes() / 2u + 1);

            // Also test iterating the ASCII mesh reader over a subset of nodes only
            std::set<unsigned> odd_indices;
            for (unsigned i=1; i<mesh_reader.GetNumNodes(); i += 2)
            {
                odd_indices.insert(i);
            }
            mesh_reader_ascii.Reset();
            mesh_reader.Reset();  //Because we aren't going to do random access for this part of the test
            count = 0u;
            for (AbstractMeshReader<3,3>::NodeIterator node_it = mesh_reader_ascii.GetNodeIteratorBegin(odd_indices);
                    node_it != mesh_reader_ascii.GetNodeIteratorEnd();
                    ++node_it)
            {
                TS_ASSERT_EQUALS(count*2+1, node_it.GetIndex());
                count++;
                ascii_location = *node_it;
                // Skip even elements
                mesh_reader.GetNextNode();
                binary_location = mesh_reader.GetNextNode();
                TS_ASSERT_DELTA(ascii_location[0],binary_location[0],1e-12);
                TS_ASSERT_DELTA(ascii_location[1],binary_location[1],1e-12);
                TS_ASSERT_DELTA(ascii_location[2],binary_location[2],1e-12);
            }
            TS_ASSERT_EQUALS(count, mesh_reader_ascii.GetNumNodes() / 2u);

            // Also test iterating the ASCII mesh reader with empty set
            std::set<unsigned> no_indices;
            mesh_reader_ascii.Reset();
            mesh_reader.Reset();  //Because we aren't going to do random access for this part of the test
            count = 0u;
            for (AbstractMeshReader<3,3>::NodeIterator node_it = mesh_reader_ascii.GetNodeIteratorBegin(no_indices);
                    node_it != mesh_reader_ascii.GetNodeIteratorEnd();
                    ++node_it)
            {
                count++;
            }
            TS_ASSERT_EQUALS(count, 0u);

        }

        /*
         * Check elements
         */
        TS_ASSERT_EQUALS(mesh_reader.GetNumElementAttributes(), 1u);
        TS_ASSERT_EQUALS(mesh_reader_ascii.GetNumElementAttributes(), 1u);
        {
            ElementData ascii_node_indices;
            ElementData binary_node_indices;
            for (unsigned i=0; i<mesh_reader.GetNumElements(); i++)
            {
                // Sequential reading in
                ascii_node_indices = mesh_reader_ascii.GetNextElementData();
                binary_node_indices = mesh_reader.GetNextElementData();
                TS_ASSERT_EQUALS(ascii_node_indices.NodeIndices[0],binary_node_indices.NodeIndices[0]);
                TS_ASSERT_EQUALS(ascii_node_indices.NodeIndices[1],binary_node_indices.NodeIndices[1]);
                TS_ASSERT_EQUALS(ascii_node_indices.NodeIndices[2],binary_node_indices.NodeIndices[2]);
                TS_ASSERT_EQUALS(ascii_node_indices.NodeIndices[3],binary_node_indices.NodeIndices[3]);
                TS_ASSERT_DELTA(ascii_node_indices.AttributeValue,binary_node_indices.AttributeValue,1e-12);
            }
            mesh_reader_ascii.Reset(); // You wouldn't believe how important this line is.
            for (unsigned i=0; i<mesh_reader.GetNumElements(); i++)
            {
                // Random access for binary mesh
                ascii_node_indices = mesh_reader_ascii.GetNextElementData();
                binary_node_indices = mesh_reader.GetElementData(i);
                TS_ASSERT_EQUALS(ascii_node_indices.NodeIndices[0],binary_node_indices.NodeIndices[0]);
                TS_ASSERT_EQUALS(ascii_node_indices.NodeIndices[1],binary_node_indices.NodeIndices[1]);
                TS_ASSERT_EQUALS(ascii_node_indices.NodeIndices[2],binary_node_indices.NodeIndices[2]);
                TS_ASSERT_EQUALS(ascii_node_indices.NodeIndices[3],binary_node_indices.NodeIndices[3]);
                TS_ASSERT_DELTA(ascii_node_indices.AttributeValue,binary_node_indices.AttributeValue,1e-12);
            }

            // Now test using the iterators for the ASCII mesh
            mesh_reader_ascii.Reset();
            mesh_reader.Reset();
            unsigned count = 0u;
            for (AbstractMeshReader<3,3>::ElementIterator elt_it = mesh_reader_ascii.GetElementIteratorBegin();
                    elt_it != mesh_reader_ascii.GetElementIteratorEnd();
                    ++elt_it)
            {
                TS_ASSERT_EQUALS(count, elt_it.GetIndex());
                count++;
                ascii_node_indices = *elt_it;
                binary_node_indices = mesh_reader.GetNextElementData();
                TS_ASSERT_EQUALS(ascii_node_indices.NodeIndices[0],binary_node_indices.NodeIndices[0]);
                TS_ASSERT_EQUALS(ascii_node_indices.NodeIndices[1],binary_node_indices.NodeIndices[1]);
                TS_ASSERT_EQUALS(ascii_node_indices.NodeIndices[2],binary_node_indices.NodeIndices[2]);
                TS_ASSERT_EQUALS(ascii_node_indices.NodeIndices[3],binary_node_indices.NodeIndices[3]);
                TS_ASSERT_DELTA(ascii_node_indices.AttributeValue,binary_node_indices.AttributeValue,1e-12);
            }
            TS_ASSERT_EQUALS(count, mesh_reader.GetNumElements());

            // Now test using the iterator for the binary mesh
            mesh_reader_ascii.Reset();
            count = 0u;
            for (AbstractMeshReader<3,3>::ElementIterator elt_it = mesh_reader.GetElementIteratorBegin();
                    elt_it != mesh_reader.GetElementIteratorEnd();
                    ++elt_it)
            {
                TS_ASSERT_EQUALS(count, elt_it.GetIndex());
                count++;
                ascii_node_indices = mesh_reader_ascii.GetNextElementData();
                binary_node_indices = *elt_it;
                TS_ASSERT_EQUALS(ascii_node_indices.NodeIndices[0],binary_node_indices.NodeIndices[0]);
                TS_ASSERT_EQUALS(ascii_node_indices.NodeIndices[1],binary_node_indices.NodeIndices[1]);
                TS_ASSERT_EQUALS(ascii_node_indices.NodeIndices[2],binary_node_indices.NodeIndices[2]);
                TS_ASSERT_EQUALS(ascii_node_indices.NodeIndices[3],binary_node_indices.NodeIndices[3]);
                TS_ASSERT_DELTA(ascii_node_indices.AttributeValue,binary_node_indices.AttributeValue,1e-12);
            }
            TS_ASSERT_EQUALS(count, mesh_reader_ascii.GetNumElements());

            // Now we test iterating the binary mesh reader over a subset of elements only
            std::set<unsigned> even_indices;
            for (unsigned i=0; i<mesh_reader.GetNumElements(); i += 2)
            {
                even_indices.insert(i);
            }
            mesh_reader_ascii.Reset();
            count = 0u;
            for (AbstractMeshReader<3,3>::ElementIterator elt_it = mesh_reader.GetElementIteratorBegin(even_indices);
                    elt_it != mesh_reader.GetElementIteratorEnd();
                    ++elt_it)
            {
                TS_ASSERT_EQUALS(count*2u, elt_it.GetIndex());
                count++;
                binary_node_indices = *elt_it;
                ascii_node_indices = mesh_reader_ascii.GetNextElementData();
                TS_ASSERT_EQUALS(ascii_node_indices.NodeIndices[0],binary_node_indices.NodeIndices[0]);
                TS_ASSERT_EQUALS(ascii_node_indices.NodeIndices[1],binary_node_indices.NodeIndices[1]);
                TS_ASSERT_EQUALS(ascii_node_indices.NodeIndices[2],binary_node_indices.NodeIndices[2]);
                TS_ASSERT_EQUALS(ascii_node_indices.NodeIndices[3],binary_node_indices.NodeIndices[3]);
                TS_ASSERT_DELTA(ascii_node_indices.AttributeValue,binary_node_indices.AttributeValue,1e-12);
                // Skip odd elements
                mesh_reader_ascii.GetNextElementData();
            }
            TS_ASSERT_EQUALS(count, mesh_reader_ascii.GetNumElements() / 2u);

            // Also test iterating the ASCII mesh reader over a subset of elements only
            std::set<unsigned> odd_indices;
            for (unsigned i=1; i<mesh_reader.GetNumElements(); i += 2)
            {
                odd_indices.insert(i);
            }
            mesh_reader_ascii.Reset();
            mesh_reader.Reset();  //Because we aren't going to do random access for this part of the test
            count = 0u;
            for (AbstractMeshReader<3,3>::ElementIterator elt_it = mesh_reader_ascii.GetElementIteratorBegin(odd_indices);
                    elt_it != mesh_reader_ascii.GetElementIteratorEnd();
                    ++elt_it)
            {
                TS_ASSERT_EQUALS(count*2+1, elt_it.GetIndex());
                count++;
                ascii_node_indices = *elt_it;
                // Skip even elements
                mesh_reader.GetNextElementData();
                binary_node_indices = mesh_reader.GetNextElementData();
                TS_ASSERT_EQUALS(ascii_node_indices.NodeIndices[0],binary_node_indices.NodeIndices[0]);
                TS_ASSERT_EQUALS(ascii_node_indices.NodeIndices[1],binary_node_indices.NodeIndices[1]);
                TS_ASSERT_EQUALS(ascii_node_indices.NodeIndices[2],binary_node_indices.NodeIndices[2]);
                TS_ASSERT_EQUALS(ascii_node_indices.NodeIndices[3],binary_node_indices.NodeIndices[3]);
                TS_ASSERT_DELTA(ascii_node_indices.AttributeValue, binary_node_indices.AttributeValue,1e-12);
            }
            TS_ASSERT_EQUALS(count, mesh_reader_ascii.GetNumElements() / 2u);

            // Also test iterating the ASCII mesh reader with empty set
            std::set<unsigned> no_indices;
            mesh_reader_ascii.Reset();
            mesh_reader.Reset();  //Because we aren't going to do random access for this part of the test
            count = 0u;
            for (AbstractMeshReader<3,3>::ElementIterator elt_it = mesh_reader_ascii.GetElementIteratorBegin(no_indices);
                    elt_it != mesh_reader_ascii.GetElementIteratorEnd();
                    ++elt_it)
            {
                count++;
            }
            TS_ASSERT_EQUALS(count, 0u);
        }

        /*
         * Check faces
         */
        ///\todo #1949 the binary writer NEVER writes face attributes
        TS_ASSERT_EQUALS(mesh_reader.GetNumFaceAttributes(), 0u);
        TS_ASSERT_EQUALS(mesh_reader_ascii.GetNumFaceAttributes(), 0u);
        {
            ElementData ascii_node_indices;
            ElementData binary_node_indices;
            for (unsigned i=0; i<mesh_reader.GetNumFaces(); i++)
            {
                // Sequential reading in
                ascii_node_indices = mesh_reader_ascii.GetNextFaceData();
                binary_node_indices = mesh_reader.GetNextFaceData();
                TS_ASSERT_EQUALS(ascii_node_indices.NodeIndices[0],binary_node_indices.NodeIndices[0]);
                TS_ASSERT_EQUALS(ascii_node_indices.NodeIndices[1],binary_node_indices.NodeIndices[1]);
                TS_ASSERT_EQUALS(ascii_node_indices.NodeIndices[2],binary_node_indices.NodeIndices[2]);
                TS_ASSERT_DELTA(ascii_node_indices.AttributeValue,binary_node_indices.AttributeValue,1e-12);
            }
            mesh_reader_ascii.Reset(); // You wouldn't believe how important this line is.
            for (unsigned i=0; i<mesh_reader.GetNumFaces(); i++)
            {
                // Sequential reading in
                ascii_node_indices = mesh_reader_ascii.GetNextFaceData();
                binary_node_indices = mesh_reader.GetFaceData(i);
                TS_ASSERT_EQUALS(ascii_node_indices.NodeIndices[0],binary_node_indices.NodeIndices[0]);
                TS_ASSERT_EQUALS(ascii_node_indices.NodeIndices[1],binary_node_indices.NodeIndices[1]);
                TS_ASSERT_EQUALS(ascii_node_indices.NodeIndices[2],binary_node_indices.NodeIndices[2]);
                TS_ASSERT_DELTA(ascii_node_indices.AttributeValue,binary_node_indices.AttributeValue,1e-12);
            }
        }

        TS_ASSERT_THROWS_THIS(mesh_reader.GetNode(9u), "Node does not exist - not enough nodes.");
        TS_ASSERT_THROWS_THIS(mesh_reader.GetFaceData(12u), "Face does not exist - not enough faces.");
        TS_ASSERT_THROWS_THIS(mesh_reader.GetElementData(12u), "Element 12 does not exist - not enough elements (only 12).");

        TS_ASSERT_THROWS_THIS(mesh_reader.GetContainingElementIndices(0u), "No NCL file available for this mesh.");
    }

    void TestReadingHexMesh()
    {
        TrianglesMeshReader<2,2> mesh_reader_2d("mesh/test/data/rectangle_608_hexa_elements");

        TS_ASSERT_EQUALS(mesh_reader_2d.GetNumNodes(), 663u);
        TS_ASSERT_EQUALS(mesh_reader_2d.GetNumElements(), 608u);
        TS_ASSERT_EQUALS(mesh_reader_2d.GetNumFaces(), 108u);

        TrianglesMeshReader<3,3> mesh_reader_3d("mesh/test/data/cuboid_140_hexa_elements");

        TS_ASSERT_EQUALS(mesh_reader_3d.GetNumNodes(), 240u);
        TS_ASSERT_EQUALS(mesh_reader_3d.GetNumElements(), 140u);
        TS_ASSERT_EQUALS(mesh_reader_3d.GetNumFaces(), 166u);
    }

    void TestReadingWithGenericReader()
    {
        std::shared_ptr<AbstractMeshReader<2,2> > p_mesh_reader = GenericMeshReader<2,2>("mesh/test/data/disk_522_elements");
        TS_ASSERT_EQUALS(p_mesh_reader->GetNumNodes(), 312u);
        TS_ASSERT_EQUALS(p_mesh_reader->GetNumElements(), 522u);
        TS_ASSERT_EQUALS(p_mesh_reader->GetNumFaceAttributes(), 1u);
        TS_ASSERT_EQUALS(p_mesh_reader->GetNumElementAttributes(), 0u);

        TS_ASSERT_THROWS_CONTAINS((GenericMeshReader<2,2>("mesh/test/data/no_such_file")),
                                  "Could not open appropriate mesh files for mesh/test/data/no_such_file");
    }

    void TestReadingNclInformation()
    {
        TrianglesMeshReader<3,3> mesh_reader_3d("mesh/test/data/simple_cube_binary");

        TS_ASSERT(mesh_reader_3d.HasNclFile());

        std::vector<unsigned> containing_element_indices = mesh_reader_3d.GetContainingElementIndices(0);
        TS_ASSERT_EQUALS(containing_element_indices.size(), 5u)
        TS_ASSERT_EQUALS(containing_element_indices[0], 0u);
        TS_ASSERT_EQUALS(containing_element_indices[1], 4u);
        TS_ASSERT_EQUALS(containing_element_indices[2], 8u);
        TS_ASSERT_EQUALS(containing_element_indices[3], 10u);
        TS_ASSERT_EQUALS(containing_element_indices[4], 11u);

        containing_element_indices = mesh_reader_3d.GetContainingElementIndices(7);
        TS_ASSERT_EQUALS(containing_element_indices.size(), 5u)
        TS_ASSERT_EQUALS(containing_element_indices[0], 3u);
        TS_ASSERT_EQUALS(containing_element_indices[1], 4u);
        TS_ASSERT_EQUALS(containing_element_indices[2], 5u);
        TS_ASSERT_EQUALS(containing_element_indices[3], 9u);
        TS_ASSERT_EQUALS(containing_element_indices[4], 10u);

        // Out of range
        TS_ASSERT_THROWS_THIS(containing_element_indices = mesh_reader_3d.GetContainingElementIndices(9), "Connectivity list does not exist - not enough nodes.");

        // Test exception if NCL file has wrong number of nodes
        TS_ASSERT_THROWS_THIS( READER_3D bad_mesh_reader("mesh/test/data/simple_cube_binary_ncl_corrupted"),
                           "NCL file does not contain the correct number of nodes for mesh");

    }

    void TestReading3dMeshWithPermutation()
    {
        const unsigned ELEMENT_DIM = 3;
        const unsigned SPACE_DIM = 3;

        TrianglesMeshReader<ELEMENT_DIM,SPACE_DIM> mesh_reader_3d_ascii("mesh/test/data/simple_cube");
        TrianglesMeshReader<ELEMENT_DIM,SPACE_DIM> mesh_reader_3d("mesh/test/data/simple_cube_binary");
        TrianglesMeshReader<ELEMENT_DIM,SPACE_DIM> mesh_reader_3d_permuted("mesh/test/data/simple_cube_binary");

        unsigned num_nodes = mesh_reader_3d.GetNumNodes();
        std::vector<unsigned> permutation(num_nodes);
        for (unsigned node_index=0; node_index < num_nodes; node_index++)
        {
            permutation[node_index] = num_nodes-node_index-1;
        }

        TS_ASSERT_THROWS_THIS(mesh_reader_3d_ascii.SetNodePermutation(permutation),
                              "Permuted read can only be used with binary files since it requires random access to the node file.");

        TS_ASSERT_EQUALS(mesh_reader_3d_permuted.HasNodePermutation(), false);
        TS_ASSERT_EQUALS(mesh_reader_3d_permuted.rGetNodePermutation().size(), 0u);

        mesh_reader_3d_permuted.SetNodePermutation(permutation);

        TS_ASSERT_EQUALS(mesh_reader_3d_permuted.HasNodePermutation(), true);
        TS_ASSERT_EQUALS(mesh_reader_3d_permuted.rGetNodePermutation().size(), 9u);
        TS_ASSERT_EQUALS(mesh_reader_3d_permuted.rGetNodePermutation()[8], 0u);
        TS_ASSERT_EQUALS(mesh_reader_3d_permuted.rGetNodePermutation()[0], 8u);

        for (unsigned node_index=0; node_index < num_nodes; node_index++)
        {
            for (unsigned dimension=0; dimension<SPACE_DIM; dimension++)
            {
                TS_ASSERT_EQUALS( mesh_reader_3d.GetNode(permutation[node_index])[dimension],
                                  mesh_reader_3d_permuted.GetNode(node_index)[dimension]);
            }
        }

        for (unsigned element_index=0; element_index < mesh_reader_3d.GetNumElements(); element_index++)
        {
            for (unsigned local_node_index = 0; local_node_index < ELEMENT_DIM+1; local_node_index++)
            {
                unsigned original_mesh_global_index = mesh_reader_3d.GetElementData(element_index).NodeIndices[local_node_index];
                unsigned permuted_mesh_global_index = mesh_reader_3d_permuted.GetElementData(element_index).NodeIndices[local_node_index];
                TS_ASSERT_EQUALS(permutation[original_mesh_global_index], permuted_mesh_global_index);
            }
        }

        for (unsigned boundary_ele_index=0; boundary_ele_index < mesh_reader_3d.GetNumElements(); boundary_ele_index++)
        {
            for (unsigned local_node_index = 0; local_node_index < ELEMENT_DIM; local_node_index++)
            {
                unsigned original_mesh_global_index = mesh_reader_3d.GetFaceData(boundary_ele_index).NodeIndices[local_node_index];
                unsigned permuted_mesh_global_index = mesh_reader_3d_permuted.GetFaceData(boundary_ele_index).NodeIndices[local_node_index];
                TS_ASSERT_EQUALS(permutation[original_mesh_global_index], permuted_mesh_global_index);
            }
        }
    }

    void TestReadingMissingAttributes()
    {
        // The reader immediately reads and caches face data so missing attributes
        // immediately cause an Exception to be thrown
        TS_ASSERT_THROWS_CONTAINS(READER_2D mesh_reader("mesh/test/data/baddata/canonical_triangle_missing_edge_attribute"),"Error in reading attribute");

        // The reader doesn't read node data until GetNextNode() is called
        READER_2D mesh_reader("mesh/test/data/baddata/canonical_triangle_missing_node_attribute");
        TS_ASSERT_EQUALS(mesh_reader.GetNumNodes(), 3u);

        // Read the data for the first node successfully
        mesh_reader.GetNextNode();
        TS_ASSERT_EQUALS(mesh_reader.GetNodeAttributes().size(), 2u);
        TS_ASSERT_DELTA(mesh_reader.GetNodeAttributes()[0], 8.234, 1e-6);
        TS_ASSERT_DELTA(mesh_reader.GetNodeAttributes()[1], 25.4,  1e-6);

        // The second node has a missing attribute
        TS_ASSERT_THROWS_CONTAINS(mesh_reader.GetNextNode(),"Error in reading attribute");
    }
};

#endif //_TESTTRIANGLESMESHREADER_HPP_
