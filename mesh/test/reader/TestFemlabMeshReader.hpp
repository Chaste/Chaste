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


#ifndef _TESTFEMLABMESHREADER_HPP_
#define _TESTFEMLABMESHREADER_HPP_

#include <cxxtest/TestSuite.h>
#include "FemlabMeshReader.hpp"

typedef FemlabMeshReader<2,2> READER_2D;
typedef FemlabMeshReader<1,1> READER_1D;

class TestFemlabMeshReaders : public CxxTest::TestSuite
{
public:

    AbstractMeshReader<2,2>* mpFemlabMeshReader;

    /**
     * Check that input files are opened correctly.
     */
    void TestFilesOpen()
    {
        TS_ASSERT_THROWS_NOTHING(mpFemlabMeshReader = new READER_2D("mesh/test/data/",
                                                                    "femlab_lshape_nodes.dat",
                                                                    "femlab_lshape_elements.dat",
                                                                    "femlab_lshape_edges.dat"));

        delete mpFemlabMeshReader;

        // Coverage test

        TS_ASSERT_THROWS_THIS(new READER_1D("mesh/test/data/",
                                            "femlab_lshape_nodes.dat",
                                            "femlab_lshape_elements.dat",
                                            "femlab_lshape_edges.dat"),
                              "SPACE_DIM  != dimension read from file");
    }

    /**
     * Check that the nodes are read correctly. Checks that the output vector
     * for a given input file is the correct length and that if the input file
     * is corrupted (missing nodes) then an exception is thrown.
     *
     */
    void TestNodesDataRead()
    {
        mpFemlabMeshReader = new FemlabMeshReader<2,2>("mesh/test/data/",
                                                       "femlab_lshape_nodes.dat",
                                                       "femlab_lshape_elements.dat",
                                                       "femlab_lshape_edges.dat");

        TS_ASSERT_EQUALS(mpFemlabMeshReader->GetNumNodes(), 151u);

        delete mpFemlabMeshReader;
    }

    /**
     * Check that the elements are read correctly. Checks that the output vector
     * for a given input file is the correct length and that if the input file
     * is corrupted (missing elements) then an exception is thrown.
     *
     */
    void TestElementsDataRead()
    {
        mpFemlabMeshReader = new FemlabMeshReader<2,2>("mesh/test/data/",
                                                       "femlab_lshape_nodes.dat",
                                                       "femlab_lshape_elements.dat",
                                                       "femlab_lshape_edges.dat");

        TS_ASSERT_EQUALS(mpFemlabMeshReader->GetNumElements(), 260u);
        TS_ASSERT_EQUALS(mpFemlabMeshReader->GetNumElementAttributes(), 0u);

        delete mpFemlabMeshReader;
    }

    /**
     * Check that the faces are read correctly. Checks that the output vector
     * for a given input file is the correct length and that if the input file
     * is corrupted (missing faces) then an exception is thrown.
     *
     */
    void TestFacesDataRead()
    {
        mpFemlabMeshReader = new FemlabMeshReader<2,2>("mesh/test/data/",
                                                       "femlab_lshape_nodes.dat",
                                                       "femlab_lshape_elements.dat",
                                                       "femlab_lshape_edges.dat");

        TS_ASSERT_EQUALS(mpFemlabMeshReader->GetNumFaces(), 40u);

        // Coverage of AbstractCachedMeshReader::GetNumEdges()
        TS_ASSERT_EQUALS(mpFemlabMeshReader->GetNumEdges(), 40u);

        delete mpFemlabMeshReader;
    }

    /**
     * Check that GetNextNode() returns the coordinates of the correct node.
     * Compares the coordinates of the first two nodes with their known
     * values, checks that no errors are thrown for the remaining nodes and
     * that an error is thrown if we try to call the function too many times.
     *
     */
    void TestGetNextNode()
    {
        mpFemlabMeshReader = new FemlabMeshReader<2,2>("mesh/test/data/",
                                                       "femlab_lshape_nodes.dat",
                                                       "femlab_lshape_elements.dat",
                                                       "femlab_lshape_edges.dat");

        std::vector<double> first_node;

        first_node = mpFemlabMeshReader->GetNextNode();

        TS_ASSERT_DELTA(first_node[0], 0.0, 1e-6);
        TS_ASSERT_DELTA(first_node[1], 1.0, 1e-6);

        std::vector<double> next_node;

        next_node = mpFemlabMeshReader->GetNextNode();

        TS_ASSERT_DELTA(next_node[0], 0.5, 1e-6);
        TS_ASSERT_DELTA(next_node[1], 1.0, 1e-6);

        for (unsigned i=2; i<mpFemlabMeshReader->GetNumNodes(); i++)
        {
            TS_ASSERT_THROWS_NOTHING(next_node = mpFemlabMeshReader->GetNextNode());
        }

        TS_ASSERT_THROWS_THIS(next_node = mpFemlabMeshReader->GetNextNode(),"All nodes already got");

        // Coverage of a default implementation of GetNodeAttributes in abstract class that returns an empty vector.
        std::vector<double> no_attributes = mpFemlabMeshReader->GetNodeAttributes();
        TS_ASSERT_EQUALS(no_attributes.size(), 0u);

        delete mpFemlabMeshReader;
    }

    /**
     * Check that GetNextElementData() works. Checks that no errors are thrown for
     * all of the elements and that an error is thrown if we try to call the
     * function too many times.
     */
    void TestGetNextElementData()
    {
        mpFemlabMeshReader = new FemlabMeshReader<2,2>("mesh/test/data/",
                                                       "femlab_lshape_nodes.dat",
                                                       "femlab_lshape_elements.dat",
                                                       "femlab_lshape_edges.dat");

        ElementData first_element_data = mpFemlabMeshReader->GetNextElementData();

        std::vector<unsigned> first_element_nodes = first_element_data.NodeIndices;

        TS_ASSERT_EQUALS(first_element_data.NodeIndices[0], 15u);
        TS_ASSERT_EQUALS(first_element_data.NodeIndices[1], 3u);
        TS_ASSERT_EQUALS(first_element_data.NodeIndices[2], 62u);
        TS_ASSERT_EQUALS(first_element_data.AttributeValue, 0u);

        ElementData next_element_data = mpFemlabMeshReader->GetNextElementData();

        TS_ASSERT_EQUALS(next_element_data.NodeIndices[0], 8u);
        TS_ASSERT_EQUALS(next_element_data.NodeIndices[1], 0u);
        TS_ASSERT_EQUALS(next_element_data.NodeIndices[2], 53u);
        TS_ASSERT_EQUALS(next_element_data.AttributeValue, 0u);

        for (unsigned i=2; i<mpFemlabMeshReader->GetNumElements(); i++)
        {
            TS_ASSERT_THROWS_NOTHING(next_element_data = mpFemlabMeshReader->GetNextElementData());
        }

        TS_ASSERT_THROWS_THIS(next_element_data = mpFemlabMeshReader->GetNextElementData(),"All elements already got");

        delete mpFemlabMeshReader;
    }

    /**
     * Check that GetNextFace() works. Checks that no errors are thrown for
     * all of the elements and that an error is thrown if we try to call the
     * function too many times.
     *
     */
    void TestGetNextFace()
    {
        mpFemlabMeshReader = new FemlabMeshReader<2,2>("mesh/test/data/",
                                                       "femlab_lshape_nodes.dat",
                                                       "femlab_lshape_elements.dat",
                                                       "femlab_lshape_edges.dat");

        std::vector<unsigned> first_face;

        first_face = mpFemlabMeshReader->GetNextFaceData().NodeIndices;

        TS_ASSERT_EQUALS(first_face[0], 0u);
        TS_ASSERT_EQUALS(first_face[1], 8u);

        std::vector<unsigned> next_face;

        next_face = mpFemlabMeshReader->GetNextFaceData().NodeIndices;

        TS_ASSERT_EQUALS(next_face[0], 8u);
        TS_ASSERT_EQUALS(next_face[1], 9u);

        for (unsigned i=2; i<mpFemlabMeshReader->GetNumFaces(); i++)
        {
            TS_ASSERT_THROWS_NOTHING(next_face = mpFemlabMeshReader->GetNextFaceData().NodeIndices);
        }

        TS_ASSERT_THROWS_THIS(next_face = mpFemlabMeshReader->GetNextFaceData().NodeIndices,
                "All faces (or edges) already got");

        delete mpFemlabMeshReader;
    }

};

#endif //_TESTFEMLABMESHREADER_HPP_
