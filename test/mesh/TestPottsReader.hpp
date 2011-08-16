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

#ifndef TESTPOTTSREADER_HPP_
#define TESTPOTTSREADER_HPP_

#include <cxxtest/TestSuite.h>

#include <string>
#include <fstream>

#include "PottsMeshReader.hpp"

class TestPottsReader : public CxxTest::TestSuite
{
public:
    /**
     * Check that input files are opened correctly.
     */
    void TestFilesOpen() throw(Exception)
    {
        PottsMeshReader<2> mesh_reader("notforrelease_cell_based/test/data/TestPottsMeshWriter/potts_mesh_2d");
    }

    /**
    * Check that the nodes are read correctly. Checks that the output vector
    * for a given input file is the correct length and that if the input file
    * is corrupted (missing nodes) then an exception is thrown.
    */
   void TestNodesDataRead() throw(Exception)
   {
       PottsMeshReader<2> mesh_reader("notforrelease_cell_based/test/data/TestPottsMeshWriter/potts_mesh_2d");

       TS_ASSERT_EQUALS(mesh_reader.GetNumNodes(), 6u);

       PottsMeshReader<2> mesh_reader2("notforrelease_cell_based/test/data/baddata/potts_mesh_bad_nodes");

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
       PottsMeshReader<2> mesh_reader("notforrelease_cell_based/test/data/TestPottsMeshWriter/potts_mesh_2d");

       TS_ASSERT_EQUALS(mesh_reader.GetNumElements(), 2u);

       // Read element 0 from file
       ElementData data = mesh_reader.GetNextElementData();

       TS_ASSERT_EQUALS(data.NodeIndices.size(), 3u);
       TS_ASSERT_EQUALS(data.NodeIndices[0], 0u);
       TS_ASSERT_EQUALS(data.NodeIndices[1], 1u);
       TS_ASSERT_EQUALS(data.NodeIndices[2], 3u);

       // Read element 1 from file
       ElementData data2 = mesh_reader.GetNextElementData();

       TS_ASSERT_EQUALS(data2.NodeIndices.size(), 3u);
       TS_ASSERT_EQUALS(data2.NodeIndices[0], 2u);
       TS_ASSERT_EQUALS(data2.NodeIndices[1], 4u);
       TS_ASSERT_EQUALS(data2.NodeIndices[2], 5u);

       TS_ASSERT_EQUALS(mesh_reader.GetNumElementAttributes(), 0u);

       mesh_reader.Reset();
       for (unsigned i=1; i<mesh_reader.GetNumElements(); i++)
       {
           ElementData data = mesh_reader.GetNextElementData();
           TS_ASSERT_EQUALS(data.AttributeValue, 0u);
       }

       PottsMeshReader<2> mesh_reader2("notforrelease_cell_based/test/data/baddata/potts_mesh_bad_elements");

       // Reads element 0 from file
       TS_ASSERT_THROWS_NOTHING(mesh_reader2.GetNextElementData());

       // Reads element 2 from file when expecting number 1
       TS_ASSERT_THROWS_THIS(mesh_reader2.GetNextElementData(), "Data for element 1 missing");
   }
};

#endif // TESTPOTTSREADER_HPP_
