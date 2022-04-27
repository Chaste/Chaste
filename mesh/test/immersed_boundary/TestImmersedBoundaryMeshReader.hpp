/*

Copyright (c) 2005-2018, University of Oxford.
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

#ifndef TESTIMMERSEDBOUNDARYMESHREADER_HPP_
#define TESTIMMERSEDBOUNDARYMESHREADER_HPP_

// Needed for test framework
#include <cxxtest/TestSuite.h>

// Includes from projects/ImmersedBoundary
#include "ImmersedBoundaryMesh.hpp"
#include "ImmersedBoundaryMeshWriter.hpp"
#include "ImmersedBoundaryMeshReader.hpp"
#include "ImmersedBoundaryPalisadeMeshGenerator.hpp"
#include "SuperellipseGenerator.hpp"

// This test is never run in parallel
#include "FakePetscSetup.hpp"

using Reader2D = ImmersedBoundaryMeshReader<2, 2>;

class TestImmersedBoundaryMeshReader : public CxxTest::TestSuite
{
public:

    void TestNonExistentFile()
    {
      TS_ASSERT_THROWS_ANYTHING(Reader2D ibReader("NONEXISTENT_FILE"));
    }
    
    void TestReadingFileHeaders2D()
    {
      Reader2D reader("mesh/test/data/ib_mesh_2d");
      TS_ASSERT_EQUALS(reader.GetNumNodes(), 22);
      TS_ASSERT_EQUALS(reader.GetNumElements(), 3);
      TS_ASSERT_EQUALS(reader.GetNumLaminas(), 2);
      TS_ASSERT_EQUALS(reader.GetNumGridPtsX(), 8);
      TS_ASSERT_EQUALS(reader.GetNumGridPtsY(), 8);
      TS_ASSERT_EQUALS(reader.GetNumElementAttributes(), 1);
      TS_ASSERT_EQUALS(reader.GetNumLaminaAttributes(), 1);
      TS_ASSERT_EQUALS(reader.GetCharacteristicNodeSpacing(), 0.0874862);
    }
    
    void TestReadingData2D()
    {
      Reader2D reader("mesh/test/data/ib_mesh_2d");
      
      // Reading node data
      std::vector<double> nodeData = reader.GetNextNode();
      std::vector<double> expectedNodeData {0, 0, 1};
      TS_ASSERT_EQUALS(nodeData, expectedNodeData);
      
      // Reading grid row
      std::vector<double> gridRowData = reader.GetNextGridRow();
      std::vector<double> expectedGridRowData {0, 0, 0, 0, 0, 0, 0, 0};
      TS_ASSERT_EQUALS(gridRowData, expectedGridRowData);
      
      ImmersedBoundaryElementData elementData = reader.GetNextImmersedBoundaryElementData();
      ImmersedBoundaryElementData expectedElementData {{0, 1, 2}, 0};
      TS_ASSERT_EQUALS(elementData.NodeIndices, expectedElementData.NodeIndices);
      TS_ASSERT_EQUALS(elementData.AttributeValue, expectedElementData.AttributeValue);
  
      ImmersedBoundaryElementData laminaData = reader.GetNextImmersedBoundaryLaminaData();
      ImmersedBoundaryElementData expectedLaminaData {{15, 16, 17, 18, 19}, 0};
      TS_ASSERT_EQUALS(laminaData.NodeIndices, expectedLaminaData.NodeIndices);
      TS_ASSERT_EQUALS(laminaData.AttributeValue, expectedLaminaData.AttributeValue);
    }
    
    void TestReset()
    {
      Reader2D reader("mesh/test/data/ib_mesh_2d");
      
      // Reading node data
      std::vector<double> nodeData = reader.GetNextNode();
      std::vector<double> expectedResult {0, 0, 1};
      TS_ASSERT_EQUALS(nodeData, expectedResult);

      reader.Reset();

      nodeData = reader.GetNextNode();
      TS_ASSERT_EQUALS(nodeData, expectedResult);
    }
    
    void TestConstructingMeshFromReader()
    {
      Reader2D reader("mesh/test/data/ib_mesh_2d");

      ImmersedBoundaryMesh<2, 2> mesh;
      mesh.ConstructFromMeshReader(reader);

      TS_ASSERT_EQUALS(mesh.GetNumNodes(), 22);
      TS_ASSERT_EQUALS(mesh.GetNumElements(), 3);
      TS_ASSERT_EQUALS(mesh.GetNumLaminas(), 2);
    }
};

#endif /*TESTIMMERSEDBOUNDARYMESHREADER_HPP_*/