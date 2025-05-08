/*

Copyright (c) 2005-2025, University of Oxford.
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
        TS_ASSERT_THROWS_ANYTHING(Reader2D ib_reader("NONEXISTENT_FILE"));
        TS_ASSERT_THROWS_ANYTHING(Reader2D missing_elem("mesh/test/data/ib_mesh/missing_elem"));
        TS_ASSERT_THROWS_ANYTHING(Reader2D missing_node("mesh/test/data/ib_mesh/missing_node"));
        TS_ASSERT_THROWS_ANYTHING(Reader2D missing_grid("mesh/test/data/ib_mesh/missing_grid"));
        TS_ASSERT_THROWS_ANYTHING(Reader2D missing_lam("mesh/test/data/ib_mesh/missing_lam"));

        Reader2D missing_elem("mesh/test/data/ib_mesh/missing_elem_data");
        TS_ASSERT_THROWS_ANYTHING(missing_elem.GetNextImmersedBoundaryElementData());

        TS_ASSERT_THROWS_ANYTHING(Reader2D missing_node("mesh/test/data/ib_mesh/missing_node_data"));

        Reader2D missing_grid("mesh/test/data/ib_mesh/missing_grid_data");
        TS_ASSERT_THROWS_ANYTHING(missing_grid.GetNextGridRow());

        Reader2D missing_lam("mesh/test/data/ib_mesh/missing_lam_data");
        TS_ASSERT_THROWS_ANYTHING(missing_lam.GetNextImmersedBoundaryLaminaData());
    }

    void TestReadingFileHeaders2D()
    {
        Reader2D reader("mesh/test/data/ib_mesh_2d");
        TS_ASSERT_EQUALS(reader.GetNumNodes(), 22u);
        TS_ASSERT_EQUALS(reader.GetNumElements(), 3u);
        TS_ASSERT_EQUALS(reader.GetNumLaminas(), 2u);
        TS_ASSERT_EQUALS(reader.GetNumGridPtsX(), 8u);
        TS_ASSERT_EQUALS(reader.GetNumGridPtsY(), 8u);
        TS_ASSERT_EQUALS(reader.GetNumElementAttributes(), 1u);
        TS_ASSERT_EQUALS(reader.GetNumLaminaAttributes(), 1u);
        TS_ASSERT_DELTA(reader.GetCharacteristicNodeSpacing(), 0.0874862, 1e-6);
    }

    void TestReadingData2D()
    {
        Reader2D reader("mesh/test/data/ib_mesh_2d");

        // Reading node data
        std::vector<double> node_data = reader.GetNextNode();
        std::vector<double> expected_node_data {0, 0, 1};
        TS_ASSERT_EQUALS(node_data, expected_node_data);

        // Reading grid row
        std::vector<double> grid_row_data = reader.GetNextGridRow();
        std::vector<double> expected_grid_row_data {0, 0, 0, 0, 0, 0, 0, 0};
        TS_ASSERT_EQUALS(grid_row_data, expected_grid_row_data);

        ImmersedBoundaryElementData element_data = reader.GetNextImmersedBoundaryElementData();
        ImmersedBoundaryElementData expected_element_data {{0, 1, 2}, 0, false, 0, {}, 0.0, false};
        TS_ASSERT_EQUALS(element_data.NodeIndices, expected_element_data.NodeIndices);
        TS_ASSERT_EQUALS(element_data.AttributeValue, expected_element_data.AttributeValue);

        ImmersedBoundaryElementData lamina_data = reader.GetNextImmersedBoundaryLaminaData();
        ImmersedBoundaryElementData expected_lamina_data {{15, 16, 17, 18, 19}, 1, false, 0, {}, 0.0, false};
        TS_ASSERT_EQUALS(lamina_data.NodeIndices, expected_lamina_data.NodeIndices);
        TS_ASSERT_EQUALS(lamina_data.AttributeValue, expected_lamina_data.AttributeValue);
    }

    void TestReset()
    {
        Reader2D reader("mesh/test/data/ib_mesh_2d");

        // Reading node data
        std::vector<double> node_data = reader.GetNextNode();
        std::vector<double> expected_result {0, 0, 1};
        TS_ASSERT_EQUALS(node_data, expected_result);

        reader.Reset();

        node_data = reader.GetNextNode();
        TS_ASSERT_EQUALS(node_data, expected_result);
    }

    void TestConstructingMeshFromReader()
    {
        Reader2D reader("mesh/test/data/ib_mesh_2d");

        ImmersedBoundaryMesh<2, 2> mesh;
        mesh.ConstructFromMeshReader(reader);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 22u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 3u);
        TS_ASSERT_EQUALS(mesh.GetNumLaminas(), 2u);
    }

    void TestNoAttributes()
    {
        Reader2D reader("mesh/test/data/ib_mesh_2d_no_attributes");

        ImmersedBoundaryMesh<2, 2> mesh;
        mesh.ConstructFromMeshReader(reader);

        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 22u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 3u);
        TS_ASSERT_EQUALS(mesh.GetNumLaminas(), 2u);
    }

    void TestMissingData()
    {
        {
            Reader2D reader("mesh/test/data/ib_mesh_2d_missing_elem");
            ImmersedBoundaryMesh<2, 2> mesh;
            TS_ASSERT_THROWS_ANYTHING(mesh.ConstructFromMeshReader(reader));
        }
        {
            Reader2D reader("mesh/test/data/ib_mesh_2d_missing_lam");
            ImmersedBoundaryMesh<2, 2> mesh;
            TS_ASSERT_THROWS_ANYTHING(mesh.ConstructFromMeshReader(reader));
        }
        {
            Reader2D reader("mesh/test/data/ib_mesh_2d_missing_node");
            ImmersedBoundaryMesh<2, 2> mesh;
            TS_ASSERT_THROWS_ANYTHING(mesh.ConstructFromMeshReader(reader));
        }
    }

    void TestOtherMethods()
    {
        Reader2D reader("mesh/test/data/ib_mesh_2d");

        // Seems to be set to return 0 by design
        TS_ASSERT_EQUALS(reader.GetNumFaces(), 0);

        auto element = reader.GetNextElementData();
        TS_ASSERT_EQUALS(element.NodeIndices.size(), 0);
        TS_ASSERT_EQUALS(element.AttributeValue, 0);

        auto face = reader.GetNextFaceData();
        TS_ASSERT_EQUALS(face.NodeIndices.size(), 0);
        TS_ASSERT_EQUALS(face.AttributeValue, 0);
    }
};

#endif /*TESTIMMERSEDBOUNDARYMESHREADER_HPP_*/