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

#ifndef TESTGMSHMESHREADER_
#define TESTGMSHMESHREADER_

#include <cxxtest/TestSuite.h>
#include <fstream>

#include "GmshMeshReader.hpp"
#include "TetrahedralMesh.hpp"
#include "GenericMeshReader.hpp"
#include "TrianglesMeshWriter.hpp"
#include "UblasVectorInclude.hpp"
#include "QuadraticMesh.hpp"

#include "PetscSetupAndFinalize.hpp"

typedef GmshMeshReader<2,2> READER_2D;
typedef GmshMeshReader<3,3> READER_3D;
typedef GmshMeshReader<2,3> READER_2D_3D;

class TestGmshMeshReader : public CxxTest::TestSuite
{
public:
    void TestFilesOpen(void)
    {
        TS_ASSERT_THROWS_NOTHING(READER_2D("mesh/test/data/square_4_elements_gmsh.msh"));
        TS_ASSERT_THROWS_THIS(READER_2D("mesh/test/data/no_file.msh"),
                              "Could not open data file: mesh/test/data/no_file.msh");
    }

    void TestCorrectVersion(void)
    {
       TS_ASSERT_THROWS_NOTHING(READER_2D("mesh/test/data/square_4_elements_gmsh.msh"));
       TS_ASSERT_THROWS_THIS(READER_2D("mesh/test/data/square_4_elements_bad_version.msh"),
                             "Only .msh version 2.2 files are supported.");
    }

    void TestErrorIfMeshContainsNodesInWeirdElements(void)
    {
        TS_ASSERT_THROWS_THIS(READER_3D reader_3d("mesh/test/data/simple_cube_gmsh_bad.msh"),
                              "Unrecognised element types present in the .msh file: check mesh generation settings in gmsh.");
    }

    void TestReadHeaders(void)
    {
        //Linear meshes
        READER_2D reader("mesh/test/data/square_4_elements_gmsh.msh");
        TS_ASSERT_EQUALS(reader.GetNumNodes(), 5u);
        TS_ASSERT_EQUALS(reader.GetNumElements(), 4u);
        TS_ASSERT_EQUALS(reader.GetNumFaces(), 4u);

        READER_3D reader_3d("mesh/test/data/simple_cube_gmsh.msh");
        TS_ASSERT_EQUALS(reader_3d.GetNumNodes(), 14u);
        TS_ASSERT_EQUALS(reader_3d.GetNumElements(), 24u);
        TS_ASSERT_EQUALS(reader_3d.GetNumFaces(), 24u);

        //Quad meshes
        READER_2D quad_reader("mesh/test/data/quad_square_4_elements_gmsh.msh");
        TS_ASSERT_EQUALS(quad_reader.GetNumNodes(), 13u);
        TS_ASSERT_EQUALS(quad_reader.GetNumElements(), 4u);
        TS_ASSERT_EQUALS(quad_reader.GetNumFaces(), 4u);

        READER_3D quad_reader_3d("mesh/test/data/quad_cube_gmsh.msh");
        TS_ASSERT_EQUALS(quad_reader_3d.GetNumNodes(), 63u);
        TS_ASSERT_EQUALS(quad_reader_3d.GetNumElements(), 24u);
        TS_ASSERT_EQUALS(quad_reader_3d.GetNumFaces(), 24u);
    }

    void TestGetNextMethods(void)
    {

        READER_2D reader("mesh/test/data/square_4_elements_gmsh.msh");

        //2D nodes
        std::vector<double> expected_coords(2);
        expected_coords[0] = 0.0; expected_coords[1] = 0.0;
        TS_ASSERT_EQUALS(reader.GetNextNode(), expected_coords);
        expected_coords[0] = 1.0; expected_coords[1] = 0.0;
        TS_ASSERT_EQUALS(reader.GetNextNode(), expected_coords);
        expected_coords[0] = 1.0; expected_coords[1] = 1.0;
        TS_ASSERT_EQUALS(reader.GetNextNode(), expected_coords);
        expected_coords[0] = 0.0; expected_coords[1] = 1.0;
        TS_ASSERT_EQUALS(reader.GetNextNode(), expected_coords);
        expected_coords[0] = 0.5; expected_coords[1] = 0.5;
        TS_ASSERT_EQUALS(reader.GetNextNode(), expected_coords);

        //2D elements
        std::vector<unsigned> expected_node_indices(3);

        expected_node_indices[0] = 0; expected_node_indices[1] = 1; expected_node_indices[2] = 4;
        ElementData data = reader.GetNextElementData();
        TS_ASSERT_EQUALS(data.NodeIndices, expected_node_indices);
        TS_ASSERT_EQUALS(data.AttributeValue, 1u);

        expected_node_indices[0] = 1; expected_node_indices[1] = 2; expected_node_indices[2] = 4;
        data = reader.GetNextElementData();
        TS_ASSERT_EQUALS(data.NodeIndices, expected_node_indices);
        TS_ASSERT_EQUALS(data.AttributeValue, 1u);

        //2D faces
        expected_node_indices.resize(2);
        expected_node_indices[0] = 0; expected_node_indices[1] = 1;
        data = reader.GetNextFaceData();
        TS_ASSERT_EQUALS(data.NodeIndices, expected_node_indices);
        TS_ASSERT_EQUALS(data.AttributeValue, 2u);

        expected_node_indices[0] = 1; expected_node_indices[1] = 2;
        data = reader.GetNextFaceData();
        TS_ASSERT_EQUALS(data.NodeIndices, expected_node_indices);
        TS_ASSERT_EQUALS(data.AttributeValue, 2u);


        READER_3D reader_3d("mesh/test/data/simple_cube_gmsh.msh");

        //3D nodes
        expected_coords.resize(3);
        expected_coords[0] = 0.0; expected_coords[1] = 0.0; expected_coords[2] = 0.0;
        TS_ASSERT_EQUALS(reader_3d.GetNextNode(), expected_coords);
        expected_coords[0] = 1.0; expected_coords[1] = 0.0; expected_coords[2] = 0.0;
        TS_ASSERT_EQUALS(reader_3d.GetNextNode(), expected_coords);
        expected_coords[0] = 1.0; expected_coords[1] = 1.0; expected_coords[2] = 0.0;
        TS_ASSERT_EQUALS(reader_3d.GetNextNode(), expected_coords);
        expected_coords[0] = 0.0; expected_coords[1] = 1.0; expected_coords[2] = 0.0;
        TS_ASSERT_EQUALS(reader_3d.GetNextNode(), expected_coords);
        expected_coords[0] = 0.0; expected_coords[1] = 0.0; expected_coords[2] = 1.0;
        TS_ASSERT_EQUALS(reader_3d.GetNextNode(), expected_coords);

        //3D elements
        expected_node_indices.resize(4);

        expected_node_indices[0] = 13; expected_node_indices[1] = 10; expected_node_indices[2] = 7; expected_node_indices[3] = 3;
        data = reader_3d.GetNextElementData();
        TS_ASSERT_EQUALS(data.NodeIndices, expected_node_indices);
        TS_ASSERT_EQUALS(data.AttributeValue, 1u);

        expected_node_indices[0] = 8; expected_node_indices[1] = 0; expected_node_indices[2] = 10; expected_node_indices[3] = 3;
        data = reader_3d.GetNextElementData();
        TS_ASSERT_EQUALS(data.NodeIndices, expected_node_indices);
        TS_ASSERT_EQUALS(data.AttributeValue, 1u);

        //3D faces
        expected_node_indices.resize(3);

        expected_node_indices[0] = 0; expected_node_indices[1] = 1; expected_node_indices[2] = 8;
        data = reader_3d.GetNextFaceData();
        TS_ASSERT_EQUALS(data.NodeIndices, expected_node_indices);
        TS_ASSERT_EQUALS(data.AttributeValue, 2u);

        expected_node_indices[0] = 0; expected_node_indices[1] = 8; expected_node_indices[2] = 3;
        data = reader_3d.GetNextFaceData();
        TS_ASSERT_EQUALS(data.NodeIndices, expected_node_indices);
        TS_ASSERT_EQUALS(data.AttributeValue, 2u);
    }

    void TestRead2dMeshes(void)
    {
        READER_2D reader("mesh/test/data/square_4_elements_gmsh.msh");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(reader);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 5u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 4u);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), 4u);
        // Check area and middle node location
        TS_ASSERT_DELTA(mesh.GetVolume(), 1.0, 1e-15);
        TS_ASSERT_DELTA(mesh.GetNode(4)->rGetLocation()[0], 0.5, 1e-15);
        TS_ASSERT_DELTA(mesh.GetNode(4)->rGetLocation()[1], 0.5, 1e-15);
    }

    void TestRead2d3dMeshes(void)
    {
        READER_2D_3D reader("mesh/test/data/square_4_elements_nonplanar_gmsh.msh");  // As previous test, but the central node is elevated
        TetrahedralMesh<2,3> mesh;
        mesh.ConstructFromMeshReader(reader);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 5u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 4u);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), 4u);
        // Check area and middle node location
        TS_ASSERT_DELTA(mesh.GetVolume(), /*1.0*/ sqrt(2), 1e-15);
        TS_ASSERT_DELTA(mesh.GetNode(4)->rGetLocation()[0], 0.5, 1e-15);
        TS_ASSERT_DELTA(mesh.GetNode(4)->rGetLocation()[1], 0.5, 1e-15);
        TS_ASSERT_DELTA(mesh.GetNode(4)->rGetLocation()[2], 0.5, 1e-15);
    }

    void TestRead2dQuadraticMeshes(void)
    {
       READER_2D reader("mesh/test/data/quad_square_4_elements_gmsh.msh",2,2);
       QuadraticMesh<2> mesh;
       mesh.ConstructFromMeshReader(reader);
       TS_ASSERT_EQUALS(mesh.GetNumNodes(), 13u);
       TS_ASSERT_EQUALS(mesh.GetNumElements(), 4u);
       TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), 4u);
    }

    void TestRead3dMeshes(void)
    {
        READER_3D reader("mesh/test/data/simple_cube_gmsh.msh");
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(reader);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 14u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 24u);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), 24u);
    }

    void TestRead3dQuadraticMeshes(void)
    {
       READER_3D reader("mesh/test/data/quad_cube_gmsh.msh",2,2);
       QuadraticMesh<3> mesh;
       mesh.ConstructFromMeshReader(reader);
       TS_ASSERT_EQUALS(mesh.GetNumNodes(), 63u);
       TS_ASSERT_EQUALS(mesh.GetNumElements(), 24u);
       TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), 24u);
    }
};

#endif /*TESTGMSHMESHREADER_*/
