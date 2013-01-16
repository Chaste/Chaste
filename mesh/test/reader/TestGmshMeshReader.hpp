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

#ifndef TESTGMSHMESHREADER_
#define TESTGMSHMESHREADER_

#include <cxxtest/TestSuite.h>
#include <fstream>

#include "GmshMeshReader.hpp"
#include "TetrahedralMesh.hpp"
#include "GenericMeshReader.hpp"
#include "TrianglesMeshWriter.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "UblasVectorInclude.hpp"

typedef GmshMeshReader<2,2> READER_2D;
typedef GmshMeshReader<3,3> READER_3D;

class TestGmshMeshReader : public CxxTest::TestSuite
{
public:
    void TestFilesOpen(void) throw(Exception)
    {
        TS_ASSERT_THROWS_NOTHING(READER_2D("mesh/test/data/square_4_elements_gmsh.msh"));
        TS_ASSERT_THROWS_ANYTHING(READER_2D("mesh/test/data/no_file.msh"));
    }

    void TestCorrectVersion(void) throw(Exception)
    {
       TS_ASSERT_THROWS_NOTHING(READER_2D("mesh/test/data/square_4_elements_gmsh.msh"));
       TS_ASSERT_THROWS_ANYTHING(READER_2D("mesh/test/data/square_4_elements_bad_version.msh"));
    }

    void TestReadHeaders(void) throw(Exception)
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


        ///\todo Check in when the above is complete
    }

};

#endif /*TESTGMSHMESHREADER_*/
