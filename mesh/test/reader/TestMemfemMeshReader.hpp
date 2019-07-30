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


#ifndef _TESTMEMFEMMESHREADER_HPP_
#define _TESTMEMFEMMESHREADER_HPP_

#include <cxxtest/TestSuite.h>
#include "MemfemMeshReader.hpp"
#include "GenericMeshReader.hpp"

#include "PetscSetupAndFinalize.hpp"

typedef MemfemMeshReader<3,3> READER_3D;
typedef MemfemMeshReader<2,2> READER_2D; // For exception coverage

class TestMemfemMeshReaders : public CxxTest::TestSuite
{
public:

    /**
     * Check that input files are opened correctly.
     *
     */
    void TestFilesOpen()
    {
        MemfemMeshReader<3,3>* pMeshReader;
        pMeshReader = new READER_3D("mesh/test/data/Memfem_slab");

        TS_ASSERT_EQUALS(pMeshReader->GetNumNodes(), 381u);
        TS_ASSERT_EQUALS(pMeshReader->GetNumElements(), 1030u);
        TS_ASSERT_EQUALS(pMeshReader->GetNumFaces(), 758u);
        TS_ASSERT_EQUALS(pMeshReader->GetNumElementAttributes(), 0u);

        std::vector<unsigned> next_face;

        next_face = pMeshReader->GetNextFaceData().NodeIndices;

        TS_ASSERT_EQUALS(next_face[0], 338u);
        TS_ASSERT_EQUALS(next_face[1], 23u);
        TS_ASSERT_EQUALS(next_face[2], 374u);

        TS_ASSERT_EQUALS(pMeshReader->GetMaxNodeIndex(), pMeshReader->GetNumNodes() - 1);

        TS_ASSERT_EQUALS(pMeshReader->GetMinNodeIndex(), 0u);

        // Coverage
        TS_ASSERT(!pMeshReader->HasNclFile());

        delete pMeshReader;
    }

    void TestGenericReader()
    {
        std::shared_ptr<AbstractMeshReader<3, 3> > p_mesh_reader = GenericMeshReader<3,3>("mesh/test/data/Memfem_slab");

        TS_ASSERT_EQUALS(p_mesh_reader->GetNumNodes(), 381u);
        TS_ASSERT_EQUALS(p_mesh_reader->GetNumElements(), 1030u);
        TS_ASSERT_EQUALS(p_mesh_reader->GetNumFaces(), 758u);
        TS_ASSERT_EQUALS(p_mesh_reader->GetNumElementAttributes(), 0u);
        TS_ASSERT_EQUALS(p_mesh_reader->GetNumFaceAttributes(), 0u);

        // The file does not exist
        TS_ASSERT_THROWS_CONTAINS((GenericMeshReader<3,3>("no_file")),
                                  "Could not open appropriate mesh files for no_file");
    }

    void TestExceptions()
    {
        // The file does not exist
        TS_ASSERT_THROWS_THIS( READER_3D mesh_reader("no_file"), "Could not open data file no_file.pts");

        // We are in the wrong dimension
        TS_ASSERT_THROWS_THIS( READER_2D reader("mesh/test/data/Memfem_slab"),
                "You have asked to read non-3D data. All Memfem data is in 3D.");

        MemfemMeshReader<3,3> mesh_reader2("mesh/test/data/Memfem_slab");
        TS_ASSERT_EQUALS(mesh_reader2.HasNodePermutation(), false);
        TS_ASSERT_THROWS_THIS(mesh_reader2.rGetNodePermutation(), "Node permutations aren't supported by this reader");
        TS_ASSERT_THROWS_THIS(mesh_reader2.GetNode(0), "Random access is only implemented in mesh readers for binary mesh files.");
        TS_ASSERT_THROWS_THIS(mesh_reader2.GetElementData(0), "Random access is only implemented in mesh readers for binary mesh files.");
        TS_ASSERT_THROWS_THIS(mesh_reader2.GetFaceData(0), "Random access is only implemented in mesh readers for binary mesh files.");
        TS_ASSERT_THROWS_THIS(mesh_reader2.GetEdgeData(0), "Random access is only implemented in mesh readers for binary mesh files.");
        TS_ASSERT_THROWS_THIS(mesh_reader2.GetContainingElementIndices(0), "Ncl files are only implemented in mesh readers for binary mesh files.");
    }
};

#endif //_TESTMEMFEMMESHREADER_HPP_
