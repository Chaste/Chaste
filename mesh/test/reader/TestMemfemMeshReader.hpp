/*

Copyright (C) University of Oxford, 2005-2012

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


#ifndef _TESTMEMFEMMESHREADER_HPP_
#define _TESTMEMFEMMESHREADER_HPP_

#include <cxxtest/TestSuite.h>
#include "MemfemMeshReader.hpp"
#include "GenericMeshReader.hpp"

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

    void TestGenericReader() throw (Exception)
    {
        std::auto_ptr<AbstractMeshReader<3, 3> > p_mesh_reader = GenericMeshReader<3,3>("mesh/test/data/Memfem_slab");

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
        TS_ASSERT_THROWS_THIS(mesh_reader2.GetNode(0), "Random access is only implemented in mesh readers for binary mesh files.");
        TS_ASSERT_THROWS_THIS(mesh_reader2.GetElementData(0), "Random access is only implemented in mesh readers for binary mesh files.");
        TS_ASSERT_THROWS_THIS(mesh_reader2.GetFaceData(0), "Random access is only implemented in mesh readers for binary mesh files.");
        TS_ASSERT_THROWS_THIS(mesh_reader2.GetEdgeData(0), "Random access is only implemented in mesh readers for binary mesh files.");
        TS_ASSERT_THROWS_THIS(mesh_reader2.GetContainingElementIndices(0), "Ncl files are only implemented in mesh readers for binary mesh files.");

    }

};

#endif //_TESTMEMFEMMESHREADER_HPP_
