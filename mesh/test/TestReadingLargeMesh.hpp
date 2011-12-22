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


#ifndef TESTREADINGLARGEMESH_HPP_
#define TESTREADINGLARGEMESH_HPP_

#include <cxxtest/TestSuite.h>
#include "TetrahedralMesh.hpp"
#include "TrianglesMeshReader.hpp"

#include <vector>

class TestReadingLargeTetrahedralMesh : public CxxTest::TestSuite
{
public:

    /**
     * This test is mainly here for performance testing, to check that loading a
     * (relatively) large mesh doesn't take too long.
     * It's a nightly test because it takes 10 hours to run under MemoryTesting!
     */
    void TestLoadingLargeMesh()
    {
        TrianglesMeshReader<3,3> meshReader("heart/test/data/heart");
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(meshReader);

        // Check we have the right number of nodes, elements and boundary elements
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 63885u);
        TS_ASSERT_EQUALS(mesh.GetNumElements(), 322267u);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), 41812u);
    }
};

#endif /*TESTREADINGLARGEMESH_HPP_*/
