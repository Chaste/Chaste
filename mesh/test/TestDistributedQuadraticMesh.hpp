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

#ifndef TESTDISTRIBUTEDQUADRATICMESH_HPP_
#define TESTDISTRIBUTEDQUADRATICMESH_HPP_

#include <cxxtest/TestSuite.h>

#include <set>
#include <vector>

#include "NodePartitioner.hpp"
#include "QuadraticMesh.hpp"
#include "TrianglesMeshReader.hpp"
#include "PetscTools.hpp"
#include "PetscSetupAndFinalize.hpp"


class TestDistributedQuadraticMesh : public CxxTest::TestSuite
{
public:
    void TestDumbMeshPartitioning() throw (Exception)
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements_quadratic",2,1, false);
        QuadraticMesh<2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        std::set<unsigned> nodes_owned;

        NodePartitioner<2, 2>::DumbPartitioning(mesh, nodes_owned);

        if(PetscTools::GetNumProcs() == 1)
        {
            TS_ASSERT_EQUALS(nodes_owned.size(), 291u);
        }
        else if(PetscTools::GetNumProcs() == 2)
        {
            if(PetscTools::GetMyRank() == 0 )
            {
                TS_ASSERT_EQUALS(nodes_owned.size(), 145u);
            }
            else //PetscTools::GetMyRank() == 1
            {
                TS_ASSERT_EQUALS(nodes_owned.size(), 144u);
            }
        }
        else if(PetscTools::GetNumProcs() == 3)
        {
            if(PetscTools::GetMyRank() == 0 )
            {
                TS_ASSERT_EQUALS(nodes_owned.size(), 97u);
            }
            else if(PetscTools::GetMyRank() == 1 )
            {
                TS_ASSERT_EQUALS(nodes_owned.size(), 96u);
            }
            else //PetscTools::GetMyRank() == 2
            {
                TS_ASSERT_EQUALS(nodes_owned.size(), 96u);
            }
        }
    }

    void TestMetisMeshPartitioning() throw (Exception)
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements_quadratic",2,1, false);
        std::vector<unsigned> nodes_permutation;
        std::set<unsigned> nodes_owned;
        std::vector<unsigned> processor_offset;

        typedef NodePartitioner<2, 2> Partioner2D;

        TS_ASSERT_THROWS_THIS(Partioner2D::MetisLibraryPartitioning(mesh_reader, nodes_permutation, nodes_owned, processor_offset),
                              "Metis cannot partition a quadratic mesh.");
    }
};

#endif // TESTDISTRIBUTEDQUADRATICMESH_HPP_
