/*

Copyright (c) 2005-2012, University of Oxford.
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
            TS_ASSERT_EQUALS(nodes_owned.size(), 289u);
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
        EXIT_IF_SEQUENTIAL //Doesn't make sense to try and partition in sequential

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
