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

#ifndef TESTHONEYCOMBVERTEXMESHGENERATOR_HPP_
#define TESTHONEYCOMBVERTEXMESHGENERATOR_HPP_

#include <cxxtest/TestSuite.h>

#include "HoneycombVertexMeshGenerator.hpp"

class TestHoneycombVertexMeshGenerator : public CxxTest::TestSuite
{
public:

    void TestSimpleMesh() throw(Exception)
    {
        HoneycombVertexMeshGenerator generator(2, 2, false, 0.1, 0.1);

        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 16u);
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 4u);
        TS_ASSERT_DELTA(p_mesh->GetCellRearrangementThreshold(), 0.1, 1e-12);
        TS_ASSERT_DELTA(p_mesh->GetT2Threshold(), 0.1, 1e-12);
    }

    void TestBoundaryNodes() throw(Exception)
    {
        HoneycombVertexMeshGenerator generator(4, 4);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 48u);
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 16u);

        unsigned num_non_boundary_nodes = 0;
        for (unsigned node_index=0; node_index<16u; node_index++)
        {
            if (!p_mesh->GetNode(node_index)->IsBoundaryNode())
            {
                num_non_boundary_nodes++;
            }
        }
        TS_ASSERT_EQUALS(num_non_boundary_nodes, 4u);
    }

    void TestLargeMesh() throw(Exception)
    {
        HoneycombVertexMeshGenerator generator(100, 100);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 20400u);
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 10000u);
    }
};

#endif /*TESTHONEYCOMBVERTEXMESHGENERATOR_HPP_*/
