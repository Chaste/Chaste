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

#ifndef TESTHONEYCOMBVERTEXMESHGENERATOR_HPP_
#define TESTHONEYCOMBVERTEXMESHGENERATOR_HPP_

#include <cxxtest/TestSuite.h>
#include "HoneycombVertexMeshGenerator.hpp"
#include "PetscSetupAndFinalize.hpp"

class TestHoneycombVertexMeshGenerator : public CxxTest::TestSuite
{
public:

    void TestSimpleMesh()
    {
        HoneycombVertexMeshGenerator generator(2, 2, false, 0.1, 0.1);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 16u);
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 4u);
        TS_ASSERT_DELTA(p_mesh->GetCellRearrangementThreshold(), 0.1, 1e-12);
        TS_ASSERT_DELTA(p_mesh->GetT2Threshold(), 0.1, 1e-12);
    }

    void TestBoundaryNodes()
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

    void TestLargeMesh()
    {
        HoneycombVertexMeshGenerator generator(100, 100);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 20400u);
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 10000u);
    }

    void TestElementArea()
    {
        HoneycombVertexMeshGenerator generator(6, 6, false, 0.01, 0.001, 2.456);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 96u);
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 36u);

        for (unsigned elem_index=0; elem_index<p_mesh->GetNumElements(); elem_index++)
        {
            TS_ASSERT_DELTA(p_mesh->GetVolumeOfElement(elem_index), 2.456, 1e-3);
        }
    }
};

#endif /*TESTHONEYCOMBVERTEXMESHGENERATOR_HPP_*/
