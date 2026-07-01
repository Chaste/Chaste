/*

Copyright (c) 2005-2026, University of Oxford.
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

#ifndef TESTSEMSINGLEELEMENTMESHGENERATOR_HPP_
#define TESTSEMSINGLEELEMENTMESHGENERATOR_HPP_

#include <cxxtest/TestSuite.h>

#include "SemSingleElementMeshGenerator.hpp"
#include "SemMesh.hpp"

// This test is never run in parallel
#include "FakePetscSetup.hpp"


class TestSemSingleElementMeshGenerator : public CxxTest::TestSuite
{
public:
    void Test1D()
    {
        SemSingleElementMeshGenerator<1> generator({5}, 2.0);
        auto p_mesh = generator.GetMesh();

        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 5u);

        TS_ASSERT_DELTA(p_mesh->GetNode(0u)->rGetLocation()[0], 0.0, 1e-6);
        TS_ASSERT_DELTA(p_mesh->GetNode(1u)->rGetLocation()[0], 0.4, 1e-6);
        TS_ASSERT_DELTA(p_mesh->GetNode(2u)->rGetLocation()[0], 0.8, 1e-6);
        TS_ASSERT_DELTA(p_mesh->GetNode(3u)->rGetLocation()[0], 1.2, 1e-6);
        TS_ASSERT_DELTA(p_mesh->GetNode(4u)->rGetLocation()[0], 1.6, 1e-6);
    }

    void Test2D()
    {
        SemSingleElementMeshGenerator<2> generator({4, 5}, 3.0);
        auto p_mesh = generator.GetMesh();

        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 20u);

        TS_ASSERT_DELTA(p_mesh->GetNode(0u)->rGetLocation()[0], 0.0, 1e-6);
        TS_ASSERT_DELTA(p_mesh->GetNode(0u)->rGetLocation()[1], 0.0, 1e-6);

        TS_ASSERT_DELTA(p_mesh->GetNode(1u)->rGetLocation()[0], 0.75, 1e-6);
        TS_ASSERT_DELTA(p_mesh->GetNode(1u)->rGetLocation()[1], 0.0, 1e-6);

        TS_ASSERT_DELTA(p_mesh->GetNode(6u)->rGetLocation()[0], 1.5, 1e-6);
        TS_ASSERT_DELTA(p_mesh->GetNode(6u)->rGetLocation()[1], 0.75, 1e-6);

        TS_ASSERT_DELTA(p_mesh->GetNode(19u)->rGetLocation()[0], 2.25, 1e-6);
        TS_ASSERT_DELTA(p_mesh->GetNode(19u)->rGetLocation()[1], 3.0, 1e-6);
    }

    void Test3D()
    {
        SemSingleElementMeshGenerator<3> generator({4, 3, 2}, 5.0);
        auto p_mesh = generator.GetMesh();

        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 24u);

        TS_ASSERT_DELTA(p_mesh->GetNode(0u)->rGetLocation()[0], 0.0, 1e-6);
        TS_ASSERT_DELTA(p_mesh->GetNode(0u)->rGetLocation()[1], 0.0, 1e-6);
        TS_ASSERT_DELTA(p_mesh->GetNode(0u)->rGetLocation()[2], 0.0, 1e-6);

        TS_ASSERT_DELTA(p_mesh->GetNode(17u)->rGetLocation()[0], 1.25, 1e-6);
        TS_ASSERT_DELTA(p_mesh->GetNode(17u)->rGetLocation()[1], 1.25, 1e-6);
        TS_ASSERT_DELTA(p_mesh->GetNode(17u)->rGetLocation()[2], 1.25, 1e-6);

        TS_ASSERT_DELTA(p_mesh->GetNode(23u)->rGetLocation()[0], 3.75, 1e-6);
        TS_ASSERT_DELTA(p_mesh->GetNode(23u)->rGetLocation()[1], 2.5, 1e-6);
        TS_ASSERT_DELTA(p_mesh->GetNode(23u)->rGetLocation()[2], 1.25, 1e-6);
    }

    void TestExceptions()
    {
        TS_ASSERT_THROWS_THIS(SemSingleElementMeshGenerator<1>({0}, 1.0),
            "SemSingleElementMeshGenerator: each entry of numNodes must be >= 1");

        TS_ASSERT_THROWS_THIS(SemSingleElementMeshGenerator<2>({2, 3}, -1.0),
            "SemSingleElementMeshGenerator: scaleFactor must be positive");
    }
};

#endif /*TESTSEMSINGLEELEMENTMESHGENERATOR_HPP_*/