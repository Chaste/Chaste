/*

Copyright (c) 2005-2018, University of Oxford.
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

#ifndef TESTIMMERSEDBOUNDARYHONEYCOMBMESHGENERATOR_HPP_
#define TESTIMMERSEDBOUNDARYHONEYCOMBMESHGENERATOR_HPP_

// Needed for test framework
#include <cxxtest/TestSuite.h>

// Includes from projects/ImmersedBoundary
#include "ImmersedBoundaryHoneycombMeshGenerator.hpp"
#include "ImmersedBoundaryMesh.hpp"

// This test is never run in parallel
#include "FakePetscSetup.hpp"

class TestImmersedBoundaryHoneycombMeshGenerator : public CxxTest::TestSuite
{
public:

    ///\todo Improve testing
    void TestBoundaryElementsAreTaggedCorrectly()
    {
        /*
         * Elements are numbered from bottom-left, upwards, then from left to right.  E.g. for a 3x3:
         *
         *      5
         *   2     8
         *      4
         *   1     7
         *      3
         *   0     6
         */

        // 3x3 mesh
        {
            ImmersedBoundaryHoneycombMeshGenerator gen(3u, 3u, 5u, 0.05, 0.2);
            auto p_mesh = gen.GetMesh();

            // Only element with idx 4 should be non-boundary
            TS_ASSERT(!p_mesh->GetElement(4u)->IsElementOnBoundary());

            // All other elements are on the boundary
            TS_ASSERT(p_mesh->GetElement(0u)->IsElementOnBoundary());
            TS_ASSERT(p_mesh->GetElement(1u)->IsElementOnBoundary());
            TS_ASSERT(p_mesh->GetElement(2u)->IsElementOnBoundary());
            TS_ASSERT(p_mesh->GetElement(3u)->IsElementOnBoundary());
            TS_ASSERT(p_mesh->GetElement(5u)->IsElementOnBoundary());
            TS_ASSERT(p_mesh->GetElement(6u)->IsElementOnBoundary());
            TS_ASSERT(p_mesh->GetElement(7u)->IsElementOnBoundary());
            TS_ASSERT(p_mesh->GetElement(8u)->IsElementOnBoundary());
        }

        // 3x4 mesh
        {
            ImmersedBoundaryHoneycombMeshGenerator gen(3u, 4u, 5u, 0.05, 0.2);
            auto p_mesh = gen.GetMesh();

            // Only elements with ids 5 and 6 should be non-boundary
            TS_ASSERT(!p_mesh->GetElement(5u)->IsElementOnBoundary());
            TS_ASSERT(!p_mesh->GetElement(6u)->IsElementOnBoundary());

            // All other elements are on the boundary
            TS_ASSERT(p_mesh->GetElement(0u)->IsElementOnBoundary());
            TS_ASSERT(p_mesh->GetElement(1u)->IsElementOnBoundary());
            TS_ASSERT(p_mesh->GetElement(2u)->IsElementOnBoundary());
            TS_ASSERT(p_mesh->GetElement(3u)->IsElementOnBoundary());
            TS_ASSERT(p_mesh->GetElement(4u)->IsElementOnBoundary());
            TS_ASSERT(p_mesh->GetElement(7u)->IsElementOnBoundary());
            TS_ASSERT(p_mesh->GetElement(8u)->IsElementOnBoundary());
            TS_ASSERT(p_mesh->GetElement(9u)->IsElementOnBoundary());
            TS_ASSERT(p_mesh->GetElement(10u)->IsElementOnBoundary());
            TS_ASSERT(p_mesh->GetElement(11u)->IsElementOnBoundary());
        }

        // 5x3 mesh
        {
            ImmersedBoundaryHoneycombMeshGenerator gen(5u, 3u, 5u, 0.05, 0.2);
            auto p_mesh = gen.GetMesh();

            // Only elements with ids 4, 7, and 10 should be non-boundary
            TS_ASSERT(!p_mesh->GetElement(4u)->IsElementOnBoundary());
            TS_ASSERT(!p_mesh->GetElement(7u)->IsElementOnBoundary());
            TS_ASSERT(!p_mesh->GetElement(10u)->IsElementOnBoundary());

            // All other elements are on the boundary
            TS_ASSERT(p_mesh->GetElement(0u)->IsElementOnBoundary());
            TS_ASSERT(p_mesh->GetElement(1u)->IsElementOnBoundary());
            TS_ASSERT(p_mesh->GetElement(2u)->IsElementOnBoundary());
            TS_ASSERT(p_mesh->GetElement(3u)->IsElementOnBoundary());
            TS_ASSERT(p_mesh->GetElement(5u)->IsElementOnBoundary());
            TS_ASSERT(p_mesh->GetElement(6u)->IsElementOnBoundary());
            TS_ASSERT(p_mesh->GetElement(8u)->IsElementOnBoundary());
            TS_ASSERT(p_mesh->GetElement(9u)->IsElementOnBoundary());
            TS_ASSERT(p_mesh->GetElement(11u)->IsElementOnBoundary());
            TS_ASSERT(p_mesh->GetElement(12u)->IsElementOnBoundary());
            TS_ASSERT(p_mesh->GetElement(13u)->IsElementOnBoundary());
            TS_ASSERT(p_mesh->GetElement(14u)->IsElementOnBoundary());
        }
    }
};

#endif /*TESTIMMERSEDBOUNDARYHONEYCOMBMESHGENERATOR_HPP_*/