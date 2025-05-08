/*

Copyright (c) 2005-2025, University of Oxford.
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
            TS_ASSERT_EQUALS(p_mesh->GetElement(4)->IsElementOnBoundary(), false);

            // All other elements are on the boundary
            TS_ASSERT_EQUALS(p_mesh->GetElement(0)->IsElementOnBoundary(), true);
            TS_ASSERT_EQUALS(p_mesh->GetElement(1)->IsElementOnBoundary(), true);
            TS_ASSERT_EQUALS(p_mesh->GetElement(2)->IsElementOnBoundary(), true);
            TS_ASSERT_EQUALS(p_mesh->GetElement(3)->IsElementOnBoundary(), true);
            TS_ASSERT_EQUALS(p_mesh->GetElement(5)->IsElementOnBoundary(), true);
            TS_ASSERT_EQUALS(p_mesh->GetElement(6)->IsElementOnBoundary(), true);
            TS_ASSERT_EQUALS(p_mesh->GetElement(7)->IsElementOnBoundary(), true);
            TS_ASSERT_EQUALS(p_mesh->GetElement(8)->IsElementOnBoundary(), true);
        }

        // 3x4 mesh
        {
            ImmersedBoundaryHoneycombMeshGenerator gen(3u, 4u, 5u, 0.05, 0.2);
            auto p_mesh = gen.GetMesh();

            // Only elements with ids 5 and 6 should be non-boundary
            TS_ASSERT_EQUALS(p_mesh->GetElement(5)->IsElementOnBoundary(), false);
            TS_ASSERT_EQUALS(p_mesh->GetElement(6)->IsElementOnBoundary(), false);

            // All other elements are on the boundary
            TS_ASSERT_EQUALS(p_mesh->GetElement(0)->IsElementOnBoundary(), true);
            TS_ASSERT_EQUALS(p_mesh->GetElement(1)->IsElementOnBoundary(), true);
            TS_ASSERT_EQUALS(p_mesh->GetElement(2)->IsElementOnBoundary(), true);
            TS_ASSERT_EQUALS(p_mesh->GetElement(3)->IsElementOnBoundary(), true);
            TS_ASSERT_EQUALS(p_mesh->GetElement(4)->IsElementOnBoundary(), true);
            TS_ASSERT_EQUALS(p_mesh->GetElement(7)->IsElementOnBoundary(), true);
            TS_ASSERT_EQUALS(p_mesh->GetElement(8)->IsElementOnBoundary(), true);
            TS_ASSERT_EQUALS(p_mesh->GetElement(9)->IsElementOnBoundary(), true);
            TS_ASSERT_EQUALS(p_mesh->GetElement(10)->IsElementOnBoundary(), true);
            TS_ASSERT_EQUALS(p_mesh->GetElement(11)->IsElementOnBoundary(), true);
        }

        // 5x3 mesh
        {
            ImmersedBoundaryHoneycombMeshGenerator gen(5u, 3u, 5u, 0.05, 0.2);
            auto p_mesh = gen.GetMesh();

            // Only elements with ids 4, 7, and 10 should be non-boundary
            TS_ASSERT_EQUALS(p_mesh->GetElement(4)->IsElementOnBoundary(), false);
            TS_ASSERT_EQUALS(p_mesh->GetElement(7)->IsElementOnBoundary(), false);
            TS_ASSERT_EQUALS(p_mesh->GetElement(10)->IsElementOnBoundary(), false);

            // All other elements are on the boundary
            TS_ASSERT_EQUALS(p_mesh->GetElement(0)->IsElementOnBoundary(), true);
            TS_ASSERT_EQUALS(p_mesh->GetElement(1)->IsElementOnBoundary(), true);
            TS_ASSERT_EQUALS(p_mesh->GetElement(2)->IsElementOnBoundary(), true);
            TS_ASSERT_EQUALS(p_mesh->GetElement(3)->IsElementOnBoundary(), true);
            TS_ASSERT_EQUALS(p_mesh->GetElement(5)->IsElementOnBoundary(), true);
            TS_ASSERT_EQUALS(p_mesh->GetElement(6)->IsElementOnBoundary(), true);
            TS_ASSERT_EQUALS(p_mesh->GetElement(8)->IsElementOnBoundary(), true);
            TS_ASSERT_EQUALS(p_mesh->GetElement(9)->IsElementOnBoundary(), true);
            TS_ASSERT_EQUALS(p_mesh->GetElement(11)->IsElementOnBoundary(), true);
            TS_ASSERT_EQUALS(p_mesh->GetElement(12)->IsElementOnBoundary(), true);
            TS_ASSERT_EQUALS(p_mesh->GetElement(13)->IsElementOnBoundary(), true);
            TS_ASSERT_EQUALS(p_mesh->GetElement(14)->IsElementOnBoundary(), true);
        }
    }

    void TestCorrectNumberNodesGenerated() {
        // 3 nodes per edge
        {
            ImmersedBoundaryHoneycombMeshGenerator gen(3u, 3u, 3u, 0.05, 0.2);
            auto p_mesh = gen.GetMesh();

            for (unsigned int i = 0; i < 9; i++) {
                TS_ASSERT_EQUALS(p_mesh->GetElement(i)->GetNumNodes(), 18);
            }
        }
        // 5 nodes per edge
        {
            ImmersedBoundaryHoneycombMeshGenerator gen(3u, 3u, 5u, 0.05, 0.2);
            auto p_mesh = gen.GetMesh();

            for (unsigned int i = 0; i < 9; i++) {
                TS_ASSERT_EQUALS(p_mesh->GetElement(i)->GetNumNodes(), 30);
            }
        }
    }

    void TestNodePositioning() {
        // 3 nodes per edge
        {
            ImmersedBoundaryHoneycombMeshGenerator gen(3u, 3u, 3u, 0.05, 0.2);
            auto p_mesh = gen.GetMesh();

            // Left hand row
            auto el0 = p_mesh->GetElement(0);
            auto el1 = p_mesh->GetElement(1);
            auto el2 = p_mesh->GetElement(2);

            // Nodes within elements 0, 1 & 1, 2 should only differ in y coordinate
            for (unsigned int i = 0; i < 18; i++) {
                TS_ASSERT_EQUALS(el0->GetNodeLocation(i)[0], el1->GetNodeLocation(i)[0]);
                TS_ASSERT_EQUALS(el0->GetNodeLocation(i)[0], el2->GetNodeLocation(i)[0]);
            }

            // y coordinates should be equally spaced for el0 to el1 & el1 to el2
            for (unsigned int i = 0; i < 18; i++) {
                auto diff1 = el0->GetNodeLocation(i)[1] - el1->GetNodeLocation(i)[1];
                auto diff2 = el1->GetNodeLocation(i)[1] - el2->GetNodeLocation(i)[1];
                TS_ASSERT_DELTA(diff1, diff2, 0.001);
            }

            // Nodes within elements 0, 6 should only differ in x coordinate
            auto el6 = p_mesh->GetElement(6);
            for (unsigned int i = 0; i < 18; i++) {
                TS_ASSERT_EQUALS(el0->GetNodeLocation(i)[1], el6->GetNodeLocation(i)[1]);
            }
        }

    }
};

#endif /*TESTIMMERSEDBOUNDARYHONEYCOMBMESHGENERATOR_HPP_*/