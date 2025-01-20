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

#ifndef TESTIMMERSEDBOUNDARYELEMENT_HPP_
#define TESTIMMERSEDBOUNDARYELEMENT_HPP_

// Needed for test framework
#include <cxxtest/TestSuite.h>

// Includes from projects/ImmersedBoundary
#include "ImmersedBoundaryElement.hpp"

// This test is never run in parallel
#include "FakePetscSetup.hpp"

class TestImmersedBoundaryElement : public CxxTest::TestSuite
{
public:

    // Note: the corner nodes and average node spacing of an element are set in the ImmersedBoundarMesh method DivideElement()

    void TestFluidSourceMethods()
    {
        // Make 4 nodes to assign to a square element
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, true, 1.0, 0.0));
        nodes.push_back(new Node<2>(2, true, 1.0, 1.0));
        nodes.push_back(new Node<2>(3, true, 0.0, 1.0));

        // Make a square element out of these nodes
        ImmersedBoundaryElement<2,2> element(0, nodes);

        TS_ASSERT_EQUALS(element.GetNumNodes(), 4u);

        // Test SetFluidSource() and GetFluidSource() work correctly
        TS_ASSERT(element.GetFluidSource() == NULL);

        auto source = std::make_shared<FluidSource<2>>(0, 0.5, 0.5);
        source->SetStrength(57.0);
        element.SetFluidSource(source);

        TS_ASSERT_EQUALS(element.GetFluidSource(), source);
        TS_ASSERT_EQUALS(element.GetFluidSource()->GetIndex(), 0u);
        TS_ASSERT_DELTA(element.GetFluidSource()->GetStrength(), 57.0, 1e-6);
        TS_ASSERT_DELTA(element.GetAverageNodeSpacing(), DOUBLE_UNSET, 1e-6);

        element.SetAverageNodeSpacing(0.123);
        TS_ASSERT_DELTA(element.GetAverageNodeSpacing(), 0.123, 1e-6);

        // Clean up
        for (auto& p_node : nodes)
        {
            delete p_node;
        }
    }

    void Test1DMethods()
    {
        std::vector<Node<1>*> nodes;
        nodes.push_back(new Node<1>(0, true, 0.0));
        nodes.push_back(new Node<1>(1, true, 1.0));

        // Fluid source
        ImmersedBoundaryElement<1, 1> element(0, nodes);
        auto source = std::make_shared<FluidSource<1>>(0, 0.5, 0.5);
        source->SetStrength(57.0);

        TS_ASSERT_THROWS_ANYTHING(element.SetFluidSource(source));
        TS_ASSERT_THROWS_ANYTHING(element.GetFluidSource());

        // Boundary methods
        element.SetIsBoundaryElement(false);
        TS_ASSERT_EQUALS(element.IsElementOnBoundary(), false);
        element.SetIsBoundaryElement(true);
        TS_ASSERT_EQUALS(element.IsElementOnBoundary(), true);
        TS_ASSERT_DIFFERS(&(element.rGetCornerNodes()), nullptr);

        for (auto& p_node : nodes)
        {
            delete p_node;
        }
    }
};

#endif /*TESTIMMERSEDBOUNDARYELEMENT_HPP_*/