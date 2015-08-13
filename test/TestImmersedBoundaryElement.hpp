/*

Copyright (c) 2005-2014, University of Oxford.
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

// Needed for test framework
#include <cxxtest/TestSuite.h>
#include "AbstractCellBasedTestSuite.hpp"

// Includes from projects/ImmersedBoundary
#include "ImmersedBoundaryElement.hpp"

// This test is never run in parallel
#include "FakePetscSetup.hpp"

class TestImmersedBoundaryElement : public CxxTest::TestSuite
{
public:

    void TestCreate2DImmersedBoundaryElement() throw(Exception)
    {
        // Make 4 nodes to assign to a square element
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, true, 1.0, 0.0));
        nodes.push_back(new Node<2>(2, true, 1.0, 1.0));
        nodes.push_back(new Node<2>(3, true, 0.0, 1.0));

        // Make a square element out of these nodes
        ImmersedBoundaryElement<2,2> element(0, nodes);

        // Check default parameters
        TS_ASSERT_DELTA(element.GetMembraneSpringConstant(), 1000.0, 1e-6);
        TS_ASSERT_DELTA(element.GetMembraneRestLength(), 0.05, 1e-6);

        TS_ASSERT_DELTA(element.GetCellCellSpringConstant(), 50.0, 1e-6);
        TS_ASSERT_DELTA(element.GetCellCellRestLength(), 0.01, 1e-6);
    }

    void TestSetAndGetMethods() throw(Exception)
    {
        // Make 4 nodes to assign to a square element
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, true, 1.0, 0.0));
        nodes.push_back(new Node<2>(2, true, 1.0, 1.0));
        nodes.push_back(new Node<2>(3, true, 0.0, 1.0));

        // Make a square element out of these nodes
        ImmersedBoundaryElement<2,2> element(0, nodes);

        // Set the settable parameters
        element.SetMembraneSpringConstant(1.23);
        element.SetMembraneRestLength(2.34);

        element.SetCellCellSpringConstant(3.45);
        element.SetCellCellRestLength(4.56);

        // Check we get the correct values
        TS_ASSERT_DELTA(element.GetMembraneSpringConstant(), 1.23, 1e-6);
        TS_ASSERT_DELTA(element.GetMembraneRestLength(), 2.34, 1e-6);

        TS_ASSERT_DELTA(element.GetCellCellSpringConstant(), 3.45, 1e-6);
        TS_ASSERT_DELTA(element.GetCellCellRestLength(), 4.56, 1e-6);

        TS_ASSERT_THROWS_ANYTHING(element.SetMembraneSpringConstant(-1.23));
    }

    void TestSetAndGetExceptions() throw(Exception)
    {
        // Make 4 nodes to assign to a square element
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, true, 1.0, 0.0));
        nodes.push_back(new Node<2>(2, true, 1.0, 1.0));
        nodes.push_back(new Node<2>(3, true, 0.0, 1.0));

        // Make a square element out of these nodes
        ImmersedBoundaryElement<2,2> element(0, nodes);

        // Test exceptions when setting bad parameter values
        TS_ASSERT_THROWS_THIS(element.SetMembraneSpringConstant(-1.23), "This parameter must be non-negative");
        TS_ASSERT_THROWS_THIS(element.SetMembraneRestLength(-2.34), "This parameter must be non-negative");

        TS_ASSERT_THROWS_THIS(element.SetCellCellSpringConstant(-3.45), "This parameter must be non-negative");
        TS_ASSERT_THROWS_THIS(element.SetCellCellRestLength(-4.56), "This parameter must be non-negative");
    }
};