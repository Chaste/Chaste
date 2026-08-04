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

#ifndef TESTSEMELEMENT_HPP_
#define TESTSEMELEMENT_HPP_

#include <cxxtest/TestSuite.h>

#include "SemElement.hpp"

// This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

class TestSemElement : public CxxTest::TestSuite
{
private:

    /**
     * Make four 2D nodes at the corners of the unit square. The caller takes ownership: SemElement
     * does not own its nodes (SemMesh does), so every test that builds nodes must delete them.
     */
    std::vector<Node<2>*> MakeUnitSquareNodes()
    {
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
        nodes.push_back(new Node<2>(2, false, 1.0, 1.0));
        nodes.push_back(new Node<2>(3, false, 0.0, 1.0));
        return nodes;
    }

public:

    void TestConstructorWithIndexOnly()
    {
        SemElement<2> element(3);

        TS_ASSERT_EQUALS(element.GetIndex(), 3u);
        TS_ASSERT_EQUALS(element.GetNumNodes(), 0u);
        TS_ASSERT_EQUALS(element.IsDeleted(), false);
        TS_ASSERT(element.rGetNodes().empty());
    }

    void TestConstructorWithNodes()
    {
        std::vector<Node<2>*> nodes = MakeUnitSquareNodes();

        {
            SemElement<2> element(7, nodes);

            TS_ASSERT_EQUALS(element.GetIndex(), 7u);
            TS_ASSERT_EQUALS(element.GetNumNodes(), 4u);
            TS_ASSERT_EQUALS(element.IsDeleted(), false);
            TS_ASSERT_EQUALS(element.rGetNodes().size(), 4u);

            for (unsigned i = 0; i < 4; ++i)
            {
                TS_ASSERT_EQUALS(element.GetNodeGlobalIndex(i), i);

                // This constructor calls RegisterWithNodes(), so membership is already in place
                std::set<unsigned> containing = nodes[i]->rGetContainingElementIndices();
                TS_ASSERT_EQUALS(containing.size(), 1u);
                TS_ASSERT_EQUALS(containing.count(7u), 1u);
            }
        }

        for (unsigned i = 0; i < nodes.size(); ++i)
        {
            delete nodes[i];
        }
    }

    void TestConstructorWithNodesIn3d()
    {
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, false, 0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, false, 1.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(2, false, 0.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(3, false, 0.0, 0.0, 1.0));

        {
            SemElement<3> element(0, nodes);

            TS_ASSERT_EQUALS(element.GetNumNodes(), 4u);
            for (unsigned i = 0; i < 4; ++i)
            {
                TS_ASSERT_EQUALS(element.GetNodeGlobalIndex(i), i);
                TS_ASSERT_EQUALS(nodes[i]->rGetContainingElementIndices().count(0u), 1u);
            }
        }

        for (unsigned i = 0; i < nodes.size(); ++i)
        {
            delete nodes[i];
        }
    }

    void TestRegisterWithNodes()
    {
        std::vector<Node<2>*> nodes = MakeUnitSquareNodes();

        {
            // Build the element via the index-only constructor, which does not register
            SemElement<2> element(2);
            for (unsigned i = 0; i < nodes.size(); ++i)
            {
                element.AddNode(nodes[i]);
            }

            TS_ASSERT_EQUALS(element.GetNumNodes(), 4u);
            for (unsigned i = 0; i < nodes.size(); ++i)
            {
                TS_ASSERT(nodes[i]->rGetContainingElementIndices().empty());
            }

            element.RegisterWithNodes();

            for (unsigned i = 0; i < nodes.size(); ++i)
            {
                TS_ASSERT_EQUALS(nodes[i]->rGetContainingElementIndices().count(2u), 1u);
            }
        }

        for (unsigned i = 0; i < nodes.size(); ++i)
        {
            delete nodes[i];
        }
    }

    void TestUpdateNode()
    {
        std::vector<Node<2>*> nodes = MakeUnitSquareNodes();
        Node<2>* p_replacement = new Node<2>(4, false, 0.5, 0.5);

        {
            SemElement<2> element(1, nodes);

            element.UpdateNode(1, p_replacement);

            // The element now holds the replacement in that slot, with its size unchanged
            TS_ASSERT_EQUALS(element.GetNumNodes(), 4u);
            TS_ASSERT_EQUALS(element.GetNodeGlobalIndex(1), 4u);
            TS_ASSERT_EQUALS(element.rGetNodes()[1], p_replacement);

            // Membership follows the swap in both directions
            TS_ASSERT(nodes[1]->rGetContainingElementIndices().empty());
            TS_ASSERT_EQUALS(p_replacement->rGetContainingElementIndices().count(1u), 1u);

            // The nodes that were not touched are unaffected
            TS_ASSERT_EQUALS(nodes[0]->rGetContainingElementIndices().count(1u), 1u);
            TS_ASSERT_EQUALS(nodes[2]->rGetContainingElementIndices().count(1u), 1u);
            TS_ASSERT_EQUALS(nodes[3]->rGetContainingElementIndices().count(1u), 1u);
        }

        for (unsigned i = 0; i < nodes.size(); ++i)
        {
            delete nodes[i];
        }
        delete p_replacement;
    }

    void TestMarkAsDeleted()
    {
        std::vector<Node<2>*> nodes = MakeUnitSquareNodes();

        {
            SemElement<2> element(5, nodes);
            TS_ASSERT_EQUALS(element.IsDeleted(), false);

            element.MarkAsDeleted();

            TS_ASSERT_EQUALS(element.IsDeleted(), true);

            // Every node must be unregistered, so that containing-element queries made by the
            // forces and by damping no longer see the deleted cell
            for (unsigned i = 0; i < nodes.size(); ++i)
            {
                TS_ASSERT(nodes[i]->rGetContainingElementIndices().empty());
            }

            // The element keeps its own node vector; only the registration is undone
            TS_ASSERT_EQUALS(element.GetNumNodes(), 4u);
        }

        for (unsigned i = 0; i < nodes.size(); ++i)
        {
            delete nodes[i];
        }
    }
};

#endif /*TESTSEMELEMENT_HPP_*/
