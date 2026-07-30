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

#ifndef TESTSEMMULTIELEMENTMESHGENERATOR_HPP_
#define TESTSEMMULTIELEMENTMESHGENERATOR_HPP_

#include <cxxtest/TestSuite.h>

#include <algorithm>
#include <cfloat>
#include <cmath>

#include "SemEnumerations.hpp"
#include "SemMultiElementMeshGenerator.hpp"
#include "SemSingleElementMeshGenerator.hpp"
#include "SemMesh.hpp"

// This test is never run in parallel
#include "FakePetscSetup.hpp"


class TestSemMultiElementMeshGenerator : public CxxTest::TestSuite
{
private:

    /**
     * @param rMesh the mesh to inspect
     * @return the smallest distance between any two nodes in the mesh
     */
    template <unsigned DIM>
    double GetMinimumNodeSeparation(SemMesh<DIM>& rMesh)
    {
        double min_separation = DBL_MAX;
        for (unsigned i = 0; i < rMesh.GetNumNodes(); ++i)
        {
            for (unsigned j = i + 1u; j < rMesh.GetNumNodes(); ++j)
            {
                const double separation = norm_2(rMesh.GetNode(i)->rGetLocation() - rMesh.GetNode(j)->rGetLocation());
                min_separation = std::min(min_separation, separation);
            }
        }
        return min_separation;
    }

public:
    void Test1DOneElem()
    {
        // Mesh should be identical to the single element generator when generating one element

        SemSingleElementMeshGenerator<1> gen_1({5}, 2.0);
        auto p_mesh_single = gen_1.GetMesh();

        SemMultiElementMeshGenerator<1> gen_2({5}, {1}, 2.0);
        auto p_mesh_multi = gen_2.GetMesh();

        TS_ASSERT_EQUALS(p_mesh_multi->GetNumNodes(), 5u);

        for (unsigned i = 0; i < p_mesh_multi->GetNumAllNodes(); ++i)
        {
            TS_ASSERT_DELTA(p_mesh_single->GetNode(i)->rGetLocation()[0], p_mesh_multi->GetNode(i)->rGetLocation()[0], 1e-6);
        }
    }

    void Test1DMultipleElements()
    {
        // A 1D mesh with more than one element previously indexed mElemSpacing (an array of
        // length DIM=1) out of bounds; check the element offsets, and hence node positions,
        // are laid out correctly.
        SemMultiElementMeshGenerator<1> generator({3}, {2}, 3.0);
        auto p_mesh = generator.GetMesh();

        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 2u);
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 6u);

        // Node spacing = scaleFactor / numNodesPerElem = 3.0 / 3 = 1.0; element spacing = 3.0,
        // so the two three-node elements sit at x in {0,1,2} and {3,4,5} respectively.
        const double expected_x[6] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
        for (unsigned i = 0; i < 6u; ++i)
        {
            TS_ASSERT_DELTA(p_mesh->GetNode(i)->rGetLocation()[0], expected_x[i], 1e-6);
        }
    }

    void Test2DOneElem()
    {
        // Mesh should be identical to the single element generator when generating one element

        SemSingleElementMeshGenerator<2> gen_1({4, 5}, 3.0);
        auto p_mesh_single = gen_1.GetMesh();

        SemMultiElementMeshGenerator<2> gen_2({4, 5}, {1, 1}, 3.0);
        auto p_mesh_multi = gen_2.GetMesh();

        TS_ASSERT_EQUALS(p_mesh_multi->GetNumNodes(), 20u);

        for (unsigned i = 0; i < p_mesh_multi->GetNumAllNodes(); ++i)
        {
            TS_ASSERT_DELTA(p_mesh_single->GetNode(i)->rGetLocation()[0], p_mesh_multi->GetNode(i)->rGetLocation()[0], 1e-6);
            TS_ASSERT_DELTA(p_mesh_single->GetNode(i)->rGetLocation()[1], p_mesh_multi->GetNode(i)->rGetLocation()[1], 1e-6);
        }
    }

    void Test3DOneElem()
    {
        // Mesh should be identical to the single element generator when generating one element

        SemSingleElementMeshGenerator<3> gen_1({4, 3, 2}, 5.0);
        auto p_mesh_single = gen_1.GetMesh();

        SemMultiElementMeshGenerator<3> gen_2({4, 3, 2}, {1, 1, 1}, 5.0);
        auto p_mesh_multi = gen_2.GetMesh();

        TS_ASSERT_EQUALS(p_mesh_multi->GetNumNodes(), 24u);

        for (unsigned i = 0; i < p_mesh_multi->GetNumAllNodes(); ++i)
        {
            TS_ASSERT_DELTA(p_mesh_single->GetNode(i)->rGetLocation()[0], p_mesh_multi->GetNode(i)->rGetLocation()[0], 1e-6);
            TS_ASSERT_DELTA(p_mesh_single->GetNode(i)->rGetLocation()[1], p_mesh_multi->GetNode(i)->rGetLocation()[1], 1e-6);
            TS_ASSERT_DELTA(p_mesh_single->GetNode(i)->rGetLocation()[2], p_mesh_multi->GetNode(i)->rGetLocation()[2], 1e-6);
        }
    }

    void TestClosePackedElementLattice2D()
    {
        // Node spacing = 2.0 / 2 = 1.0, so each element is a unit square of four nodes
        SemMultiElementMeshGenerator<2> generator({2, 2}, {2, 2}, 2.0,
            SEM_LATTICE_CUBIC, SEM_LATTICE_CLOSE_PACKED);
        auto p_mesh = generator.GetMesh();

        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 4u);
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 16u);

        // The limiting direction is the diagonal to a staggered neighbour, along which a unit
        // square is (0.5 + sqrt(3)/2) wide; adding one node spacing of clearance gives the spacing
        const double elem_spacing = 1.5 + 0.5 * sqrt(3.0);

        // The first node of each element sits on the element offset, as the first node of an
        // element is at the element's local origin
        const double expected_x[4] = {0.0, elem_spacing, 0.5 * elem_spacing, 1.5 * elem_spacing};
        const double expected_y[4] = {0.0, 0.0, 0.5 * sqrt(3.0) * elem_spacing, 0.5 * sqrt(3.0) * elem_spacing};

        for (unsigned elem = 0; elem < 4u; ++elem)
        {
            const unsigned first_node = 4u * elem;
            TS_ASSERT_DELTA(p_mesh->GetNode(first_node)->rGetLocation()[0], expected_x[elem], 1e-6);
            TS_ASSERT_DELTA(p_mesh->GetNode(first_node)->rGetLocation()[1], expected_y[elem], 1e-6);
        }
    }

    void TestClosePackedNodeLatticeTightensElementSpacing()
    {
        // Node spacing = 3.0 / 3 = 1.0. Close packing the nodes brings the three rows of an
        // element from two node spacings apart to sqrt(3) apart, so the elements above and below
        // may move closer by the same amount.
        SemMultiElementMeshGenerator<2> cubic_gen({3, 3}, {1, 2}, 3.0, SEM_LATTICE_CUBIC);
        SemMultiElementMeshGenerator<2> packed_gen({3, 3}, {1, 2}, 3.0, SEM_LATTICE_CLOSE_PACKED);

        auto p_cubic_mesh = cubic_gen.GetMesh();
        auto p_packed_mesh = packed_gen.GetMesh();

        // The first node of the second element is offset by the element spacing in y
        TS_ASSERT_DELTA(p_cubic_mesh->GetNode(9u)->rGetLocation()[1], 3.0, 1e-6);
        TS_ASSERT_DELTA(p_packed_mesh->GetNode(9u)->rGetLocation()[1], sqrt(3.0) + 1.0, 1e-6);
    }

    void TestElementsNeverOverlapForAnyLatticeCombination()
    {
        // Node spacing = 3.0 / 3 = 1.0. However the nodes and the elements are laid out, nodes in
        // different elements must be no closer than nodes within an element, so that the mesh is
        // usable as an initial condition without the intercellular forces blowing up.
        const SemLatticeType lattices[2] = {SEM_LATTICE_CUBIC, SEM_LATTICE_CLOSE_PACKED};

        for (const SemLatticeType node_lattice : lattices)
        {
            for (const SemLatticeType element_lattice : lattices)
            {
                SemMultiElementMeshGenerator<3> generator({3, 3, 3}, {2, 2, 2}, 3.0,
                    node_lattice, element_lattice);
                auto p_mesh = generator.GetMesh();

                TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 8u);
                TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 216u);
                TS_ASSERT_DELTA(GetMinimumNodeSeparation<3>(*p_mesh), 1.0, 1e-6);
            }
        }
    }

    void TestCubicLatticesAreTheDefault()
    {
        // Omitting the lattice arguments must reproduce the original axis-aligned layout
        SemMultiElementMeshGenerator<3> default_gen({2, 3, 2}, {2, 1, 2}, 4.0);
        SemMultiElementMeshGenerator<3> cubic_gen({2, 3, 2}, {2, 1, 2}, 4.0,
            SEM_LATTICE_CUBIC, SEM_LATTICE_CUBIC);

        auto p_default_mesh = default_gen.GetMesh();
        auto p_cubic_mesh = cubic_gen.GetMesh();

        TS_ASSERT_EQUALS(p_default_mesh->GetNumNodes(), p_cubic_mesh->GetNumNodes());
        for (unsigned i = 0; i < p_cubic_mesh->GetNumNodes(); ++i)
        {
            TS_ASSERT_DELTA(norm_2(p_default_mesh->GetNode(i)->rGetLocation()
                                   - p_cubic_mesh->GetNode(i)->rGetLocation()), 0.0, 1e-6);
        }
    }

    void TestExceptions()
    {
        TS_ASSERT_THROWS_THIS(SemMultiElementMeshGenerator<1>({0}, {1}, 2.0),
            "SemMultiElementMeshGenerator: each entry of numNodesPerElem must be >= 1");

        TS_ASSERT_THROWS_THIS(SemMultiElementMeshGenerator<2>({1, 3}, {1, 0}, 3.0),
            "SemMultiElementMeshGenerator: each entry of numElems must be >= 1");

        TS_ASSERT_THROWS_THIS(SemMultiElementMeshGenerator<3>({2, 3, 4}, {1, 2, 3}, -1.0),
            "SemMultiElementMeshGenerator: scaleFactor must be positive");
    }
};

#endif /*TESTSEMMULTIELEMENTMESHGENERATOR_HPP_*/