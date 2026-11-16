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

#include <algorithm>
#include <cfloat>
#include <cmath>

#include "SemEnumerations.hpp"
#include "SemSingleElementMeshGenerator.hpp"
#include "SemMesh.hpp"

// This test is never run in parallel
#include "FakePetscSetup.hpp"


class TestSemSingleElementMeshGenerator : public CxxTest::TestSuite
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

    /**
     * @param rMesh the mesh to inspect
     * @param nodeIndex the index of the node whose neighbours are counted
     * @param separation the separation at which a node counts as a nearest neighbour
     * @return the number of nodes exactly one separation away from the given node
     */
    template <unsigned DIM>
    unsigned CountNearestNeighbours(SemMesh<DIM>& rMesh, unsigned nodeIndex, double separation)
    {
        unsigned num_neighbours = 0u;
        for (unsigned i = 0; i < rMesh.GetNumNodes(); ++i)
        {
            if (i == nodeIndex)
            {
                continue;
            }
            const double distance = norm_2(rMesh.GetNode(i)->rGetLocation() - rMesh.GetNode(nodeIndex)->rGetLocation());
            if (fabs(distance - separation) < 1e-6)
            {
                num_neighbours++;
            }
        }
        return num_neighbours;
    }

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

    void TestClosePacked1DMatchesCubic()
    {
        // There is only one way to space points along a line, so the lattice type is immaterial in 1D
        SemSingleElementMeshGenerator<1> cubic_gen({5}, 2.0, SemLatticeType::SEM_LATTICE_CUBIC);
        SemSingleElementMeshGenerator<1> packed_gen({5}, 2.0, SemLatticeType::SEM_LATTICE_CLOSE_PACKED);

        auto p_cubic_mesh = cubic_gen.GetMesh();
        auto p_packed_mesh = packed_gen.GetMesh();

        TS_ASSERT_EQUALS(p_packed_mesh->GetNumNodes(), p_cubic_mesh->GetNumNodes());
        for (unsigned i = 0; i < p_cubic_mesh->GetNumNodes(); ++i)
        {
            TS_ASSERT_DELTA(p_packed_mesh->GetNode(i)->rGetLocation()[0],
                            p_cubic_mesh->GetNode(i)->rGetLocation()[0], 1e-6);
        }
    }

    void TestClosePacked2D()
    {
        // Node spacing = scaleFactor / numNodes[0] = 4.0 / 4 = 1.0
        SemSingleElementMeshGenerator<2> generator({4, 3}, 4.0, SemLatticeType::SEM_LATTICE_CLOSE_PACKED);
        auto p_mesh = generator.GetMesh();

        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 12u);

        // Odd-numbered rows are offset half a spacing in x, and rows are sqrt(3)/2 spacings apart in y
        const double row_offset = 0.5;
        const double row_spacing = 0.5 * sqrt(3.0);

        for (unsigned j = 0; j < 3u; ++j)
        {
            for (unsigned i = 0; i < 4u; ++i)
            {
                const unsigned node_index = i + 4u * j;
                const double expected_x = static_cast<double>(i) + (j % 2u == 1u ? row_offset : 0.0);
                const double expected_y = static_cast<double>(j) * row_spacing;

                TS_ASSERT_DELTA(p_mesh->GetNode(node_index)->rGetLocation()[0], expected_x, 1e-6);
                TS_ASSERT_DELTA(p_mesh->GetNode(node_index)->rGetLocation()[1], expected_y, 1e-6);
            }
        }

        // Staggering the rows must not change which nodes are on the boundary of the node grid
        TS_ASSERT_EQUALS(p_mesh->GetNode(0u)->GetRegion(), SEM_BOUNDARY_REGION);
        TS_ASSERT_EQUALS(p_mesh->GetNode(4u)->GetRegion(), SEM_BOUNDARY_REGION);
        TS_ASSERT_EQUALS(p_mesh->GetNode(11u)->GetRegion(), SEM_BOUNDARY_REGION);
        TS_ASSERT_EQUALS(p_mesh->GetNode(5u)->GetRegion(), SEM_INTERIOR_REGION);
        TS_ASSERT_EQUALS(p_mesh->GetNode(6u)->GetRegion(), SEM_INTERIOR_REGION);
    }

    void TestClosePacked2DIsTriangularLattice()
    {
        // Node spacing = 5.0 / 5 = 1.0
        SemSingleElementMeshGenerator<2> generator({5, 5}, 5.0, SemLatticeType::SEM_LATTICE_CLOSE_PACKED);
        auto p_mesh = generator.GetMesh();

        // No pair of nodes may be closer than one spacing, and an interior node of a triangular
        // lattice has six nearest neighbours
        TS_ASSERT_DELTA(GetMinimumNodeSeparation<2>(*p_mesh), 1.0, 1e-6);
        TS_ASSERT_EQUALS(CountNearestNeighbours<2>(*p_mesh, 2u + 5u * 2u, 1.0), 6u);
    }

    void TestClosePacked3DIsHexagonalClosePacked()
    {
        // Node spacing = 5.0 / 5 = 1.0
        SemSingleElementMeshGenerator<3> generator({5, 5, 5}, 5.0, SemLatticeType::SEM_LATTICE_CLOSE_PACKED);
        auto p_mesh = generator.GetMesh();

        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 125u);

        // Layers are sqrt(2/3) spacings apart, which is closer than the cubic lattice's one
        // spacing, so the same number of nodes occupies a shorter box in z
        TS_ASSERT_DELTA(p_mesh->GetNode(100u)->rGetLocation()[2], 4.0 * sqrt(2.0 / 3.0), 1e-6);

        // The defining property of a close packing: no pair of nodes is closer than one spacing,
        // and a node in the interior has twelve nearest neighbours
        TS_ASSERT_DELTA(GetMinimumNodeSeparation<3>(*p_mesh), 1.0, 1e-6);
        TS_ASSERT_EQUALS(CountNearestNeighbours<3>(*p_mesh, 2u + 5u * (2u + 5u * 2u), 1.0), 12u);
    }

    void TestCubicLatticeIsTheDefault()
    {
        // Omitting the lattice argument must reproduce the original axis-aligned grid
        SemSingleElementMeshGenerator<3> default_gen({3, 3, 3}, 3.0);
        SemSingleElementMeshGenerator<3> cubic_gen({3, 3, 3}, 3.0, SemLatticeType::SEM_LATTICE_CUBIC);

        auto p_default_mesh = default_gen.GetMesh();
        auto p_cubic_mesh = cubic_gen.GetMesh();

        for (unsigned i = 0; i < p_cubic_mesh->GetNumNodes(); ++i)
        {
            TS_ASSERT_DELTA(norm_2(p_default_mesh->GetNode(i)->rGetLocation()
                                   - p_cubic_mesh->GetNode(i)->rGetLocation()), 0.0, 1e-6);
        }

        // A cubic lattice has 2*DIM nearest neighbours in the interior
        TS_ASSERT_EQUALS(CountNearestNeighbours<3>(*p_cubic_mesh, 1u + 3u * (1u + 3u * 1u), 1.0), 6u);
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