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

#ifndef TESTSEMSPHERICALELEMENTMESHGENERATOR_HPP_
#define TESTSEMSPHERICALELEMENTMESHGENERATOR_HPP_

#include <cxxtest/TestSuite.h>

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <vector>

#include "SemEnumerations.hpp"
#include "SemMesh.hpp"
#include "SemSphericalElementMeshGenerator.hpp"

#include "FakePetscSetup.hpp"

class TestSemSphericalElementMeshGenerator : public CxxTest::TestSuite
{
private:

    /**
     * @param rMesh the mesh to measure
     * @return the smallest distance between any two nodes in the mesh
     */
    template <unsigned DIM>
    double GetMinimumNodeSeparation(SemMesh<DIM>& rMesh)
    {
        double min_separation = DBL_MAX;

        for (unsigned i = 0; i < rMesh.GetNumNodes(); ++i)
        {
            for (unsigned j = i + 1; j < rMesh.GetNumNodes(); ++j)
            {
                const double separation
                    = norm_2(rMesh.GetNode(i)->rGetLocation() - rMesh.GetNode(j)->rGetLocation());
                min_separation = std::min(min_separation, separation);
            }
        }

        return min_separation;
    }

    /**
     * @param rMesh the mesh to measure
     * @return the distance from the origin to the outermost node
     */
    template <unsigned DIM>
    double GetMaximumNodeRadius(SemMesh<DIM>& rMesh)
    {
        double max_radius = 0.0;

        for (unsigned i = 0; i < rMesh.GetNumNodes(); ++i)
        {
            max_radius = std::max(max_radius, norm_2(rMesh.GetNode(i)->rGetLocation()));
        }

        return max_radius;
    }

    /**
     * @param rMesh the mesh to measure
     * @param region the region label to count
     * @return the number of nodes in the mesh carrying the given region label
     */
    template <unsigned DIM>
    unsigned CountNodesInRegion(SemMesh<DIM>& rMesh, unsigned region)
    {
        unsigned num_nodes = 0u;

        for (unsigned i = 0; i < rMesh.GetNumNodes(); ++i)
        {
            if (rMesh.GetNode(i)->GetRegion() == region)
            {
                num_nodes++;
            }
        }

        return num_nodes;
    }

    /**
     * The number of other nodes within a shade over one node spacing, i.e. the number of nodes
     * this node is in direct contact with in the packing. A node in the interior of a close
     * packing has a full shell of these: six in 2D and twelve in 3D.
     *
     * @param rMesh the mesh to measure
     * @param nodeIndex the node whose neighbours are counted
     * @param spacing the node spacing of the mesh
     * @return the number of nearest neighbours of the given node
     */
    template <unsigned DIM>
    unsigned CountNearestNeighbours(SemMesh<DIM>& rMesh, unsigned nodeIndex, double spacing)
    {
        unsigned num_neighbours = 0u;

        for (unsigned i = 0; i < rMesh.GetNumNodes(); ++i)
        {
            if (i == nodeIndex)
            {
                continue;
            }

            const double separation
                = norm_2(rMesh.GetNode(i)->rGetLocation() - rMesh.GetNode(nodeIndex)->rGetLocation());
            if (separation < 1.01 * spacing)
            {
                num_neighbours++;
            }
        }

        return num_neighbours;
    }

public:

    /**
     * The generator is asked for a node count and a radius, and must deliver exactly those: the
     * requested number of nodes, with the outermost sitting exactly at the requested radius. The
     * cell is carved from a close packing, so the closest pair of nodes anywhere in it are exactly
     * one node spacing apart.
     */
    void TestNodeCountRadiusAndSpacingAreExact()
    {
        for (unsigned num_nodes : { 2u, 7u, 10u, 57u, 200u })
        {
            SemSphericalElementMeshGenerator<3> generator(num_nodes, 0.5);
            auto p_mesh = generator.GetMesh();

            TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), num_nodes);
            TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 1u);
            TS_ASSERT_EQUALS(p_mesh->GetElement(0)->GetNumNodes(), num_nodes);

            TS_ASSERT_DELTA(GetMaximumNodeRadius<3>(*p_mesh), 0.5, 1e-12);
            TS_ASSERT_DELTA(GetMinimumNodeSeparation<3>(*p_mesh), generator.GetNodeSpacing(), 1e-12);
        }

        // A different radius scales the whole cell, and with it the node spacing
        SemSphericalElementMeshGenerator<3> small_generator(200, 0.5);
        SemSphericalElementMeshGenerator<3> large_generator(200, 1.5);

        TS_ASSERT_DELTA(GetMaximumNodeRadius<3>(*large_generator.GetMesh()), 1.5, 1e-12);
        TS_ASSERT_DELTA(large_generator.GetNodeSpacing(), 3.0 * small_generator.GetNodeSpacing(), 1e-12);
    }

    /**
     * Sandersius & Newman (2008) Section 1.2 relate the equilibrium element separation to the cell
     * radius by r_eq = 2 * R_cell * (p / N)^(1/3), where p is the packing density of spheres. The
     * generator does not use this relation - the spacing falls out of the geometry instead - so
     * agreement with it is a check that the packing really is dense, and that the paper's
     * parameterisation can be used to choose N for a target r_eq.
     */
    void TestNodeSpacingMatchesTheClosePackedRelation()
    {
        // The density of a hexagonal close packing of spheres
        const double packing_density = 0.7405;
        const double cell_radius = 0.5;

        for (unsigned num_nodes : { 200u, 1000u })
        {
            SemSphericalElementMeshGenerator<3> generator(num_nodes, cell_radius);

            const double paper_spacing = 2.0 * cell_radius
                * std::pow(packing_density / static_cast<double>(num_nodes), 1.0 / 3.0);

            TS_ASSERT_DELTA(generator.GetNodeSpacing(), paper_spacing, 0.02 * paper_spacing);
        }
    }

    /**
     * The key property of the radial boundary rule: everything it labels interior really is
     * interior, in the sense of having a complete shell of nearest neighbours. This checks the
     * labelling against an independent, purely local criterion.
     */
    void TestInteriorNodesHaveAFullNeighbourShell()
    {
        SemSphericalElementMeshGenerator<3> generator_3d(500, 0.5);
        auto p_mesh_3d = generator_3d.GetMesh();

        unsigned num_interior_3d = 0u;
        for (unsigned i = 0; i < p_mesh_3d->GetNumNodes(); ++i)
        {
            if (p_mesh_3d->GetNode(i)->GetRegion() == SEM_INTERIOR_REGION)
            {
                TS_ASSERT_EQUALS(CountNearestNeighbours<3>(*p_mesh_3d, i, generator_3d.GetNodeSpacing()), 12u);
                num_interior_3d++;
            }
        }
        TS_ASSERT_LESS_THAN(0u, num_interior_3d);

        SemSphericalElementMeshGenerator<2> generator_2d(200, 0.5);
        auto p_mesh_2d = generator_2d.GetMesh();

        unsigned num_interior_2d = 0u;
        for (unsigned i = 0; i < p_mesh_2d->GetNumNodes(); ++i)
        {
            if (p_mesh_2d->GetNode(i)->GetRegion() == SEM_INTERIOR_REGION)
            {
                TS_ASSERT_EQUALS(CountNearestNeighbours<2>(*p_mesh_2d, i, generator_2d.GetNodeSpacing()), 6u);
                num_interior_2d++;
            }
        }
        TS_ASSERT_LESS_THAN(0u, num_interior_2d);
    }

    /**
     * The converse check: boundary nodes are, on average, less well connected than interior ones,
     * as surface elements of a real aggregate are. Individual boundary nodes may still have a full
     * shell, since the region is a shell of finite thickness rather than only the outermost layer.
     */
    void TestBoundaryNodesAreLessWellConnectedThanInteriorNodes()
    {
        SemSphericalElementMeshGenerator<3> generator(500, 0.5);
        auto p_mesh = generator.GetMesh();

        double boundary_total = 0.0;
        double interior_total = 0.0;
        unsigned num_boundary = 0u;
        unsigned num_interior = 0u;

        for (unsigned i = 0; i < p_mesh->GetNumNodes(); ++i)
        {
            const unsigned coordination
                = CountNearestNeighbours<3>(*p_mesh, i, generator.GetNodeSpacing());

            if (p_mesh->GetNode(i)->GetRegion() == SEM_BOUNDARY_REGION)
            {
                boundary_total += static_cast<double>(coordination);
                num_boundary++;
            }
            else
            {
                interior_total += static_cast<double>(coordination);
                num_interior++;
            }
        }

        TS_ASSERT_LESS_THAN(0u, num_boundary);
        TS_ASSERT_LESS_THAN(0u, num_interior);
        TS_ASSERT_LESS_THAN(boundary_total / num_boundary, interior_total / num_interior);
    }

    /**
     * Every node is labelled either interior or boundary, and the boundary region is the outer
     * shell: no boundary node lies further from the surface than any interior node.
     */
    void TestRegionsPartitionTheNodesIntoAnInteriorAndAShell()
    {
        SemSphericalElementMeshGenerator<3> generator(500, 0.5);
        auto p_mesh = generator.GetMesh();

        const unsigned num_boundary = CountNodesInRegion<3>(*p_mesh, SEM_BOUNDARY_REGION);
        const unsigned num_interior = CountNodesInRegion<3>(*p_mesh, SEM_INTERIOR_REGION);
        TS_ASSERT_EQUALS(num_boundary + num_interior, p_mesh->GetNumNodes());

        double innermost_boundary_radius = DBL_MAX;
        double outermost_interior_radius = 0.0;

        for (unsigned i = 0; i < p_mesh->GetNumNodes(); ++i)
        {
            const double radius = norm_2(p_mesh->GetNode(i)->rGetLocation());

            if (p_mesh->GetNode(i)->GetRegion() == SEM_BOUNDARY_REGION)
            {
                innermost_boundary_radius = std::min(innermost_boundary_radius, radius);
            }
            else
            {
                outermost_interior_radius = std::max(outermost_interior_radius, radius);
            }
        }

        TS_ASSERT_LESS_THAN(outermost_interior_radius, innermost_boundary_radius);
    }

    /**
     * The cortex thickness is measured in node spacings, so it interpolates between labelling only
     * the outermost shell of nodes and labelling every node in the cell.
     */
    void TestBoundaryThicknessControlsTheCortexDepth()
    {
        const unsigned num_nodes = 1000u;

        std::vector<unsigned> boundary_counts;
        for (double thickness : { 0.0, 0.5, 1.0, 2.0 })
        {
            SemSphericalElementMeshGenerator<3> generator(num_nodes, 0.5, thickness);
            boundary_counts.push_back(CountNodesInRegion<3>(*generator.GetMesh(), SEM_BOUNDARY_REGION));
        }

        // A thickness of zero labels exactly the nodes lying at the cell radius, of which there is
        // at least one; each increase in thickness then draws in strictly more nodes
        TS_ASSERT_LESS_THAN(0u, boundary_counts[0]);
        for (unsigned i = 1; i < boundary_counts.size(); ++i)
        {
            TS_ASSERT_LESS_THAN(boundary_counts[i - 1], boundary_counts[i]);
        }

        // A thickness exceeding the cell radius leaves no interior
        SemSphericalElementMeshGenerator<3> thick_generator(num_nodes, 0.5, 100.0);
        TS_ASSERT_EQUALS(CountNodesInRegion<3>(*thick_generator.GetMesh(), SEM_BOUNDARY_REGION), num_nodes);
    }

    /**
     * Lattice sites are frequently exactly equidistant from the centre of the cell, so which of
     * them are kept depends on how ties are broken. The generator must therefore produce the same
     * mesh every time rather than one that depends on the ordering the sort happens to see.
     */
    void TestGeneratorIsDeterministic()
    {
        SemSphericalElementMeshGenerator<3> first_generator(300, 0.5);
        SemSphericalElementMeshGenerator<3> second_generator(300, 0.5);

        auto p_first_mesh = first_generator.GetMesh();
        auto p_second_mesh = second_generator.GetMesh();

        TS_ASSERT_EQUALS(p_first_mesh->GetNumNodes(), p_second_mesh->GetNumNodes());

        for (unsigned i = 0; i < p_first_mesh->GetNumNodes(); ++i)
        {
            for (unsigned dim = 0; dim < 3; ++dim)
            {
                TS_ASSERT_DELTA(p_first_mesh->GetNode(i)->rGetLocation()[dim],
                                p_second_mesh->GetNode(i)->rGetLocation()[dim], 0.0);
            }
            TS_ASSERT_EQUALS(p_first_mesh->GetNode(i)->GetRegion(),
                             p_second_mesh->GetNode(i)->GetRegion());
        }
    }

    /**
     * In 2D the cell is a disc carved from a triangular lattice.
     */
    void TestSphericalElementIn2d()
    {
        const unsigned num_nodes = 200u;
        SemSphericalElementMeshGenerator<2> generator(num_nodes, 0.5);
        auto p_mesh = generator.GetMesh();

        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), num_nodes);
        TS_ASSERT_DELTA(GetMaximumNodeRadius<2>(*p_mesh), 0.5, 1e-12);
        TS_ASSERT_DELTA(GetMinimumNodeSeparation<2>(*p_mesh), generator.GetNodeSpacing(), 1e-12);
    }

    /**
     * In 1D the cell is a line segment, and the nodes nearest the centre are an evenly spaced run
     * either side of it. With an even number of nodes the run cannot be symmetric, and the tie is
     * broken towards the negative end.
     */
    void TestSphericalElementIn1d()
    {
        SemSphericalElementMeshGenerator<1> odd_generator(5, 0.5);
        auto p_odd_mesh = odd_generator.GetMesh();

        TS_ASSERT_EQUALS(p_odd_mesh->GetNumNodes(), 5u);
        TS_ASSERT_DELTA(odd_generator.GetNodeSpacing(), 0.25, 1e-12);

        const std::vector<double> expected_positions = { 0.0, -0.25, 0.25, -0.5, 0.5 };
        for (unsigned i = 0; i < expected_positions.size(); ++i)
        {
            TS_ASSERT_DELTA(p_odd_mesh->GetNode(i)->rGetLocation()[0], expected_positions[i], 1e-12);
        }

        // Only the node at the centre lies more than one spacing from an end
        TS_ASSERT_EQUALS(CountNodesInRegion<1>(*p_odd_mesh, SEM_INTERIOR_REGION), 1u);
        TS_ASSERT_EQUALS(CountNodesInRegion<1>(*p_odd_mesh, SEM_BOUNDARY_REGION), 4u);

        SemSphericalElementMeshGenerator<1> even_generator(2, 0.5);
        auto p_even_mesh = even_generator.GetMesh();

        TS_ASSERT_EQUALS(p_even_mesh->GetNumNodes(), 2u);
        TS_ASSERT_DELTA(p_even_mesh->GetNode(0)->rGetLocation()[0], 0.0, 1e-12);
        TS_ASSERT_DELTA(p_even_mesh->GetNode(1)->rGetLocation()[0], -0.5, 1e-12);
    }

    void TestExceptions()
    {
        TS_ASSERT_THROWS_THIS(SemSphericalElementMeshGenerator<3> generator(0),
            "SemSphericalElementMeshGenerator: numNodes must be >= 2");
        TS_ASSERT_THROWS_THIS(SemSphericalElementMeshGenerator<3> generator(1),
            "SemSphericalElementMeshGenerator: numNodes must be >= 2");

        TS_ASSERT_THROWS_THIS(SemSphericalElementMeshGenerator<3> generator(100, 0.0),
            "SemSphericalElementMeshGenerator: cellRadius must be positive");
        TS_ASSERT_THROWS_THIS(SemSphericalElementMeshGenerator<3> generator(100, -1.0),
            "SemSphericalElementMeshGenerator: cellRadius must be positive");

        TS_ASSERT_THROWS_THIS(SemSphericalElementMeshGenerator<3> generator(100, 0.5, -1e-6),
            "SemSphericalElementMeshGenerator: boundaryThickness must be non-negative");

        // A thickness of exactly zero is permitted, and labels the outermost shell of nodes
        TS_ASSERT_THROWS_NOTHING(SemSphericalElementMeshGenerator<3> generator(100, 0.5, 0.0));
    }
};

#endif /*TESTSEMSPHERICALELEMENTMESHGENERATOR_HPP_*/
