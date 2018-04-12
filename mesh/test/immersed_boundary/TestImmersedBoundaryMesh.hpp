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

#ifndef TESTIMMERSEDBOUNDARYMESH_HPP_
#define TESTIMMERSEDBOUNDARYMESH_HPP_

// Needed for test framework
#include <cxxtest/TestSuite.h>

#include "ImmersedBoundaryEnumerations.hpp"
#include "ImmersedBoundaryHoneycombMeshGenerator.hpp"
#include "ImmersedBoundaryPalisadeMeshGenerator.hpp"
#include "ImmersedBoundaryMesh.hpp"
#include "Node.hpp"
#include "RandomNumberGenerator.hpp"
#include "UblasVectorInclude.hpp"

// This test is never run in parallel
#include "FakePetscSetup.hpp"


class TestImmersedBoundaryMesh : public CxxTest::TestSuite
{
public:
    void TestSolveNodeAndElementMapping()
    {
    }

    void TestSetupFluidVelocityGrids()
    {
    }

    void TestArchiving()
    {
    }

    void TestImmersedBoundaryElementAndLaminaIterators()
    {
        // Empty mesh object with no elements (we only have a != operator available)
        {
            ImmersedBoundaryMesh<2, 2> ib_mesh;
            TS_ASSERT_EQUALS(ib_mesh.GetElementIteratorBegin() != ib_mesh.GetElementIteratorEnd(), false);
            TS_ASSERT_EQUALS(ib_mesh.GetLaminaIteratorBegin() != ib_mesh.GetLaminaIteratorEnd(), false);
        }

        // Mesh with elements and laminas
        {
            // Make a few nodes
            std::vector<Node<2>*> nodes;
            nodes.push_back(new Node<2>(0, true, 0.1, 0.1));
            nodes.push_back(new Node<2>(1, true, 0.2, 0.1));
            nodes.push_back(new Node<2>(2, true, 0.3, 0.2));

            // Make two elements
            std::vector<ImmersedBoundaryElement<2, 2>*> elements;
            elements.push_back(new ImmersedBoundaryElement<2, 2>(0, nodes));
            elements.push_back(new ImmersedBoundaryElement<2, 2>(1, nodes));

            // Make four laminas
            std::vector<ImmersedBoundaryElement<1, 2>*> lams;
            lams.push_back(new ImmersedBoundaryElement<1, 2>(0, nodes));
            lams.push_back(new ImmersedBoundaryElement<1, 2>(1, nodes));
            lams.push_back(new ImmersedBoundaryElement<1, 2>(2, nodes));
            lams.push_back(new ImmersedBoundaryElement<1, 2>(3, nodes));

            // Make a mesh
            ImmersedBoundaryMesh<2, 2> ib_mesh(nodes, elements, lams);

            unsigned elem_count = 0u;
            for (auto elem_it = ib_mesh.GetElementIteratorBegin(); elem_it != ib_mesh.GetElementIteratorEnd(); ++elem_it)
            {
                TS_ASSERT_EQUALS(elem_it->GetIndex(), elem_count);
                elem_count++;
            }

            unsigned lam_count = 0u;
            for (auto lam_it = ib_mesh.GetLaminaIteratorBegin(); lam_it != ib_mesh.GetLaminaIteratorBegin(); ++lam_it)
            {
                TS_ASSERT_EQUALS(lam_it->GetIndex(), lam_count);
                lam_count++;
            }
        }
    }

    void TestSetAndGetMethods()
    {
        // Default-construct a mesh object
        ImmersedBoundaryMesh<2, 2> ib_mesh;

        ib_mesh.SetNeighbourDist(1.23);
        ib_mesh.SetCharacteristicNodeSpacing(2.34);
        ib_mesh.SetElementDivisionSpacing(3.45);
        ib_mesh.SetNumGridPtsXAndY(678u);

        TS_ASSERT_DELTA(ib_mesh.GetNeighbourDist(), 1.23, 1e-12);
        TS_ASSERT_DELTA(ib_mesh.GetCharacteristicNodeSpacing(), 2.34, 1e-12);
        TS_ASSERT_DELTA(ib_mesh.GetElementDivisionSpacing(), 3.45, 1e-12);
        TS_ASSERT_EQUALS(ib_mesh.GetNumGridPtsX(), 678u);
        TS_ASSERT_EQUALS(ib_mesh.GetNumGridPtsY(), 678u);


        ib_mesh.SetNumGridPtsX(456u);
        ib_mesh.SetNumGridPtsY(567u);
        TS_ASSERT_EQUALS(ib_mesh.GetNumGridPtsX(), 456u);
        TS_ASSERT_EQUALS(ib_mesh.GetNumGridPtsY(), 567u);


    }

    void TestGetVectorFromAtoB()
    {
        // Create a small mesh
        ImmersedBoundaryHoneycombMeshGenerator gen(1, 1, 3, 0.1, 0.3);
        ImmersedBoundaryMesh<2, 2>* p_mesh = gen.GetMesh();

        // Two helper vectors
        c_vector<double, 2> x_unit = unit_vector<double>(2, 0);
        c_vector<double, 2> y_unit = unit_vector<double>(2, 1);

        c_vector<double, 2> point_a;
        c_vector<double, 2> point_b;
        c_vector<double, 2> vec_a2b;

        // Immersed boundary meshes are always doubly-periodic on the square [0, 1] x [0, 1]

        // Regular cases where periodicity plays no part
        point_a = 0.0 * x_unit;
        point_b = 0.3 * x_unit + 0.4 * y_unit;
        vec_a2b = p_mesh->GetVectorFromAtoB(point_a, point_b);

        TS_ASSERT_DELTA(vec_a2b[0], point_b[0], 1e-6);
        TS_ASSERT_DELTA(vec_a2b[1], point_b[1], 1e-6);
        TS_ASSERT_DELTA(norm_2(vec_a2b), 0.5, 1e-6);

        // x-periodicity
        point_a = 0.2 * x_unit;
        point_b = 0.8 * x_unit;
        vec_a2b = p_mesh->GetVectorFromAtoB(point_a, point_b);

        TS_ASSERT_DELTA(vec_a2b[0], -0.4, 1e-6);
        TS_ASSERT_DELTA(vec_a2b[1], 0.0, 1e-6);

        // y-periodicity
        point_a = 0.2 * y_unit;
        point_b = 0.8 * y_unit;
        vec_a2b = p_mesh->GetVectorFromAtoB(point_a, point_b);

        TS_ASSERT_DELTA(vec_a2b[0], 0.0, 1e-6);
        TS_ASSERT_DELTA(vec_a2b[1], -0.4, 1e-6);

        // x and y-periodicity
        point_a = 0.1 * x_unit + 0.1 * y_unit;
        point_b = 0.9 * x_unit + 0.9 * y_unit;
        vec_a2b = p_mesh->GetVectorFromAtoB(point_a, point_b);

        TS_ASSERT_DELTA(vec_a2b[0], -0.2, 1e-6);
        TS_ASSERT_DELTA(vec_a2b[1], -0.2, 1e-6);
    }

    void TestGetSkewnessOfElementMassDistributionAboutAxis()
    {
        // A square should have no skewness about any axis
        {
            std::vector<Node<2>*> nodes;
            nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
            nodes.push_back(new Node<2>(1, true, 0.1, 0.0));
            nodes.push_back(new Node<2>(2, true, 0.1, 0.1));
            nodes.push_back(new Node<2>(3, true, 0.0, 0.1));

            std::vector<ImmersedBoundaryElement<2, 2>*> elems;
            elems.push_back(new ImmersedBoundaryElement<2, 2>(0, nodes));

            ImmersedBoundaryMesh<2,2> mesh(nodes, elems);

            for (unsigned i = 0; i < 16; i++)
            {
                double theta = 2.0 * M_PI * (double)i / 17.0;

                c_vector<double, 2> axis;
                axis[0] = cos(theta);
                axis[1] = sin(theta);

                TS_ASSERT_DELTA(mesh.GetSkewnessOfElementMassDistributionAboutAxis(0, axis), 0.0, 1e-12);
            }
        }

        // A triangle should have skewness
        {
            std::vector<Node<2>*> nodes;
            nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
            nodes.push_back(new Node<2>(1, true, 0.1, 0.0));
            nodes.push_back(new Node<2>(2, true, 0.1, 0.1));

            std::vector<ImmersedBoundaryElement<2, 2>*> elems;
            elems.push_back(new ImmersedBoundaryElement<2, 2>(0, nodes));

            ImmersedBoundaryMesh<2,2> mesh(nodes, elems);

            c_vector<double, 2> axis;
            axis[0] = 0.0;
            axis[1] = 1.0;

            double hand_calculated_skewness = -0.5656854249;

            // Test that the skewness is equal to the hand calculated value
            TS_ASSERT_DELTA(mesh.GetSkewnessOfElementMassDistributionAboutAxis(0, axis) - hand_calculated_skewness, 0.0, 1e-9);

            // If we flip the axis, the skewness should be minus what it was before
            axis[1] = -1.0;
            TS_ASSERT_DELTA(mesh.GetSkewnessOfElementMassDistributionAboutAxis(0, axis) + hand_calculated_skewness, 0.0, 1e-9);
        }
    }

    void TestReMesh()
    {
        /*
         * In this test, we generate a mesh with multiple elements and lamina, in which nodes are already evenly spaced.
         * We simply check that none of the node locations are altered in a ReMesh.
         *
         * The specific ReMeshElements and ReMeshLaminas methods are tested separately.
         */

        std::vector<Node<2>*> nodes;

        std::vector<Node<2>*> nodes_elem1;
        std::vector<Node<2>*> nodes_elem2;
        std::vector<Node<2>*> nodes_elem3;

        std::vector<Node<2>*> nodes_lam1;
        std::vector<Node<2>*> nodes_lam2;

        nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, true, 0.1, 0.0));
        nodes.push_back(new Node<2>(2, true, 0.05, 0.05 * sqrt(3)));

        nodes.push_back(new Node<2>(3, true, 0.0, 0.0));
        nodes.push_back(new Node<2>(4, true, 0.1, 0.0));
        nodes.push_back(new Node<2>(5, true, 0.1, 0.1));
        nodes.push_back(new Node<2>(6, true, 0.0, 0.1));

        nodes.push_back(new Node<2>(7, true, 0.5 + 0.1 * cos(0.0 * M_PI / 8.0), 0.5 + 0.1 * sin(0.0 * M_PI / 8.0)));
        nodes.push_back(new Node<2>(8, true, 0.5 + 0.1 * cos(2.0 * M_PI / 8.0), 0.5 + 0.1 * sin(2.0 * M_PI / 8.0)));
        nodes.push_back(new Node<2>(9, true, 0.5 + 0.1 * cos(4.0 * M_PI / 8.0), 0.5 + 0.1 * sin(4.0 * M_PI / 8.0)));
        nodes.push_back(new Node<2>(10, true, 0.5 + 0.1 * cos(6.0 * M_PI / 8.0), 0.5 + 0.1 * sin(6.0 * M_PI / 8.0)));
        nodes.push_back(new Node<2>(11, true, 0.5 + 0.1 * cos(8.0 * M_PI / 8.0), 0.5 + 0.1 * sin(8.0 * M_PI / 8.0)));
        nodes.push_back(new Node<2>(12, true, 0.5 + 0.1 * cos(10.0 * M_PI / 8.0), 0.5 + 0.1 * sin(10.0 * M_PI / 8.0)));
        nodes.push_back(new Node<2>(13, true, 0.5 + 0.1 * cos(12.0 * M_PI / 8.0), 0.5 + 0.1 * sin(12.0 * M_PI / 8.0)));
        nodes.push_back(new Node<2>(14, true, 0.5 + 0.1 * cos(14.0 * M_PI / 8.0), 0.5 + 0.1 * sin(14.0 * M_PI / 8.0)));

        nodes.push_back(new Node<2>(15, true, 0.1, 0.3));
        nodes.push_back(new Node<2>(16, true, 0.3, 0.3));
        nodes.push_back(new Node<2>(17, true, 0.5, 0.3));
        nodes.push_back(new Node<2>(18, true, 0.7, 0.3));
        nodes.push_back(new Node<2>(19, true, 0.9, 0.3));

        nodes.push_back(new Node<2>(20, true, 0.2, 0.2));
        nodes.push_back(new Node<2>(21, true, 0.7, 0.7));

        // Triangle
        nodes_elem1.push_back(nodes[0]);
        nodes_elem1.push_back(nodes[1]);
        nodes_elem1.push_back(nodes[2]);

        // Square
        nodes_elem2.push_back(nodes[3]);
        nodes_elem2.push_back(nodes[4]);
        nodes_elem2.push_back(nodes[5]);
        nodes_elem2.push_back(nodes[6]);

        // Octagon
        nodes_elem3.push_back(nodes[7]);
        nodes_elem3.push_back(nodes[8]);
        nodes_elem3.push_back(nodes[9]);
        nodes_elem3.push_back(nodes[10]);
        nodes_elem3.push_back(nodes[11]);
        nodes_elem3.push_back(nodes[12]);
        nodes_elem3.push_back(nodes[13]);
        nodes_elem3.push_back(nodes[14]);

        // Lam 1 (x varies)
        nodes_lam1.push_back(nodes[15]);
        nodes_lam1.push_back(nodes[16]);
        nodes_lam1.push_back(nodes[17]);
        nodes_lam1.push_back(nodes[18]);
        nodes_lam1.push_back(nodes[19]);

        // Lam 2 (y varies)
        nodes_lam2.push_back(nodes[20]);
        nodes_lam2.push_back(nodes[21]);

        std::vector<ImmersedBoundaryElement<2, 2>*> elems;
        elems.push_back(new ImmersedBoundaryElement<2, 2>(0, nodes_elem1));
        elems.push_back(new ImmersedBoundaryElement<2, 2>(1, nodes_elem2));
        elems.push_back(new ImmersedBoundaryElement<2, 2>(2, nodes_elem3));

        std::vector<ImmersedBoundaryElement<1, 2>*> lams;
        lams.push_back(new ImmersedBoundaryElement<1, 2>(0, nodes_lam1));
        lams.push_back(new ImmersedBoundaryElement<1, 2>(1, nodes_lam2));
        
        ImmersedBoundaryMesh<2,2> mesh(nodes, elems, lams);

        // Make a copy of the node locations prior to ReMesh
        std::vector<c_vector<double, 2> > old_locations;
        for (unsigned node_idx = 0; node_idx < mesh.GetNumNodes(); node_idx++)
        {
            old_locations.push_back(mesh.GetNode(node_idx)->rGetLocation());
        }

        // Second ReMesh - node locations should remain unchanged
        mesh.ReMesh();

        // Verify everything's still the same
        for (unsigned node_idx = 0; node_idx < mesh.GetNumNodes(); node_idx++)
        {
            TS_ASSERT_DELTA(old_locations[node_idx][0], mesh.GetNode(node_idx)->rGetLocation()[0], 1e-12);
            TS_ASSERT_DELTA(old_locations[node_idx][1], mesh.GetNode(node_idx)->rGetLocation()[1], 1e-12);
        }
    }

    void TestReMeshElement()
    {
        // ReMeshElement where nothing should change
        {
            std::vector<Node<2>*> nodes;
            nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
            nodes.push_back(new Node<2>(1, true, 0.1, 0.0));
            nodes.push_back(new Node<2>(2, true, 0.1, 0.1));
            nodes.push_back(new Node<2>(3, true, 0.0, 0.1));

            std::vector<ImmersedBoundaryElement<2, 2>*> elems;
            elems.push_back(new ImmersedBoundaryElement<2, 2>(0, nodes));

            ImmersedBoundaryMesh<2,2> mesh(nodes, elems);

            unsigned num_nodes = mesh.GetElement(0)->GetNumNodes();

            // Get locations before ReMesh
            std::vector<double> old_pos;
            for (unsigned node_idx = 0; node_idx < num_nodes; node_idx++)
            {
                old_pos.push_back(norm_2(mesh.GetElement(0)->GetNode(node_idx)->rGetLocation()));
            }

            // Remesh fixing the first location
            mesh.ReMesh(false);

            for (unsigned node_idx = 0; node_idx < num_nodes; node_idx++)
            {
                TS_ASSERT_DELTA(old_pos[node_idx], norm_2(mesh.GetElement(0)->GetNode(node_idx)->rGetLocation()), 1e-12);
            }

            // Remesh starting from a random location
            mesh.ReMesh(true);

            for (unsigned node_idx = 0; node_idx < num_nodes; node_idx++)
            {
                TS_ASSERT_DELTA(old_pos[node_idx], norm_2(mesh.GetElement(0)->GetNode(node_idx)->rGetLocation()), 1e-12);
            }
        }

        // ReMeshElement where nothing should change, coping with an overlap
        {
            std::vector<Node<2>*> nodes;
            nodes.push_back(new Node<2>(0, true, 0.95, 0.1));
            nodes.push_back(new Node<2>(1, true, 0.05, 0.1));
            nodes.push_back(new Node<2>(2, true, 0.05, 0.2));
            nodes.push_back(new Node<2>(3, true, 0.95, 0.2));

            std::vector<ImmersedBoundaryElement<2, 2>*> elems;
            elems.push_back(new ImmersedBoundaryElement<2, 2>(0, nodes));

            ImmersedBoundaryMesh<2,2> mesh(nodes, elems);

            unsigned num_nodes = mesh.GetElement(0)->GetNumNodes();

            // Get locations before ReMesh
            std::vector<c_vector<double, 2>> old_pos;
            for (unsigned node_idx = 0; node_idx < num_nodes; node_idx++)
            {
                old_pos.push_back(mesh.GetElement(0)->GetNode(node_idx)->rGetLocation());
            }

            // Remesh with random location
            mesh.ReMesh(true);

            for (unsigned node_idx = 0; node_idx < num_nodes; node_idx++)
            {
                TS_ASSERT_DELTA(old_pos[node_idx][0], mesh.GetElement(0)->GetNode(node_idx)->rGetLocation()[0], 1e-12);
                TS_ASSERT_DELTA(old_pos[node_idx][1], mesh.GetElement(0)->GetNode(node_idx)->rGetLocation()[1], 1e-12);
            }
        }

        // ReMeshElement where positions should change
        {
            std::vector<Node<2>*> nodes;
            nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
            nodes.push_back(new Node<2>(1, true, 0.1, 0.0));
            nodes.push_back(new Node<2>(2, true, 0.1, 0.1));

            std::vector<ImmersedBoundaryElement<2, 2>*> elems;
            elems.push_back(new ImmersedBoundaryElement<2, 2>(0, nodes));

            ImmersedBoundaryMesh<2,2> mesh(nodes, elems);

            double cumulative_dist = 0.2 + sqrt(0.02);
            double node_spacing = cumulative_dist / 3.0;
            double epsilon = node_spacing - 0.1;

            // Remesh fixing the first location
            mesh.ReMesh(false);

            // First node should not move
            TS_ASSERT_DELTA(mesh.GetElement(0)->GetNode(0)->rGetLocation()[0], 0.0, 1e-12);
            TS_ASSERT_DELTA(mesh.GetElement(0)->GetNode(0)->rGetLocation()[1], 0.0, 1e-12);

            // Second node up back edge slightly
            TS_ASSERT_DELTA(mesh.GetElement(0)->GetNode(1)->rGetLocation()[0], 0.1, 1e-12);
            TS_ASSERT_DELTA(mesh.GetElement(0)->GetNode(1)->rGetLocation()[1], epsilon, 1e-12);

            // Third node down diagonal towards (0,0)
            double down_diagonal = 0.1 - 0.2 * epsilon / sqrt(0.02);
            TS_ASSERT_DELTA(mesh.GetElement(0)->GetNode(2)->rGetLocation()[0], down_diagonal, 1e-12);
            TS_ASSERT_DELTA(mesh.GetElement(0)->GetNode(2)->rGetLocation()[1], down_diagonal, 1e-12);
        }

        // ReMeshElement where positions should change
        {
            std::vector<Node<2>*> nodes;
            nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
            nodes.push_back(new Node<2>(1, true, 0.1, 0.0));
            nodes.push_back(new Node<2>(2, true, 0.1, 0.1));

            std::vector<ImmersedBoundaryElement<2, 2>*> elems;
            elems.push_back(new ImmersedBoundaryElement<2, 2>(0, nodes));

            ImmersedBoundaryMesh<2,2> mesh(nodes, elems);

            double cumulative_dist = 0.2 + sqrt(0.02);
            double node_spacing = cumulative_dist / 3.0;

            // Remesh using a random location (fixes the second index)
            RandomNumberGenerator::Instance()->Reseed(0u);
            mesh.ReMesh(true);

            // Third node should not move
            TS_ASSERT_DELTA(mesh.GetElement(0)->GetNode(2)->rGetLocation()[0], 0.1, 1e-12);
            TS_ASSERT_DELTA(mesh.GetElement(0)->GetNode(2)->rGetLocation()[1], 0.1, 1e-12);

            // First node up diagonal slightly
            double up_diagonal = (sqrt(0.02) - node_spacing) / sqrt(2.0);
            TS_ASSERT_DELTA(mesh.GetElement(0)->GetNode(0)->rGetLocation()[0], up_diagonal, 1e-12);
            TS_ASSERT_DELTA(mesh.GetElement(0)->GetNode(0)->rGetLocation()[1], up_diagonal, 1e-12);

            // Second node should be along the first edge, near (0.1, 0)
            double along_edge = 2.0 * node_spacing - sqrt(0.02);
            TS_ASSERT_DELTA(mesh.GetElement(0)->GetNode(1)->rGetLocation()[0], along_edge, 1e-12);
            TS_ASSERT_DELTA(mesh.GetElement(0)->GetNode(1)->rGetLocation()[1], 0.0, 1e-12);
        }
    }

    void TestReMeshLamina() throw(Exception)
    {
        // ReMeshLamina where nothing should change
        {
            std::vector<Node<2>*> nodes;
            nodes.push_back(new Node<2>(0, true, 0.1, 0.3));
            nodes.push_back(new Node<2>(1, true, 0.3, 0.3));
            nodes.push_back(new Node<2>(2, true, 0.5, 0.3));
            nodes.push_back(new Node<2>(3, true, 0.7, 0.3));
            nodes.push_back(new Node<2>(4, true, 0.9, 0.3));

            std::vector<ImmersedBoundaryElement<2, 2>*> elems; // empty, for constructor
            std::vector<ImmersedBoundaryElement<1, 2>*> lams;
            lams.push_back(new ImmersedBoundaryElement<1, 2>(0, nodes));

            ImmersedBoundaryMesh<2,2> mesh(nodes, elems, lams);

            unsigned num_nodes = mesh.GetLamina(0)->GetNumNodes();

            // Get locations before ReMesh
            std::vector<c_vector<double, 2>> old_pos;
            for (unsigned node_idx = 0; node_idx < num_nodes; node_idx++)
            {
                old_pos.push_back(mesh.GetLamina(0)->GetNode(node_idx)->rGetLocation());
            }

            // Should make no different at all which node is selected
            mesh.ReMesh(true);

            for (unsigned node_idx = 0; node_idx < num_nodes; node_idx++)
            {
                TS_ASSERT_DELTA(old_pos[node_idx][0], mesh.GetLamina(0)->GetNode(node_idx)->rGetLocation()[0], 1e-12);
                TS_ASSERT_DELTA(old_pos[node_idx][1], mesh.GetLamina(0)->GetNode(node_idx)->rGetLocation()[1], 1e-12);
            }
        }

        // ReMeshLamina where the points should become evenly spaced
        {
            std::vector<Node<2>*> nodes;
            nodes.push_back(new Node<2>(0, true, 0.1, 0.3));
            nodes.push_back(new Node<2>(1, true, 0.32, 0.3));
            nodes.push_back(new Node<2>(2, true, 0.53, 0.3));
            nodes.push_back(new Node<2>(3, true, 0.64, 0.3));
            nodes.push_back(new Node<2>(4, true, 0.88, 0.3));

            std::vector<ImmersedBoundaryElement<2, 2>*> elems; // empty, for constructor
            std::vector<ImmersedBoundaryElement<1, 2>*> lams;
            lams.push_back(new ImmersedBoundaryElement<1, 2>(0, nodes));

            ImmersedBoundaryMesh<2,2> mesh(nodes, elems, lams);

            mesh.ReMesh(false);

            TS_ASSERT_DELTA(mesh.GetLamina(0)->GetNode(0)->rGetLocation()[0], 0.1, 1e-12);
            TS_ASSERT_DELTA(mesh.GetLamina(0)->GetNode(0)->rGetLocation()[1], 0.3, 1e-12);

            TS_ASSERT_DELTA(mesh.GetLamina(0)->GetNode(1)->rGetLocation()[0], 0.3, 1e-12);
            TS_ASSERT_DELTA(mesh.GetLamina(0)->GetNode(1)->rGetLocation()[1], 0.3, 1e-12);

            TS_ASSERT_DELTA(mesh.GetLamina(0)->GetNode(2)->rGetLocation()[0], 0.5, 1e-12);
            TS_ASSERT_DELTA(mesh.GetLamina(0)->GetNode(2)->rGetLocation()[1], 0.3, 1e-12);

            TS_ASSERT_DELTA(mesh.GetLamina(0)->GetNode(3)->rGetLocation()[0], 0.7, 1e-12);
            TS_ASSERT_DELTA(mesh.GetLamina(0)->GetNode(3)->rGetLocation()[1], 0.3, 1e-12);

            TS_ASSERT_DELTA(mesh.GetLamina(0)->GetNode(4)->rGetLocation()[0], 0.9, 1e-12);
            TS_ASSERT_DELTA(mesh.GetLamina(0)->GetNode(4)->rGetLocation()[1], 0.3, 1e-12);
        }
    }

    void TestNodesInDifferentElementOrLamina()
    {
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, true, 0.1, 0.1));
        nodes.push_back(new Node<2>(1, true, 0.1, 0.3));
        nodes.push_back(new Node<2>(2, true, 0.1, 0.5));

        nodes.push_back(new Node<2>(3, true, 0.1, 0.1));
        nodes.push_back(new Node<2>(4, true, 0.3, 0.1));
        nodes.push_back(new Node<2>(5, true, 0.5, 0.1));

        nodes.push_back(new Node<2>(6, true, 0.1, 0.1));
        nodes.push_back(new Node<2>(7, true, 0.2, 0.1));
        nodes.push_back(new Node<2>(8, true, 0.2, 0.2));

        nodes.push_back(new Node<2>(9, true, 0.1, 0.1));
        nodes.push_back(new Node<2>(10, true, 0.2, 0.1));
        nodes.push_back(new Node<2>(11, true, 0.2, 0.2));

        for (unsigned i = 0; i < 5; i++)
        {
            nodes[i]->SetRegion(LAMINA_REGION);
        }

        std::vector<Node<2>*> nodes_lam1;
        nodes_lam1.push_back(nodes[0]);
        nodes_lam1.push_back(nodes[1]);
        nodes_lam1.push_back(nodes[2]);

        std::vector<Node<2>*> nodes_lam2;
        nodes_lam2.push_back(nodes[3]);
        nodes_lam2.push_back(nodes[4]);
        nodes_lam2.push_back(nodes[5]);

        std::vector<Node<2>*> nodes_elem1;
        nodes_elem1.push_back(nodes[6]);
        nodes_elem1.push_back(nodes[7]);
        nodes_elem1.push_back(nodes[8]);

        std::vector<Node<2>*> nodes_elem2;
        nodes_elem2.push_back(nodes[9]);
        nodes_elem2.push_back(nodes[10]);
        nodes_elem2.push_back(nodes[11]);


        std::vector<ImmersedBoundaryElement<1, 2>*> lams;
        lams.push_back(new ImmersedBoundaryElement<1, 2>(0, nodes_lam1));
        lams.push_back(new ImmersedBoundaryElement<1, 2>(1, nodes_lam2));

        std::vector<ImmersedBoundaryElement<2, 2>*> elems;
        elems.push_back(new ImmersedBoundaryElement<2, 2>(0, nodes_elem1));
        elems.push_back(new ImmersedBoundaryElement<2, 2>(1, nodes_elem2));

        ImmersedBoundaryMesh<2,2> mesh(nodes, elems, lams);

        // Lamina 0 and Elem 0
        TS_ASSERT(mesh.NodesInDifferentElementOrLamina(mesh.GetNode(0), mesh.GetNode(6)));

        // Lamina 0 and Elem 1
        TS_ASSERT(mesh.NodesInDifferentElementOrLamina(mesh.GetNode(0), mesh.GetNode(9)));

        // Lamina 1 and Elem 0
        TS_ASSERT(mesh.NodesInDifferentElementOrLamina(mesh.GetNode(3), mesh.GetNode(6)));

        // Lamina 1 and Elem 1
        TS_ASSERT(mesh.NodesInDifferentElementOrLamina(mesh.GetNode(3), mesh.GetNode(9)));

        // Lamina 0 and Lamina 0
        TS_ASSERT(mesh.NodesInDifferentElementOrLamina(mesh.GetNode(0), mesh.GetNode(1)));

        // Lamina 1 and Lamina 1
        TS_ASSERT(mesh.NodesInDifferentElementOrLamina(mesh.GetNode(3), mesh.GetNode(4)));

        // Elem 0 and Elem 0
        TS_ASSERT(!mesh.NodesInDifferentElementOrLamina(mesh.GetNode(6), mesh.GetNode(7)));

        // Elem 1 and Elem 1
        TS_ASSERT(!mesh.NodesInDifferentElementOrLamina(mesh.GetNode(9), mesh.GetNode(10)));
    }

    void TestGeometricMethods()
    {
        // Make six nodes
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, true, 0.1, 0.1));
        nodes.push_back(new Node<2>(1, true, 0.2, 0.1));
        nodes.push_back(new Node<2>(2, true, 0.3, 0.2));
        nodes.push_back(new Node<2>(3, true, 0.3, 0.3));
        nodes.push_back(new Node<2>(4, true, 0.1, 0.2));

        // Make one element out of these nodes
        std::vector<Node<2>*> nodes_elem;
        for (unsigned i=0; i<5; i++)
        {
            nodes_elem.push_back(nodes[i]);
        }

        std::vector<ImmersedBoundaryElement<2,2>*> elements;
        elements.push_back(new ImmersedBoundaryElement<2,2>(0, nodes_elem));

        // Make a mesh
        ImmersedBoundaryMesh<2,2> mesh(nodes, elements);

        TS_ASSERT_EQUALS(mesh.GetNumElements(), 1u);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), 5u);

        // Test that the centroid, moments and short axis of the element are calculated correctly
        // (i.e. agreee with Matlab and pen-and-paper calculations)
        c_vector<double,2> centroid = mesh.GetCentroidOfElement(0);
        TS_ASSERT_DELTA(centroid[0], 0.2000, 1e-4);
        TS_ASSERT_DELTA(centroid[1], 0.1866, 1e-4);

        c_vector<double,3> moments = mesh.CalculateMomentsOfElement(0);
        TS_ASSERT_DELTA(moments[0], 5.388e-5, 1e-8);
        TS_ASSERT_DELTA(moments[1], 7.500e-5, 1e-8);
        TS_ASSERT_DELTA(moments[2], 3.750e-5, 1e-8);

        c_vector<double,2> short_axis = mesh.GetShortAxisOfElement(0);
        TS_ASSERT_DELTA(short_axis[0],  0.6037, 1e-4);
        TS_ASSERT_DELTA(short_axis[1], -0.7971, 1e-4);
    }

    void TestGetNeighbouringElementIndices()
    {
        /* This 3x3 honeycomb will look like:
         *    5
         * 2     8
         *    4
         * 1     7
         *    3
         * 0     6
         */
        ImmersedBoundaryHoneycombMeshGenerator gen(3u, 3u, 5u, 0.05, 0.2);

        auto p_mesh = gen.GetMesh();
        p_mesh->SetNeighbourDist(0.1);

        TS_ASSERT_EQUALS(p_mesh->GetNeighbouringElementIndices(0u), std::set<unsigned>({1u, 3u}));
        TS_ASSERT_EQUALS(p_mesh->GetNeighbouringElementIndices(1u), std::set<unsigned>({0u, 2u, 3u, 4u}));
        TS_ASSERT_EQUALS(p_mesh->GetNeighbouringElementIndices(2u), std::set<unsigned>({1u, 4u, 5u}));
        TS_ASSERT_EQUALS(p_mesh->GetNeighbouringElementIndices(3u), std::set<unsigned>({0u, 1u, 4u, 6u, 7u}));
        TS_ASSERT_EQUALS(p_mesh->GetNeighbouringElementIndices(4u), std::set<unsigned>({1u, 2u, 3u, 5u, 7u, 8u}));
        TS_ASSERT_EQUALS(p_mesh->GetNeighbouringElementIndices(5u), std::set<unsigned>({2u, 4u, 8u}));
        TS_ASSERT_EQUALS(p_mesh->GetNeighbouringElementIndices(6u), std::set<unsigned>({3u, 7u}));
        TS_ASSERT_EQUALS(p_mesh->GetNeighbouringElementIndices(7u), std::set<unsigned>({3u, 4u, 6u, 8u}));
        TS_ASSERT_EQUALS(p_mesh->GetNeighbouringElementIndices(8u), std::set<unsigned>({4u, 5u, 7u}));
    }

    void TestGetPolygonDistribution()
    {
        // A 3x3 honeycomb will just have a single non-boundary element
        ImmersedBoundaryHoneycombMeshGenerator gen(3u, 3u, 5u, 0.05, 0.2);

        auto p_mesh = gen.GetMesh();
        p_mesh->SetNeighbourDist(0.1);

        std::array<unsigned, 13> polygon_dist = p_mesh->GetPolygonDistribution();
        TS_ASSERT_EQUALS(polygon_dist[0], 0u);
        TS_ASSERT_EQUALS(polygon_dist[1], 0u);
        TS_ASSERT_EQUALS(polygon_dist[2], 0u);
        TS_ASSERT_EQUALS(polygon_dist[3], 0u);
        TS_ASSERT_EQUALS(polygon_dist[4], 0u);
        TS_ASSERT_EQUALS(polygon_dist[5], 0u);
        TS_ASSERT_EQUALS(polygon_dist[6], 1u);
        TS_ASSERT_EQUALS(polygon_dist[7], 0u);
        TS_ASSERT_EQUALS(polygon_dist[8], 0u);
        TS_ASSERT_EQUALS(polygon_dist[9], 0u);
        TS_ASSERT_EQUALS(polygon_dist[10], 0u);
        TS_ASSERT_EQUALS(polygon_dist[11], 0u);
        TS_ASSERT_EQUALS(polygon_dist[12], 0u);

        // If we pretend none of the elements are on the boundary, the answer should changes
        for (unsigned elem_idx = 0; elem_idx < p_mesh->GetNumElements(); ++elem_idx)
        {
            p_mesh->GetElement(elem_idx)->SetIsBoundaryElement(false);
        }

        std::array<unsigned, 13> new_polygon_dist = p_mesh->GetPolygonDistribution();
        TS_ASSERT_EQUALS(new_polygon_dist[0], 0u);
        TS_ASSERT_EQUALS(new_polygon_dist[1], 0u);
        TS_ASSERT_EQUALS(new_polygon_dist[2], 2u);
        TS_ASSERT_EQUALS(new_polygon_dist[3], 3u);
        TS_ASSERT_EQUALS(new_polygon_dist[4], 2u);
        TS_ASSERT_EQUALS(new_polygon_dist[5], 1u);
        TS_ASSERT_EQUALS(new_polygon_dist[6], 1u);
        TS_ASSERT_EQUALS(new_polygon_dist[7], 0u);
        TS_ASSERT_EQUALS(new_polygon_dist[8], 0u);
        TS_ASSERT_EQUALS(new_polygon_dist[9], 0u);
        TS_ASSERT_EQUALS(new_polygon_dist[10], 0u);
        TS_ASSERT_EQUALS(new_polygon_dist[11], 0u);
        TS_ASSERT_EQUALS(new_polygon_dist[12], 0u);
    }

    void TestGetMaxIndexMethods()
    {
        // Get default mesh with no nodes, elements, or laminas
        {
            ImmersedBoundaryMesh<2, 2> ib_mesh;

            TS_ASSERT_EQUALS(ib_mesh.GetMaxNodeIndex(), UINT_MAX);
            TS_ASSERT_EQUALS(ib_mesh.GetMaxElementIndex(), UINT_MAX);
            TS_ASSERT_EQUALS(ib_mesh.GetMaxLaminaIndex(), UINT_MAX);
        }

        // Get real mesh with nodes, elements, and laminas
        {
            // Make a few nodes
            std::vector<Node<2>*> nodes;
            nodes.push_back(new Node<2>(0, true, 0.1, 0.1));
            nodes.push_back(new Node<2>(1, true, 0.2, 0.1));
            nodes.push_back(new Node<2>(2, true, 0.3, 0.2));

            // Make two elements
            std::vector<ImmersedBoundaryElement<2,2>*> elements;
            elements.push_back(new ImmersedBoundaryElement<2,2>(0, nodes));
            elements.push_back(new ImmersedBoundaryElement<2,2>(1, nodes));

            // Make four laminas
            std::vector<ImmersedBoundaryElement<1,2>*> lams;
            lams.push_back(new ImmersedBoundaryElement<1,2>(0, nodes));
            lams.push_back(new ImmersedBoundaryElement<1,2>(1, nodes));
            lams.push_back(new ImmersedBoundaryElement<1,2>(2, nodes));
            lams.push_back(new ImmersedBoundaryElement<1,2>(3, nodes));

            // Make a mesh
            ImmersedBoundaryMesh<2,2> ib_mesh(nodes, elements, lams);

            TS_ASSERT_EQUALS(ib_mesh.GetMaxNodeIndex(), 2);
            TS_ASSERT_EQUALS(ib_mesh.GetMaxElementIndex(), 1);
            TS_ASSERT_EQUALS(ib_mesh.GetMaxLaminaIndex(), 3);
        }
    }


};

#endif /*TESTIMMERSEDBOUNDARYMESH_HPP_*/
