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

#ifndef TESTSEMREGIONALFORCE_HPP_
#define TESTSEMREGIONALFORCE_HPP_

#include <cxxtest/TestSuite.h>

#include "CheckpointArchiveTypes.hpp"

#include <fstream>
#include <string>
#include <vector>

#include "OutputFileHandler.hpp"
#include "SemRegionalForce.hpp"
#include "SemBasedCellPopulation.hpp"
#include "SemMesh.hpp"
#include "SemElement.hpp"
#include "NodesOnlyMesh.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "CellsGenerator.hpp"
#include "NoCellCycleModel.hpp"
#include "AbstractCellBasedTestSuite.hpp"

// This test is never run in parallel
#include "FakePetscSetup.hpp"


class TestSemRegionalForce : public AbstractCellBasedTestSuite
{
private:

    /** Build a c_vector<double,2> from two coordinates. */
    static c_vector<double, 2> Vec(double x, double y)
    {
        c_vector<double, 2> v;
        v[0] = x;
        v[1] = y;
        return v;
    }

    /**
     * A minimal SEM population: one SemElement holding all supplied nodes (one cell), with each
     * node's region set. The SemMesh is held by value, so the fixture must outlive any population
     * built on it.
     */
    struct SingleElementFixture
    {
        std::vector<Node<2>*> nodes;
        std::vector<SemElement<2>*> elements;
        SemMesh<2> mesh;

        SingleElementFixture(const std::vector<c_vector<double, 2> >& rLocations,
                             const std::vector<unsigned>& rRegions)
            : nodes(MakeNodes(rLocations, rRegions)),
              elements(MakeElements(nodes)),
              mesh(nodes, elements)
        {
        }

        static std::vector<Node<2>*> MakeNodes(const std::vector<c_vector<double, 2> >& rLocations,
                                               const std::vector<unsigned>& rRegions)
        {
            std::vector<Node<2>*> nodes;
            for (unsigned i = 0; i < rLocations.size(); ++i)
            {
                Node<2>* p_node = new Node<2>(i, rLocations[i], false);
                p_node->SetRegion(rRegions[i]);
                nodes.push_back(p_node);
            }
            return nodes;
        }

        static std::vector<SemElement<2>*> MakeElements(std::vector<Node<2>*>& rNodes)
        {
            std::vector<SemElement<2>*> elements;
            elements.push_back(new SemElement<2>(0u, rNodes));
            return elements;
        }
    };

    static std::vector<CellPtr> CreateCells(unsigned numCells)
    {
        std::vector<CellPtr> cells;
        CellsGenerator<NoCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, numCells);
        return cells;
    }

public:

    void TestForceMatchesHarmonicLawWithDefaults()
    {
        // Two nodes 0.3 apart, in regions 0 (interior) and 1 (boundary)
        SingleElementFixture fixture({Vec(0.0, 0.0), Vec(0.0, 0.3)}, {0u, 1u});
        std::vector<CellPtr> cells = CreateCells(1u);
        std::vector<unsigned> location_indices = {0u};
        SemBasedCellPopulation<2> population(fixture.mesh, cells, false, true, location_indices);

        SemRegionalForce<2> force;

        // Defaults: spring constants {1.0, 2.0}, rest lengths {0.2, 0.15}
        //   k = 0.5*(1.0 + 2.0) = 1.5, L = 0.5*(0.2 + 0.15) = 0.175, r = 0.3
        //   coeff = k*(1 - L/r) = 1.5*(1 - 0.175/0.3) = 0.625
        //   F_A = coeff * (r_B - r_A) = 0.625 * (0, 0.3) = (0, 0.1875)
        c_vector<double, 2> force_on_a = force.CalculateForceBetweenNodes(0u, 1u, population);
        TS_ASSERT_DELTA(force_on_a[0], 0.0, 1e-9);
        TS_ASSERT_DELTA(force_on_a[1], 0.1875, 1e-9);
    }

    void TestArbitraryNumberOfRegions()
    {
        // Three regions; interacting nodes in regions 0 and 2
        SingleElementFixture fixture({Vec(0.0, 0.0), Vec(0.0, 0.3)}, {0u, 2u});
        std::vector<CellPtr> cells = CreateCells(1u);
        std::vector<unsigned> location_indices = {0u};
        SemBasedCellPopulation<2> population(fixture.mesh, cells, false, true, location_indices);

        SemRegionalForce<2> force;
        force.SetSpringConstants({1.0, 2.0, 3.0});
        force.SetRestLengths({0.2, 0.15, 0.1});

        //   k = 0.5*(1.0 + 3.0) = 2.0, L = 0.5*(0.2 + 0.1) = 0.15, r = 0.3
        //   coeff = 2.0*(1 - 0.15/0.3) = 1.0, F_A = 1.0 * (0, 0.3) = (0, 0.3)
        c_vector<double, 2> force_on_a = force.CalculateForceBetweenNodes(0u, 1u, population);
        TS_ASSERT_DELTA(force_on_a[0], 0.0, 1e-9);
        TS_ASSERT_DELTA(force_on_a[1], 0.3, 1e-9);
    }

    void TestCutOffDistance()
    {
        // Nodes 0.6 apart, beyond the default 0.5 cut-off
        SingleElementFixture fixture({Vec(0.0, 0.0), Vec(0.0, 0.6)}, {0u, 0u});
        std::vector<CellPtr> cells = CreateCells(1u);
        std::vector<unsigned> location_indices = {0u};
        SemBasedCellPopulation<2> population(fixture.mesh, cells, false, true, location_indices);

        SemRegionalForce<2> force;
        c_vector<double, 2> force_beyond_cutoff = force.CalculateForceBetweenNodes(0u, 1u, population);
        TS_ASSERT_DELTA(force_beyond_cutoff[0], 0.0, 1e-12);
        TS_ASSERT_DELTA(force_beyond_cutoff[1], 0.0, 1e-12);

        // Raise the cut-off so the pair now interacts; both nodes region 0 => k=1.0, L=0.2, r=0.6
        //   coeff = 1.0*(1 - 0.2/0.6) = 2/3, F_A = (0, (2/3)*0.6) = (0, 0.4)
        force.SetCutOffDistance(0.7);
        TS_ASSERT_DELTA(force.GetCutOffDistance(), 0.7, 1e-12);
        c_vector<double, 2> force_within_cutoff = force.CalculateForceBetweenNodes(0u, 1u, population);
        TS_ASSERT_DELTA(force_within_cutoff[1], 0.4, 1e-9);
    }

    void TestCoincidentNodesReturnZero()
    {
        SingleElementFixture fixture({Vec(0.0, 0.0), Vec(0.0, 0.0)}, {0u, 0u});
        std::vector<CellPtr> cells = CreateCells(1u);
        std::vector<unsigned> location_indices = {0u};
        SemBasedCellPopulation<2> population(fixture.mesh, cells, false, true, location_indices);

        SemRegionalForce<2> force;
        c_vector<double, 2> force_on_a = force.CalculateForceBetweenNodes(0u, 1u, population);
        TS_ASSERT_DELTA(force_on_a[0], 0.0, 1e-12);
        TS_ASSERT_DELTA(force_on_a[1], 0.0, 1e-12);
    }

    void TestForceIsAntisymmetric()
    {
        SingleElementFixture fixture({Vec(0.1, 0.2), Vec(0.3, 0.45)}, {0u, 1u});
        std::vector<CellPtr> cells = CreateCells(1u);
        std::vector<unsigned> location_indices = {0u};
        SemBasedCellPopulation<2> population(fixture.mesh, cells, false, true, location_indices);

        SemRegionalForce<2> force;
        c_vector<double, 2> force_on_a = force.CalculateForceBetweenNodes(0u, 1u, population);
        c_vector<double, 2> force_on_b = force.CalculateForceBetweenNodes(1u, 0u, population);

        TS_ASSERT_DELTA(force_on_a[0], -force_on_b[0], 1e-12);
        TS_ASSERT_DELTA(force_on_a[1], -force_on_b[1], 1e-12);
    }

    void TestSettersAndGetters()
    {
        SemRegionalForce<2> force;

        force.SetSpringConstants({5.0, 6.0, 7.0});
        force.SetRestLengths({0.5, 0.4, 0.3});
        force.SetCutOffDistance(1.25);

        TS_ASSERT_EQUALS(force.GetSpringConstants().size(), 3u);
        TS_ASSERT_DELTA(force.GetSpringConstants()[0], 5.0, 1e-12);
        TS_ASSERT_DELTA(force.GetSpringConstants()[2], 7.0, 1e-12);
        TS_ASSERT_DELTA(force.GetRestLengths()[1], 0.4, 1e-12);
        TS_ASSERT_DELTA(force.GetCutOffDistance(), 1.25, 1e-12);
    }

    void TestArchiving()
    {
        EXIT_IF_PARALLEL;
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "SemRegionalForce.arch";

        {
            SemRegionalForce<2> force;
            force.SetSpringConstants({5.0, 6.0, 7.0});
            force.SetRestLengths({0.5, 0.4, 0.3});
            force.SetCutOffDistance(1.25);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Serialize via pointer to most abstract class possible
            AbstractForce<2>* const p_force = &force;
            output_arch << p_force;
        }

        {
            AbstractForce<2>* p_force;

            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> p_force;

            SemRegionalForce<2>* p_loaded = static_cast<SemRegionalForce<2>*>(p_force);
            TS_ASSERT_EQUALS(p_loaded->GetSpringConstants().size(), 3u);
            TS_ASSERT_DELTA(p_loaded->GetSpringConstants()[0], 5.0, 1e-12);
            TS_ASSERT_DELTA(p_loaded->GetSpringConstants()[2], 7.0, 1e-12);
            TS_ASSERT_EQUALS(p_loaded->GetRestLengths().size(), 3u);
            TS_ASSERT_DELTA(p_loaded->GetRestLengths()[1], 0.4, 1e-12);
            TS_ASSERT_DELTA(p_loaded->GetCutOffDistance(), 1.25, 1e-12);

            delete p_force;
        }
    }

    void TestExceptions()
    {
        // Mismatched parameter-array lengths are rejected at force-application time
        {
            SingleElementFixture fixture({Vec(0.0, 0.0), Vec(0.0, 0.3)}, {0u, 1u});
            std::vector<CellPtr> cells = CreateCells(1u);
            std::vector<unsigned> location_indices = {0u};
            SemBasedCellPopulation<2> population(fixture.mesh, cells, false, true, location_indices);

            SemRegionalForce<2> force;
            force.SetSpringConstants({1.0, 2.0, 3.0});
            force.SetRestLengths({0.2, 0.15});
            TS_ASSERT_THROWS_CONTAINS(force.AddForceContribution(population),
                                      "must be non-empty and of equal length");
        }

        // The force may only be used with a SemBasedCellPopulation
        {
            std::vector<Node<2>*> nodes;
            nodes.push_back(new Node<2>(0u, false, 0.0, 0.0));
            nodes.push_back(new Node<2>(1u, false, 0.0, 0.3));
            NodesOnlyMesh<2> mesh;
            mesh.ConstructNodesWithoutMesh(nodes, 1.5);
            std::vector<CellPtr> cells = CreateCells(mesh.GetNumNodes());
            NodeBasedCellPopulation<2> population(mesh, cells);

            SemRegionalForce<2> force;
            TS_ASSERT_THROWS_THIS(force.AddForceContribution(population),
                                  "SemRegionalForce is to be used with a SemBasedCellPopulation only");

            for (unsigned i = 0; i < nodes.size(); ++i)
            {
                delete nodes[i];
            }
        }
    }
};

#endif /*TESTSEMREGIONALFORCE_HPP_*/
