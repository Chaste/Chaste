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
#ifndef TESTPERIODICFORCES_HPP_
#define TESTPERIODICFORCES_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "CellsGenerator.hpp"
#include "OffLatticeSimulation.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "GeneralisedPeriodicLinearSpringForce.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "WildTypeCellMutationState.hpp"
#include "MutableMesh.hpp"
#include "SmartPointers.hpp"
#include "FileComparison.hpp"
//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

class TestPeriodicForces : public AbstractCellBasedTestSuite
{

private:

    MutableMesh<3,3>* Make3dMesh(unsigned width=3, unsigned height=3, unsigned depth=3)
    {
        /*          _ _ _ _ _
         *        /        /|
         *       /        / |
         *      /_ _ _ _ /  | depth
         *     |         |  |
         *     |         |  |
         *     |         |  /
         *     |         | / height
         *     |_ _ _ _ _|/
         *        width
         */

         MutableMesh<3,3>* p_mesh = new MutableMesh<3,3>;
         p_mesh->ConstructCuboid(width, height, depth);
         TrianglesMeshWriter<3,3> mesh_writer("", "3dSpringMesh");
         mesh_writer.WriteFilesUsingMesh(*p_mesh);

         return p_mesh;
    }

public:

    void TestSetDomainDimensions() throw (Exception)
    {
        GeneralisedPeriodicLinearSpringForce<2> linear_force2d;
        TS_ASSERT_DELTA(linear_force2d.GetPeriodicDomainWidth(), DOUBLE_UNSET, 1e-3);
        linear_force2d.SetPeriodicDomainWidth(1.576);
        TS_ASSERT_DELTA(linear_force2d.GetPeriodicDomainWidth(), 1.576, 1e-3);

        GeneralisedPeriodicLinearSpringForce<3> linear_force3d;
        TS_ASSERT_DELTA(linear_force3d.GetPeriodicDomainWidth(), DOUBLE_UNSET, 1e-3);
        TS_ASSERT_DELTA(linear_force3d.GetPeriodicDomainDepth(), DOUBLE_UNSET, 1e-3);
        linear_force3d.SetPeriodicDomainWidth(1.576);
        linear_force3d.SetPeriodicDomainDepth(1.865);
        TS_ASSERT_DELTA(linear_force3d.GetPeriodicDomainWidth(), 1.576, 1e-3);
        TS_ASSERT_DELTA(linear_force3d.GetPeriodicDomainDepth(), 1.865, 1e-3);
    }

    void TestPeriodicForceOnHoneycombMesh() throw (Exception)
    {
        EXIT_IF_PARALLEL; //HoneycombMeshGenerator doesn't work in parallel

        // Create a simple mesh
        HoneycombMeshGenerator generator(2, 2, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumNodes());

        // Create cell population
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Create force
        GeneralisedPeriodicLinearSpringForce<2> linear_force;

        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
            cell_population.GetNode(i)->ClearAppliedForce();
        }

        linear_force.AddForceContribution(cell_population);

        TS_ASSERT_DELTA(cell_population.GetNode(0)->rGetAppliedForce()[0], 7.5, 1e-3);
        TS_ASSERT_DELTA(cell_population.GetNode(0)->rGetAppliedForce()[1], -2.00962, 1e-3);
        TS_ASSERT_DELTA(cell_population.GetNode(1)->rGetAppliedForce()[0], -7.5, 1e-3);
        TS_ASSERT_DELTA(cell_population.GetNode(1)->rGetAppliedForce()[1], -2.88444e-15, 1e-3);
        TS_ASSERT_DELTA(cell_population.GetNode(2)->rGetAppliedForce()[0], 7.5, 1e-3);
        TS_ASSERT_DELTA(cell_population.GetNode(2)->rGetAppliedForce()[1], 2.88444e-15, 1e-3);
        TS_ASSERT_DELTA(cell_population.GetNode(3)->rGetAppliedForce()[0], -7.5, 1e-3);
        TS_ASSERT_DELTA(cell_population.GetNode(3)->rGetAppliedForce()[1], 2.00962, 1e-3);
    }

    void TestPeriodicForceOnNonUniformMesh() throw(Exception)
    {
        // Create 2D mesh
        std::vector<Node<2> *> nodes;
        nodes.push_back(new Node<2>(0, true, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, true, 4.0, 0.0));
        nodes.push_back(new Node<2>(2, true, 1.5, 0.5));
        nodes.push_back(new Node<2>(3, true, 3.0, 0.5));
        nodes.push_back(new Node<2>(4, false, 4.5, 0.5));
        nodes.push_back(new Node<2>(5, false, 1.0, 1.0));
        nodes.push_back(new Node<2>(6, false, 3.5, 1.0));
        MutableMesh<2,2> mesh(nodes);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Create cell population
        MeshBasedCellPopulation<2> cell_population(mesh, cells);

        // Create Voronoi tessellation
        cell_population.CreateVoronoiTessellation();

        // Create force
        GeneralisedPeriodicLinearSpringForce<2> linear_force;
        linear_force.SetPeriodicDomainWidth(4.5);

        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
            cell_population.GetNode(i)->ClearAppliedForce();
        }

        linear_force.AddForceContribution(cell_population);

        // Close enough
        TS_ASSERT_DELTA(cell_population.GetNode(0)->rGetAppliedForce()[0], 60.7698, 1e-3);
        TS_ASSERT_DELTA(cell_population.GetNode(0)->rGetAppliedForce()[1], -4.74342, 1e-3);
        TS_ASSERT_DELTA(cell_population.GetNode(1)->rGetAppliedForce()[0], -80.7733, 1e-3);
        TS_ASSERT_DELTA(cell_population.GetNode(1)->rGetAppliedForce()[1], 3.82704, 1e-3);
        TS_ASSERT_DELTA(cell_population.GetNode(2)->rGetAppliedForce()[0], 17.6281, 1e-3);
        TS_ASSERT_DELTA(cell_population.GetNode(2)->rGetAppliedForce()[1], -10.4214, 1e-3);
        TS_ASSERT_DELTA(cell_population.GetNode(3)->rGetAppliedForce()[0], -24.4709, 1e-3);
        TS_ASSERT_DELTA(cell_population.GetNode(3)->rGetAppliedForce()[1], -0.0364322, 1e-3);
        TS_ASSERT_DELTA(cell_population.GetNode(4)->rGetAppliedForce()[0], 10.6066, 1e-3);
        TS_ASSERT_DELTA(cell_population.GetNode(4)->rGetAppliedForce()[1], 12.1902, 1e-3);
        TS_ASSERT_DELTA(cell_population.GetNode(5)->rGetAppliedForce()[0], 18.2577, 1e-3);
        TS_ASSERT_DELTA(cell_population.GetNode(5)->rGetAppliedForce()[1], -1.54716, 1e-3);
        TS_ASSERT_DELTA(cell_population.GetNode(6)->rGetAppliedForce()[0], -2.01801, 1e-3);
        TS_ASSERT_DELTA(cell_population.GetNode(6)->rGetAppliedForce()[1], 0.731214, 1e-3);
    }

    void TestPeriodicForceOnHoneycombMeshWithGhostNodes() throw (Exception)
    {
        EXIT_IF_PARALLEL; // HoneycombMeshGenerator doesn't work in parallel

        HoneycombMeshGenerator generator(2, 2, 1);

        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Set up cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel,2> cells_generator;
        cells_generator.GenerateGivenLocationIndices(cells, location_indices);

        // Create a cell population
        MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);

        // Create force
        GeneralisedPeriodicLinearSpringForce<2> linear_force;
        linear_force.SetPeriodicDomainWidth(2.0);

        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
            cell_population.GetNode(i)->ClearAppliedForce();
        }

        linear_force.AddForceContribution(cell_population);

        TS_ASSERT_DELTA(cell_population.GetNode(5)->rGetAppliedForce()[0], 0.0, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(5)->rGetAppliedForce()[1], 0.0, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(6)->rGetAppliedForce()[0], 0.0, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(6)->rGetAppliedForce()[1], 0.0, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(9)->rGetAppliedForce()[0], 0.0, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(9)->rGetAppliedForce()[1], 0.0, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(10)->rGetAppliedForce()[0], 0.0, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(10)->rGetAppliedForce()[1], 0.0, 1e-4);
    }

    void TestSimpleNonZeroForces() throw(Exception)
    {
                // Create 2D mesh
        std::vector<Node<2> *> nodes;
        nodes.push_back(new Node<2>(0, true, 0.1, 0.0));
        nodes.push_back(new Node<2>(1, true, 0.9, 0.0));
        nodes.push_back(new Node<2>(2, true, 0.1, 1.0));
        nodes.push_back(new Node<2>(3, true, 0.9, 1.0));
        MutableMesh<2,2> mesh(nodes);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Create cell population
        MeshBasedCellPopulation<2> cell_population(mesh, cells);

        // Create Voronoi tessellation
        cell_population.CreateVoronoiTessellation();

        // Create force
        GeneralisedPeriodicLinearSpringForce<2> linear_force;
        linear_force.SetPeriodicDomainWidth(1.0);
        linear_force.SetCutOffLength(0.5);

        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
            cell_population.GetNode(i)->ClearAppliedForce();
        }

        linear_force.AddForceContribution(cell_population);

        TS_ASSERT_DELTA(cell_population.GetNode(0)->rGetAppliedForce()[0], 12.0, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(0)->rGetAppliedForce()[1], 0.0, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(1)->rGetAppliedForce()[0], -12.0, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(1)->rGetAppliedForce()[1], 0.0, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(2)->rGetAppliedForce()[0], 12.0, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(2)->rGetAppliedForce()[1], 0.0, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(3)->rGetAppliedForce()[0], -12.0, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(3)->rGetAppliedForce()[1], 0.0, 1e-4);
    }

    void TestSimpleNonZeroForcesWithGhostNodes() throw(Exception)
    {
        // Create 2D mesh
        std::vector<Node<2> *> nodes;
        nodes.push_back(new Node<2>(0, true, 0.1, 0.0));
        nodes.push_back(new Node<2>(1, true, 0.9, 0.0));
        nodes.push_back(new Node<2>(2, true, 0.1, 1.0));
        nodes.push_back(new Node<2>(3, true, 0.9, 1.0));
        nodes.push_back(new Node<2>(4, true, 0.5, -0.2));
        nodes.push_back(new Node<2>(5, true, 0.5, 1.1));
        MutableMesh<2,2> mesh(nodes);

        std::vector<unsigned> location_indices;
        location_indices.push_back(0);
        location_indices.push_back(1);
        location_indices.push_back(2);
        location_indices.push_back(3);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateGivenLocationIndices(cells, location_indices);

        // Create cell population
        MeshBasedCellPopulationWithGhostNodes<2> cell_population(mesh, cells, location_indices);

        // Create Voronoi tessellation
        cell_population.CreateVoronoiTessellation();

        // Create force
        GeneralisedPeriodicLinearSpringForce<2> linear_force;
        linear_force.SetPeriodicDomainWidth(1.0);
        linear_force.SetCutOffLength(0.5);

        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
            cell_population.GetNode(i)->ClearAppliedForce();
        }

        linear_force.AddForceContribution(cell_population);

        TS_ASSERT_DELTA(cell_population.GetNode(0)->rGetAppliedForce()[0], 12.0, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(0)->rGetAppliedForce()[1], 0.0, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(1)->rGetAppliedForce()[0], -12.0, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(1)->rGetAppliedForce()[1], 0.0, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(2)->rGetAppliedForce()[0], 12.0, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(2)->rGetAppliedForce()[1], 0.0, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(3)->rGetAppliedForce()[0], -12.0, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(3)->rGetAppliedForce()[1], 0.0, 1e-4);
    }

    void TestPeriodicSpringForces3d() throw(Exception)
    {
        unsigned width = 2;
        unsigned height = 2;
        unsigned depth = 1;

        MutableMesh<3,3>* p_mesh = Make3dMesh(width, height, depth);

        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumNodes());

        // Test Save with a MeshBasedCellPopulationWithGhostNodes
        MeshBasedCellPopulation<3> cell_population(*p_mesh, cells);

        OffLatticeSimulation<3> simulator(cell_population);

        // Create periodic force law
        MAKE_PTR(GeneralisedPeriodicLinearSpringForce<3>, p_periodic_force); // Variable spring strengths
        p_periodic_force->SetPeriodicDomainWidth(3.0);
        simulator.AddForce(p_periodic_force);

        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
             cell_population.GetNode(i)->ClearAppliedForce();
        }

        SimulationTime::Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(0.01, 1);

        p_periodic_force->AddForceContribution(cell_population);

        SimulationTime::Destroy();

        TS_ASSERT_DELTA(cell_population.GetNode(0)->rGetAppliedForce()[0], 10.7331, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(0)->rGetAppliedForce()[1], 10.7331, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(0)->rGetAppliedForce()[2], 36.5928, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(1)->rGetAppliedForce()[0], 25.8597, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(1)->rGetAppliedForce()[1], 0.0, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(1)->rGetAppliedForce()[2], 25.8597, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(2)->rGetAppliedForce()[0], 15.1265, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(2)->rGetAppliedForce()[1], 10.7331, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(2)->rGetAppliedForce()[2], 19.5199, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(3)->rGetAppliedForce()[0], 0.0, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(3)->rGetAppliedForce()[1], 25.8597, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(3)->rGetAppliedForce()[2], 25.8597, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(4)->rGetAppliedForce()[0], 10.7331, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(4)->rGetAppliedForce()[1], 10.7331, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(4)->rGetAppliedForce()[2], 15.1265, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(5)->rGetAppliedForce()[0], 10.7331, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(5)->rGetAppliedForce()[1], 10.7331, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(5)->rGetAppliedForce()[2], 15.1265, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(6)->rGetAppliedForce()[0], 10.7331, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(6)->rGetAppliedForce()[1], 15.1265, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(6)->rGetAppliedForce()[2], 19.5199, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(7)->rGetAppliedForce()[0], 10.7331, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(7)->rGetAppliedForce()[1], 10.7331, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(7)->rGetAppliedForce()[2], 15.1265, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(8)->rGetAppliedForce()[0], 0.0, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(8)->rGetAppliedForce()[1], 0.0, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(8)->rGetAppliedForce()[2], 8.7868, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(9)->rGetAppliedForce()[0], 4.3934, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(9)->rGetAppliedForce()[1], 4.3934, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(9)->rGetAppliedForce()[2], 0.0, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(10)->rGetAppliedForce()[0], 4.3934, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(10)->rGetAppliedForce()[1], 0.0, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(10)->rGetAppliedForce()[2], -4.3934, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(11)->rGetAppliedForce()[0], 0.0, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(11)->rGetAppliedForce()[1], 4.3934, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(11)->rGetAppliedForce()[2], -4.3934, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(12)->rGetAppliedForce()[0], 0.0, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(12)->rGetAppliedForce()[1], 4.3934, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(12)->rGetAppliedForce()[2], -4.3934, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(13)->rGetAppliedForce()[0], -10.7331, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(13)->rGetAppliedForce()[1], -10.7331, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(13)->rGetAppliedForce()[2], -15.1265, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(14)->rGetAppliedForce()[0], -10.7331, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(14)->rGetAppliedForce()[1], -10.7331, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(14)->rGetAppliedForce()[2], -15.1265, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(15)->rGetAppliedForce()[0], 4.3934, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(15)->rGetAppliedForce()[1], 0.0, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(15)->rGetAppliedForce()[2], -4.3934, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(16)->rGetAppliedForce()[0], -10.7331, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(16)->rGetAppliedForce()[1], -10.7331, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(16)->rGetAppliedForce()[2], -15.1265, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(17)->rGetAppliedForce()[0], -15.1265, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(17)->rGetAppliedForce()[1], -15.1265, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetNode(17)->rGetAppliedForce()[2], -15.1265, 1e-4);

        // Tidy up
        delete p_mesh;
    }

    void TestGeneralisedPeriodicLinearSpringForceArchiving() throw (Exception)
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "GeneralisedPeriodicLinearSpringForce.arch";

        {
            // Create force object and set member variables
            GeneralisedPeriodicLinearSpringForce<2> force;

            force.SetPeriodicDomainWidth(7.3);
            force.SetMeinekeSpringStiffness(13.0);
            force.SetMeinekeDivisionRestingSpringLength(0.6);
            force.SetMeinekeSpringGrowthDuration(2.0);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Serialize via pointer to most abstract class possible
            AbstractForce<2>* const p_force = &force;

            output_arch << p_force;
        }

        {
            AbstractForce<2>* p_force;

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // Restore from the archive
            input_arch >> p_force;

            // Check member variables have been correctly archived
            TS_ASSERT_DELTA(dynamic_cast<GeneralisedPeriodicLinearSpringForce<2>*>(p_force)->GetPeriodicDomainWidth(), 7.3, 1e-6);
            TS_ASSERT_DELTA(dynamic_cast<GeneralisedPeriodicLinearSpringForce<2>*>(p_force)->GetMeinekeSpringStiffness(), 13.0, 1e-6);
            TS_ASSERT_DELTA(dynamic_cast<GeneralisedPeriodicLinearSpringForce<2>*>(p_force)->GetMeinekeDivisionRestingSpringLength(), 0.6, 1e-6);
            TS_ASSERT_DELTA(dynamic_cast<GeneralisedPeriodicLinearSpringForce<2>*>(p_force)->GetMeinekeSpringGrowthDuration(), 2.0, 1e-6);

            // Tidy up
            delete p_force;
        }
    }

    void TestForceOutputParameters()
    {
        std::string output_directory = "TestPeriodicForcesOutputParameters";
        OutputFileHandler output_file_handler(output_directory, false);

        // Test with GeneralisedPeriodicLinearSpringForce
        GeneralisedPeriodicLinearSpringForce<2> force;

        force.SetPeriodicDomainWidth(7.3);
        force.SetMeinekeSpringStiffness(13.0);
        force.SetMeinekeDivisionRestingSpringLength(0.6);
        force.SetMeinekeSpringGrowthDuration(2.0);

        TS_ASSERT_EQUALS(force.GetIdentifier(), "GeneralisedPeriodicLinearSpringForce-2");

        out_stream force_parameter_file = output_file_handler.OpenOutputFile("periodic_results.parameters");
        force.OutputForceParameters(force_parameter_file);
        force_parameter_file->close();

        FileFinder generated_file = output_file_handler.FindFile("periodic_results.parameters");
        FileFinder reference_file("notforrelease_cell_based/test/data/TestPeriodicForces/periodic_results.parameters",RelativeTo::ChasteSourceRoot);
        FileComparison comparer(generated_file,reference_file);
        TS_ASSERT(comparer.CompareFiles());
    }
};

#endif /* TESTPERIODICFORCES_HPP_ */
