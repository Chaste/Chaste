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

#ifndef TESTOFFLATTICESIMULATIONWITHNODEBASEDCELLPOPULATION_HPP_
#define TESTOFFLATTICESIMULATIONWITHNODEBASEDCELLPOPULATION_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "CellsGenerator.hpp"
#include "OffLatticeSimulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "RandomCellKiller.hpp"
#include "PlaneBasedCellKiller.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "LogFile.hpp"
#include "WildTypeCellMutationState.hpp"
#include "PlaneBoundaryCondition.hpp"
#include "SmartPointers.hpp"
#include "CellVolumesWriter.hpp"

// Cell population writers
#include "CellMutationStatesCountWriter.hpp"
#include "CellProliferativePhasesCountWriter.hpp"
#include "CellProliferativeTypesCountWriter.hpp"
#include "CellPopulationAreaWriter.hpp"
#include "NodeVelocityWriter.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestOffLatticeSimulationWithNodeBasedCellPopulation : public AbstractCellBasedTestSuite
{
private:

    double mLastStartTime;
    void setUp()
    {
        mLastStartTime = (double) std::clock();
        AbstractCellBasedTestSuite::setUp();
    }
    void tearDown()
    {
        double time = (double) std::clock();
        double elapsed_time = (time - mLastStartTime)/(CLOCKS_PER_SEC);
        std::cout << "Elapsed time: " << elapsed_time << std::endl;
        AbstractCellBasedTestSuite::tearDown();
    }

public:

    /**
     * Create a simulation of a NodeBasedCellPopulation with a NodeBasedCellPopulationMechanicsSystem.
     * Test that no exceptions are thrown, and write the results to file.
     */
    void TestSimpleMonolayer() throw (Exception)
    {
        EXIT_IF_PARALLEL;    // HoneycombMeshGenereator does not work in parallel.

        // Create a simple mesh
        unsigned num_cells_depth = 5;
        unsigned num_cells_width = 5;
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0);
        TetrahedralMesh<2,2>* p_generating_mesh = generator.GetMesh();

        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh, 1.5);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes());

        // Create a node-based cell population
        NodeBasedCellPopulation<2> node_based_cell_population(mesh, cells);

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(node_based_cell_population);
        simulator.SetOutputDirectory("TestOffLatticeSimulationWithNodeBasedCellPopulation");

        // No need to go for long, don't want any birth or regular grid will be disrupted.
        simulator.SetEndTime(0.5);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetCutOffLength(1.5);
        simulator.AddForce(p_linear_force);

        simulator.Solve();

        // Check that nothing's gone badly wrong by testing that nodes aren't too close together
        double min_distance_between_cells = 1.0;

        for (unsigned i=0; i<simulator.rGetCellPopulation().GetNumNodes(); i++)
        {
            for (unsigned j=i+1; j<simulator.rGetCellPopulation().GetNumNodes(); j++)
            {
                double distance = norm_2(simulator.rGetCellPopulation().GetNode(i)->rGetLocation()-simulator.rGetCellPopulation().GetNode(j)->rGetLocation());

                if (distance < min_distance_between_cells)
                {
                    min_distance_between_cells = distance;
                }
            }
        }

        TS_ASSERT(min_distance_between_cells > 0.999);
    }

    /**
     * Create a simulation of a NodeBasedCellPopulation with different cell radii
     */
    void TestSimpleMonolayerWithDifferentRadii() throw (Exception)
    {
        // Creates nodes and mesh
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0,  false,  0.0, 0.0));
        nodes.push_back(new Node<2>(1,  false,  1.0, 0.0));
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 5.0);    // Large cut off as larger cells.

        // Modify the radii of the cells
        if (PetscTools::AmMaster())
        {
            mesh.GetNode(0)->SetRadius(1.0);
            mesh.GetNode(PetscTools::GetNumProcs())->SetRadius(2.0);
        }

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_transit_type);

        // Create a node-based cell population
        NodeBasedCellPopulation<2> node_based_cell_population(mesh, cells);
        node_based_cell_population.AddCellWriter<CellVolumesWriter>();
        node_based_cell_population.Update();

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(node_based_cell_population);
        simulator.SetOutputDirectory("TestOffLatticeSimulationWithNodeBasedCellPopulationAndDifferentRadi");
        simulator.SetSamplingTimestepMultiple(12);
        simulator.SetEndTime(12.0);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetCutOffLength(5.0); // Different as bigger cells
        simulator.AddForce(p_linear_force);

//        MAKE_PTR_ARGS(PlaneBasedCellKiller<2>, p_cell_killer,(&node_based_cell_population, 5.0*unit_vector<double>(2,0), unit_vector<double>(2,1)));
//        simulator.AddCellKiller(p_cell_killer);

        simulator.Solve();

        // Check that the radii of all the cells are correct
        // (cell 0 divided into 0 and 3 and cell 1 divided into 1 and 2)
        if (PetscTools::AmMaster())
        {
            TS_ASSERT_DELTA(mesh.GetNode(0)->GetRadius(), 1.0, 1e-6);
            TS_ASSERT_DELTA(mesh.GetNode(PetscTools::GetNumProcs())->GetRadius(), 2.0, 1e-6);
            TS_ASSERT_DELTA(mesh.GetNode(2*PetscTools::GetNumProcs())->GetRadius(), 2.0, 1e-6);
            TS_ASSERT_DELTA(mesh.GetNode(3*PetscTools::GetNumProcs())->GetRadius(), 1.0, 1e-6);

            // Check the separation of some node pairs
            TS_ASSERT_DELTA(norm_2(simulator.rGetCellPopulation().GetNode(0)->rGetLocation()-simulator.rGetCellPopulation().GetNode(PetscTools::GetNumProcs())->rGetLocation()), 2.9710, 1e-1);
            TS_ASSERT_DELTA(norm_2(simulator.rGetCellPopulation().GetNode(0)->rGetLocation()-simulator.rGetCellPopulation().GetNode(2*PetscTools::GetNumProcs())->rGetLocation()), 4.7067, 1e-1);
            TS_ASSERT_DELTA(norm_2(simulator.rGetCellPopulation().GetNode(0)->rGetLocation()-simulator.rGetCellPopulation().GetNode(3*PetscTools::GetNumProcs())->rGetLocation()), 2.0, 1e-1);
        }

        // Clean up memory
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

    /**
     * Create a simulation of a NodeBasedCellPopulation with variable cell radii
     */
    void TestSimpleMonolayerWithVariableRadii() throw (Exception)
    {
        // Creates nodes and mesh
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0,  false,  0.0, 0.0));
        nodes.push_back(new Node<2>(1,  false,  1.0, 0.0));
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 5.0);    // Larger cut off as bigger cells.

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_transit_type);

        // Store the radius of the cells in Cell Data
        if (PetscTools::AmMaster())
        {
            cells[0]->GetCellData()->SetItem("Radius", 1.0);
            cells[1]->GetCellData()->SetItem("Radius", 2.0);
        }

        // Create a node-based cell population
        NodeBasedCellPopulation<2> node_based_cell_population(mesh, cells);
        node_based_cell_population.SetUseVariableRadii(true);
        node_based_cell_population.AddCellWriter<CellVolumesWriter>();
        node_based_cell_population.Update();

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(node_based_cell_population);
        simulator.SetOutputDirectory("TestOffLatticeSimulationWithNodeBasedCellPopulationAndVariableRadii");
        simulator.SetSamplingTimestepMultiple(12);
        simulator.SetEndTime(10.0);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetCutOffLength(5.0); // Different as bigger cells
        simulator.AddForce(p_linear_force);

        simulator.Solve();

        // Check the Radii of all the cells are correct cell 0 divided into 0 and 3 and cell 1 divided into 1 and 2.
        // This testing is designed for sequential code.
        if (PetscTools::IsSequential())
        {
            TS_ASSERT_DELTA(mesh.GetNode(0)->GetRadius(), 1.0, 1e-6);
            TS_ASSERT_DELTA(mesh.GetNode(1)->GetRadius(), 2.0, 1e-6);
            TS_ASSERT_DELTA(mesh.GetNode(2)->GetRadius(), 2.0, 1e-6);
            TS_ASSERT_DELTA(mesh.GetNode(3)->GetRadius(), 1.0, 1e-6);

            // Check the separation of some node pairs
            TS_ASSERT_DELTA(norm_2(simulator.rGetCellPopulation().GetNode(0)->rGetLocation()-simulator.rGetCellPopulation().GetNode(1)->rGetLocation()), 3.0, 1e-1);
            TS_ASSERT_DELTA(norm_2(simulator.rGetCellPopulation().GetNode(0)->rGetLocation()-simulator.rGetCellPopulation().GetNode(2)->rGetLocation()), 4.70670, 1e-1);
            TS_ASSERT_DELTA(norm_2(simulator.rGetCellPopulation().GetNode(0)->rGetLocation()-simulator.rGetCellPopulation().GetNode(3)->rGetLocation()), 2.0, 1e-1);

            // Now set all the Radii to 2.0 Note this could be done inside a cell cycle model.
            for (AbstractCellPopulation<2>::Iterator cell_iter = simulator.rGetCellPopulation().Begin();
                 cell_iter != simulator.rGetCellPopulation().End();
                 ++cell_iter)
            {
                cell_iter->GetCellData()->SetItem("Radius",2.0);
            }

            simulator.SetEndTime(12.0);
            simulator.Solve();

            for (unsigned i=0; i<simulator.rGetCellPopulation().GetNumNodes(); i++)
            {
                TS_ASSERT_DELTA(mesh.GetNode(i)->GetRadius(), 2.0, 1e-6);
            }

            // Check the separation of some node pairs
            TS_ASSERT_DELTA(norm_2(simulator.rGetCellPopulation().GetNode(0)->rGetLocation()-simulator.rGetCellPopulation().GetNode(1)->rGetLocation()), 4.0, 1e-3);
            TS_ASSERT_DELTA(norm_2(simulator.rGetCellPopulation().GetNode(0)->rGetLocation()-simulator.rGetCellPopulation().GetNode(2)->rGetLocation()), 6.9282, 1e-3);
            TS_ASSERT_DELTA(norm_2(simulator.rGetCellPopulation().GetNode(0)->rGetLocation()-simulator.rGetCellPopulation().GetNode(3)->rGetLocation()), 4.0, 1e-3);
        }

        // Clean up memory
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

    void TestSimulationWithBoxes() throw (Exception)
    {
        EXIT_IF_PARALLEL;    // HoneycombMeshGenereator does not work in parallel.

        // Create a simple mesh
        int num_cells_depth = 5;
        int num_cells_width = 5;
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0);
        TetrahedralMesh<2,2>* p_generating_mesh = generator.GetMesh();

        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh, 1.5);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes());

        // Create a node-based cell population
        NodeBasedCellPopulation<2> node_based_cell_population(mesh, cells);

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(node_based_cell_population);
        simulator.SetOutputDirectory("TestOffLatticeSimulationWithNodeBasedCellPopulation");
        simulator.SetEndTime(1.0);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetCutOffLength(1.5);
        simulator.AddForce(p_linear_force);

        TS_ASSERT_THROWS_NOTHING(simulator.Solve());

        // Check that nothing's gone badly wrong by testing that nodes aren't too close together
        double min_distance_between_cells = 1.0;

        for (unsigned i=0; i<simulator.rGetCellPopulation().GetNumNodes(); i++)
        {
            for (unsigned j=i+1; j<simulator.rGetCellPopulation().GetNumNodes(); j++)
            {
                double distance = norm_2(simulator.rGetCellPopulation().GetNode(i)->rGetLocation()-simulator.rGetCellPopulation().GetNode(j)->rGetLocation());
                if (distance < min_distance_between_cells)
                {
                    min_distance_between_cells = distance;
                }
            }
        }

        TS_ASSERT(min_distance_between_cells > 1e-3);
    }

    /**
     * Create a simulation of a NodeBasedCellPopulation with a NodeBasedCellPopulationMechanicsSystem
     * and a CellKiller. Test that no exceptions are thrown, and write the results to file.
     */
    void TestCellDeath() throw (Exception)
    {
        EXIT_IF_PARALLEL;    // HoneycombMeshGenereator does not work in parallel.

        // Create a simple mesh
        int num_cells_depth = 5;
        int num_cells_width = 5;
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0);
        TetrahedralMesh<2,2>* p_generating_mesh = generator.GetMesh();

        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh, 1.5);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes());

        // Create a node based cell population
        NodeBasedCellPopulation<2> node_based_cell_population(mesh, cells);

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(node_based_cell_population);
        simulator.SetOutputDirectory("TestOffLatticeSimulationWithNodeBasedCellPopulationCellPtrDeath");
        simulator.SetEndTime(0.5);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetCutOffLength(1.5);
        simulator.AddForce(p_linear_force);

        // Add cell killer
        MAKE_PTR_ARGS(RandomCellKiller<2>, p_killer, (&node_based_cell_population, 0.997877574));
        simulator.AddCellKiller(p_killer);

        // Solve
        simulator.Solve();

        // Check some results
        TS_ASSERT_EQUALS(simulator.GetNumBirths(), 0u);
        TS_ASSERT_EQUALS(simulator.GetNumDeaths(), 20u);

        std::vector<double> node_8_location = simulator.GetNodeLocation(8);
        TS_ASSERT_DELTA(node_8_location[0], 3.4729, 1e-4);
        TS_ASSERT_DELTA(node_8_location[1], 1.0051, 1e-4);

        std::vector<double> node_3_location = simulator.GetNodeLocation(3);
        TS_ASSERT_DELTA(node_3_location[0], 2.9895, 1e-4);
        TS_ASSERT_DELTA(node_3_location[1], 0.3105, 1e-4);
    }

    double mNode3x, mNode4x, mNode3y, mNode4y; // To preserve locations between the below test and test load.

    void TestStandardResultForArchivingTestsBelow() throw (Exception)
    {
        EXIT_IF_PARALLEL;    // HoneycombMeshGenereator does not work in parallel.

        // Create a simple mesh
        int num_cells_depth = 5;
        int num_cells_width = 5;
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0);
        TetrahedralMesh<2,2>* p_generating_mesh = generator.GetMesh();

        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh, 1.5);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes());

        // Create a node based cell population
        NodeBasedCellPopulation<2> node_based_cell_population(mesh, cells);

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(node_based_cell_population);
        simulator.SetOutputDirectory("TestOffLatticeSimulationWithNodeBasedCellPopulationStandardResult");
        simulator.SetEndTime(2.5);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetCutOffLength(1.5);
        simulator.AddForce(p_linear_force);

        // Create some boundary conditions and pass them to the simulation
        c_vector<double,2> normal = zero_vector<double>(2);
        normal(1) =-1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc, (&node_based_cell_population, zero_vector<double>(2), normal)); // y>0
        simulator.AddCellPopulationBoundaryCondition(p_bc);

        // Solve
        simulator.Solve();

        // Check some results
        mNode3x = 3.0454;
        mNode3y = 0.0000;
        mNode4x = 4.0468;
        mNode4y = 0.0101;

        std::vector<double> node_3_location = simulator.GetNodeLocation(3);
        TS_ASSERT_DELTA(node_3_location[0], mNode3x, 1e-4);
        TS_ASSERT_DELTA(node_3_location[1], mNode3y, 1e-4);

        std::vector<double> node_4_location = simulator.GetNodeLocation(4);
        TS_ASSERT_DELTA(node_4_location[0], mNode4x, 1e-4);
        TS_ASSERT_DELTA(node_4_location[1], mNode4y, 1e-4);
    }

    // Testing Save
    void TestSave() throw (Exception)
    {
        EXIT_IF_PARALLEL; // HoneycombMeshGenereator does not work in parallel

        // Create a simple mesh
        int num_cells_depth = 5;
        int num_cells_width = 5;
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0);
        TetrahedralMesh<2,2>* p_generating_mesh = generator.GetMesh();

        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh, 1.5);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes());

        // Create a node based cell population
        NodeBasedCellPopulation<2> node_based_cell_population(mesh, cells);

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(node_based_cell_population);
        simulator.SetOutputDirectory("TestOffLatticeSimulationWithNodeBasedCellPopulationSaveAndLoad");
        simulator.SetEndTime(0.1);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetCutOffLength(1.5);
        simulator.AddForce(p_linear_force);

        // Create some boundary conditions and pass them to the simulation
        c_vector<double,2> normal = zero_vector<double>(2);
        normal(1) =-1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc, (&node_based_cell_population, zero_vector<double>(2), normal)); // y>0
        simulator.AddCellPopulationBoundaryCondition(p_bc);

        // Solve
        simulator.Solve();

        // Save the results
        CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Save(&simulator);
    }

    // Testing Load (based on previous two tests)
    void TestLoad() throw (Exception)
    {
        EXIT_IF_PARALLEL;    // Cell based archiving doesn't work in parallel.

        // Load the simulation from the TestSave method above and
        // run it from 0.1 to 1.0
        OffLatticeSimulation<2>* p_simulator1;
        p_simulator1 = CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Load("TestOffLatticeSimulationWithNodeBasedCellPopulationSaveAndLoad", 0.1);

        p_simulator1->SetEndTime(1.0);
        p_simulator1->Solve();

        // Save, then reload and run from 1.0 to 2.5
        CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Save(p_simulator1);
        OffLatticeSimulation<2>* p_simulator2
            = CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Load("TestOffLatticeSimulationWithNodeBasedCellPopulationSaveAndLoad", 1.0);

        p_simulator2->SetEndTime(2.5);
        p_simulator2->Solve();

        // These results are from time 2.5 in TestStandardResultForArchivingTestBelow() (above!)
        std::vector<double> node_3_location = p_simulator2->GetNodeLocation(3);
        TS_ASSERT_DELTA(node_3_location[0], mNode3x, 1e-4);
        TS_ASSERT_DELTA(node_3_location[1], mNode3y, 1e-4);

        std::vector<double> node_4_location = p_simulator2->GetNodeLocation(4);
        TS_ASSERT_DELTA(node_4_location[0], mNode4x, 1e-4);
        TS_ASSERT_DELTA(node_4_location[1], mNode4y, 1e-4);

        // Tidy up
        delete p_simulator1;
        delete p_simulator2;
    }

    /**
     * Create a simulation of a NodeBasedCellPopulation to test movement threshold.
     */
    void TestMovementThreshold() throw (Exception)
    {
        EXIT_IF_PARALLEL;   // This test doesn't work in parallel because only one process will throw.

        // Creates nodes and mesh
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0,  false,  0.0, 0.0));
        nodes.push_back(new Node<2>(0,  false,  0.0, 0.3));
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes());

        // Create a node based cell population
        NodeBasedCellPopulation<2> node_based_cell_population(mesh, cells);
        node_based_cell_population.SetAbsoluteMovementThreshold(1e-6);

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(node_based_cell_population);
        simulator.SetEndTime(0.1);
        simulator.SetOutputDirectory("TestOffLatticeSimulationWithNodeBasedCellPopulationThreshold");

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetCutOffLength(1.5);
        simulator.AddForce(p_linear_force);

        // Solve
        TS_ASSERT_THROWS_CONTAINS(simulator.Solve(),
                "which is more than the AbsoluteMovementThreshold:");

        // Avoid memory leak
        delete nodes[0];
        delete nodes[1];
    }

    void TestUpdateCellLocationsAndTopologyWithNoForce()
    {
        // Creates nodes and mesh
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0,  false,  0.0, 0.0));
        nodes.push_back(new Node<2>(0,  false,  0.0, 0.3));
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes());

        // Create a node based cell population
        NodeBasedCellPopulation<2> node_based_cell_population(mesh, cells);
        node_based_cell_population.SetAbsoluteMovementThreshold(1e-6);

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(node_based_cell_population);
        simulator.SetEndTime(0.1);
        simulator.SetOutputDirectory("TestOffLatticeSimulationUpdateCellLocationsAndTopologyWithNoForce");

        simulator.UpdateCellLocationsAndTopology();

        if (PetscTools::AmMaster())
        {
            for (AbstractMesh<2,2>::NodeIterator node_iter = mesh.GetNodeIteratorBegin();
                    node_iter != mesh.GetNodeIteratorEnd();
                    ++node_iter)
            {
                for (unsigned d=0; d<2; d++)
                {
                    TS_ASSERT_DELTA(node_iter->rGetAppliedForce()[d], 0.0, 1e-15);
                }
            }
        }

        // Avoid memory leak
        delete nodes[0];
        delete nodes[1];
    }
};

#endif /*TESTOFFLATTICESIMULATIONWITHNODEBASEDCELLPOPULATION_HPP_*/
