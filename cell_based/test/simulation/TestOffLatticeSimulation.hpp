/*

Copyright (c) 2005-2012, University of Oxford.
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

#ifndef TESTOFFLATTICESIMULATION_HPP_
#define TESTOFFLATTICESIMULATION_HPP_

#include <cxxtest/TestSuite.h>

#include <cstdio>
#include <ctime>
#include <cmath>

#include "CheckpointArchiveTypes.hpp"
#include "OffLatticeSimulation.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "CylindricalHoneycombMeshGenerator.hpp"
#include "CellsGenerator.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "CellwiseData.hpp"
#include "ChemotacticForce.hpp"
#include "RandomCellKiller.hpp"
#include "PlaneBoundaryCondition.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "NumericFileComparison.hpp"
#include "CellBasedEventHandler.hpp"
#include "WildTypeCellMutationState.hpp"
#include "OffLatticeSimulationWithMyStoppingEvent.hpp"
#include "SmartPointers.hpp"


class TestOffLatticeSimulation : public AbstractCellBasedTestSuite
{
private:

    double mLastStartTime;
    void setUp()
    {
        mLastStartTime = std::clock();
        AbstractCellBasedTestSuite::setUp();
    }
    void tearDown()
    {
        double time = std::clock();
        double elapsed_time = (time - mLastStartTime)/(CLOCKS_PER_SEC);
        std::cout << "Elapsed time: " << elapsed_time << std::endl;
        AbstractCellBasedTestSuite::tearDown();
    }

public:

    void TestOutputNodeVelocities() throw(Exception)
    {
        // Create a simple mesh
        HoneycombMeshGenerator generator(5, 5, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes(), DIFFERENTIATED);

        // Create a cell population
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Output Voronoi Data
        cell_population.SetOutputVoronoiData(true);

        // Set up simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestOutputNodeVelocities");
        simulator.SetEndTime(0.5);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        p_force->SetCutOffLength(1.5);
        simulator.AddForce(p_force);

        // Record node velocities
        TS_ASSERT_EQUALS(simulator.GetOutputNodeVelocities(), false);
        simulator.SetOutputNodeVelocities(true);

        // Run simulation
        simulator.Solve();

        // Check node velocities file
        // The velocities should all be zero(ish), as the cell population is in mechanical equilibrium
        OutputFileHandler handler("TestOutputNodeVelocities", false);

        std::string node_velocities_file = handler.GetOutputDirectoryFullPath() + "results_from_time_0/nodevelocities.dat";
        NumericFileComparison node_velocities(node_velocities_file, "cell_based/test/data/TestOutputNodeVelocities/nodevelocities.dat");
        TS_ASSERT(node_velocities.CompareFiles(1e-2));

        // Test vtk files exist.
#ifdef CHASTE_VTK
        std::string results_dir = handler.GetOutputDirectoryFullPath();

        // Initial condition file
        FileFinder vtk_file(results_dir + "results_from_time_0/results_0.vtu", RelativeTo::Absolute);
        TS_ASSERT(vtk_file.Exists());

        // Final file
        FileFinder vtk_file2(results_dir + "results_from_time_0/results_60.vtu", RelativeTo::Absolute);
        TS_ASSERT(vtk_file2.Exists());
#endif //CHASTE_VTK
    }

    void TestOutputNodeVelocitiesWithGhostNodes() throw(Exception)
    {
        // Create a simple mesh with a surrounding layer of ghost nodes
        HoneycombMeshGenerator generator(3, 3, 1);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Get location indices corresponding to real cells
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Create a differentiated cell for each non-ghost node
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, location_indices.size(), location_indices, DIFFERENTIATED);

        for (unsigned i=0; i<cells.size(); i++)
        {
            cells[i]->SetBirthTime(-1.0);
        }

        // Create cell population
        MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);

        //Output Voronoi data
        cell_population.SetOutputVoronoiData(true);

        // Set up simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestOutputNodeVelocitiesWithGhostNodes");
        simulator.SetEndTime(0.5);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        p_force->SetCutOffLength(1.5);
        simulator.AddForce(p_force);

        // Record node velocities
        simulator.SetOutputNodeVelocities(true);

        // Run simulation
        simulator.Solve();

        // Test Node Velocities File.
        std::string output_directory = "TestOutputNodeVelocitiesWithGhostNodes";
        OutputFileHandler output_file_handler(output_directory, false);

        std::string node_velocities_file = output_file_handler.GetOutputDirectoryFullPath() + "results_from_time_0/nodevelocities.dat";
        NumericFileComparison node_velocities(node_velocities_file, "cell_based/test/data/TestOutputNodeVelocitiesWithGhostNodes/nodevelocities.dat");
        TS_ASSERT(node_velocities.CompareFiles(1e-2));

        // Test vtk files exist.
#ifdef CHASTE_VTK

        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        // Initial condition file
        FileFinder vtk_file(results_dir + "results_from_time_0/results_0.vtu", RelativeTo::Absolute);
        TS_ASSERT(vtk_file.Exists());

        // Final file
        FileFinder vtk_file2(results_dir + "results_from_time_0/results_60.vtu", RelativeTo::Absolute);
        TS_ASSERT(vtk_file2.Exists());
#endif //CHASTE_VTK
    }

    /**
     * Test a cell-based simulation with a cell killer.
     *
     * In this test, we solve a cell-based simulation without ghost nodes and
     * check that the numbers of nodes and cells match at the end of the
     * simulation.
     */
    void TestOffLatticeSimulationWithCellDeath() throw (Exception)
    {
        // Create a simple mesh
        int num_cells_depth = 5;
        int num_cells_width = 5;
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes());

        // Create a cell population
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestOffLatticeSimulationWithCellDeath");
        simulator.SetEndTime(0.5);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        p_force->SetCutOffLength(1.5);
        simulator.AddForce(p_force);

        // Add cell killer
        MAKE_PTR_ARGS(RandomCellKiller<2>, p_killer, (&cell_population, 0.997877574));
        simulator.AddCellKiller(p_killer);

        // For coverage of an exception.
        simulator.SetUpdateCellPopulationRule(false);
        TS_ASSERT_THROWS_THIS(simulator.Solve(),"CellPopulation has had births or deaths but mUpdateCellPopulation is set to false, please set it to true.");
        CellBasedEventHandler::Reset(); // Otherwise logging has been started but not stopped due to exception above.

        simulator.SetUpdateCellPopulationRule(true);
        simulator.Solve();

        // Check that the number of nodes is equal to the number of cells
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumNodes(), simulator.rGetCellPopulation().GetNumRealCells());

        // For coverage of these 'Get' functions
        TS_ASSERT_EQUALS(simulator.GetNumBirths(), 0u);
        TS_ASSERT_EQUALS(simulator.GetNumDeaths(), 17u);

        // Note that this test used to take an extra time step
        TS_ASSERT_EQUALS(SimulationTime::Instance()->GetTime(), 0.5);
    }

    /**
     * Test a cell-based simulation with multiple forces.
     */
    void TestOffLatticeSimulationWithMultipleForces() throw (Exception)
    {
        // Create a simple mesh
        int num_cells_depth = 5;
        int num_cells_width = 5;
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes());

        // Create a cell population
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestOffLatticeSimulationWithMultipleForces");
        simulator.SetEndTime(0.5);

        // Create some force laws and pass them to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        simulator.AddForce(p_linear_force);

        // Need to set this up for the chemotactic force.
        CellwiseData<2>* p_data = CellwiseData<2>::Instance();
        p_data->SetNumCellsAndVars(cell_population.GetNumRealCells(), 1);
        p_data->SetCellPopulation(&cell_population);
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            double x = p_mesh->GetNode(i)->rGetLocation()[0];
            p_data->SetValue(x/50.0, p_mesh->GetNode(i)->GetIndex());
        }

        MAKE_PTR(ChemotacticForce<2>, p_chemotactic_force);
        simulator.AddForce(p_chemotactic_force);

        simulator.Solve();

        // Check that the number of nodes is equal to the number of cells
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumNodes(), simulator.rGetCellPopulation().GetNumRealCells());

        // For coverage of these 'Get' functions
        TS_ASSERT_EQUALS(simulator.GetNumBirths(), 0u);
        TS_ASSERT_EQUALS(simulator.GetNumDeaths(), 0u);
    }


    /**
     * Test a cell-based simulation with a periodic mesh.
     */
    void TestOffLatticeSimulationWithPeriodicMesh() throw (Exception)
    {
        // Create a simple mesh
        int cells_up = 6;
        int cells_across = 6;

        CylindricalHoneycombMeshGenerator generator(cells_across, cells_up, 0);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes());

        // Create a cell population
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestOffLatticeSimulationWithPeriodicMesh");
        simulator.SetEndTime(0.5);

        // Create some force laws and pass them to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        simulator.AddForce(p_linear_force);

        simulator.Solve();

        // Check that the number of nodes is equal to the number of cells
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumNodes(), simulator.rGetCellPopulation().GetNumRealCells());

        // Check that the setup file is written correctly
        std::string output_directory = "TestOffLatticeSimulationWithPeriodicMesh/results_from_time_0/";
        OutputFileHandler output_file_handler(output_directory, false);
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "results.vizsetup  cell_based/test/data/TestOffLatticeSimulationWithPeriodicMesh/results.vizsetup").c_str()), 0);
    }

    /**
     * Test a cell-based simulation with multiple boundary conditions. y<2 and y>0
     */
    void TestOffLatticeSimulationWithMultipleCellBoundaryConditions() throw (Exception)
    {
        // Create a simple mesh
        int num_cells_depth = 5;
        int num_cells_width = 5;
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes());

        // Create a cell population
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestOffLatticeSimulationWithMultipleCellBoundaryConditions");
        simulator.SetEndTime(0.5);

        // Create some force laws and pass them to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        simulator.AddForce(p_linear_force);

        // Create some boundary conditions and pass them to the simulation
        c_vector<double,2> point = zero_vector<double>(2);
        c_vector<double,2> normal = zero_vector<double>(2);
        normal(1) = -1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc1, (&cell_population, point, normal)); // y>0
        simulator.AddCellPopulationBoundaryCondition(p_bc1);
        point(1) = 2.0;
        normal(1) = 1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc2, (&cell_population, point, normal)); // y<2
        simulator.AddCellPopulationBoundaryCondition(p_bc2);

        simulator.Solve();

        // Check that the number of nodes is equal to the number of cells
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumNodes(), simulator.rGetCellPopulation().GetNumRealCells());

        /**
         * Test a cell-based simulation with contradicting boundary conditions. y<2, y>0 and y>3
         */
        point(1) = 3.0;
        normal(1) = -1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc3, (&cell_population, point, normal)); // y>3
        simulator.AddCellPopulationBoundaryCondition(p_bc3);

        simulator.SetEndTime(1.0);
        TS_ASSERT_THROWS_THIS(simulator.Solve(),"The cell population boundary conditions are incompatible.");
        CellBasedEventHandler::Reset(); // Otherwise logging has been started but not stopped due to exception above.
    }

    void TestOffLatticeSimulationWithStoppingEvent() throw (Exception)
    {
        HoneycombMeshGenerator generator(2, 2, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        CellPropertyRegistry::Instance()->Clear();
        RandomNumberGenerator* p_random_num_gen = RandomNumberGenerator::Instance();

        // Set up cells
        std::vector<CellPtr> cells;
        cells.clear();
        unsigned num_cells = location_indices.empty() ? p_mesh->GetNumNodes() : location_indices.size();
        cells.reserve(num_cells);

        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            CellProliferativeType cell_type;
            unsigned generation;
            double y = 0.0;

            if (std::find(location_indices.begin(), location_indices.end(), i) != location_indices.end())
            {
                y = p_mesh->GetNode(i)->GetPoint().rGetLocation()[1];
            }

            FixedDurationGenerationBasedCellCycleModel* p_cell_cycle_model = new FixedDurationGenerationBasedCellCycleModel;
            p_cell_cycle_model->SetDimension(2);

            double typical_transit_cycle_time = p_cell_cycle_model->GetAverageTransitCellCycleTime();
            double typical_stem_cycle_time = p_cell_cycle_model->GetAverageStemCellCycleTime();

            double birth_time = 0.0;
            birth_time = -p_random_num_gen->ranf();

            if (y <= 0.3)
            {
                cell_type = STEM;
                generation = 0;
                birth_time *= typical_stem_cycle_time; // hours
            }
            else if (y < 2.0)
            {
                cell_type = TRANSIT;
                generation = 1;
                birth_time *= typical_transit_cycle_time; // hours
            }
            else if (y < 3.0)
            {
                cell_type = TRANSIT;
                generation = 2;
                birth_time *= typical_transit_cycle_time; // hours
            }
            else if (y < 4.0)
            {
                cell_type = TRANSIT;
                generation = 3;
                birth_time *= typical_transit_cycle_time; // hours
            }
            else
            {
                cell_type = p_cell_cycle_model->CanCellTerminallyDifferentiate() ? DIFFERENTIATED : TRANSIT;
                generation = 4;
                birth_time *= typical_transit_cycle_time; // hours
            }

            p_cell_cycle_model->SetGeneration(generation);
            p_cell_cycle_model->SetCellProliferativeType(cell_type);

            boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());

            CellPtr p_cell(new Cell(p_state, p_cell_cycle_model));
            p_cell->SetBirthTime(birth_time);

            if (std::find(location_indices.begin(), location_indices.end(), i) != location_indices.end())
            {
                cells.push_back(p_cell);
            }
        }

        // Create a cell population
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up cell-based simulation WITH the stopping event
        OffLatticeSimulationWithMyStoppingEvent simulator(cell_population);
        simulator.SetOutputDirectory("TestCellPopulationSimWithStoppingEvent");

        // Set the end time to 10.0 - the stopping event is, however, t>3.1415.
        simulator.SetEndTime(10.0);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetCutOffLength(1.5);
        simulator.AddForce(p_linear_force);

        // Run cell-based simulation
        simulator.Solve();

        double time = SimulationTime::Instance()->GetTime();
        TS_ASSERT_DELTA(time, 3.1415, 1e-1); // big tol, doesn't matter, just want t~3.14 and t!=10
        // t should be strictly greater than the 3.1415
        TS_ASSERT_LESS_THAN(3.1415, time);
    }

    void TestApoptosisSpringLengths() throw (Exception)
    {
        unsigned num_cells_depth = 2;
        unsigned num_cells_width = 2;

        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 2);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, location_indices.size(), TRANSIT);

        // Create cell population
        MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);

        // Create crypt simulation from cell population and force law
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("2dSpheroidApoptosis");
        simulator.SetEndTime(1.0);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetCutOffLength(1.5);
        simulator.AddForce(p_linear_force);

        cell_population.GetCellUsingLocationIndex(14)->SetApoptosisTime(2.0);
        cell_population.GetCellUsingLocationIndex(15)->SetApoptosisTime(2.0);
        cell_population.GetCellUsingLocationIndex(14)->StartApoptosis();
        cell_population.GetCellUsingLocationIndex(15)->StartApoptosis();
        simulator.SetNoBirth(true);

        // Run cell-based simulation
        simulator.Solve();

        /*
         * We track the locations of two dying cells (a and b) and two
         * live cells adjacent to them (c and d)
         *
         * All cells begin distance 1 apart.
         *
         * a and b move together to leave a gap of 0.
         * a and c (and b and d) move to a distance of 0.5 apart.
         */

        c_vector<double, 2> a_location = cell_population.rGetMesh().GetNode(14)->rGetLocation();
        c_vector<double, 2> b_location = cell_population.rGetMesh().GetNode(15)->rGetLocation();
        c_vector<double, 2> c_location = cell_population.rGetMesh().GetNode(20)->rGetLocation();
        c_vector<double, 2> d_location = cell_population.rGetMesh().GetNode(21)->rGetLocation();

        double a_b_separation = sqrt((a_location[0]-b_location[0])*(a_location[0]-b_location[0]) +
                                (a_location[1]-b_location[1])*(a_location[1]-b_location[1]));
        double a_c_separation = sqrt((a_location[0]-c_location[0])*(a_location[0]-c_location[0]) +
                                (a_location[1]-c_location[1])*(a_location[1]-c_location[1]));
        double c_d_separation = sqrt((d_location[0]-c_location[0])*(d_location[0]-c_location[0]) +
                                (d_location[1]-c_location[1])*(d_location[1]-c_location[1]));

        TS_ASSERT_DELTA(a_b_separation, 0.5, 1e-1);
        TS_ASSERT_DELTA(a_c_separation, 0.75, 1e-1);
        TS_ASSERT_DELTA(c_d_separation, 1.0, 1e-1);

        // Reset end time and run cell-based simulation
        simulator.SetEndTime(1.99);
        simulator.Solve();

        a_location = cell_population.rGetMesh().GetNode(14)->rGetLocation();
        b_location = cell_population.rGetMesh().GetNode(15)->rGetLocation();
        c_location = cell_population.rGetMesh().GetNode(20)->rGetLocation();
        d_location = cell_population.rGetMesh().GetNode(21)->rGetLocation();

        a_b_separation = sqrt((a_location[0]-b_location[0])*(a_location[0]-b_location[0]) +
                         (a_location[1]-b_location[1])*(a_location[1]-b_location[1]));
        a_c_separation = sqrt((a_location[0]-c_location[0])*(a_location[0]-c_location[0]) +
                         (a_location[1]-c_location[1])*(a_location[1]-c_location[1]));
        c_d_separation = sqrt((d_location[0]-c_location[0])*(d_location[0]-c_location[0]) +
                         (d_location[1]-c_location[1])*(d_location[1]-c_location[1]));

        TS_ASSERT_DELTA(a_b_separation, 0.01, 1e-1);
        TS_ASSERT_DELTA(a_c_separation, 0.5, 1e-1);
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 4u);

        // Reset end time and run cell-based simulation
        simulator.SetEndTime(2.01);
        simulator.Solve();

        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 2u);
    }

    void TestOffLatticeSimulationParameterOutputMethods() throw (Exception)
    {
        // Create a simple mesh
        HoneycombMeshGenerator generator(5, 5, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Set up cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumNodes());

        // Create a cell population
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(cell_population);
        TS_ASSERT_EQUALS(simulator.GetIdentifier(), "OffLatticeSimulation-2");

        //#1453 should have forces and cell killer included here to make it a better test.

        std::string output_directory = "TestOffLatticeSimulationOutputParameters";
        OutputFileHandler output_file_handler(output_directory, false);
        out_stream parameter_file = output_file_handler.OpenOutputFile("cell_based_sim_results.parameters");
        simulator.OutputSimulationParameters(parameter_file);
        parameter_file->close();

        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "cell_based_sim_results.parameters  cell_based/test/data/TestOffLatticeSimulationOutputParameters/cell_based_sim_results.parameters").c_str()), 0);

        simulator.SetOutputDirectory("TestOffLatticeSimulationOutputParameters");
        simulator.OutputSimulationSetup();
        ///\todo #1453 This is to do with the pre Boost 1.37 problem ---
        //TS_ASSERT_EQUALS(system(("diff " + results_dir + "results.parameters  cell_based/test/data/TestOffLatticeSimulationOutputParameters/results.parameters").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff --ignore-matching-lines=\"CellPopulation\" " + results_dir + "results.parameters  cell_based/test/data/TestOffLatticeSimulationOutputParameters/results.parameters").c_str()), 0);

        // Check that the files which we don't want to compare actually exist
        std::ifstream machine_file;
        std::string command = results_dir + "/system_info_0.txt";
        machine_file.open(command.c_str());
        TS_ASSERT(machine_file.is_open());
        machine_file.close();

        std::ifstream info_file;
        command = results_dir + "/build.info";
        info_file.open(command.c_str());
        TS_ASSERT(info_file.is_open());
        info_file.close();
    }

    /**
     * The purpose of this test is to check that it is possible to construct and run
     * a short 1D cell-based simulation without throwing any exceptions.
     */
    void Test1dOffLatticeSimulation() throw (Exception)
    {
        // Create mesh
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        MutableMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Set up cells so that cell 10 divides at time t=0.5, cell 9 at time t=1.5, etc
        for (unsigned i=0; i<cells.size(); i++)
        {
            cells[i]->SetBirthTime(-13.5 - i);
        }

        // Create a cell population
        MeshBasedCellPopulation<1> cell_population(mesh, cells);

        // Coverage
        cell_population.SetOutputCellIdData(true);

        // Set up cell-based simulation
        OffLatticeSimulation<1> simulator(cell_population);
        simulator.SetOutputDirectory("Test1DOffLatticeSimulation");
        simulator.SetEndTime(0.6);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<1>, p_linear_force);
        simulator.AddForce(p_linear_force);

        unsigned initial_num_cells = simulator.rGetCellPopulation().GetNumRealCells();
        unsigned initial_num_nodes = simulator.rGetCellPopulation().GetNumNodes();
        unsigned initial_num_elements = (static_cast<MeshBasedCellPopulation<1>* >(&(simulator.rGetCellPopulation())))->rGetMesh().GetNumElements();

        // Run simulation for a short time
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());

        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), initial_num_cells + 1);
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumNodes(), initial_num_nodes + 1);
        TS_ASSERT_EQUALS((static_cast<MeshBasedCellPopulation<1>* >(&(simulator.rGetCellPopulation())))->rGetMesh().GetNumElements(), initial_num_elements + 1);
    }

    void TestSettingEndTimeIssue() throw(Exception)
    {
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(0.1, 1);

        HoneycombMeshGenerator generator(2, 2, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Create some cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes(), DIFFERENTIATED);

        // Create a mesh-based cell population
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        OffLatticeSimulation<2> simulator(cell_population);

        simulator.SetOutputDirectory("TestSettingEndTimeIssue");
        simulator.SetEndTime(1.0);
        TS_ASSERT_THROWS_THIS(simulator.Solve(),
                "End time and number of timesteps already setup. You should not use SimulationTime::SetEndTimeAndNumberOfTimeSteps in cell-based tests.");
    }


};

#endif /*TESTOFFLATTICESIMULATION_HPP_*/
