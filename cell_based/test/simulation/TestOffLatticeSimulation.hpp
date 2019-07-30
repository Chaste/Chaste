/*

Copyright (c) 2005-2019, University of Oxford.
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
#include <cmath>

#include "CheckpointArchiveTypes.hpp"
#include "OffLatticeSimulation.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "CylindricalHoneycombMeshGenerator.hpp"
#include "CellsGenerator.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "UniformCellCycleModel.hpp"
#include "NoCellCycleModel.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "ChemotacticForce.hpp"
#include "RandomCellKiller.hpp"
#include "PlaneBasedCellKiller.hpp"
#include "PlaneBoundaryCondition.hpp"
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "NumericFileComparison.hpp"
#include "CellBasedEventHandler.hpp"
#include "WildTypeCellMutationState.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "OffLatticeSimulationWithMyStoppingEvent.hpp"
#include "TransitCellProliferativeType.hpp"
#include "SmartPointers.hpp"
#include "FileComparison.hpp"
#include "CellIdWriter.hpp"
#include "VolumeTrackingModifier.hpp"
#include "CellBasedSimulationArchiver.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"
#include "DefaultCellProliferativeType.hpp"
#include "ForwardEulerNumericalMethod.hpp"

// Cell population writers
#include "CellMutationStatesCountWriter.hpp"
#include "CellProliferativeTypesCountWriter.hpp"
#include "NodeVelocityWriter.hpp"
#include "VoronoiDataWriter.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestOffLatticeSimulation : public AbstractCellBasedWithTimingsTestSuite
{
public:

    void TestOutputNodeAndCellVelocitiesAndDivisionLocations()
    {
        EXIT_IF_PARALLEL;    // HoneycombMeshGenerator does not work in parallel

        // Create a simple 2D MeshBasedCellPopulation
        HoneycombMeshGenerator generator(5, 5, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes(), p_transit_type);

        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Output Voronoi data
        cell_population.AddPopulationWriter<VoronoiDataWriter>();

        // Set up simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestOutputNodeAndCellVelocitiesAndDivisionLocations");
        simulator.SetEndTime(1.0);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        p_force->SetCutOffLength(1.5);
        simulator.AddForce(p_force);

        // Record node velocities
        cell_population.AddPopulationWriter<NodeVelocityWriter>();

        // Record division locations
        TS_ASSERT_EQUALS(simulator.GetOutputDivisionLocations(), false);
        simulator.SetOutputDivisionLocations(true);

        // Record cell velocities
        TS_ASSERT_EQUALS(simulator.GetOutputCellVelocities(), false);
        simulator.SetOutputCellVelocities(true);

        // Run simulation
        simulator.Solve();

        // Check node velocities file
        OutputFileHandler handler("TestOutputNodeAndCellVelocitiesAndDivisionLocations", false);

        std::string node_velocities_file = handler.GetOutputDirectoryFullPath() + "results_from_time_0/nodevelocities.dat";
        NumericFileComparison node_velocities(node_velocities_file, "cell_based/test/data/TestOutputNodeAndCellVelocitiesAndDivisionLocations/nodevelocities.dat");
        TS_ASSERT(node_velocities.CompareFiles(1e-2));

        // Check division locations file
        std::string division_locations_file = handler.GetOutputDirectoryFullPath() + "results_from_time_0/divisions.dat";
        NumericFileComparison division_locations(division_locations_file, "cell_based/test/data/TestOutputNodeAndCellVelocitiesAndDivisionLocations/divisions.dat");
        TS_ASSERT(division_locations.CompareFiles(1e-2));

        // Check cell velocities file
        std::string cell_velocities_file = handler.GetOutputDirectoryFullPath() + "results_from_time_0/cellvelocities.dat";
        NumericFileComparison cell_velocities(cell_velocities_file, "cell_based/test/data/TestOutputNodeAndCellVelocitiesAndDivisionLocations/cellvelocities.dat");
        TS_ASSERT(cell_velocities.CompareFiles(1e-2));

        // Test vtk files exist
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

    void TestOutputNodeVelocitiesWithGhostNodes()
    {
        EXIT_IF_PARALLEL;    // HoneycombMeshGenerator does not work in parallel

        // Create a simple mesh with a surrounding layer of ghost nodes
        HoneycombMeshGenerator generator(3, 3, 1);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Get location indices corresponding to real cells
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Create a differentiated cell for each non-ghost node
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, location_indices.size(), location_indices, p_diff_type);

        for (unsigned i=0; i<cells.size(); i++)
        {
            cells[i]->SetBirthTime(-1.0);
        }

        // Create cell population
        MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);

        // Output Voronoi data
        cell_population.AddPopulationWriter<VoronoiDataWriter>();

        // Record node velocities
        cell_population.AddPopulationWriter<NodeVelocityWriter>();

        // Set up simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestOutputNodeVelocitiesWithGhostNodes");
        simulator.SetEndTime(0.5);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        p_force->SetCutOffLength(1.5);
        simulator.AddForce(p_force);

        // Run simulation
        simulator.Solve();

        // Test node velocities file
        std::string output_directory = "TestOutputNodeVelocitiesWithGhostNodes";
        OutputFileHandler output_file_handler(output_directory, false);

        std::string node_velocities_file = output_file_handler.GetOutputDirectoryFullPath() + "results_from_time_0/nodevelocities.dat";
        NumericFileComparison node_velocities(node_velocities_file, "cell_based/test/data/TestOutputNodeVelocitiesWithGhostNodes/nodevelocities.dat");
        TS_ASSERT(node_velocities.CompareFiles(1e-2));

        // Test vtk files exist
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
    void TestOffLatticeSimulationWithCellDeath()
    {
        EXIT_IF_PARALLEL;    // HoneycombMeshGenerator does not work in parallel

        // Create a simple 2D MeshBasedCellPopulation
        int num_cells_depth = 5;
        int num_cells_width = 5;
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes());

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
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumNodes(),
                         simulator.rGetCellPopulation().GetNumRealCells());

        // For coverage of these 'Get' functions
        TS_ASSERT_EQUALS(simulator.GetNumBirths(), 0u);
        TS_ASSERT_EQUALS(simulator.GetNumDeaths(), 20u);

        // Note that this test used to take an extra time step
        TS_ASSERT_EQUALS(SimulationTime::Instance()->GetTime(), 0.5);
    }

    /**
     * Test a cell-based simulation with multiple cell killers.
     */
    void TestOffLatticeSimulationWithMultipleCellKillers()
    {
        EXIT_IF_PARALLEL;    // HoneycombMeshGenerator does not work in parallel

        // Create a simple 2D MeshBasedCellPopulation
        int num_cells_depth = 5;
        int num_cells_width = 5;
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes(), p_transit_type);

        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestOffLatticeSimulationWithMultipleCellKillers");
        simulator.SetEndTime(0.5);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        p_force->SetCutOffLength(1.5);
        simulator.AddForce(p_force);

        // Add cell killers
        c_vector<double,2> y_normal = zero_vector<double>(2);
        y_normal[1] = -1.0;
        MAKE_PTR_ARGS(PlaneBasedCellKiller<2>, p_killer_1, (&cell_population, zero_vector<double>(2), y_normal) ); // y<0
        simulator.AddCellKiller(p_killer_1);

        // Add cell killers
        c_vector<double,2> x_normal = zero_vector<double>(2);
        x_normal[0] = -1.0;
        MAKE_PTR_ARGS(PlaneBasedCellKiller<2>, p_killer_2, (&cell_population, zero_vector<double>(2), x_normal) ); // x<0
        simulator.AddCellKiller(p_killer_2);

        simulator.SetUpdateCellPopulationRule(true);
        simulator.Solve();

        // Check that the number of nodes is equal to the number of cells
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumNodes(), simulator.rGetCellPopulation().GetNumRealCells());

        //Check that the correct number of cells are killed
        unsigned num_deaths = 6u; // sensitive to changes in random number generation
        TS_ASSERT_EQUALS(simulator.GetNumDeaths(), num_deaths);

        // Now remove the killers and check no more cells are killed
        simulator.RemoveAllCellKillers();
        simulator.SetEndTime(1.0);

        simulator.Solve();
        TS_ASSERT_EQUALS(simulator.GetNumDeaths(), num_deaths);
    }

    /**
     * Test a cell-based simulation with multiple forces.
     */
    void TestOffLatticeSimulationWithMultipleForces()
    {
        EXIT_IF_PARALLEL;    // HoneycombMeshGenerator does not work in parallel

        // Create a simple 2D MeshBasedCellPopulation
        int num_cells_depth = 5;
        int num_cells_width = 5;
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes());

        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestOffLatticeSimulationWithMultipleForces");
        simulator.SetEndTime(0.5);

        // Create some force laws and pass them to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        simulator.AddForce(p_linear_force);

        // Need to set this up for the chemotactic force.
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            double x = p_mesh->GetNode(i)->rGetLocation()[0];
            CellPtr p_cell = cell_population.GetCellUsingLocationIndex(p_mesh->GetNode(i)->GetIndex());
            p_cell->GetCellData()->SetItem("nutrient", x/50.0);
        }

        MAKE_PTR(ChemotacticForce<2>, p_chemotactic_force);
        simulator.AddForce(p_chemotactic_force);

        simulator.Solve();

        // Check that the number of nodes is equal to the number of cells
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumNodes(), simulator.rGetCellPopulation().GetNumRealCells());

        // For coverage of these 'Get' functions
        TS_ASSERT_EQUALS(simulator.GetNumBirths(), 0u);
        TS_ASSERT_EQUALS(simulator.GetNumDeaths(), 0u);

        // Now remove these forces and only have the spring force
        simulator.RemoveAllForces();
        simulator.AddForce(p_linear_force);
        simulator.SetEndTime(1.0);

        simulator.Solve();

        // Check that the number of nodes is equal to the number of cells
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumNodes(), simulator.rGetCellPopulation().GetNumRealCells());
    }

    /**
     * Test a cell-based simulation with variabe rest lengths forces.
     */
    void TestOffLatticeSimulationWithVariableRestLengths()
    {
        EXIT_IF_PARALLEL;    // Cell population output doesn't work in parallel

        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        mesh.Scale(8,8);

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_diff_type);

        // Create a cell population
        MeshBasedCellPopulation<2> cell_population(mesh, cells);
        // Calculate the rest lengths of the springs assuming that the current configuration is in equilibrium.
        cell_population.CalculateRestLengths();

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestOffLatticeSimulationWithVariableRestLengths");
        simulator.SetEndTime(0.5);
        // Turn off remeshing so we only have the same mesh connectivity over time this is needed to use the variable rest length
        simulator.SetUpdateCellPopulationRule(false);

        // Create some force laws and pass them to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        simulator.AddForce(p_linear_force);

        simulator.Solve();

        // Check nothing has moved
        TS_ASSERT_DELTA(simulator.GetNodeLocation(0)[0], 0.0, 1e-6);
        TS_ASSERT_DELTA(simulator.GetNodeLocation(0)[0], 0.0, 1e-6);

        TS_ASSERT_DELTA(simulator.GetNodeLocation(1)[0], 8.0, 1e-6);
        TS_ASSERT_DELTA(simulator.GetNodeLocation(1)[1], 0.0, 1e-6);

        TS_ASSERT_DELTA(simulator.GetNodeLocation(2)[0], 8.0, 1e-6);
        TS_ASSERT_DELTA(simulator.GetNodeLocation(2)[1], 8.0, 1e-6);

        TS_ASSERT_DELTA(simulator.GetNodeLocation(3)[0], 0.0, 1e-6);
        TS_ASSERT_DELTA(simulator.GetNodeLocation(3)[1], 8.0, 1e-6);

        // Create some boundary conditions and pass them to the simulation
        c_vector<double,2> point = zero_vector<double>(2);
        c_vector<double,2> normal = zero_vector<double>(2);
        normal(1) = -1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc1, (&cell_population, point, normal)); // y>0
        simulator.AddCellPopulationBoundaryCondition(p_bc1);
        point(1) = 7.5;
        normal(1) = 1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc2, (&cell_population, point, normal)); // y<7.5
        simulator.AddCellPopulationBoundaryCondition(p_bc2);

        simulator.SetEndTime(1.0);
        simulator.Solve();

        // Check top has shifted down
        TS_ASSERT_DELTA(simulator.GetNodeLocation(0)[0], -0.0207, 1e-3);
        TS_ASSERT_DELTA(simulator.GetNodeLocation(0)[1], 0.0, 1e-6);

        TS_ASSERT_DELTA(simulator.GetNodeLocation(1)[0], 8.0207, 1e-3);
        TS_ASSERT_DELTA(simulator.GetNodeLocation(1)[1], 0.0, 1e-6);

        TS_ASSERT_DELTA(simulator.GetNodeLocation(2)[0], 8.1097, 1e-3);
        TS_ASSERT_DELTA(simulator.GetNodeLocation(2)[1], 7.5, 1e-6);

        TS_ASSERT_DELTA(simulator.GetNodeLocation(3)[0], -0.1097, 1e-3);
        TS_ASSERT_DELTA(simulator.GetNodeLocation(3)[1], 7.5, 1e-6);

        simulator.RemoveAllCellPopulationBoundaryConditions();
        simulator.AddCellPopulationBoundaryCondition(p_bc1);
        simulator.SetEndTime(10.0);
        simulator.Solve();

        // Check has relaxed back to original shape
        TS_ASSERT_DELTA(simulator.GetNodeLocation(0)[0], 0.0, 1e-3);
        TS_ASSERT_DELTA(simulator.GetNodeLocation(0)[1], 0.0, 1e-3);

        TS_ASSERT_DELTA(simulator.GetNodeLocation(1)[0], 8.0, 1e-3);
        TS_ASSERT_DELTA(simulator.GetNodeLocation(1)[1], 0.0, 1e-3);

        TS_ASSERT_DELTA(simulator.GetNodeLocation(2)[0], 8.0, 1e-3);
        TS_ASSERT_DELTA(simulator.GetNodeLocation(2)[1], 8.0, 1e-2);

        TS_ASSERT_DELTA(simulator.GetNodeLocation(3)[0], 0.0, 1e-3);
        TS_ASSERT_DELTA(simulator.GetNodeLocation(3)[1], 8.0, 1e-2);
    }

    /**
     * This tests running a simulation in 3d with a 2d mesh (see #2112)
     */
    void TestOffLatticeSimulationWith2dMeshIn3d()
    {
        EXIT_IF_PARALLEL;    // Cell population output doesn't work in parallel

        // Load mesh
        TrianglesMeshReader<2,3> mesh_reader("cell_based/test/data/Square2dMeshIn3d/Square2dMeshIn3d");
        MutableMesh<2,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(),p_diff_type);

        // Create a cell population
        MeshBasedCellPopulation<2,3> cell_population(mesh, cells);
        cell_population.SetWriteVtkAsPoints(true);

        // Set up cell-based simulation
        OffLatticeSimulation<2,3> simulator(cell_population);
        simulator.SetOutputDirectory("TestOffLatticeSimulationWith2dMeshIn3d");
        simulator.SetEndTime(1.0);

        // Create a force law and pass it to the simulation
        typedef GeneralisedLinearSpringForce<2,3> Force;
        MAKE_PTR(Force, p_force);
        p_force->SetCutOffLength(1.5);
        simulator.AddForce(p_force);

        // Stops remeshing as not possible for 2d in 3d meshes
        simulator.SetUpdateCellPopulationRule(false);

        simulator.Solve();

        // Check that nodes are all sat at resting length (1.0) apart
        TS_ASSERT_DELTA(norm_2(simulator.rGetCellPopulation().GetNode(0)->rGetLocation()-simulator.rGetCellPopulation().GetNode(1)->rGetLocation()),1.0,1e-5);
        TS_ASSERT_DELTA(norm_2(simulator.rGetCellPopulation().GetNode(1)->rGetLocation()-simulator.rGetCellPopulation().GetNode(2)->rGetLocation()),1.0,1e-5);
        TS_ASSERT_DELTA(norm_2(simulator.rGetCellPopulation().GetNode(2)->rGetLocation()-simulator.rGetCellPopulation().GetNode(0)->rGetLocation()),1.0,1e-5);
    }

    /**
     * This tests archiving a simulation in 3d with a 2d mesh (see #2112)
     */
    void TestArchivingOffLatticeSimulationWith2dMeshIn3d()
    {
        EXIT_IF_PARALLEL;    // Cell population output doesn't work in parallel

        // Load mesh
        TrianglesMeshReader<2,3> mesh_reader("cell_based/test/data/Square2dMeshIn3d/Square2dMeshIn3d");
        MutableMesh<2,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Create cells (we use a NoCellCycleModel here for simplicity, since there is no proliferation)
        std::vector<CellPtr> cells;
        CellsGenerator<NoCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes());

        // Create a cell population
        MeshBasedCellPopulation<2,3> cell_population(mesh, cells);
        cell_population.SetWriteVtkAsPoints(true);

        // Add a writer
        cell_population.AddCellWriter<CellIdWriter>();

        // Set up cell-based simulation
        OffLatticeSimulation<2,3> simulator(cell_population);
        std::string output_dir = "TestOffLatticeSimulationWith2dMeshIn3dArchive";
        simulator.SetOutputDirectory(output_dir);
        simulator.SetEndTime(1.0);

        // Create a force law and pass it to the simulation
        typedef GeneralisedLinearSpringForce<2,3> Force;
        MAKE_PTR(Force, p_force);
        p_force->SetCutOffLength(1.5);
        simulator.AddForce(p_force);

        // Stop remeshing, as this is not possible for 2d in 3d meshes
        simulator.SetUpdateCellPopulationRule(false);

        CellBasedSimulationArchiver<2,OffLatticeSimulation<2,3>, 3>::Save(&simulator);

        double start_time = 0.;
        OffLatticeSimulation<2,3>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2,3>, 3 >::Load(output_dir, start_time);

        p_simulator->Solve();

        // Check that nodes are all sat at resting length  (1.0) apart.
        TS_ASSERT_DELTA(norm_2(p_simulator->rGetCellPopulation().GetNode(0)->rGetLocation()-p_simulator->rGetCellPopulation().GetNode(1)->rGetLocation()),1.0,1e-5);
        TS_ASSERT_DELTA(norm_2(p_simulator->rGetCellPopulation().GetNode(1)->rGetLocation()-p_simulator->rGetCellPopulation().GetNode(2)->rGetLocation()),1.0,1e-5);
        TS_ASSERT_DELTA(norm_2(p_simulator->rGetCellPopulation().GetNode(2)->rGetLocation()-p_simulator->rGetCellPopulation().GetNode(0)->rGetLocation()),1.0,1e-5);

        // Avoid memory leak
        delete p_simulator;
    }

    /**
     * Test a cell-based simulation with a periodic mesh.
     */
    void TestOffLatticeSimulationWithPeriodicMesh()
    {
        EXIT_IF_PARALLEL;    // HoneycombMeshGenerator does not work in parallel

        // Create a simple mesh
        int cells_up = 6;
        int cells_across = 6;

        CylindricalHoneycombMeshGenerator generator(cells_across, cells_up, 0);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
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
        FileComparison( results_dir + "results.vizsetup", "cell_based/test/data/TestOffLatticeSimulationWithPeriodicMesh/results.vizsetup").CompareFiles();
    }

    /**
     * Test a cell-based simulation with multiple boundary conditions. y<2 and y>0
     */
    void TestOffLatticeSimulationWithMultipleCellBoundaryConditions()
    {
        EXIT_IF_PARALLEL;    // HoneycombMeshGenerator does not work in parallel

        // Create a simple mesh
        int num_cells_depth = 5;
        int num_cells_width = 5;
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
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

    void TestOffLatticeSimulationWithStoppingEvent()
    {
        EXIT_IF_PARALLEL;    // HoneycombMeshGenerator does not work in parallel

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
            unsigned generation;
            double y = 0.0;

            if (std::find(location_indices.begin(), location_indices.end(), i) != location_indices.end())
            {
                y = p_mesh->GetNode(i)->GetPoint().rGetLocation()[1];
            }

            FixedG1GenerationalCellCycleModel* p_cell_cycle_model = new FixedG1GenerationalCellCycleModel;
            p_cell_cycle_model->SetDimension(2);

            double typical_transit_cycle_time = p_cell_cycle_model->GetAverageTransitCellCycleTime();
            double typical_stem_cycle_time = p_cell_cycle_model->GetAverageStemCellCycleTime();

            boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
            boost::shared_ptr<AbstractCellProperty> p_stem_type(CellPropertyRegistry::Instance()->Get<StemCellProliferativeType>());
            boost::shared_ptr<AbstractCellProperty> p_transit_type(CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>());
            boost::shared_ptr<AbstractCellProperty> p_diff_type(CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>());

            CellPtr p_cell(new Cell(p_state, p_cell_cycle_model));

            double birth_time = -p_random_num_gen->ranf();

            if (y <= 0.3)
            {
                p_cell->SetCellProliferativeType(p_stem_type);
                generation = 0;
                birth_time *= typical_stem_cycle_time; // hours
            }
            else if (y < 2.0)
            {
                p_cell->SetCellProliferativeType(p_transit_type);
                generation = 1;
                birth_time *= typical_transit_cycle_time; // hours
            }
            else if (y < 3.0)
            {
                p_cell->SetCellProliferativeType(p_transit_type);
                generation = 2;
                birth_time *= typical_transit_cycle_time; // hours
            }
            else if (y < 4.0)
            {
                p_cell->SetCellProliferativeType(p_transit_type);
                generation = 3;
                birth_time *= typical_transit_cycle_time; // hours
            }
            else
            {
                if (p_cell_cycle_model->CanCellTerminallyDifferentiate())
                {
                    p_cell->SetCellProliferativeType(p_diff_type);
                }
                else
                {
                    p_cell->SetCellProliferativeType(p_transit_type);
                }
                generation = 4;
                birth_time *= typical_transit_cycle_time; // hours
            }

            p_cell_cycle_model->SetGeneration(generation);
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

    void TestApoptosisSpringLengths()
    {
        EXIT_IF_PARALLEL;    // HoneycombMeshGenerator does not work in parallel

        unsigned num_cells_depth = 2;
        unsigned num_cells_width = 2;

        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 2);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, location_indices.size(), p_transit_type);

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

    void TestOffLatticeSimulationParameterOutputMethods()
    {
        EXIT_IF_PARALLEL;    // HoneycombMeshGenerator does not work in parallel

        // Create a simple mesh
        HoneycombMeshGenerator generator(5, 5, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Set up cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumNodes());

        // Create a cell population
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetEndTime(0.5);

        // We call SetupSolve() here to set up the default numerical method, which is otherwise NULL
        simulator.SetupSolve();

        TS_ASSERT_EQUALS(simulator.GetIdentifier(), "OffLatticeSimulation-2-2");

        ///\todo (#1453) should have forces and cell killer included here to make it a better test

        std::string output_directory = "TestOffLatticeSimulationOutputParameters";
        OutputFileHandler output_file_handler(output_directory, false);
        out_stream parameter_file = output_file_handler.OpenOutputFile("cell_based_sim_results.parameters");
        simulator.OutputSimulationParameters(parameter_file);
        parameter_file->close();

        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();
        FileComparison(results_dir + "cell_based_sim_results.parameters", "cell_based/test/data/TestOffLatticeSimulationOutputParameters/cell_based_sim_results.parameters").CompareFiles();

        simulator.SetOutputDirectory("TestOffLatticeSimulationOutputParameters");
        simulator.OutputSimulationSetup();

        FileFinder generated(results_dir + "results.parameters", RelativeTo::Absolute);
        FileFinder reference("cell_based/test/data/TestOffLatticeSimulationOutputParameters/results.parameters", RelativeTo::ChasteSourceRoot);
        FileComparison comparer(generated,reference);
        comparer.IgnoreLinesContaining("CellPopulation");
        TS_ASSERT(comparer.CompareFiles());

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
    void Test1dOffLatticeSimulation()
    {
        // Create a 1D MeshBasedCellPopulation
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements");
        MutableMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Set up cells so that cell 10 divides at time t=0.5, cell 9 at time t=1.5, etc
        for (unsigned i=0; i<cells.size(); i++)
        {
            cells[i]->SetBirthTime(-13.5 - i);
        }

        MeshBasedCellPopulation<1> cell_population(mesh, cells);

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
        simulator.Solve();

        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), initial_num_cells + 1);
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumNodes(), initial_num_nodes + 1);
        TS_ASSERT_EQUALS((static_cast<MeshBasedCellPopulation<1>* >(&(simulator.rGetCellPopulation())))->rGetMesh().GetNumElements(), initial_num_elements + 1);
    }

    void TestSettingEndTimeIssue()
    {
        EXIT_IF_PARALLEL;    // HoneycombMeshGenerator does not work in parallel

        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(0.1, 1);

        HoneycombMeshGenerator generator(2, 2, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Create some cells
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes(), p_diff_type);

        // Create a mesh-based cell population
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestSettingEndTimeIssue");

        TS_ASSERT_THROWS_THIS(simulator.Solve(), "SetEndTime has not yet been called.");

        simulator.SetEndTime(1.0);
        TS_ASSERT_THROWS_THIS(simulator.Solve(),
            "End time and number of timesteps already setup. You should not use SimulationTime::SetEndTimeAndNumberOfTimeSteps in cell-based tests.");
    }

    void TestCellProliferativeTypeCounts()
    {
        EXIT_IF_PARALLEL;    // HoneycombMeshGenerator does not work in parallel

        // Create a simple 2D MeshBasedCellPopulation
        HoneycombMeshGenerator generator(5, 5, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes());

        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();
        cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();

        // An ordering must be specified for cell mutation states and cell proliferative types
        cell_population.SetDefaultCellMutationStateAndProliferativeTypeOrdering();

        // Test we have the correct cell mutation state counts
        CellPropertyRegistry* p_registry_before_solve = cell_population.Begin()->rGetCellPropertyCollection().GetCellPropertyRegistry();
        TS_ASSERT_EQUALS(p_registry_before_solve->Get<WildTypeCellMutationState>()->GetCellCount(), 25u);
        TS_ASSERT_EQUALS(p_registry_before_solve->Get<ApcOneHitCellMutationState>()->GetCellCount(), 0u);
        TS_ASSERT_EQUALS(p_registry_before_solve->Get<ApcTwoHitCellMutationState>()->GetCellCount(), 0u);
        TS_ASSERT_EQUALS(p_registry_before_solve->Get<BetaCateninOneHitCellMutationState>()->GetCellCount(), 0u);

        std::vector<unsigned> mutation_state_count_before_solve = cell_population.GetCellMutationStateCount();
        TS_ASSERT_EQUALS(mutation_state_count_before_solve.size(), 4u);
        TS_ASSERT_EQUALS(mutation_state_count_before_solve[0], 25u);
        TS_ASSERT_EQUALS(mutation_state_count_before_solve[1], 0u);
        TS_ASSERT_EQUALS(mutation_state_count_before_solve[2], 0u);
        TS_ASSERT_EQUALS(mutation_state_count_before_solve[3], 0u);

        // Test we have the correct cell proliferative type counts
        TS_ASSERT_EQUALS(p_registry_before_solve->Get<StemCellProliferativeType>()->GetCellCount(), 25u);
        TS_ASSERT_EQUALS(p_registry_before_solve->Get<TransitCellProliferativeType>()->GetCellCount(), 0u);
        TS_ASSERT_EQUALS(p_registry_before_solve->Get<DifferentiatedCellProliferativeType>()->GetCellCount(), 0u);
        TS_ASSERT_EQUALS(p_registry_before_solve->Get<DefaultCellProliferativeType>()->GetCellCount(), 0u);

        std::vector<unsigned> prolif_type_count_before_solve = cell_population.GetCellProliferativeTypeCount();
        TS_ASSERT_EQUALS(prolif_type_count_before_solve.size(), 4u);
        TS_ASSERT_EQUALS(prolif_type_count_before_solve[0], 25u);
        TS_ASSERT_EQUALS(prolif_type_count_before_solve[1], 0u);
        TS_ASSERT_EQUALS(prolif_type_count_before_solve[2], 0u);
        TS_ASSERT_EQUALS(prolif_type_count_before_solve[3], 0u);

        // Set up simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestCellProliferativeTypeCounts");
        simulator.SetEndTime(0.5);

        // Run simulation
        simulator.Solve();

        // Test we still have the correct cell mutation state counts
        CellPropertyRegistry* p_registry_after_solve = cell_population.Begin()->rGetCellPropertyCollection().GetCellPropertyRegistry();
        TS_ASSERT_EQUALS(p_registry_after_solve->Get<WildTypeCellMutationState>()->GetCellCount(), 25u);
        TS_ASSERT_EQUALS(p_registry_after_solve->Get<ApcOneHitCellMutationState>()->GetCellCount(), 0u);
        TS_ASSERT_EQUALS(p_registry_after_solve->Get<ApcTwoHitCellMutationState>()->GetCellCount(), 0u);
        TS_ASSERT_EQUALS(p_registry_after_solve->Get<BetaCateninOneHitCellMutationState>()->GetCellCount(), 0u);

        std::vector<unsigned> mutation_state_count_after_solve = simulator.rGetCellPopulation().GetCellMutationStateCount();
        TS_ASSERT_EQUALS(mutation_state_count_after_solve.size(), 4u);
        TS_ASSERT_EQUALS(mutation_state_count_after_solve[0], 25u);
        TS_ASSERT_EQUALS(mutation_state_count_after_solve[1], 0u);
        TS_ASSERT_EQUALS(mutation_state_count_after_solve[2], 0u);
        TS_ASSERT_EQUALS(mutation_state_count_after_solve[3], 0u);

        // Test we still have the correct cell proliferative type counts
        TS_ASSERT_EQUALS(p_registry_after_solve->Get<StemCellProliferativeType>()->GetCellCount(), 25u);
        TS_ASSERT_EQUALS(p_registry_after_solve->Get<TransitCellProliferativeType>()->GetCellCount(), 0u);
        TS_ASSERT_EQUALS(p_registry_after_solve->Get<DifferentiatedCellProliferativeType>()->GetCellCount(), 0u);
        TS_ASSERT_EQUALS(p_registry_after_solve->Get<DefaultCellProliferativeType>()->GetCellCount(), 0u);

        std::vector<unsigned> prolif_type_count_after_solve = simulator.rGetCellPopulation().GetCellProliferativeTypeCount();
        TS_ASSERT_EQUALS(prolif_type_count_after_solve.size(), 4u);
        TS_ASSERT_EQUALS(prolif_type_count_after_solve[0], 25u);
        TS_ASSERT_EQUALS(prolif_type_count_after_solve[1], 0u);
        TS_ASSERT_EQUALS(prolif_type_count_after_solve[2], 0u);
        TS_ASSERT_EQUALS(prolif_type_count_after_solve[3], 0u);
    }

    /**
     * Test a cell-based simulation with modifiers.
     *
     * Test to check we can add and remove them
     */
    void TestOffLatticeSimulationWithModifiers()
    {
        EXIT_IF_PARALLEL;    // HoneycombMeshGenerator does not work in parallel

        // Create a simple 2D MeshBasedCellPopulation
        int num_cells_depth = 5;
        int num_cells_width = 5;
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes());

        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestOffLatticeSimulationWithModifiers");
        simulator.SetEndTime(0.5);

        // Create a force law and pass them to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        simulator.AddForce(p_linear_force);

        // Add a modifier
        MAKE_PTR(VolumeTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        // Check there is a modifier and it's the correct type
        std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<2> > >::iterator iter = simulator.GetSimulationModifiers()->begin();
        TS_ASSERT(boost::static_pointer_cast<VolumeTrackingModifier<2> >(*iter));
        TS_ASSERT_EQUALS(simulator.GetSimulationModifiers()->size(), 1u);

        simulator.Solve();

        // Check that the number of nodes is equal to the number of cells
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumNodes(), simulator.rGetCellPopulation().GetNumRealCells());
    }

    void TestWriterIteratorsWithCellDeath()
    {
        /*
         * This test was added because of the issue raised by #2422 and #2689,
         * where cell mutation state labels appear to flip between cells when
         * the output of a cell-based simulation is visualized.
         *
         * This issue occurred because, for a MeshBasedCellPopulation, the
         * NodeLocationWriter wrote each line of its file using the NodeIterator
         * on pCellPopulation->rGetMesh(), while the CellProliferativeTypesWriter
         * writes each line using the AbstractCellPopulation iterator.
         */
        EXIT_IF_PARALLEL;

        HoneycombMeshGenerator generator(10,10);
        MutableMesh<2,2>* p_mesh = generator.GetCircularMesh(5);

        MAKE_PTR(ApcOneHitCellMutationState, p_mute_state);
        MAKE_PTR(WildTypeCellMutationState, p_wild_state);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);

        std::vector<CellPtr> cells;

        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            UniformCellCycleModel* p_model = new UniformCellCycleModel();
            p_model->SetDimension(2);
            CellPtr p_cell(new Cell(p_wild_state, p_model));
            p_cell->SetCellProliferativeType(p_transit_type);
            cells.push_back(p_cell);
        }

        cells[27]->SetMutationState(p_mute_state);
        cells[34]->SetMutationState(p_mute_state);
        cells[55]->SetMutationState(p_mute_state);

        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestWriterIteratorsWithCellDeath");
        simulator.SetSamplingTimestepMultiple(12);
        simulator.SetEndTime(15);

        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetCutOffLength(3);
        simulator.AddForce(p_linear_force);

        MAKE_PTR_ARGS(RandomCellKiller<2>, p_killer, (&cell_population, 0.01));
        simulator.AddCellKiller(p_killer);

        simulator.Solve();

        FileFinder generated_type_file("TestWriterIteratorsWithCellDeath/results_from_time_0/results.vizcelltypes", RelativeTo::ChasteTestOutput);
        FileFinder generated_node_file("TestWriterIteratorsWithCellDeath/results_from_time_0/results.viznodes", RelativeTo::ChasteTestOutput);

        FileFinder reference_type_file("cell_based/test/data/TestWriterIteratorsWithCellDeath/results.vizcelltypes",RelativeTo::ChasteSourceRoot);
        FileFinder reference_node_file("cell_based/test/data/TestWriterIteratorsWithCellDeath/results.viznodes",RelativeTo::ChasteSourceRoot);

        FileComparison type_files(generated_type_file,reference_type_file);
        FileComparison node_files(generated_node_file,reference_node_file);

        TS_ASSERT(type_files.CompareFiles());
        TS_ASSERT(node_files.CompareFiles());
    }

    void TestOffLatticeSimulationWithAdaptiveTimestep()
    {
        EXIT_IF_PARALLEL;    // HoneycombMeshGenerator does not work in parallel

        // Create a simple 2D MeshBasedCellPopulation
        HoneycombMeshGenerator generator(3, 3, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes(), p_transit_type);

        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Use smaller movement threshold to maintain smooth motion
        cell_population.SetAbsoluteMovementThreshold(0.1);

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestOffLatticeSimulationWithAdaptivity");
        simulator.SetEndTime(5.0);
        simulator.SetDt(0.1);

        // Pass an adaptive numerical method to the simulation
        boost::shared_ptr<AbstractNumericalMethod<2,2> > p_method(new ForwardEulerNumericalMethod<2,2>());
        p_method->SetUseAdaptiveTimestep(true);
        simulator.SetNumericalMethod(p_method);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        p_force->SetCutOffLength(1.5);
        simulator.AddForce(p_force);

        // Record division locations
        simulator.SetOutputDivisionLocations(true);

        // Run simulation
        simulator.Solve();

        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), 14u);

        // Check cells have moved to the correct location
        TS_ASSERT_DELTA(simulator.rGetCellPopulation().rGetMesh().GetNode(0)->rGetLocation()[0], 0.3906,1e-4);
        TS_ASSERT_DELTA(simulator.rGetCellPopulation().rGetMesh().GetNode(0)->rGetLocation()[1], -0.1782,1e-4);
    }
};

#endif /*TESTOFFLATTICESIMULATION_HPP_*/
