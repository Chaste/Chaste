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

#ifndef TESTONLATTICESIMULATIONWITHPOTTSBASEDCELLPOPULATION_HPP_
#define TESTONLATTICESIMULATIONWITHPOTTSBASEDCELLPOPULATION_HPP_

#include <cxxtest/TestSuite.h>

#include <algorithm> //For max_element of vector

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "HoneycombMeshGenerator.hpp"
#include "PottsMeshGenerator.hpp"
#include "CellsGenerator.hpp"
#include "CellLabel.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "NoCellCycleModel.hpp"
#include "WildTypeCellMutationState.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "VolumeConstraintPottsUpdateRule.hpp"
#include "AdhesionPottsUpdateRule.hpp"
#include "DifferentialAdhesionPottsUpdateRule.hpp"
#include "ChemotaxisPottsUpdateRule.hpp"
#include "PlaneBasedCellKiller.hpp"
#include "OnLatticeSimulation.hpp"
#include "OffLatticeSimulation.hpp"
#include "Warnings.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "SmartPointers.hpp"
#include "CellMutationStatesCountWriter.hpp"
#include "PetscSetupAndFinalize.hpp"

class TestOnLatticeSimulationWithPottsBasedCellPopulation : public AbstractCellBasedWithTimingsTestSuite
{
private:
    void RandomlyLabelCells(std::vector<CellPtr>& rCells, boost::shared_ptr<AbstractCellProperty> pLabel, double labelledRatio)
    {
        for (unsigned i = 0; i<rCells.size(); i++)
        {
            if (RandomNumberGenerator::Instance()->ranf() < labelledRatio)
            {
                rCells[i]->AddCellProperty(pLabel);
            }
        }
    }

public:

    void TestOnLatticeSimulationExceptions()
    {
        EXIT_IF_PARALLEL;    // HoneycombMeshGenerator doesn't work in parallel.

        // Create a simple tetrahedral mesh
        HoneycombMeshGenerator generator(3, 3, 0);
        TetrahedralMesh<2,2>* p_generating_mesh = generator.GetMesh();

        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh, 1.5);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes());

        // Create a node-based cell population
        NodeBasedCellPopulation<2> node_based_cell_population(mesh, cells);

        TS_ASSERT_THROWS_THIS(OnLatticeSimulation<2> simulator(node_based_cell_population),
            "OnLatticeSimulations require a subclass of AbstractOnLatticeCellPopulation.");
    }

    void TestMoreOnLatticeSimulationExceptions()
    {
        // Create a simple 2D PottsMesh
        PottsMeshGenerator<2> generator(6, 2, 2, 6, 2, 2);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_diff_type);

        // Create cell population
        PottsBasedCellPopulation<2> potts_based_cell_population(*p_mesh, cells);

        // Try to set up off lattice simulation
        TS_ASSERT_THROWS_THIS(OffLatticeSimulation<2> simulator(potts_based_cell_population),
            "OffLatticeSimulations require a subclass of AbstractOffLatticeCellPopulation.");
    }

    void TestPottsMonolayerWithNoBirthOrDeath()
    {
        EXIT_IF_PARALLEL;    // Potts simulations don't work in parallel because they depend on NodesOnlyMesh for writing.

        // Create a simple 2D PottsMesh
        PottsMeshGenerator<2> generator(6, 2, 2, 6, 2, 2);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells (we use a NoCellCycleModel here for simplicity, since there is no proliferation)
        std::vector<CellPtr> cells;
        CellsGenerator<NoCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());

        // Create cell population
        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up cell-based simulation
        OnLatticeSimulation<2> simulator(cell_population);

        // Test that the member mDt has been initialised correctly
        TS_ASSERT_DELTA(simulator.GetDt(), 0.1, 1e-6);

        simulator.SetOutputDirectory("TestSimplePottsMonolayer");
        simulator.SetEndTime(0.1);

        // Create update rules and pass to the simulation
        MAKE_PTR(VolumeConstraintPottsUpdateRule<2>, p_volume_constraint_update_rule);
        simulator.AddUpdateRule(p_volume_constraint_update_rule);
        MAKE_PTR(AdhesionPottsUpdateRule<2>, p_adhesion_update_rule);
        simulator.AddUpdateRule(p_adhesion_update_rule);

        // Run simulation
        simulator.Solve();

        // Check that the same number of cells
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), 4u);

        // Test no births or deaths
        TS_ASSERT_EQUALS(simulator.GetNumBirths(), 0u);
        TS_ASSERT_EQUALS(simulator.GetNumDeaths(), 0u);

        // Now remove the update rules and check nothing happens when the simulator runs again
        simulator.RemoveAllUpdateRules();
        simulator.SetEndTime(0.2);
        simulator.Solve();

        // Check that the same number of cells
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), 4u);
    }

    void TestPottsMonolayerWithNonRandomSweep()
    {
        EXIT_IF_PARALLEL;    // Potts simulations don't work in parallel because they depend on NodesOnlyMesh for writing.

        // Create a simple 2D PottsMesh
        PottsMeshGenerator<2> generator(6, 2, 2, 6, 2, 2);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_diff_type);

        // Create cell population
        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.SetUpdateNodesInRandomOrder(false);

        // Set up cell-based simulation
        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestSimplePottsMonolayerWithRandomSweep");
        simulator.SetEndTime(0.1);

        // Create update rules and pass to the simulation
        MAKE_PTR(VolumeConstraintPottsUpdateRule<2>, p_volume_constraint_update_rule);
        simulator.AddUpdateRule(p_volume_constraint_update_rule);
        MAKE_PTR(AdhesionPottsUpdateRule<2>, p_adhesion_update_rule);
        simulator.AddUpdateRule(p_adhesion_update_rule);

        // Run simulation
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());
    }

    void TestPottsMonolayerWithDeath()
    {
        EXIT_IF_PARALLEL;    // Potts simulations don't work in parallel because they depend on NodesOnlyMesh for writing.

        // Create a simple 2D PottsMesh
        PottsMeshGenerator<2> generator(16, 4, 4, 24, 8, 2);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_diff_type);

        // Create cell population
        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up cell-based simulation
        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestPottsMonolayerWithDeath");
        simulator.SetDt(0.1);
        simulator.SetEndTime(1.0);

        // Create update rules and pass to the simulation
        MAKE_PTR(VolumeConstraintPottsUpdateRule<2>, p_volume_constraint_update_rule);
        simulator.AddUpdateRule(p_volume_constraint_update_rule);
        MAKE_PTR(AdhesionPottsUpdateRule<2>, p_adhesion_update_rule);
        simulator.AddUpdateRule(p_adhesion_update_rule);

        // Create cell killer and pass in to simulation
        c_vector<double,2> point = zero_vector<double>(2);
        point(1) = 16.0;
        c_vector<double,2> normal = zero_vector<double>(2);
        normal(1) = 1.0;
        MAKE_PTR_ARGS(PlaneBasedCellKiller<2>, p_killer, (&cell_population, point, normal));
        simulator.AddCellKiller(p_killer);

        // Run simulation
        simulator.Solve();

        // Check the number of cells
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), 18u);

        // Test no births or deaths
        TS_ASSERT_EQUALS(simulator.GetNumBirths(), 0u);
        TS_ASSERT_EQUALS(simulator.GetNumDeaths(), 14u);
    }

    void TestPottsMonolayerWithBirth()
    {
        EXIT_IF_PARALLEL;    // Potts simulations don't work in parallel because they depend on NodesOnlyMesh for writing.

        // Create a simple 2D PottsMesh
        PottsMeshGenerator<2> generator(8, 1, 4, 10, 1, 4);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        CellsGenerator<UniformG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_stem_type);

        // Create cell population
        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up cell-based simulation
        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestPottsMonolayerWithBirth");
        simulator.SetDt(0.1);
        simulator.SetEndTime(20);
        simulator.SetSamplingTimestepMultiple(20);

        // Create update rules and pass to the simulation
        MAKE_PTR(VolumeConstraintPottsUpdateRule<2>, p_volume_constraint_update_rule);
        simulator.AddUpdateRule(p_volume_constraint_update_rule);
        MAKE_PTR(AdhesionPottsUpdateRule<2>, p_adhesion_update_rule);
        simulator.AddUpdateRule(p_adhesion_update_rule);

        // Run simulation
        simulator.Solve();

        // Check the number of cells
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), 2u);

        // Test no deaths and some births
        TS_ASSERT_EQUALS(simulator.GetNumBirths(), 1u);
        TS_ASSERT_EQUALS(simulator.GetNumDeaths(), 0u);
    }

    void TestPottsMonolayerCellSorting()
    {
        EXIT_IF_PARALLEL;    // Potts simulations don't work in parallel because they depend on NodesOnlyMesh for writing.

        // Create a simple 2D PottsMesh
        PottsMeshGenerator<2> generator(30, 4, 4, 30, 4, 4);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_diff_type);

        // Randomly label some cells
        boost::shared_ptr<AbstractCellProperty> p_label(CellPropertyRegistry::Instance()->Get<CellLabel>());
        RandomlyLabelCells(cells, p_label, 0.5);

        // Create cell population
        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>(); // So outputs the labelled cells

        // Set up cell-based simulation
        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestPottsCellSorting");
        simulator.SetDt(0.1);
        simulator.SetEndTime(1.0);

        // Create update rules and pass to the simulation
        MAKE_PTR(VolumeConstraintPottsUpdateRule<2>, p_volume_constraint_update_rule);
        p_volume_constraint_update_rule->SetMatureCellTargetVolume(16);
        p_volume_constraint_update_rule->SetDeformationEnergyParameter(0.2);
        simulator.AddUpdateRule(p_volume_constraint_update_rule);

        MAKE_PTR(DifferentialAdhesionPottsUpdateRule<2>, p_differential_adhesion_update_rule);
        p_differential_adhesion_update_rule->SetLabelledCellLabelledCellAdhesionEnergyParameter(0.16);
        p_differential_adhesion_update_rule->SetLabelledCellCellAdhesionEnergyParameter(0.11);
        p_differential_adhesion_update_rule->SetCellCellAdhesionEnergyParameter(0.02);
        p_differential_adhesion_update_rule->SetLabelledCellBoundaryAdhesionEnergyParameter(0.16);
        p_differential_adhesion_update_rule->SetCellBoundaryAdhesionEnergyParameter(0.16);
        simulator.AddUpdateRule(p_differential_adhesion_update_rule);

        // Run simulation
        simulator.Solve();

        // Check that the same number of cells
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), 16u);

        // Test no births or deaths
        TS_ASSERT_EQUALS(simulator.GetNumBirths(), 0u);
        TS_ASSERT_EQUALS(simulator.GetNumDeaths(), 0u);
#ifdef CHASTE_VTK
        //Test that VTK writer has produced some files
        OutputFileHandler handler("TestPottsCellSorting", false);
        std::string results_dir = handler.GetOutputDirectoryFullPath();
        // Initial condition file
        FileFinder vtk_file_1(results_dir + "results_from_time_0/results_0.vtu", RelativeTo::Absolute);
        TS_ASSERT(vtk_file_1.Exists());
        FileFinder vtk_file_2(results_dir + "results_from_time_0/outlines_0.vtu", RelativeTo::Absolute);
        TS_ASSERT(vtk_file_2.Exists());

        // Final file
        FileFinder vtk_file_3(results_dir + "results_from_time_0/results_10.vtu", RelativeTo::Absolute);
        TS_ASSERT(vtk_file_3.Exists());
        FileFinder vtk_file_4(results_dir + "results_from_time_0/outlines_10.vtu", RelativeTo::Absolute);
        TS_ASSERT(vtk_file_4.Exists());

 #endif //CHASTE_VTK
    }

    void TestPottsMonolayerCellSortingPeriodic()
    {
        EXIT_IF_PARALLEL;    // Potts simulations don't work in parallel because they depend on NodesOnlyMesh for writing.

        // Create a simple 2D PottsMesh
        PottsMeshGenerator<2> generator(20, 5, 4, 20, 5, 4, 1, 1, 1, false, true, true); // Periodic in x and y
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_diff_type);

        // Randomly label some cells
        boost::shared_ptr<AbstractCellProperty> p_label(CellPropertyRegistry::Instance()->Get<CellLabel>());
        RandomlyLabelCells(cells, p_label, 0.5);

        // Create cell population
        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>(); // So outputs the labelled cells

        // Set up cell-based simulation
        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestPottsCellSortingPeriodic");
        simulator.SetDt(0.1);
        simulator.SetEndTime(1.0);

        // Create update rules and pass to the simulation
        MAKE_PTR(VolumeConstraintPottsUpdateRule<2>, p_volume_constraint_update_rule);
        p_volume_constraint_update_rule->SetMatureCellTargetVolume(16);
        p_volume_constraint_update_rule->SetDeformationEnergyParameter(0.2);
        simulator.AddUpdateRule(p_volume_constraint_update_rule);

        MAKE_PTR(DifferentialAdhesionPottsUpdateRule<2>, p_differential_adhesion_update_rule);
        p_differential_adhesion_update_rule->SetLabelledCellLabelledCellAdhesionEnergyParameter(0.16);
        p_differential_adhesion_update_rule->SetLabelledCellCellAdhesionEnergyParameter(0.11);
        p_differential_adhesion_update_rule->SetCellCellAdhesionEnergyParameter(0.02);
        p_differential_adhesion_update_rule->SetLabelledCellBoundaryAdhesionEnergyParameter(0.16);
        p_differential_adhesion_update_rule->SetCellBoundaryAdhesionEnergyParameter(0.16);
        simulator.AddUpdateRule(p_differential_adhesion_update_rule);

        // Run simulation
        simulator.Solve();

        // Check that the same number of cells
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), 25u);

        // Test no births or deaths
        TS_ASSERT_EQUALS(simulator.GetNumBirths(), 0u);
        TS_ASSERT_EQUALS(simulator.GetNumDeaths(), 0u);
#ifdef CHASTE_VTK
        //Test that VTK writer has produced some files
        OutputFileHandler handler("TestPottsCellSortingPeriodic", false);
        std::string results_dir = handler.GetOutputDirectoryFullPath();
        // Initial condition file
        FileFinder vtk_file_1(results_dir + "results_from_time_0/results_0.vtu", RelativeTo::Absolute);
        TS_ASSERT(vtk_file_1.Exists());
        FileFinder vtk_file_2(results_dir + "results_from_time_0/outlines_0.vtu", RelativeTo::Absolute);
        TS_ASSERT(vtk_file_2.Exists());

        // Final file
        FileFinder vtk_file_3(results_dir + "results_from_time_0/results_10.vtu", RelativeTo::Absolute);
        TS_ASSERT(vtk_file_3.Exists());
        FileFinder vtk_file_4(results_dir + "results_from_time_0/outlines_10.vtu", RelativeTo::Absolute);
        TS_ASSERT(vtk_file_4.Exists());
 #endif //CHASTE_VTK
    }

    void TestPottsSpheroidWithNoBirthOrDeath()
    {
        EXIT_IF_PARALLEL;    // Potts simulations don't work in parallel because they depend on NodesOnlyMesh for writing.

        // Create a simple 3D PottsMesh
        PottsMeshGenerator<3> generator(10, 2, 2, 10, 2, 2, 10, 2, 2);
        PottsMesh<3>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<FixedG1GenerationalCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_diff_type);

        // Create cell population
        PottsBasedCellPopulation<3> cell_population(*p_mesh, cells);

        // Set up cell-based simulation
        OnLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory("TestSimplePottsSpheroid");
        simulator.SetEndTime(1.0);

        // Create update rules and pass to the simulation
        MAKE_PTR(VolumeConstraintPottsUpdateRule<3>, p_volume_constraint_update_rule);
        p_volume_constraint_update_rule->SetMatureCellTargetVolume(8);
        simulator.AddUpdateRule(p_volume_constraint_update_rule);
        MAKE_PTR(AdhesionPottsUpdateRule<3>, p_adhesion_update_rule);
        simulator.AddUpdateRule(p_adhesion_update_rule);

        // Run simulation
        simulator.Solve();

        // Check that the same number of cells
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), 8u);

        // Test no births or deaths
        TS_ASSERT_EQUALS(simulator.GetNumBirths(), 0u);
        TS_ASSERT_EQUALS(simulator.GetNumDeaths(), 0u);
    }

    void TestPottsChemotaxis()
    {
        EXIT_IF_PARALLEL;    // Potts simulations don't work in parallel because they depend on NodesOnlyMesh for writing.

        // Create a simple 3D PottsMesh
        PottsMeshGenerator<3> generator(12, 1, 2, 6, 1, 2, 6, 1, 2, false, true, true, true);
        PottsMesh<3>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<FixedG1GenerationalCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_diff_type);

        // Create cell population
        PottsBasedCellPopulation<3> cell_population(*p_mesh, cells);

        // Set up cell-based simulation
        OnLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory("TestPottsChemotaxis");
        simulator.SetSamplingTimestepMultiple(10);
        simulator.SetEndTime(50.0);

        // Create update rules and pass to the simulation
        MAKE_PTR(VolumeConstraintPottsUpdateRule<3>, p_volume_constraint_update_rule);
        p_volume_constraint_update_rule->SetMatureCellTargetVolume(8);
        simulator.AddUpdateRule(p_volume_constraint_update_rule);
        MAKE_PTR(AdhesionPottsUpdateRule<3>, p_adhesion_update_rule);
        simulator.AddUpdateRule(p_adhesion_update_rule);
        MAKE_PTR(ChemotaxisPottsUpdateRule<3>, p_chemotaxis_update_rule);
        simulator.AddUpdateRule(p_chemotaxis_update_rule);

        // Run simulation
        simulator.Solve();

        // Check that the same number of cells
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), 1u);

        // Test no births or deaths
        TS_ASSERT_EQUALS(simulator.GetNumBirths(), 0u);
        TS_ASSERT_EQUALS(simulator.GetNumDeaths(), 0u);
    }

    void TestRandomIterationOverUpdateRules()
    {
        EXIT_IF_PARALLEL;    // Potts simulations don't work in parallel because they depend on NodesOnlyMesh for writing.

        // Create a simple 2D PottsMesh
        PottsMeshGenerator<2> generator(8, 1, 4, 10, 1, 4);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        CellsGenerator<UniformG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_stem_type);

        // Create cell population
        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.SetIterateRandomlyOverUpdateRuleCollection(true);

        // Set up cell-based simulation
        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestPottsMonolayerWithBirth");
        simulator.SetDt(0.1);
        simulator.SetEndTime(20);
        simulator.SetSamplingTimestepMultiple(20);

        // Create update rules and pass to the simulation
        MAKE_PTR(VolumeConstraintPottsUpdateRule<2>, p_volume_constraint_update_rule);
        simulator.AddUpdateRule(p_volume_constraint_update_rule);
        MAKE_PTR(AdhesionPottsUpdateRule<2>, p_adhesion_update_rule);
        simulator.AddUpdateRule(p_adhesion_update_rule);

        // Run simulation
        simulator.Solve();

        // Check the number of cells
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), 2u);

        // Test no deaths and some births
        TS_ASSERT_EQUALS(simulator.GetNumBirths(), 1u);
        TS_ASSERT_EQUALS(simulator.GetNumDeaths(), 0u);
    }

    void TestPottsSpheroidCellSorting()
    {
        EXIT_IF_PARALLEL;    // Potts simulations don't work in parallel because they depend on NodesOnlyMesh for writing.

        // Create a simple 3D PottsMesh
        unsigned domain_size = 10;
        unsigned element_number = 4;
        unsigned element_size = 2;

        PottsMeshGenerator<3> generator(domain_size, element_number, element_size, domain_size, element_number, element_size, domain_size, element_number, element_size);
        PottsMesh<3>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<FixedG1GenerationalCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_diff_type);

        // Randomly label some cells
        boost::shared_ptr<AbstractCellProperty> p_label(CellPropertyRegistry::Instance()->Get<CellLabel>());
        RandomlyLabelCells(cells, p_label, 0.5);

        // Create cell population
        PottsBasedCellPopulation<3> cell_population(*p_mesh, cells);
        cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>(); // So outputs the labelled cells

        // Set up cell-based simulation
        OnLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory("TestPotts3DCellSorting");
        simulator.SetDt(0.1);
        simulator.SetEndTime(1.0);

        // Create update rules and pass to the simulation
        MAKE_PTR(VolumeConstraintPottsUpdateRule<3>, p_volume_constraint_update_rule);
        p_volume_constraint_update_rule->SetMatureCellTargetVolume(element_size*element_size*element_size);
        p_volume_constraint_update_rule->SetDeformationEnergyParameter(0.2);
        simulator.AddUpdateRule(p_volume_constraint_update_rule);

        MAKE_PTR(DifferentialAdhesionPottsUpdateRule<3>, p_differential_adhesion_update_rule);
        p_differential_adhesion_update_rule->SetLabelledCellLabelledCellAdhesionEnergyParameter(0.16);
        p_differential_adhesion_update_rule->SetLabelledCellCellAdhesionEnergyParameter(0.11);
        p_differential_adhesion_update_rule->SetCellCellAdhesionEnergyParameter(0.02);
        p_differential_adhesion_update_rule->SetLabelledCellBoundaryAdhesionEnergyParameter(0.16);
        p_differential_adhesion_update_rule->SetCellBoundaryAdhesionEnergyParameter(0.16);
        simulator.AddUpdateRule(p_differential_adhesion_update_rule);

        // Run simulation
        simulator.Solve();

        // Check that the same number of cells
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), 64u);

        // Test no births or deaths
        TS_ASSERT_EQUALS(simulator.GetNumBirths(), 0u);
        TS_ASSERT_EQUALS(simulator.GetNumDeaths(), 0u);
#ifdef CHASTE_VTK
        // Test that VTK writer has produced some files
        OutputFileHandler handler("TestPotts3DCellSorting", false);
        std::string results_dir = handler.GetOutputDirectoryFullPath();

        // Initial condition file
        FileFinder vtk_file(results_dir + "results_from_time_0/results_0.vtu", RelativeTo::Absolute);
        TS_ASSERT(vtk_file.Exists());

        // Final file
        FileFinder vtk_file2(results_dir + "results_from_time_0/results_10.vtu", RelativeTo::Absolute);
        TS_ASSERT(vtk_file2.Exists());

        // Check that the second VTK file for Potts specific data (element IDs)
        VtkMeshReader<3,3> vtk_reader(results_dir + "results_from_time_0/results_10.vtu");
        std::vector<double> cell_types;
        vtk_reader.GetPointData("Cell types", cell_types);
        TS_ASSERT_EQUALS(cell_types.size(), 1000u);

        // The cell types are between -1 and 5. Check the maximum
        TS_ASSERT_DELTA(*max_element(cell_types.begin(), cell_types.end()), 5.0, 1e-12);

        std::vector<double> cell_ids;
        vtk_reader.GetPointData("Cell IDs", cell_ids);
        // Prior to release 3.2 this was called "Element index".
        // It is changed to cell_id as this is preferable for VTK output.
        TS_ASSERT_EQUALS(cell_ids.size(), 1000u);
        TS_ASSERT_DELTA(*max_element(cell_ids.begin(), cell_ids.end()), 63.0, 1e-12);

 #endif //CHASTE_VTK
    }

    c_vector<unsigned, 6> mNodes; // TO check after save and load.

    void TestStandardResultForArchivingTestsBelow()
    {
        EXIT_IF_PARALLEL;    // Potts simulations don't work in parallel because they depend on NodesOnlyMesh for writing.

        // Create a simple 2D PottsMesh
        PottsMeshGenerator<2> generator(10, 1, 4, 10, 1, 4);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        CellsGenerator<UniformG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_stem_type);

        // Create cell population
        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up cell-based simulation
        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestOnLatticeSimulationWithPottsBasedCellPopulationStandardResult");
        simulator.SetDt(0.1);
        simulator.SetEndTime(20);
        simulator.SetSamplingTimestepMultiple(10);

        // Create update rules and pass to the simulation
        MAKE_PTR(VolumeConstraintPottsUpdateRule<2>, p_volume_constraint_update_rule);
        simulator.AddUpdateRule(p_volume_constraint_update_rule);
        MAKE_PTR(AdhesionPottsUpdateRule<2>, p_adhesion_update_rule);
        simulator.AddUpdateRule(p_adhesion_update_rule);

        // Run simulation
        simulator.Solve();

        // Check some results
        PottsElement<2>* element_0 = static_cast <PottsBasedCellPopulation<2>*>(&simulator.rGetCellPopulation())->GetElement(0u);
        TS_ASSERT_EQUALS(element_0->GetNumNodes(), 16u);
        mNodes[0] = 33;
        mNodes[1] = 36;
        mNodes[2] = 63;
        TS_ASSERT_EQUALS(element_0->GetNode(0)->GetIndex(), mNodes[0]);
        TS_ASSERT_EQUALS(element_0->GetNode(8)->GetIndex(), mNodes[1]);
        TS_ASSERT_EQUALS(element_0->GetNode(15)->GetIndex(), mNodes[2]);

        mNodes[3] = 55;
        mNodes[4] = 83;
        mNodes[5] = 94;
        PottsElement<2>* element_1 = static_cast <PottsBasedCellPopulation<2>*>(&simulator.rGetCellPopulation())->GetElement(1u);
        TS_ASSERT_EQUALS(element_1->GetNumNodes(), 16u);
        TS_ASSERT_EQUALS(element_1->GetNode(0)->GetIndex(), mNodes[3]);
        TS_ASSERT_EQUALS(element_1->GetNode(8)->GetIndex(), mNodes[4]);
        TS_ASSERT_EQUALS(element_1->GetNode(15)->GetIndex(), mNodes[5]);
    }

    void TestSave()
    {
        EXIT_IF_PARALLEL;    // Potts simulations don't work in parallel because they depend on NodesOnlyMesh for writing.

        // Create a simple 2D PottsMesh
        PottsMeshGenerator<2> generator(10, 1, 4, 10, 1, 4);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        CellsGenerator<UniformG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_stem_type);

        // Create cell population
        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up cell-based simulation
        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestOnLatticeSimulationWithPottsBasedCellPopulationSaveAndLoad");
        simulator.SetDt(0.1);
        simulator.SetEndTime(10);
        simulator.SetSamplingTimestepMultiple(10);

        // Create update rules and pass to the simulation
        MAKE_PTR(VolumeConstraintPottsUpdateRule<2>, p_volume_constraint_update_rule);
        simulator.AddUpdateRule(p_volume_constraint_update_rule);
        MAKE_PTR(AdhesionPottsUpdateRule<2>, p_adhesion_update_rule);
        simulator.AddUpdateRule(p_adhesion_update_rule);

        // Run simulation
        simulator.Solve();

        // Save the results
        CellBasedSimulationArchiver<2, OnLatticeSimulation<2> >::Save(&simulator);
    }

    void TestLoad()
    {
        EXIT_IF_PARALLEL;    // Potts simulations don't work in parallel because they depend on NodesOnlyMesh for writing.

        // Load the simulation from the TestSave() method above and run it from 10.0 to 15.0
        OnLatticeSimulation<2>* p_simulator1;
        p_simulator1 = CellBasedSimulationArchiver<2, OnLatticeSimulation<2> >::Load("TestOnLatticeSimulationWithPottsBasedCellPopulationSaveAndLoad", 10.0);

        p_simulator1->SetEndTime(15.0);
        p_simulator1->Solve();

        // Save, then reload and run from 15.0 to 20.0
        CellBasedSimulationArchiver<2, OnLatticeSimulation<2> >::Save(p_simulator1);
        OnLatticeSimulation<2>* p_simulator2
            = CellBasedSimulationArchiver<2, OnLatticeSimulation<2> >::Load("TestOnLatticeSimulationWithPottsBasedCellPopulationSaveAndLoad", 15.0);

        p_simulator2->SetEndTime(20.0);
        p_simulator2->Solve();

        // These results are from time 20.0 in TestStandardResultForArchivingTestsBelow()
        PottsElement<2>* element_0 = static_cast <PottsBasedCellPopulation<2>*>(&p_simulator2->rGetCellPopulation())->GetElement(0u);
        TS_ASSERT_EQUALS(element_0->GetNumNodes(), 16u);
        TS_ASSERT_EQUALS(element_0->GetNode(0)->GetIndex(), mNodes[0]);
        TS_ASSERT_EQUALS(element_0->GetNode(8)->GetIndex(), mNodes[1]);
        TS_ASSERT_EQUALS(element_0->GetNode(15)->GetIndex(), mNodes[2]);

        PottsElement<2>* element_1 = static_cast <PottsBasedCellPopulation<2>*>(&p_simulator2->rGetCellPopulation())->GetElement(1u);
        TS_ASSERT_EQUALS(element_1->GetNumNodes(), 16u);
        TS_ASSERT_EQUALS(element_1->GetNode(0)->GetIndex(), mNodes[3]);
        TS_ASSERT_EQUALS(element_1->GetNode(8)->GetIndex(), mNodes[4]);
        TS_ASSERT_EQUALS(element_1->GetNode(15)->GetIndex(), mNodes[5]);

        // Tidy up
        delete p_simulator1;
        delete p_simulator2;
    }
};

#endif /*TESTONLATTICESIMULATIONWITHPOTTSBASEDCELLPOPULATION_HPP_*/
