/*

Copyright (C) University of Oxford, 2005-2011

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/

#ifndef TESTONLATTICESIMULATIONWITHPOTTSBASEDCELLPOPULATION_HPP_
#define TESTONLATTICESIMULATIONWITHPOTTSBASEDCELLPOPULATION_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "HoneycombMeshGenerator.hpp"
#include "PottsMeshGenerator.hpp"
#include "CellsGenerator.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "StochasticDurationGenerationBasedCellCycleModel.hpp"
#include "WntConcentration.hpp"
#include "SimpleWntCellCycleModel.hpp"
#include "WildTypeCellMutationState.hpp"
#include "CryptCellsGenerator.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "VolumeConstraintPottsUpdateRule.hpp"
#include "AdhesionPottsUpdateRule.hpp"
#include "DifferentialAdhesionPottsUpdateRule.hpp"
#include "SloughingCellKiller.hpp"
#include "OnLatticeSimulation.hpp"
#include "OffLatticeSimulation.hpp"
#include "Warnings.hpp"
#include "LogFile.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "SmartPointers.hpp"

class TestOnLatticeSimulationWithPottsBasedCellPopulation : public AbstractCellBasedTestSuite
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

    void RandomlyLabelCells(std::vector<CellPtr>& rCells, boost::shared_ptr<AbstractCellProperty> pLabel, double labeledRatio)
    {
        for (unsigned i = 0; i<rCells.size(); i++)
        {
            if (RandomNumberGenerator::Instance()->ranf() < labeledRatio)
            {
                rCells[i]->AddCellProperty(pLabel);
            }
        }
    }

public:

    void TestOnLatticeSimulationExceptions()
    {
        // Create a simple tetrahedral mesh
        unsigned num_cells_depth = 5;
        unsigned num_cells_width = 5;
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0);
        TetrahedralMesh<2,2>* p_generating_mesh = generator.GetMesh();

        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes());

        // Create a node-based cell population
        NodeBasedCellPopulation<2> node_based_cell_population(mesh, cells);
        node_based_cell_population.SetMechanicsCutOffLength(1.5);

        TS_ASSERT_THROWS_THIS(OnLatticeSimulation<2> simulator(node_based_cell_population),"OnLatticeSimulations require a PottsBasedCellPopulation.");
    }

    void TestMoreOnLatticeSimulationExceptions()
    {
        // Create a simple 2D PottsMesh
        PottsMeshGenerator<2> generator(16, 4, 4, 18, 4, 4);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), DIFFERENTIATED);

        // Create cell population
        PottsBasedCellPopulation<2> potts_based_cell_population(*p_mesh, cells);

        // Try to set up off lattice simulation
        TS_ASSERT_THROWS_THIS(OffLatticeSimulation<2> simulator(potts_based_cell_population),"OffLatticeSimulations require a VertexBasedCellPopulation or a subclass of AbstractCentreBasedCellPopulation.");
    }

    void TestPottsMonolayerWithNoBirthOrDeath() throw (Exception)
    {
        // Create a simple 2D PottsMesh
        PottsMeshGenerator<2> generator(16, 4, 4, 18, 4, 4);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), DIFFERENTIATED);

        // Create cell population
        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up cell-based simulation
        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestSimplePottsMonolayer");
        simulator.SetEndTime(0.1);

        // Create update rules and pass to the simulation
        MAKE_PTR(VolumeConstraintPottsUpdateRule<2>, p_volume_constraint_update_rule,);
        simulator.AddUpdateRule(p_volume_constraint_update_rule);
        MAKE_PTR(AdhesionPottsUpdateRule<2>, p_adhesion_update_rule,);
        simulator.AddUpdateRule(p_adhesion_update_rule);

        // Run simulation
        simulator.Solve();

        // Check that the same number of cells
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), 16u);

        // Test no births or deaths
        TS_ASSERT_EQUALS(simulator.GetNumBirths(), 0u);
        TS_ASSERT_EQUALS(simulator.GetNumDeaths(), 0u);
    }
    
    void TestPottsMonolayerWithNonRandomSweep() throw (Exception)
    {
        // Create a simple 2D PottsMesh
        PottsMeshGenerator<2> generator(16, 4, 4, 18, 4, 4);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), DIFFERENTIATED);

        // Create cell population
        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.SetUpdateNodesInRandomOrder(false);

        // Set up cell-based simulation
        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestSimplePottsMonolayerWithRandomSweep");
        simulator.SetEndTime(0.1);

        // Create update rules and pass to the simulation
        MAKE_PTR(VolumeConstraintPottsUpdateRule<2>, p_volume_constraint_update_rule,);
        simulator.AddUpdateRule(p_volume_constraint_update_rule);
        MAKE_PTR(AdhesionPottsUpdateRule<2>, p_adhesion_update_rule,);
        simulator.AddUpdateRule(p_adhesion_update_rule);

        // Run simulation
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());
    }

    void TestPottsMonolayerWithDeath() throw (Exception)
    {
        // Create a simple 2D PottsMesh
        PottsMeshGenerator<2> generator(16, 4, 4, 24, 8, 2);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), DIFFERENTIATED);

        // Create cell population
        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up cell-based simulation
        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestPottsMonolayerWithDeath");
        simulator.SetDt(0.1);
        simulator.SetEndTime(1.0);

        // Create update rules and pass to the simulation
        MAKE_PTR(VolumeConstraintPottsUpdateRule<2>, p_volume_constraint_update_rule,);
        simulator.AddUpdateRule(p_volume_constraint_update_rule);
        MAKE_PTR(AdhesionPottsUpdateRule<2>, p_adhesion_update_rule,);
        simulator.AddUpdateRule(p_adhesion_update_rule);

        // Create cell killer and pass in to simulation
        SloughingCellKiller<2> sloughing_cell_killer(&cell_population,16u);
        simulator.AddCellKiller(&sloughing_cell_killer);

        // Run simulation
        simulator.Solve();

        // Check the number of cells
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), 20u);

        // Test no births or deaths
        TS_ASSERT_EQUALS(simulator.GetNumBirths(), 0u);
        TS_ASSERT_EQUALS(simulator.GetNumDeaths(), 12u);
    }

    void TestPottsMonolayerWithBirth() throw (Exception)
    {
        // Create a simple 2D PottsMesh
        PottsMeshGenerator<2> generator(8, 1, 4, 10, 1, 4);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<StochasticDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), STEM);

        // Create cell population
        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up cell-based simulation
        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestPottsMonolayerWithBirth");
        simulator.SetDt(0.1);
        simulator.SetEndTime(20);
        simulator.SetSamplingTimestepMultiple(20);

        // Create update rules and pass to the simulation
        MAKE_PTR(VolumeConstraintPottsUpdateRule<2>, p_volume_constraint_update_rule,);
        simulator.AddUpdateRule(p_volume_constraint_update_rule);
        MAKE_PTR(AdhesionPottsUpdateRule<2>, p_adhesion_update_rule,);
        simulator.AddUpdateRule(p_adhesion_update_rule);

        // Run simulation
        simulator.Solve();

        // Check the number of cells
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), 4u);

        // Test no deaths and some births
        TS_ASSERT_EQUALS(simulator.GetNumBirths(), 3u);
        TS_ASSERT_EQUALS(simulator.GetNumDeaths(), 0u);
    }

    void TestPottsMonolayerCellSorting() throw (Exception)
    {
        // Create a simple 2D PottsMesh
        PottsMeshGenerator<2> generator(30, 4, 4, 30, 4, 4);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), DIFFERENTIATED);

        // Make this pointer first as if we move it after creating the cell population the label numbers aren't tracked
        boost::shared_ptr<AbstractCellProperty> p_label(new CellLabel);

        RandomlyLabelCells(cells,p_label,0.5);

        // Create cell population
        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.SetOutputCellMutationStates(true); // So outputs the labelled cells

        // Set up cell-based simulation
        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestPottsCellSorting");
        simulator.SetDt(0.1);
        simulator.SetEndTime(1.0);

        // Create update rules and pass to the simulation
        MAKE_PTR(VolumeConstraintPottsUpdateRule<2>, p_volume_constraint_update_rule,);
        p_volume_constraint_update_rule->SetMatureCellTargetVolume(16);
        p_volume_constraint_update_rule->SetDeformationEnergyParameter(0.2);
        simulator.AddUpdateRule(p_volume_constraint_update_rule);

        MAKE_PTR(DifferentialAdhesionPottsUpdateRule<2>, p_differential_adhesion_update_rule,);
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
        FileFinder vtk_file(results_dir + "results_from_time_0/results_0.vtu", RelativeTo::Absolute);
        TS_ASSERT(vtk_file.Exists());

        // Final file
        FileFinder vtk_file2(results_dir + "results_from_time_0/results_10.vtu", RelativeTo::Absolute);
        TS_ASSERT(vtk_file2.Exists());
 #endif //CHASTE_VTK
    }

    void TestPottsSpheroidWithNoBirthOrDeath() throw (Exception)
    {
        // Create a simple 3D PottsMesh
        PottsMeshGenerator<3> generator(10, 2, 2, 10, 2, 2, 10, 2, 2);
        PottsMesh<3>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), DIFFERENTIATED);

        // Create cell population
        PottsBasedCellPopulation<3> cell_population(*p_mesh, cells);

        // Set up cell-based simulation
        OnLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory("TestSimplePottsSpheroid");
        simulator.SetDt(0.1);
        simulator.SetEndTime(1.0);

        // Create update rules and pass to the simulation
        MAKE_PTR(VolumeConstraintPottsUpdateRule<3>, p_volume_constraint_update_rule,);
        p_volume_constraint_update_rule->SetMatureCellTargetVolume(16);
        simulator.AddUpdateRule(p_volume_constraint_update_rule);
        MAKE_PTR(AdhesionPottsUpdateRule<3>, p_adhesion_update_rule,);
        simulator.AddUpdateRule(p_adhesion_update_rule);


        // Run simulation
        simulator.Solve();

        // Check that the same number of cells
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), 8u);

        // Test no births or deaths
        TS_ASSERT_EQUALS(simulator.GetNumBirths(), 0u);
        TS_ASSERT_EQUALS(simulator.GetNumDeaths(), 0u);
    }

    void TestPottsSpheroidCellSorting() throw (Exception)
    {
        // Create a simple 3D PottsMesh
        unsigned domain_size = 10;
        unsigned element_number = 4;
        unsigned element_size = 2;

        PottsMeshGenerator<3> generator(domain_size, element_number, element_size, domain_size, element_number, element_size, domain_size, element_number, element_size);
        PottsMesh<3>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), DIFFERENTIATED);

        // Make this pointer first as if we move it after creating the cell population the label numbers aren't tracked
        boost::shared_ptr<AbstractCellProperty> p_label(new CellLabel);

        RandomlyLabelCells(cells,p_label,0.5);

        // Create cell population
        PottsBasedCellPopulation<3> cell_population(*p_mesh, cells);
        cell_population.SetOutputCellMutationStates(true); // So outputs the labelled cells

        // Set up cell-based simulation
        OnLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory("TestPotts3DCellSorting");
        simulator.SetDt(0.1);
        simulator.SetEndTime(1.0);

        // Create update rules and pass to the simulation
        MAKE_PTR(VolumeConstraintPottsUpdateRule<3>, p_volume_constraint_update_rule,);
        p_volume_constraint_update_rule->SetMatureCellTargetVolume(element_size*element_size*element_size);
        p_volume_constraint_update_rule->SetDeformationEnergyParameter(0.2);
        simulator.AddUpdateRule(p_volume_constraint_update_rule);

        MAKE_PTR(DifferentialAdhesionPottsUpdateRule<3>, p_differential_adhesion_update_rule,);
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
        //Test that VTK writer has produced some files
        OutputFileHandler handler("TestPotts3DCellSorting", false);
        std::string results_dir = handler.GetOutputDirectoryFullPath();
        // Initial condition file
        FileFinder vtk_file(results_dir + "results_from_time_0/results_0.vtu", RelativeTo::Absolute);
        TS_ASSERT(vtk_file.Exists());

        // Final file
        FileFinder vtk_file2(results_dir + "results_from_time_0/results_10.vtu", RelativeTo::Absolute);
        TS_ASSERT(vtk_file2.Exists());
 #endif //CHASTE_VTK
    }

    void TestStandardResultForArchivingTestsBelow() throw (Exception)
    {
        // Create a simple 2D PottsMesh
        PottsMeshGenerator<2> generator(10, 1, 4, 10, 1, 4);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<StochasticDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), STEM);

        // Create cell population
        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up cell-based simulation
        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestOnLatticeSimulationWithPottsBasedCellPopulationStandardResult");
        simulator.SetDt(0.1);
        simulator.SetEndTime(20);
        simulator.SetSamplingTimestepMultiple(10);

        // Create update rules and pass to the simulation
        MAKE_PTR(VolumeConstraintPottsUpdateRule<2>, p_volume_constraint_update_rule,);
        simulator.AddUpdateRule(p_volume_constraint_update_rule);
        MAKE_PTR(AdhesionPottsUpdateRule<2>, p_adhesion_update_rule,);
        simulator.AddUpdateRule(p_adhesion_update_rule);

        // Run simulation
        simulator.Solve();

        // Check some results
        PottsElement<2>* element_0 = static_cast <PottsBasedCellPopulation<2>*>(&simulator.rGetCellPopulation())->GetElement(0u);
        TS_ASSERT_EQUALS(element_0->GetNumNodes(), 16u);
        TS_ASSERT_EQUALS(element_0->GetNode(0)->GetIndex(), 33u);
        TS_ASSERT_EQUALS(element_0->GetNode(8)->GetIndex(), 57u);
        TS_ASSERT_EQUALS(element_0->GetNode(15)->GetIndex(), 78u);

        PottsElement<2>* element_1 = static_cast <PottsBasedCellPopulation<2>*>(&simulator.rGetCellPopulation())->GetElement(1u);
        TS_ASSERT_EQUALS(element_1->GetNumNodes(), 16u);
        TS_ASSERT_EQUALS(element_1->GetNode(0)->GetIndex(), 80u);
        TS_ASSERT_EQUALS(element_1->GetNode(8)->GetIndex(), 72u);
        TS_ASSERT_EQUALS(element_1->GetNode(15)->GetIndex(), 74u);
    }

    // Testing Save
    void TestSave() throw (Exception)
    {
        // Create a simple 2D PottsMesh
        PottsMeshGenerator<2> generator(10, 1, 4, 10, 1, 4);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<StochasticDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), STEM);

        // Create cell population
        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up cell-based simulation
        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestOnLatticeSimulationWithPottsBasedCellPopulationSaveAndLoad");
        simulator.SetDt(0.1);
        simulator.SetEndTime(10);
        simulator.SetSamplingTimestepMultiple(10);

        // Create update rules and pass to the simulation
        MAKE_PTR(VolumeConstraintPottsUpdateRule<2>, p_volume_constraint_update_rule,);
        simulator.AddUpdateRule(p_volume_constraint_update_rule);
        MAKE_PTR(AdhesionPottsUpdateRule<2>, p_adhesion_update_rule,);
        simulator.AddUpdateRule(p_adhesion_update_rule);

        // Run simulation
        simulator.Solve();

        // Save the results
        CellBasedSimulationArchiver<2, OnLatticeSimulation<2> >::Save(&simulator);
    }

    // Testing Load (based on previous two tests)
    void TestLoad() throw (Exception)
    {
        // Load the simulation from the TestSave method above and
        // run it from 10.0 to 15.0
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

        // These results are from time 20.0 in TestStandardResultForArchivingTestBelow()
        PottsElement<2>* element_0 = static_cast <PottsBasedCellPopulation<2>*>(&p_simulator2->rGetCellPopulation())->GetElement(0u);
        TS_ASSERT_EQUALS(element_0->GetNumNodes(), 16u);
        TS_ASSERT_EQUALS(element_0->GetNode(0)->GetIndex(), 33u);
        TS_ASSERT_EQUALS(element_0->GetNode(8)->GetIndex(), 57u);
        TS_ASSERT_EQUALS(element_0->GetNode(15)->GetIndex(), 78u);

        PottsElement<2>* element_1 = static_cast <PottsBasedCellPopulation<2>*>(&p_simulator2->rGetCellPopulation())->GetElement(1u);
        TS_ASSERT_EQUALS(element_1->GetNumNodes(), 16u);
        TS_ASSERT_EQUALS(element_1->GetNode(0)->GetIndex(), 80u);
        TS_ASSERT_EQUALS(element_1->GetNode(8)->GetIndex(), 72u);
        TS_ASSERT_EQUALS(element_1->GetNode(15)->GetIndex(), 74u);

        // Tidy up
        delete p_simulator1;
        delete p_simulator2;
    }

    void TestPottsCrypt() throw (Exception)
    {

        double crypt_length = 80;


        // Create a simple 2D PottsMesh
        //PottsMeshGenerator<2> generator(40, 10, 4, 90, 20, 4, 1, 1, 1, true);
        PottsMeshGenerator<2> generator(20, 5, 4, 45, 10, 4, 1, 1, 1, true);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<SimpleWntCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), TRANSIT);

        // Create cell population
        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.SetOutputCellVolumes(true);

        // Create an instance of a Wnt concentration
        WntConcentration<2>::Instance()->SetType(LINEAR);
        WntConcentration<2>::Instance()->SetCellPopulation(cell_population);
        WntConcentration<2>::Instance()->SetCryptLength(crypt_length);

        // Set up cell-based simulation
        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestPottsCrypt");
        simulator.SetDt(0.1);
        simulator.SetSamplingTimestepMultiple(1);
        simulator.SetEndTime(1.0);
        simulator.SetOutputCellVelocities(true);

        // Create cell killer and pass in to crypt simulation
        SloughingCellKiller<2> sloughing_cell_killer(&cell_population,crypt_length);
        simulator.AddCellKiller(&sloughing_cell_killer);


        // Create update rules and pass to the simulation
        MAKE_PTR(VolumeConstraintPottsUpdateRule<2>, p_volume_constraint_update_rule,);
        simulator.AddUpdateRule(p_volume_constraint_update_rule);
        MAKE_PTR(AdhesionPottsUpdateRule<2>, p_adhesion_update_rule,);
        simulator.AddUpdateRule(p_adhesion_update_rule);

        // Run simulation
        simulator.Solve();

        // Check the number of cells
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), 56u);

        // Test number of births or deaths
        TS_ASSERT_EQUALS(simulator.GetNumBirths(), 6u);
        TS_ASSERT_EQUALS(simulator.GetNumDeaths(), 0u);
    }
};

#endif /*TESTONLATTICESIMULATIONWITHPOTTSBASEDCELLPOPULATION_HPP_*/
