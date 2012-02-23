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

#ifndef TESTOFFLATTICESIMULATIONWITHPDES_HPP_
#define TESTOFFLATTICESIMULATIONWITHPDES_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "OffLatticeSimulation.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "OxygenBasedCellKiller.hpp"
#include "SimpleOxygenBasedCellCycleModel.hpp"
#include "CellsGenerator.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "SimpleUniformSourcePde.hpp"
#include "CellwiseSourcePde.hpp"
#include "ConstBoundaryCondition.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "ReplicatableVector.hpp"
#include "NumericFileComparison.hpp"
#include "WildTypeCellMutationState.hpp"
#include "FunctionalBoundaryCondition.hpp"
#include "AveragedSourcePde.hpp"
#include "VolumeDependentAveragedSourcePde.hpp"
#include "SmartPointers.hpp"

class SimplePdeForTesting : public AbstractLinearEllipticPde<2,2>
{
public:
    double ComputeConstantInUSourceTerm(const ChastePoint<2>&, Element<2,2>* pElement)
    {
        return -1.0;
    }

    double ComputeLinearInUCoeffInSourceTerm(const ChastePoint<2>&, Element<2,2>*)
    {
        return 0.0;
    }

    c_matrix<double,2,2> ComputeDiffusionTerm(const ChastePoint<2>& )
    {
        return identity_matrix<double>(2);
    }
};

/**
 * For use in TestOffLatticeSimulationWithPdes::TestWithBoundaryConditionVaryingInTime.
 */
double bc_func(const ChastePoint<2>& p)
{
    double value = SimulationTime::Instance()->GetTime();
    return value;
}

class TestOffLatticeSimulationWithPdes : public AbstractCellBasedTestSuite
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

    /*
     * A two-part test for the UpdateAtEndOfTimeStep() method.
     *
     * Firstly, test the PDE solver using the problem del squared C = 1
     * on the unit disc, with boundary condition C=1 on r=1, which has
     * analytic solution C = 1-0.25*(1-r^2).
     *
     * Secondly, test that cells' hypoxic durations are correctly updated when a
     * nutrient distribution is prescribed.
     */
    void TestUpdateAtEndOfTimeStep() throw(Exception)
    {
        EXIT_IF_PARALLEL;

        // Set up mesh
        MutableMesh<2,2> mesh;
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_522_elements");
        mesh.ConstructFromMeshReader(mesh_reader);

        // Set up cells
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            SimpleOxygenBasedCellCycleModel* p_model = new SimpleOxygenBasedCellCycleModel();
            p_model->SetDimension(2);
            p_model->SetCellProliferativeType(STEM);
            p_model->SetHypoxicConcentration(0.9);
            p_model->SetQuiescentConcentration(0.9);

            // Use non-default G1 durations
            p_model->SetStemCellG1Duration(8.0);
            p_model->SetTransitCellG1Duration(8.0);

            CellPtr p_cell(new Cell(p_state, p_model));

            double birth_time = -RandomNumberGenerator::Instance()->ranf()*18.0;
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        // Set up cell population
        MeshBasedCellPopulation<2> cell_population(mesh, cells);

        // Set up CellwiseData and associate it with the cell population
        CellwiseData<2>* p_data = CellwiseData<2>::Instance();
        p_data->SetNumCellsAndVars(cell_population.GetNumRealCells(), 1);
        p_data->SetCellPopulation(&cell_population);

        // Since values are first passed in to CellwiseData before it is updated in UpdateAtEndOfTimeStep(),
        // we must initialise it here to avoid memory errors
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            p_data->SetValue(1.0, mesh.GetNode(i)->GetIndex());
        }

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestUpdateAtEndOfTimeStepMethod");
        simulator.SetEndTime(2.0/120.0);

        // Set up PDE and pass to simulation via handler
        CellBasedPdeHandler<2> pde_handler(&cell_population);
        pde_handler.SetImposeBcsOnCoarseBoundary(false);

        SimplePdeForTesting pde;
        ConstBoundaryCondition<2> bc(1.0);
        PdeAndBoundaryConditions<2> pde_and_bc(&pde, &bc, false);

        pde_handler.AddPdeAndBc(&pde_and_bc);
        simulator.SetCellBasedPdeHandler(&pde_handler);

        /*
         * Create a force law and pass it to the simulation. Use an extremely small
         * cutoff so that no cells interact - this is to ensure that in the Solve()
         * method, the cells don't move (we need to call Solve() to set up the
         * .vizpdesolution file).
         */
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetCutOffLength(0.0001);
        simulator.AddForce(p_linear_force);

        // Set up cell killer and pass into simulation
        MAKE_PTR_ARGS(OxygenBasedCellKiller<2>, p_killer, (&cell_population));
        simulator.AddCellKiller(p_killer);

        // Run cell-based simulation
        simulator.Solve();

        // Check the correct solution was obtained
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            double radius = norm_2(cell_population.GetLocationOfCellCentre(*cell_iter));
            double analytic_solution = 1 - 0.25*(1 - pow(radius,2.0));

            // Get cell model
            AbstractCellCycleModel* p_abstract_model = cell_iter->GetCellCycleModel();
            SimpleOxygenBasedCellCycleModel* p_oxygen_model = static_cast<SimpleOxygenBasedCellCycleModel*> (p_abstract_model);

            // First part of test - check that PDE solver is working correctly
            TS_ASSERT_DELTA(p_data->GetValue(*cell_iter), analytic_solution, 1e-2);

            // Second part of test - check that each cell's hypoxic duration is correctly updated
            if (p_data->GetValue(*cell_iter) >= p_oxygen_model->GetHypoxicConcentration())
            {
                TS_ASSERT_DELTA(p_oxygen_model->GetCurrentHypoxicDuration(), 0.0, 1e-5);
            }
            else
            {
                TS_ASSERT_DELTA(p_oxygen_model->GetCurrentHypoxicDuration(), 2/120.0, 1e-5);
            }
        }

        // Tidy up
        CellwiseData<2>::Destroy();
    }

    void TestWithOxygen() throw(Exception)
    {
        EXIT_IF_PARALLEL;

        // Set up mesh
        HoneycombMeshGenerator generator(5, 5, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Set up cells
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            SimpleOxygenBasedCellCycleModel* p_model = new SimpleOxygenBasedCellCycleModel();
            p_model->SetDimension(2);
            p_model->SetCellProliferativeType(STEM);

            // Use non-default G1 durations
            p_model->SetStemCellG1Duration(8.0);
            p_model->SetTransitCellG1Duration(8.0);

            CellPtr p_cell(new Cell(p_state, p_model));

            double birth_time = -1.0 - ((double) i/p_mesh->GetNumNodes())*18.0;
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        // Set up cell population
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.SetOutputCellPopulationVolumes(true); // record the spheroid radius and apoptotic radius

        // Set up CellwiseData and associate it with the cell population
        CellwiseData<2>* p_data = CellwiseData<2>::Instance();
        p_data->SetNumCellsAndVars(cell_population.GetNumRealCells(), 1);
        p_data->SetCellPopulation(&cell_population);
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            p_data->SetValue(1.0, p_mesh->GetNode(i)->GetIndex());
        }

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(cell_population);
        TS_ASSERT_EQUALS(simulator.GetIdentifier(), "OffLatticeSimulation-2");

        // Set up PDE and pass to simulation via handler
        SimpleUniformSourcePde<2> pde(-0.1);
        ConstBoundaryCondition<2> bc(1.0);
        PdeAndBoundaryConditions<2> pde_and_bc(&pde, &bc, false);

        CellBasedPdeHandler<2> pde_handler(&cell_population);
        pde_handler.AddPdeAndBc(&pde_and_bc);
        simulator.SetCellBasedPdeHandler(&pde_handler);

        simulator.SetOutputDirectory("OffLatticeSimulationWithOxygen");
        simulator.SetEndTime(0.5);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetCutOffLength(1.5);
        simulator.AddForce(p_linear_force);

        // Set up cell killer and pass into simulation
        MAKE_PTR_ARGS(OxygenBasedCellKiller<2>, p_killer, (&cell_population));
        simulator.AddCellKiller(p_killer);

        // Run cell-based simulation
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());

        // Tidy up
        CellwiseData<2>::Destroy();
    }

    /*
     * This test compares the visualizer output from the previous test
     * with a known file.
     *
     * Note: if the previous test is changed we need to update the file
     * this test refers to.
     */
    void TestVisualizerOutput() throw (Exception)
    {
        EXIT_IF_PARALLEL;

        // Work out where one of the previous tests wrote its files
        OutputFileHandler handler("OffLatticeSimulationWithOxygen", false);
        std::string results_dir = handler.GetOutputDirectoryFullPath() + "results_from_time_0";

        NumericFileComparison comp_nut(results_dir + "/results.vizpdesolution", "cell_based/test/data/OffLatticeSimulationWithOxygen/results.vizpdesolution");
        TS_ASSERT(comp_nut.CompareFiles());
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "/results.vizpdesolution cell_based/test/data/OffLatticeSimulationWithOxygen/results.vizpdesolution").c_str()), 0);

        NumericFileComparison comp_ele(results_dir + "/results.vizelements", "cell_based/test/data/OffLatticeSimulationWithOxygen/results.vizelements");
        TS_ASSERT(comp_ele.CompareFiles());
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "/results.vizelements cell_based/test/data/OffLatticeSimulationWithOxygen/results.vizelements").c_str()), 0);

        NumericFileComparison comp_nodes(results_dir + "/results.viznodes", "cell_based/test/data/OffLatticeSimulationWithOxygen/results.viznodes");
        TS_ASSERT(comp_nodes.CompareFiles(1e-15));

        NumericFileComparison comp_celltypes(results_dir + "/results.vizcelltypes", "cell_based/test/data/OffLatticeSimulationWithOxygen/results.vizcelltypes");
        TS_ASSERT(comp_celltypes.CompareFiles(1e-15));

        TS_ASSERT_EQUALS(system(("diff " + results_dir + "/results.vizsetup cell_based/test/data/OffLatticeSimulationWithOxygen/results.vizsetup").c_str()), 0);
    }

    void TestWithPointwiseSource() throw(Exception)
    {
        EXIT_IF_PARALLEL;

        // Set up mesh
        HoneycombMeshGenerator generator(5, 5, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Set up cells
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(ApoptoticCellProperty, p_apoptotic_state);

        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            SimpleOxygenBasedCellCycleModel* p_model = new SimpleOxygenBasedCellCycleModel();
            p_model->SetDimension(2);
            p_model->SetCellProliferativeType(STEM);
            CellPtr p_cell(new Cell(p_state, p_model));

            // Use non-default G1 durations
            p_model->SetStemCellG1Duration(8.0);
            p_model->SetTransitCellG1Duration(8.0);

            double birth_time = -1.0 - ((double) i/p_mesh->GetNumNodes())*18.0;
            p_cell->SetBirthTime(birth_time);

            // Make the cell apoptotic if near the centre
            double x = p_mesh->GetNode(i)->rGetLocation()[0];
            double y = p_mesh->GetNode(i)->rGetLocation()[1];
            double dist_from_centre = sqrt( (x-2.5)*(x-2.5) + (y-2.5)*(y-2.5) );

            if (dist_from_centre < 1.5)
            {
                p_cell->AddCellProperty(p_apoptotic_state);
            }

            cells.push_back(p_cell);
        }

        // Set up cell population
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.SetOutputCellPopulationVolumes(true); // record the spheroid radius and apoptotic radius

        // Set up CellwiseData and associate it with the cell population
        CellwiseData<2>* p_data = CellwiseData<2>::Instance();
        p_data->SetNumCellsAndVars(cell_population.GetNumRealCells(), 1);
        p_data->SetCellPopulation(&cell_population);
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            p_data->SetValue(1.0, p_mesh->GetNode(i)->GetIndex());
        }

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("OffLatticeSimulationWithPdes");
        simulator.SetEndTime(0.5);

        // Set up PDE and pass to simulation via handler
        CellwiseSourcePde<2> pde(cell_population, -0.1);
        ConstBoundaryCondition<2> bc(1.0);
        PdeAndBoundaryConditions<2> pde_and_bc(&pde, &bc, false);

        CellBasedPdeHandler<2> pde_handler(&cell_population);
        pde_handler.AddPdeAndBc(&pde_and_bc);
        pde_handler.SetImposeBcsOnCoarseBoundary(false);

        simulator.SetCellBasedPdeHandler(&pde_handler);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetCutOffLength(1.5);
        simulator.AddForce(p_linear_force);

        // Set up cell killer and pass into simulation
        MAKE_PTR_ARGS(OxygenBasedCellKiller<2>, p_killer, (&cell_population));
        simulator.AddCellKiller(p_killer);

        // Run cell-based simulation
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());

        // A few hardcoded tests to check nothing has changed
        std::vector<double> node_5_location = simulator.GetNodeLocation(5);
        TS_ASSERT_DELTA(node_5_location[0], 0.6576, 1e-4);
        TS_ASSERT_DELTA(node_5_location[1], 1.1358, 1e-4);
        TS_ASSERT_DELTA(p_data->GetValue(simulator.rGetCellPopulation().GetCellUsingLocationIndex(5)), 0.9702, 1e-4);

        // Tidy up
        CellwiseData<2>::Destroy();
    }

    void TestWithPointwiseTwoSource() throw(Exception)
    {
        EXIT_IF_PARALLEL;

        // Set up mesh
        HoneycombMeshGenerator generator(5, 5, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Set up cells
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(ApoptoticCellProperty, p_apoptotic_state);

        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            SimpleOxygenBasedCellCycleModel* p_model = new SimpleOxygenBasedCellCycleModel();
            p_model->SetDimension(2);
            p_model->SetCellProliferativeType(STEM);

            // Use non-default G1 durations
            p_model->SetStemCellG1Duration(8.0);
            p_model->SetTransitCellG1Duration(8.0);

            CellPtr p_cell(new Cell(p_state, p_model));
            double birth_time = -1.0 - ((double) i/p_mesh->GetNumNodes())*18.0;
            p_cell->SetBirthTime(birth_time);

            // Make the cell apoptotic if near the centre
            double x = p_mesh->GetNode(i)->rGetLocation()[0];
            double y = p_mesh->GetNode(i)->rGetLocation()[1];
            double dist_from_centre = sqrt( (x-2.5)*(x-2.5) + (y-2.5)*(y-2.5) );

            if (dist_from_centre < 1.5)
            {
                p_cell->AddCellProperty(p_apoptotic_state);
            }

            cells.push_back(p_cell);
        }

        // Set up cell population
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.SetOutputCellPopulationVolumes(true); // record the spheroid radius and apoptotic radius

        // Set up CellwiseData and associate it with the cell population
        CellwiseData<2>* p_data = CellwiseData<2>::Instance();
        p_data->SetNumCellsAndVars(cell_population.GetNumRealCells(), 2);
        p_data->SetCellPopulation(&cell_population);
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            p_data->SetValue(1.0, p_mesh->GetNode(i)->GetIndex(), 0);
            p_data->SetValue(1.0, p_mesh->GetNode(i)->GetIndex(), 1);
        }

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("OffLatticeSimulationWithPointwiseSource");
        simulator.SetEndTime(0.5);

        CellBasedPdeHandler<2> pde_handler(&cell_population);

        // Set up first PDE and pass to handler
        CellwiseSourcePde<2> pde(cell_population, -0.1);
        ConstBoundaryCondition<2> bc(1.0);
        PdeAndBoundaryConditions<2> pde_and_bc(&pde, &bc, false);
        pde_handler.AddPdeAndBc(&pde_and_bc);

        // Set up second PDE and pass to handler
        CellwiseSourcePde<2> pde2(cell_population, -0.8);
        ConstBoundaryCondition<2> bc2(0.0);
        PdeAndBoundaryConditions<2> pde_and_bc2(&pde2, &bc2, true);
        pde_handler.AddPdeAndBc(&pde_and_bc2);

        // Pass PDE handler to simulation
        pde_handler.SetImposeBcsOnCoarseBoundary(false);
        simulator.SetCellBasedPdeHandler(&pde_handler);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetCutOffLength(1.5);
        simulator.AddForce(p_linear_force);

        // Set up cell killer and pass into simulation
        MAKE_PTR_ARGS(OxygenBasedCellKiller<2>, p_killer, (&cell_population));
        simulator.AddCellKiller(p_killer);

        // Run cell-based simulation
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());

        // A few hardcoded tests to check nothing has changed
        std::vector<double> node_5_location = simulator.GetNodeLocation(5);
        TS_ASSERT_DELTA(node_5_location[0], 0.6576, 1e-4);
        TS_ASSERT_DELTA(node_5_location[1], 1.1358, 1e-4);
        TS_ASSERT_DELTA(p_data->GetValue(simulator.rGetCellPopulation().GetCellUsingLocationIndex(5),0), 0.9702, 1e-4);
        TS_ASSERT_DELTA(p_data->GetValue(simulator.rGetCellPopulation().GetCellUsingLocationIndex(5),1), 0.0000, 1e-4);
        TS_ASSERT_LESS_THAN(p_data->GetValue(simulator.rGetCellPopulation().GetCellUsingLocationIndex(5),1),
                            p_data->GetValue(simulator.rGetCellPopulation().GetCellUsingLocationIndex(5),0));

        // Tidy up
        CellwiseData<2>::Destroy();
    }

    void TestSpheroidStatistics() throw (Exception)
    {
        EXIT_IF_PARALLEL;

        // Set up mesh
        HoneycombMeshGenerator generator(5, 5, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Set up cells
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(ApoptoticCellProperty, p_apoptotic_state);

        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            SimpleOxygenBasedCellCycleModel* p_model = new SimpleOxygenBasedCellCycleModel();
            p_model->SetDimension(2);
            p_model->SetCellProliferativeType(STEM);

            // Use non-default G1 durations
            p_model->SetStemCellG1Duration(8.0);
            p_model->SetTransitCellG1Duration(8.0);

            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetBirthTime(-0.1);

            // Label three neighbouring cells as apoptotic
            if (i==12 || i==13 || i==17)
            {
                p_cell->AddCellProperty(p_apoptotic_state);
            }
            cells.push_back(p_cell);
        }

        // Set up cell population
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.SetOutputCellPopulationVolumes(true); // record the spheroid radius and apoptotic radius

        // Set up CellwiseData and associate it with the cell population
        CellwiseData<2>* p_data = CellwiseData<2>::Instance();
        p_data->SetNumCellsAndVars(cell_population.GetNumRealCells(), 1);
        p_data->SetCellPopulation(&cell_population);

        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            p_data->SetValue(1.0, p_mesh->GetNode(i)->GetIndex());
        }

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestSpheroidStatistics");
        simulator.SetEndTime(1.0/120.0);

        // Set up PDE and pass to simulation via handler
        SimpleUniformSourcePde<2> pde(-0.1);
        ConstBoundaryCondition<2> bc(1.0);
        PdeAndBoundaryConditions<2> pde_and_bc(&pde, &bc, false);

        CellBasedPdeHandler<2> pde_handler(&cell_population);
        pde_handler.AddPdeAndBc(&pde_and_bc);
        pde_handler.SetWriteAverageRadialPdeSolution(5);
        pde_handler.SetImposeBcsOnCoarseBoundary(false);

        simulator.SetCellBasedPdeHandler(&pde_handler);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetCutOffLength(1.5);
        simulator.AddForce(p_linear_force);

        // Add an oxygen-dependent cell killer to the cell-based simulation
        MAKE_PTR_ARGS(OxygenBasedCellKiller<2>, p_killer, (&cell_population));
        simulator.AddCellKiller(p_killer);

        // Run the cell-based simulation for one timestep
        simulator.Solve();

        // Just check that we do indeed have three apoptotic cells
        unsigned num_apoptotic_cells = 0;
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            if (cell_iter->HasCellProperty<ApoptoticCellProperty>())
            {
                num_apoptotic_cells++;
            }
        }
        TS_ASSERT_EQUALS(num_apoptotic_cells, 3u);

        /**
         * We have 25 cells. Adding up the boundary cell areas, we
         * should have the equivalent area of 16 full regular hexagonal
         * cells.
         *
         * The area of a single hexagonal cell is sqrt(3)/2, so the
         * correct spheroid radius is given by sqrt((16*sqrt(3)/2)/pi).
         *
         * Since there are 3 apoptotic cells, the correct apoptotic radius
         * is given by sqrt((3*sqrt(3)/2)/pi).
         */

        // Work out where the previous test wrote its files
        OutputFileHandler handler("TestSpheroidStatistics", false);
        std::string areas_results_file = handler.GetOutputDirectoryFullPath() + "results_from_time_0/cellpopulationareas.dat";
        TS_ASSERT_EQUALS(system(("diff " + areas_results_file + " cell_based/test/data/TestSpheroidStatistics/cellpopulationareas.dat").c_str()), 0);

        std::string dist_results_file = handler.GetOutputDirectoryFullPath() + "results_from_time_0/radial_dist.dat";
        TS_ASSERT_EQUALS(system(("diff " + dist_results_file + " cell_based/test/data/TestSpheroidStatistics/radial_dist.dat").c_str()), 0);

        // Coverage
        TS_ASSERT_THROWS_NOTHING(pde_handler.WriteAverageRadialPdeSolution(SimulationTime::Instance()->GetTime()));

        // Tidy up
        CellwiseData<2>::Destroy();
    }

    void TestCoarseSourceMesh() throw(Exception)
    {
        EXIT_IF_PARALLEL;

        // Set up mesh
        HoneycombMeshGenerator generator(5, 5, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Set up cells
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            SimpleOxygenBasedCellCycleModel* p_model = new SimpleOxygenBasedCellCycleModel();
            p_model->SetDimension(2);
            p_model->SetCellProliferativeType(STEM);

            // Use non-default G1 durations
            p_model->SetStemCellG1Duration(8.0);
            p_model->SetTransitCellG1Duration(8.0);

            CellPtr p_cell(new Cell(p_state, p_model));
            double birth_time = -RandomNumberGenerator::Instance()->ranf()*18.0;
            p_cell->SetBirthTime(birth_time);

            cells.push_back(p_cell);
        }

        // Set up cell population
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up CellwiseData and associate it with the cell population
        CellwiseData<2>* p_data = CellwiseData<2>::Instance();
        p_data->SetNumCellsAndVars(cell_population.GetNumRealCells(), 2);
        p_data->SetCellPopulation(&cell_population);

        // Since values are first passed in to CellwiseData before it is updated in UpdateAtEndOfTimeStep(),
        // we need to pass it some initial conditions to avoid memory errors
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            p_data->SetValue(1.0, p_mesh->GetNode(i)->GetIndex(),0);
            p_data->SetValue(1.0, p_mesh->GetNode(i)->GetIndex(),1);
        }

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestCoarseSourceMesh");
        simulator.SetEndTime(0.05);

        CellBasedPdeHandler<2> pde_handler(&cell_population);

        // Set up PDE and pass to handler
        AveragedSourcePde<2> pde(cell_population, -0.1);
        ConstBoundaryCondition<2> bc(1.0);
        PdeAndBoundaryConditions<2> pde_and_bc(&pde, &bc, false);
        pde_handler.AddPdeAndBc(&pde_and_bc);

        // Set up second PDE and pass to handler
        AveragedSourcePde<2> pde2(cell_population, -0.5);
        PdeAndBoundaryConditions<2> pde_and_bc2(&pde2, &bc, false);
        pde_handler.AddPdeAndBc(&pde_and_bc2);

        // Pass PDE handler to simulation
        pde_handler.UseCoarsePdeMesh(10.0, 50.0);
        pde_handler.SetImposeBcsOnCoarseBoundary(false);
        simulator.SetCellBasedPdeHandler(&pde_handler);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetCutOffLength(1.5);
        simulator.AddForce(p_linear_force);

        // Set up cell killer and pass into simulation
        MAKE_PTR_ARGS(OxygenBasedCellKiller<2>, p_killer, (&cell_population));
        simulator.AddCellKiller(p_killer);

        // Find centre of cell population
        c_vector<double,2> centre_of_cell_population = cell_population.GetCentroidOfCellPopulation();

        // Find centre of coarse PDE mesh
        c_vector<double,2> centre_of_coarse_pde_mesh = zero_vector<double>(2);

        TetrahedralMesh<2,2>* p_coarse_mesh = simulator.GetCellBasedPdeHandler()->GetCoarsePdeMesh();

        for (unsigned i=0; i<p_coarse_mesh->GetNumNodes(); i++)
        {
            centre_of_coarse_pde_mesh += p_coarse_mesh->GetNode(i)->rGetLocation();
        }
        centre_of_coarse_pde_mesh /= p_coarse_mesh->GetNumNodes();

        // Test that the two centres match
        c_vector<double,2> centre_diff = centre_of_cell_population - centre_of_coarse_pde_mesh;
        TS_ASSERT_DELTA(norm_2(centre_diff), 0.0, 1e-4);

        // Test FindCoarseElementContainingCell() and initialisation of mCellPdeElementMap
        simulator.GetCellBasedPdeHandler()->InitialiseCellPdeElementMap();
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            unsigned containing_element_index = simulator.GetCellBasedPdeHandler()->mCellPdeElementMap[*cell_iter];
            TS_ASSERT_LESS_THAN(containing_element_index, p_coarse_mesh->GetNumElements());
            TS_ASSERT_EQUALS(containing_element_index, simulator.GetCellBasedPdeHandler()->FindCoarseElementContainingCell(*cell_iter));
        }

        // Run cell-based simulation
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());

        TS_ASSERT(p_coarse_mesh != NULL);

        ReplicatableVector pde_solution0(simulator.GetCellBasedPdeHandler()->GetPdeSolution(0));
        ReplicatableVector pde_solution1(simulator.GetCellBasedPdeHandler()->GetPdeSolution(1));

        TS_ASSERT_EQUALS(pde_solution0.GetSize(),pde_solution1.GetSize());

        // Test the nutrient concentration is equal to 1.0 at each coarse mesh node far from the cells
        for (unsigned i=0; i<pde_solution0.GetSize(); i++)
        {
            c_vector<double,2> centre;
            centre(0) = 2.5; // assuming 5 by 5 honeycomb mesh
            centre(1) = 2.5;
            c_vector<double,2> posn = p_coarse_mesh->GetNode(i)->rGetLocation();
            double dist = norm_2(centre - posn);
            double u0 = pde_solution0[i];
            double u1 = pde_solution1[i];

            if (dist > 4.0)
            {
                TS_ASSERT_DELTA(u0, 1.0, 1e-5);
                TS_ASSERT_DELTA(u1, 1.0, 1e-5);
            }
        }

        /*
         * Loop over cells, find the coarse mesh element containing it, then
         * check the interpolated PDE solution is between the min and max of
         * the PDE solution on the nodes of that element.
         */
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
            cell_iter != cell_population.End();
            ++cell_iter)
        {
            unsigned elem_index = simulator.GetCellBasedPdeHandler()->GetCoarsePdeMesh()->GetContainingElementIndex(cell_population.GetLocationOfCellCentre(*cell_iter));
            Element<2,2>* p_element = p_coarse_mesh->GetElement(elem_index);

            double max0 = std::max(pde_solution0[p_element->GetNodeGlobalIndex(0)], pde_solution0[p_element->GetNodeGlobalIndex(1)]);
            max0 = std::max(max0, pde_solution0[p_element->GetNodeGlobalIndex(2)]);

            double max1 = std::max(pde_solution1[p_element->GetNodeGlobalIndex(0)], pde_solution1[p_element->GetNodeGlobalIndex(1)]);
            max1 = std::max(max1, pde_solution1[p_element->GetNodeGlobalIndex(2)]);

            double min0 = std::min(pde_solution0[p_element->GetNodeGlobalIndex(0)], pde_solution0[p_element->GetNodeGlobalIndex(1)]);
            min0 = std::min(min0, pde_solution0[p_element->GetNodeGlobalIndex(2)]);

            double min1 = std::min(pde_solution1[p_element->GetNodeGlobalIndex(0)], pde_solution1[p_element->GetNodeGlobalIndex(1)]);
            min1 = std::min(min1, pde_solution1[p_element->GetNodeGlobalIndex(2)]);

            double value0_at_cell = CellwiseData<2>::Instance()->GetValue(*cell_iter, 0);
            double value1_at_cell = CellwiseData<2>::Instance()->GetValue(*cell_iter, 1);

            TS_ASSERT_LESS_THAN_EQUALS(value1_at_cell, value0_at_cell);
            TS_ASSERT_LESS_THAN_EQUALS(min0, value0_at_cell + DBL_EPSILON);
            TS_ASSERT_LESS_THAN_EQUALS(value0_at_cell, max0 + DBL_EPSILON);
            TS_ASSERT_LESS_THAN_EQUALS(min1, value1_at_cell + DBL_EPSILON);
            TS_ASSERT_LESS_THAN_EQUALS(value1_at_cell, max1 + DBL_EPSILON);
        }

        // Tidy up
        CellwiseData<2>::Destroy();
    }

    void TestArchivingWithSimplePde() throw (Exception)
    {
        EXIT_IF_PARALLEL;

        // Set up mesh
        HoneycombMeshGenerator generator(5, 5, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Set up cells
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            SimpleOxygenBasedCellCycleModel* p_model = new SimpleOxygenBasedCellCycleModel();
            p_model->SetDimension(2);
            p_model->SetCellProliferativeType(STEM);

            // Use non-default G1 durations
            p_model->SetStemCellG1Duration(8.0);
            p_model->SetTransitCellG1Duration(8.0);

            CellPtr p_cell(new Cell(p_state, p_model));
            double birth_time = -1.0 - ((double) i/p_mesh->GetNumNodes())*18.0;
            p_cell->SetBirthTime(birth_time);

            cells.push_back(p_cell);
        }

        // Set up cell population
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up CellwiseData and associate it with the cell population
        CellwiseData<2>* p_data = CellwiseData<2>::Instance();
        p_data->SetNumCellsAndVars(cell_population.GetNumRealCells(), 1);
        p_data->SetCellPopulation(&cell_population);

        // Since values are first passed in to CellwiseData before it is updated in UpdateAtEndOfTimeStep(),
        // we need to pass it some initial conditions to avoid memory errors
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            p_data->SetValue(1.0, p_mesh->GetNode(i)->GetIndex());
        }

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("OffLatticeSimulationWithPdesSaveAndLoad");
        simulator.SetEndTime(0.2);

        // Set up PDE and pass to simulation via handler
        SimpleUniformSourcePde<2> pde(-0.1);
        ConstBoundaryCondition<2> bc(1.0);
        PdeAndBoundaryConditions<2> pde_and_bc(&pde, &bc, false);

        CellBasedPdeHandler<2> pde_handler(&cell_population);
        pde_handler.AddPdeAndBc(&pde_and_bc);

        pde_handler.SetImposeBcsOnCoarseBoundary(false);
        simulator.SetCellBasedPdeHandler(&pde_handler);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetCutOffLength(1.5);
        simulator.AddForce(p_linear_force);

        // Set up cell killer and pass into simulation
        MAKE_PTR_ARGS(OxygenBasedCellKiller<2>, p_killer, (&cell_population));
        simulator.AddCellKiller(p_killer);

        // Run cell-based simulation
        simulator.Solve();

        // Save cell-based simulation
        CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Save(&simulator);

        OffLatticeSimulation<2>* p_simulator
            = CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Load("OffLatticeSimulationWithPdesSaveAndLoad", 0.2);

        p_simulator->SetEndTime(0.5);
        p_simulator->Solve();

        // These results are from time 0.5 in TestWithOxygen.
        std::vector<double> node_5_location = p_simulator->GetNodeLocation(5);
        TS_ASSERT_DELTA(node_5_location[0], 0.4968, 1e-4);
        TS_ASSERT_DELTA(node_5_location[1], 0.8635, 1e-4);

        std::vector<double> node_15_location = p_simulator->GetNodeLocation(15);
        TS_ASSERT_DELTA(node_15_location[0], 0.4976, 1e-4);
        TS_ASSERT_DELTA(node_15_location[1], 2.5977, 1e-4);

        // Test CellwiseData was set up correctly
        TS_ASSERT_EQUALS(CellwiseData<2>::Instance()->IsSetUp(), true);

        // Test the CellwiseData result
        TS_ASSERT_DELTA(p_data->GetValue(p_simulator->rGetCellPopulation().GetCellUsingLocationIndex(5)), 0.9604, 1e-4);
        TS_ASSERT_DELTA(p_data->GetValue(p_simulator->rGetCellPopulation().GetCellUsingLocationIndex(15)), 0.9584, 1e-4);

        // Run cell-based simulation
        delete p_simulator;
        CellwiseData<2>::Destroy();
    }

    /**
     * This test demonstrates how to archive a OffLatticeSimulation
     * in the case where the PDE has the cell population as a member variable.
     */
    void TestArchivingWithCellwisePde() throw (Exception)
    {
        EXIT_IF_PARALLEL;

        std::string output_directory = "TestArchivingWithCellwisePde";
        double end_time = 0.1;

        // Set up mesh
        HoneycombMeshGenerator generator(5, 5, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Set up cells
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            SimpleOxygenBasedCellCycleModel* p_model = new SimpleOxygenBasedCellCycleModel();
            p_model->SetDimension(2);
            p_model->SetCellProliferativeType(STEM);

            // Use non-default G1 durations
            p_model->SetStemCellG1Duration(8.0);
            p_model->SetTransitCellG1Duration(8.0);

            CellPtr p_cell(new Cell(p_state, p_model));
            double birth_time = -1.0 - ((double) i/p_mesh->GetNumNodes())*18.0;
            p_cell->SetBirthTime(birth_time);

            cells.push_back(p_cell);
        }

        // Set up cell population
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up CellwiseData and associate it with the cell population
        CellwiseData<2>* p_data = CellwiseData<2>::Instance();
        p_data->SetNumCellsAndVars(cell_population.GetNumRealCells(), 1);
        p_data->SetCellPopulation(&cell_population);
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            p_data->SetValue(1.0, p_mesh->GetNode(i)->GetIndex());
        }

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory(output_directory);
        simulator.SetEndTime(end_time);

        // Set up PDE and pass to simulation via handler
        CellwiseSourcePde<2> pde(cell_population, -0.03);
        ConstBoundaryCondition<2> bc(1.0);
        PdeAndBoundaryConditions<2> pde_and_bc(&pde, &bc, false);

        CellBasedPdeHandler<2> pde_handler(&cell_population);
        pde_handler.AddPdeAndBc(&pde_and_bc);
        pde_handler.SetImposeBcsOnCoarseBoundary(false);

        simulator.SetCellBasedPdeHandler(&pde_handler);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetCutOffLength(3.0);
        simulator.AddForce(p_linear_force);

        // Run cell-based simulation
        simulator.Solve();

        // Save cell-based simulation
        CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Save(&simulator);

        // Load simulation
        OffLatticeSimulation<2>* p_simulator
            = CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Load(output_directory, end_time);

        p_simulator->SetEndTime(2.0*end_time);

        // Run cell-based simulation
        TS_ASSERT_THROWS_NOTHING(p_simulator->Solve());

        // Tidy up
        delete p_simulator;
        CellwiseData<2>::Destroy();
    }

    void Test3DOffLatticeSimulationWithPdes() throw(Exception)
    {
        EXIT_IF_PARALLEL;

        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_136_elements");
        MutableMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        TrianglesMeshWriter<3,3> mesh_writer("TestSolveMethodSpheroidSimulation3DMesh", "StartMesh");
        mesh_writer.WriteFilesUsingMesh(mesh);

        // Set up cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 3> generator;
        generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Set some model parameters for the cell-cycle model
        for (unsigned index=0; index < cells.size(); index++)
        {
            cells[index]->GetCellCycleModel()->SetTransitCellG1Duration(8.0);
            cells[index]->GetCellCycleModel()->SetStemCellG1Duration(8.0);
        }

        // Set up cell population
        MeshBasedCellPopulation<3> cell_population(mesh, cells);

        // Set up CellwiseData and associate it with the cell population
        CellwiseData<3>* p_data = CellwiseData<3>::Instance();
        p_data->SetNumCellsAndVars(cell_population.GetNumRealCells(), 1);
        p_data->SetCellPopulation(&cell_population);
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            p_data->SetValue(1.0, mesh.GetNode(i)->GetIndex());
        }

        // Set up cell-based simulation
        OffLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory("OffLatticeSimulationWithOxygen3d");
        simulator.SetSamplingTimestepMultiple(12);
        simulator.SetEndTime(0.1);

        // Set up PDE and pass to simulation via handler
        SimpleUniformSourcePde<3> pde(-0.03);
        ConstBoundaryCondition<3> bc(1.0);
        PdeAndBoundaryConditions<3> pde_and_bc(&pde, &bc, false);

        CellBasedPdeHandler<3> pde_handler(&cell_population);
        pde_handler.AddPdeAndBc(&pde_and_bc);
        pde_handler.SetImposeBcsOnCoarseBoundary(false);

        simulator.SetCellBasedPdeHandler(&pde_handler);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<3>, p_linear_force);
        p_linear_force->SetCutOffLength(1.5);
        simulator.AddForce(p_linear_force);

        // Set up cell killer and pass into simulation
        MAKE_PTR_ARGS(OxygenBasedCellKiller<3>, p_killer, (&cell_population));
        simulator.AddCellKiller(p_killer);

        // Run cell-based simulation
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());

        // Tidy up
        CellwiseData<3>::Destroy();
    }

    /**
     * Test that the simulation runs correctly when a PDE with time-dependent
     * boundary condition is used. Here we test the PDE solver using the problem
     * del squared C = 1 on the unit disc, with boundary condition C=t on r=1,
     * which has analytic solution C = t - 0.25*(1-r^2).
     */
    void TestWithBoundaryConditionVaryingInTime() throw(Exception)
    {
        EXIT_IF_PARALLEL;

        // Set up mesh
        MutableMesh<2,2> mesh;
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_522_elements");
        mesh.ConstructFromMeshReader(mesh_reader);

        // Set up cells
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            SimpleOxygenBasedCellCycleModel* p_model = new SimpleOxygenBasedCellCycleModel();
            p_model->SetDimension(2);
            p_model->SetCellProliferativeType(DIFFERENTIATED);
            p_model->SetHypoxicConcentration(0.9);
            p_model->SetQuiescentConcentration(0.9);

            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetBirthTime(-1.0);
            cells.push_back(p_cell);
        }

        // Set up cell population
        MeshBasedCellPopulation<2> cell_population(mesh, cells);

        // Set up CellwiseData and associate it with the cell population
        CellwiseData<2>* p_data = CellwiseData<2>::Instance();
        p_data->SetNumCellsAndVars(cell_population.GetNumRealCells(), 1);
        p_data->SetCellPopulation(&cell_population);
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            p_data->SetValue(1.0, mesh.GetNode(i)->GetIndex());
        }

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestUpdateAtEndOfTimeStepMethod");

        // Create PDE and pass to simulation via handler
        SimplePdeForTesting pde;
        FunctionalBoundaryCondition<2> functional_bc(&bc_func);
        PdeAndBoundaryConditions<2> pde_and_bc(&pde, &functional_bc, false);

        CellBasedPdeHandler<2> pde_handler(&cell_population);
        pde_handler.AddPdeAndBc(&pde_and_bc);
        pde_handler.SetImposeBcsOnCoarseBoundary(false);

        simulator.SetCellBasedPdeHandler(&pde_handler);

        double end_time = 0.5;
        simulator.SetEndTime(end_time);

        // Do not pass in a force law so that no cells interact mechanically

        // Run cell-based simulation
        simulator.SetSamplingTimestepMultiple(12);
        simulator.Solve();

        // Check the correct solution was obtained
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            double radius = norm_2(cell_population.GetLocationOfCellCentre(*cell_iter));
            double analytic_solution = end_time - 0.25*(1 - pow(radius,2.0));

            // Test that PDE solver is working correctly
            TS_ASSERT_DELTA(p_data->GetValue(*cell_iter), analytic_solution, 0.02);
        }

        // Tidy up
        CellwiseData<2>::Destroy();
    }

    void TestOffLatticeSimulationWithPdesParameterOutputMethods() throw (Exception)
    {
        EXIT_IF_PARALLEL;

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

        // Set up PDE and pass to simulation via handler
        CellwiseSourcePde<2> pde(cell_population, -0.03);
        ConstBoundaryCondition<2> bc(1.0);
        PdeAndBoundaryConditions<2> pde_and_bc(&pde, &bc, false);

        CellBasedPdeHandler<2> pde_handler(&cell_population);
        pde_handler.AddPdeAndBc(&pde_and_bc);
        simulator.SetCellBasedPdeHandler(&pde_handler);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetCutOffLength(1.5);
        simulator.AddForce(p_linear_force);

        std::string output_directory = "TestOffLatticeSimulationOutputParameters";
        OutputFileHandler output_file_handler(output_directory, false);
        out_stream parameter_file = output_file_handler.OpenOutputFile("cell_based_sim_with_pde_results.parameters");
        simulator.OutputSimulationParameters(parameter_file);
        parameter_file->close();

        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "cell_based_sim_with_pde_results.parameters  cell_based/test/data/TestOffLatticeSimulationOutputParameters/cell_based_sim_with_pde_results.parameters").c_str()), 0);

        ///\todo check output of simulator.OutputSimulationSetup();
    }

    void TestNodeBasedWithoutCoarseMeshThrowsException() throw(Exception)
    {
        EXIT_IF_PARALLEL;

        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_522_elements");
        TetrahedralMesh<2,2> temp_mesh;
        temp_mesh.ConstructFromMeshReader(mesh_reader);
        temp_mesh.Scale(5.0,1.0);

        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(temp_mesh);

        // Set up cells
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
            p_model->SetDimension(2);
            p_model->SetCellProliferativeType(DIFFERENTIATED);

            CellPtr p_cell(new Cell(p_state, p_model));
            double birth_time = -RandomNumberGenerator::Instance()->ranf()*18.0;
            p_cell->SetBirthTime(birth_time);

            cells.push_back(p_cell);
        }

        // Set up cell population
        NodeBasedCellPopulation<2> cell_population(mesh, cells);
        cell_population.SetMechanicsCutOffLength(1.5);

        // Set up CellwiseData and associate it with the cell population
        CellwiseData<2>* p_data = CellwiseData<2>::Instance();
        p_data->SetNumCellsAndVars(cell_population.GetNumRealCells(), 1);
        p_data->SetCellPopulation(&cell_population);

        /*
         * Since values are first passed in to CellwiseData before it is updated in UpdateAtEndOfTimeStep(),
         * we need to pass it some initial conditions to avoid memory errors.
         */
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            p_data->SetValue(5.0, mesh.GetNode(i)->GetIndex(), 0);
        }

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestNodeBasedCellPopulationWithpoutCoarseMeshThrowsException");
        simulator.SetEndTime(0.01);

        // Set up PDE and pass to simulation via handler
        AveragedSourcePde<2> pde(cell_population, -1.0);
        ConstBoundaryCondition<2> bc(1.0);
        PdeAndBoundaryConditions<2> pde_and_bc(&pde, &bc, false);

        CellBasedPdeHandler<2> pde_handler(&cell_population);
        pde_handler.AddPdeAndBc(&pde_and_bc);
        pde_handler.SetImposeBcsOnCoarseBoundary(false);
        simulator.SetCellBasedPdeHandler(&pde_handler);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetCutOffLength(1.5);
        simulator.AddForce(p_linear_force);

        TS_ASSERT_THROWS_THIS(simulator.Solve(), "Trying to solve a PDE on a cell population that doesn't have a mesh. Try calling UseCoarsePdeMesh().");

        // Tidy up
        CellwiseData<2>::Destroy();
    }

    void TestNodeBasedWithCoarseMesh() throw(Exception)
    {
        EXIT_IF_PARALLEL;

        // Create a simple mesh - larger to use coarse graining.
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_522_elements");
        TetrahedralMesh<2,2> temp_mesh;
        temp_mesh.ConstructFromMeshReader(mesh_reader);
        temp_mesh.Scale(20.0,20.0);

        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(temp_mesh);

        // Set up cells
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
            p_model->SetDimension(2);
            p_model->SetCellProliferativeType(DIFFERENTIATED);

            CellPtr p_cell(new Cell(p_state, p_model));
            double birth_time = -RandomNumberGenerator::Instance()->ranf()*18.0;
            p_cell->SetBirthTime(birth_time);

            cells.push_back(p_cell);
        }

        // Set up cell population
        NodeBasedCellPopulation<2> cell_population(mesh, cells);
        cell_population.SetMechanicsCutOffLength(1.5);

        // Set up CellwiseData and associate it with the cell population
        CellwiseData<2>* p_data = CellwiseData<2>::Instance();
        p_data->SetNumCellsAndVars(cell_population.GetNumRealCells(), 1);
        p_data->SetCellPopulation(&cell_population);

        /*
         * Since values are first passed in to CellwiseData before it is updated in UpdateAtEndOfTimeStep(),
         * we need to pass it some initial conditions to avoid memory errors.
         */
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            p_data->SetValue(1.0, mesh.GetNode(i)->GetIndex(), 0);
        }

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestNodeBasedCellPopulationWithCoarseMeshPDE");
        simulator.SetEndTime(0.01);

        // Set up PDE and pass to simulation via handler (zero uptake to check analytic solution)
        AveragedSourcePde<2> pde(cell_population, 0.0);
        ConstBoundaryCondition<2> bc(1.0);
        PdeAndBoundaryConditions<2> pde_and_bc(&pde, &bc, false);

        CellBasedPdeHandler<2> pde_handler(&cell_population);
        pde_handler.AddPdeAndBc(&pde_and_bc);
        pde_handler.UseCoarsePdeMesh(10.0, 50.0);
        pde_handler.SetImposeBcsOnCoarseBoundary(false);

        simulator.SetCellBasedPdeHandler(&pde_handler);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetCutOffLength(1.5);
        simulator.AddForce(p_linear_force);

        // Solve the system
        simulator.Solve();

        // Test solution is constant
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {

            double analytic_solution = 1.0;
            // Test that PDE solver is working correctly
            TS_ASSERT_DELTA(p_data->GetValue(*cell_iter), analytic_solution, 1e-2);
        }

        // Find centre of cell population
        c_vector<double,2> centre_of_cell_population = zero_vector<double>(2);

        for (unsigned i=0; i<simulator.rGetCellPopulation().GetNumNodes(); i++)
        {
            centre_of_cell_population += simulator.rGetCellPopulation().GetNode(i)->rGetLocation();
        }
        centre_of_cell_population /= simulator.rGetCellPopulation().GetNumNodes();

        // Find centre of coarse PDE mesh
        c_vector<double,2> centre_of_coarse_pde_mesh = zero_vector<double>(2);
        TetrahedralMesh<2,2>* p_coarse_mesh = simulator.GetCellBasedPdeHandler()->GetCoarsePdeMesh();
        for (unsigned i=0; i<p_coarse_mesh->GetNumNodes(); i++)
        {
            centre_of_coarse_pde_mesh += p_coarse_mesh->GetNode(i)->rGetLocation();
        }
        centre_of_coarse_pde_mesh /= p_coarse_mesh->GetNumNodes();

        // Test that the two centres match
        c_vector<double,2> centre_diff = centre_of_cell_population - centre_of_coarse_pde_mesh;
        TS_ASSERT_DELTA(norm_2(centre_diff), 0.0, 1e-4);

        // Test FindCoarseElementContainingCell() and initialisation of mCellPdeElementMap
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
            cell_iter != cell_population.End();
            ++cell_iter)
        {
            unsigned containing_element_index = simulator.GetCellBasedPdeHandler()->mCellPdeElementMap[*cell_iter];
            TS_ASSERT_LESS_THAN(containing_element_index, p_coarse_mesh->GetNumElements());
            TS_ASSERT_EQUALS(containing_element_index, simulator.GetCellBasedPdeHandler()->FindCoarseElementContainingCell(*cell_iter));
        }

        // Tidy up
        CellwiseData<2>::Destroy();
    }

    void TestVolumeDependentAveragedPde() throw(Exception)
    {
        EXIT_IF_PARALLEL;

        // Create a simple mesh
        TetrahedralMesh<2,2> temp_mesh;
        temp_mesh.ConstructRegularSlabMesh(2.0,2,2);

        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(temp_mesh);

        // Set some cell radius
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            mesh.SetCellRadius(i, 0.9);
        }

        // Set up cells
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
            p_model->SetDimension(2);
            p_model->SetCellProliferativeType(DIFFERENTIATED);

            CellPtr p_cell(new Cell(p_state, p_model));
            double birth_time = -RandomNumberGenerator::Instance()->ranf()*18.0;
            p_cell->SetBirthTime(birth_time);

            cells.push_back(p_cell);
        }

        // Set up cell population
        NodeBasedCellPopulation<2> cell_population(mesh, cells);
        cell_population.SetMechanicsCutOffLength(1.5);

        // Set up CellwiseData and associate it with the cell population
        CellwiseData<2>* p_data = CellwiseData<2>::Instance();
        p_data->SetNumCellsAndVars(cell_population.GetNumRealCells(), 1);
        p_data->SetCellPopulation(&cell_population);

        c_vector<double,2> centre_of_mesh;
        centre_of_mesh[0] = 5.0;
        centre_of_mesh[1] = 5.0;

        /*
         * Since values are first passed in to CellwiseData before it is updated in UpdateAtEndOfTimeStep(),
         * we need to pass it some initial conditions to avoid memory errors.
         */
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            c_vector<double,2> location = mesh.GetNode(i)->rGetLocation();

            double initial_condition = 0.0;
            if (norm_2(location-centre_of_mesh) >= 1.0)
            {
                initial_condition = 1.0;
            }

            p_data->SetValue(initial_condition, mesh.GetNode(i)->GetIndex(),0);
        }

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestCoarsePdeSolutionOnNodeBased");
        simulator.SetEndTime(0.01);

        // Set up PDE and pass to simulation via handler (uniform uptake at each cell)
        VolumeDependentAveragedSourcePde<2> pde(cell_population, -0.01);
        ConstBoundaryCondition<2> bc(1.0);
        PdeAndBoundaryConditions<2> pde_and_bc(&pde, &bc, false);

        CellBasedPdeHandler<2> pde_handler(&cell_population);
        pde_handler.AddPdeAndBc(&pde_and_bc);

        pde_handler.UseCoarsePdeMesh(10.0, 50.0);
        pde_handler.SetImposeBcsOnCoarseBoundary(false);
        simulator.SetCellBasedPdeHandler(&pde_handler);

        // Solve the system
        simulator.Solve();

        TetrahedralMesh<2,2>* p_coarse_mesh = simulator.GetCellBasedPdeHandler()->GetCoarsePdeMesh();

        // Check the correct cell density is in each element
        unsigned num_elements_in_coarse_mesh = p_coarse_mesh->GetNumElements();
        std::vector<unsigned> cells_in_each_coarse_element(num_elements_in_coarse_mesh);
        for (unsigned i=0; i<num_elements_in_coarse_mesh; i++)
        {
            cells_in_each_coarse_element[i] = 0;
        }

        // Find out how many cells lie in each element
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
            cell_iter != cell_population.End();
            ++cell_iter)
        {
            // Get containing element
            unsigned containing_element_index = simulator.GetCellBasedPdeHandler()->mCellPdeElementMap[*cell_iter];
            cells_in_each_coarse_element[containing_element_index] += 1;
        }

        for (unsigned i=0; i<num_elements_in_coarse_mesh; i++)
        {
            c_matrix<double, 2, 2> jacobian;
            double det;
            p_coarse_mesh->GetElement(i)->CalculateJacobian(jacobian, det);
            double element_volume = p_coarse_mesh->GetElement(i)->GetVolume(det);
            double uptake_rate = pde.GetUptakeRateForElement(i);
            double expected_uptake = 0.9*0.9*(cells_in_each_coarse_element[i]/element_volume);

            TS_ASSERT_DELTA(uptake_rate, expected_uptake, 1e-4);
        }

        // Tidy up
        CellwiseData<2>::Destroy();
    }

    void TestCoarsePdeSolutionOnNodeBased1d() throw(Exception)
    {
        EXIT_IF_PARALLEL;

        // Create mesh
        std::vector<Node<1>*> nodes;
        nodes.push_back(new Node<1>(0, true,  0.0));

        NodesOnlyMesh<1> mesh;
        mesh.ConstructNodesWithoutMesh(nodes);

        // Set up differentiated cells
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
            p_model->SetDimension(1);
            p_model->SetCellProliferativeType(DIFFERENTIATED);

            CellPtr p_cell(new Cell(p_state, p_model));
            double birth_time = -RandomNumberGenerator::Instance()->ranf()*18.0;
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        // Set up cell population
        NodeBasedCellPopulation<1> cell_population(mesh, cells);
        cell_population.SetMechanicsCutOffLength(1.5);

        // Set up CellwiseData and associate it with the cell population
        CellwiseData<1>* p_data = CellwiseData<1>::Instance();
        p_data->SetNumCellsAndVars(cell_population.GetNumRealCells(), 1);
        p_data->SetCellPopulation(&cell_population);

        // Since values are first passed in to CellwiseData before it is updated in UpdateAtEndOfTimeStep(),
        // we need to pass it some initial conditions to avoid memory errors
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            p_data->SetValue(1.0, mesh.GetNode(i)->GetIndex(),0);
        }

        // Set up cell-based simulation
        OffLatticeSimulation<1> simulator(cell_population);
        simulator.SetEndTime(0.01);

        // Set up PDE and pass to simulation via handler
        AveragedSourcePde<1> pde(cell_population, -1.0);
        ConstBoundaryCondition<1> bc(0.0);
        PdeAndBoundaryConditions<1> pde_and_bc(&pde, &bc, false);

        CellBasedPdeHandler<1> pde_handler(&cell_population);
        pde_handler.AddPdeAndBc(&pde_and_bc);

        unsigned num_nodes = pow(2,3);
        double mesh_size = 2.0/(double)(num_nodes-1);
        pde_handler.UseCoarsePdeMesh(mesh_size, 2.0);

        simulator.SetCellBasedPdeHandler(&pde_handler);

        simulator.SetOutputDirectory("TestCoarsePdeSolutionOnNodeBased1d");

        // Solve the system
        simulator.Solve();

        // Tidy up
        CellwiseData<1>::Destroy();

        // When the mesh goes out of scope, then it's a different set of nodes that get destroyed
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

    void TestCoarsePdeSolutionOnNodeBased2d() throw(Exception)
    {
        EXIT_IF_PARALLEL;

        // Create a simple mesh
        TetrahedralMesh<2,2> temp_mesh;
        temp_mesh.ConstructRectangularMesh(10,10);

        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(temp_mesh);

        // Set up cells
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
            p_model->SetDimension(2);
            p_model->SetCellProliferativeType(DIFFERENTIATED);

            CellPtr p_cell(new Cell(p_state, p_model));
            double birth_time = -RandomNumberGenerator::Instance()->ranf()*18.0;
            p_cell->SetBirthTime(birth_time);

            cells.push_back(p_cell);
        }

        // Set up cell population
        NodeBasedCellPopulation<2> cell_population(mesh, cells);
        cell_population.SetMechanicsCutOffLength(1.5);

        // Set up CellwiseData and associate it with the cell population
        CellwiseData<2>* p_data = CellwiseData<2>::Instance();
        p_data->SetNumCellsAndVars(cell_population.GetNumRealCells(), 1);
        p_data->SetCellPopulation(&cell_population);

        c_vector<double,2> centre_of_mesh;
        centre_of_mesh[0] = 5.0;
        centre_of_mesh[1] = 5.0;

        /*
         * Since values are first passed in to CellwiseData before it is updated in UpdateAtEndOfTimeStep(),
         * we need to pass it some initial conditions to avoid memory errors.
         */
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            c_vector<double,2> location = mesh.GetNode(i)->rGetLocation();

            double initial_condition = 0.0;
            if (norm_2(location-centre_of_mesh) >= 1.0)
            {
                initial_condition = 1.0;
            }

            p_data->SetValue(initial_condition, mesh.GetNode(i)->GetIndex(),0);
        }

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestCoarsePdeSolutionOnNodeBased2d");
        simulator.SetEndTime(0.01);

        // Set up PDE and pass to simulation via handler (uniform uptake at each cell)
        AveragedSourcePde<2> pde(cell_population, -0.01);
        ConstBoundaryCondition<2> bc(1.0);
        PdeAndBoundaryConditions<2> pde_and_bc(&pde, &bc, false);

        CellBasedPdeHandler<2> pde_handler(&cell_population);
        pde_handler.AddPdeAndBc(&pde_and_bc);
        pde_handler.UseCoarsePdeMesh(10.0, 50.0);
        pde_handler.SetImposeBcsOnCoarseBoundary(false);

        simulator.SetCellBasedPdeHandler(&pde_handler);

        // Solve the system
        simulator.Solve();

        // Test solution is unchanged at first two cells
        TS_ASSERT_DELTA(p_data->GetValue(*(cell_population.Begin())), 0.9543, 1e-4);
        TS_ASSERT_DELTA(p_data->GetValue(*(++cell_population.Begin())), 0.9589, 1e-4);

        // Tidy up
        CellwiseData<2>::Destroy();
    }

    void TestCoarsePdeSolutionOnNodeBased3d() throw(Exception)
    {
        EXIT_IF_PARALLEL;

        // Create a simple mesh
        TetrahedralMesh<3,3> temp_mesh;
        temp_mesh.ConstructCuboid(5,5,5);

        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(temp_mesh);

        // Set up cells
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
            p_model->SetDimension(3);
            p_model->SetCellProliferativeType(DIFFERENTIATED);

            CellPtr p_cell(new Cell(p_state, p_model));
            double birth_time = -RandomNumberGenerator::Instance()->ranf()*18.0;
            p_cell->SetBirthTime(birth_time);

            cells.push_back(p_cell);
        }

        // Set up cell population
        NodeBasedCellPopulation<3> cell_population(mesh, cells);
        cell_population.SetMechanicsCutOffLength(1.5);

        // Set up CellwiseData and associate it with the cell population
        CellwiseData<3>* p_data = CellwiseData<3>::Instance();
        p_data->SetNumCellsAndVars(cell_population.GetNumRealCells(), 1);
        p_data->SetCellPopulation(&cell_population);

        c_vector<double,3> centre_of_mesh;
        centre_of_mesh[0] = 5.0;
        centre_of_mesh[1] = 5.0;
        centre_of_mesh[2] = 5.0;

        /*
         * Since values are first passed in to CellwiseData before it is updated in UpdateAtEndOfTimeStep(),
         * we need to pass it some initial conditions to avoid memory errors.
         */
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            c_vector<double,3> location = mesh.GetNode(i)->rGetLocation();

            double initial_condition = 0.0;
            if (norm_2(location-centre_of_mesh) >= 1.0)
            {
                initial_condition = 1.0;
            }

            p_data->SetValue(initial_condition, mesh.GetNode(i)->GetIndex(),0);
        }

        // Set up cell-based simulation
        OffLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory("TestCoarsePdeSolutionOnNodeBased3d");
        simulator.SetEndTime(0.01);

        // Set up PDE and pass to simulation via handler (uniform uptake at each cell)
        AveragedSourcePde<3> pde(cell_population, -0.01);
        ConstBoundaryCondition<3> bc(1.0);
        PdeAndBoundaryConditions<3> pde_and_bc(&pde, &bc, false);

        CellBasedPdeHandler<3> pde_handler(&cell_population);
        pde_handler.AddPdeAndBc(&pde_and_bc);
        pde_handler.UseCoarsePdeMesh(10.0, 50.0);
        pde_handler.SetImposeBcsOnCoarseBoundary(false);
        simulator.SetCellBasedPdeHandler(&pde_handler);

        // Solve the system
        simulator.Solve();

        // Test solution is unchanged at first two cells
        TS_ASSERT_DELTA(p_data->GetValue(*(cell_population.Begin())), 1.000, 1e-4);
        TS_ASSERT_DELTA(p_data->GetValue(*(++cell_population.Begin())), 1.000, 1e-4);

        // Tidy up
        CellwiseData<3>::Destroy();
    }
};

#endif /*TESTOFFLATTICESIMULATIONWITHPDES_HPP_*/
