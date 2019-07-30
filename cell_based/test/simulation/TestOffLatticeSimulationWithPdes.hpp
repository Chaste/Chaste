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

#ifndef TESTOFFLATTICESIMULATIONWITHPDES_HPP_
#define TESTOFFLATTICESIMULATIONWITHPDES_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "NagaiHondaForce.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "MutableVertexMesh.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "OffLatticeSimulation.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "ApoptoticCellKiller.hpp"
#include "SimpleOxygenBasedCellCycleModel.hpp"
#include "CellsGenerator.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "ApoptoticCellProperty.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "UniformSourceEllipticPde.hpp"
#include "CellwiseSourceEllipticPde.hpp"
#include "ConstBoundaryCondition.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "ReplicatableVector.hpp"
#include "NumericFileComparison.hpp"
#include "WildTypeCellMutationState.hpp"
#include "FunctionalBoundaryCondition.hpp"
#include "AveragedSourceEllipticPde.hpp"
#include "VolumeDependentAveragedSourceEllipticPde.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "SmartPointers.hpp"
#include "FileComparison.hpp"
#include "CellPopulationAreaWriter.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "EllipticBoxDomainPdeModifier.hpp"
#include "EllipticGrowingDomainPdeModifier.hpp"
#include "RadialCellDataDistributionWriter.hpp"

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

///\todo move into cell_based/test/cell_based_pde
///\todo merge content into TestSimulationsWith*DomainPdeModifier and remove this test suite
class TestOffLatticeSimulationWithPdes : public AbstractCellBasedWithTimingsTestSuite
{
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
    void TestUpdateAtEndOfTimeStep()
    {
        EXIT_IF_PARALLEL;

        // Set up mesh
        MutableMesh<2,2> mesh;
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_522_elements");
        mesh.ConstructFromMeshReader(mesh_reader);

        // Set up cells
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            SimpleOxygenBasedCellCycleModel* p_model = new SimpleOxygenBasedCellCycleModel();
            p_model->SetDimension(2);
            p_model->SetHypoxicConcentration(0.9);
            p_model->SetQuiescentConcentration(0.9);

            // Use non-default G1 durations
            p_model->SetStemCellG1Duration(8.0);
            p_model->SetTransitCellG1Duration(8.0);

            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetCellProliferativeType(p_stem_type);

            double birth_time = -RandomNumberGenerator::Instance()->ranf()*18.0;
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        // Set up cell population
        MeshBasedCellPopulation<2> cell_population(mesh, cells);

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestUpdateAtEndOfTimeStepMethod");
        simulator.SetEndTime(2.0/120.0);

        // Create PDE and boundary condition objects
        MAKE_PTR(SimplePdeForTesting, p_pde);
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(EllipticGrowingDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false));
        p_pde_modifier->SetDependentVariableName("oxygen");

        simulator.AddSimulationModifier(p_pde_modifier);

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
        MAKE_PTR_ARGS(ApoptoticCellKiller<2>, p_killer, (&cell_population));
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
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("oxygen"), analytic_solution, 1e-2);

            // Second part of test - check that each cell's hypoxic duration is correctly updated
            if (cell_iter->GetCellData()->GetItem("oxygen") >= p_oxygen_model->GetHypoxicConcentration())
            {
                TS_ASSERT_DELTA(p_oxygen_model->GetCurrentHypoxicDuration(), 0.0, 1e-5);
            }
            else
            {
                TS_ASSERT_DELTA(p_oxygen_model->GetCurrentHypoxicDuration(), 2/120.0, 1e-5);
            }
        }
    }

    void TestWithOxygen()
    {
        EXIT_IF_PARALLEL;

        // Set up mesh
        HoneycombMeshGenerator generator(5, 5, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Set up cells
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            SimpleOxygenBasedCellCycleModel* p_model = new SimpleOxygenBasedCellCycleModel();
            p_model->SetDimension(2);

            // Use non-default G1 durations
            p_model->SetStemCellG1Duration(8.0);
            p_model->SetTransitCellG1Duration(8.0);

            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetCellProliferativeType(p_stem_type);

            double birth_time = -1.0 - ((double) i/p_mesh->GetNumNodes())*18.0;
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        // Set up cell population
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.AddPopulationWriter<CellPopulationAreaWriter>(); // record the spheroid radius and apoptotic radius

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(cell_population);
        TS_ASSERT_EQUALS(simulator.GetIdentifier(), "OffLatticeSimulation-2-2");

        // Create PDE and boundary condition objects
        MAKE_PTR_ARGS(UniformSourceEllipticPde<2>, p_pde, (-0.1));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(EllipticGrowingDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false));
        p_pde_modifier->SetDependentVariableName("oxygen");

        // Output the PDE solution at each time step
        p_pde_modifier->SetOutputSolutionAtPdeNodes(true);

        simulator.AddSimulationModifier(p_pde_modifier);

        simulator.SetOutputDirectory("OffLatticeSimulationWithOxygen");
        simulator.SetEndTime(0.5);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetCutOffLength(1.5);
        simulator.AddForce(p_linear_force);

        // Set up cell killer and pass into simulation
        MAKE_PTR_ARGS(ApoptoticCellKiller<2>, p_killer, (&cell_population));
        simulator.AddCellKiller(p_killer);

        // Run cell-based simulation
        simulator.Solve();
    }

    /*
     * This test compares the visualizer output from the previous test
     * with a known file.
     *
     * Note: if the previous test is changed we need to update the file
     * this test refers to.
     */
    void TestVisualizerOutput()
    {
        EXIT_IF_PARALLEL;

        // Work out where one of the previous tests wrote its files
        OutputFileHandler handler("OffLatticeSimulationWithOxygen", false);
        std::string results_dir = handler.GetOutputDirectoryFullPath() + "results_from_time_0";

        NumericFileComparison comp_nut(results_dir + "/results.vizpdesolution", "cell_based/test/data/OffLatticeSimulationWithOxygen/results.vizpdesolution");
        TS_ASSERT(comp_nut.CompareFiles());

        NumericFileComparison comp_ele(results_dir + "/results.vizelements", "cell_based/test/data/OffLatticeSimulationWithOxygen/results.vizelements");
        TS_ASSERT(comp_ele.CompareFiles());

        NumericFileComparison comp_nodes(results_dir + "/results.viznodes", "cell_based/test/data/OffLatticeSimulationWithOxygen/results.viznodes");
        TS_ASSERT(comp_nodes.CompareFiles(1e-15));

        NumericFileComparison comp_celltypes(results_dir + "/results.vizcelltypes", "cell_based/test/data/OffLatticeSimulationWithOxygen/results.vizcelltypes");
        TS_ASSERT(comp_celltypes.CompareFiles(1e-15));

        FileComparison(results_dir + "/results.vizsetup", "cell_based/test/data/OffLatticeSimulationWithOxygen/results.vizsetup").CompareFiles();
    }

    void TestWithPointwiseSource()
    {
        EXIT_IF_PARALLEL;

        // Set up mesh
        HoneycombMeshGenerator generator(5, 5, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Set up cells
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(ApoptoticCellProperty, p_apoptotic_state);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);

        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            SimpleOxygenBasedCellCycleModel* p_model = new SimpleOxygenBasedCellCycleModel();
            p_model->SetDimension(2);
            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetCellProliferativeType(p_stem_type);

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
        cell_population.AddPopulationWriter<CellPopulationAreaWriter>(); // record the spheroid radius and apoptotic radius

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("OffLatticeSimulationWithPdes");
        simulator.SetEndTime(0.5);

        // Create PDE and boundary condition objects
        MAKE_PTR_ARGS(CellwiseSourceEllipticPde<2>, p_pde, (cell_population, -0.1));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(EllipticGrowingDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false));
        p_pde_modifier->SetDependentVariableName("oxygen");

        simulator.AddSimulationModifier(p_pde_modifier);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetCutOffLength(1.5);
        simulator.AddForce(p_linear_force);

        // Set up cell killer and pass into simulation
        MAKE_PTR_ARGS(ApoptoticCellKiller<2>, p_killer, (&cell_population));
        simulator.AddCellKiller(p_killer);

        // Run cell-based simulation
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());

        // A few hardcoded tests to check nothing has changed
        std::vector<double> node_5_location = simulator.GetNodeLocation(5);
        TS_ASSERT_DELTA(node_5_location[0], 0.6605, 1e-4);
        TS_ASSERT_DELTA(node_5_location[1], 1.1422, 1e-4);
        TS_ASSERT_DELTA( (simulator.rGetCellPopulation().GetCellUsingLocationIndex(5))->GetCellData()->GetItem("oxygen"), 0.9704, 1e-4);
    }

    void TestWithPointwiseTwoSource()
    {
        EXIT_IF_PARALLEL;

        // Set up mesh
        HoneycombMeshGenerator generator(5, 5, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Set up cells
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(ApoptoticCellProperty, p_apoptotic_state);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);

        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            SimpleOxygenBasedCellCycleModel* p_model = new SimpleOxygenBasedCellCycleModel();
            p_model->SetDimension(2);

            // Use non-default G1 durations
            p_model->SetStemCellG1Duration(8.0);
            p_model->SetTransitCellG1Duration(8.0);

            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetCellProliferativeType(p_stem_type);
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
        cell_population.AddPopulationWriter<CellPopulationAreaWriter>(); // record the spheroid radius and apoptotic radius

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("OffLatticeSimulationWithPointwiseSource");
        simulator.SetEndTime(0.5);

        // Create PDE and boundary condition objects
        MAKE_PTR_ARGS(CellwiseSourceEllipticPde<2>, p_pde, (cell_population, -0.1));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(EllipticGrowingDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false));
        p_pde_modifier->SetDependentVariableName("oxygen");

        simulator.AddSimulationModifier(p_pde_modifier);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetCutOffLength(1.5);
        simulator.AddForce(p_linear_force);

        // Set up cell killer and pass into simulation
        MAKE_PTR_ARGS(ApoptoticCellKiller<2>, p_killer, (&cell_population));
        simulator.AddCellKiller(p_killer);

        // Run cell-based simulation
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());

        // A few hardcoded tests to check nothing has changed
        std::vector<double> node_5_location = simulator.GetNodeLocation(5);
        TS_ASSERT_DELTA(node_5_location[0], 0.6605, 1e-4);
        TS_ASSERT_DELTA(node_5_location[1], 1.1422, 1e-4);
        CellPtr p_cell_at_5 = simulator.rGetCellPopulation().GetCellUsingLocationIndex(5);
        TS_ASSERT_DELTA(p_cell_at_5->GetCellData()->GetItem("oxygen"), 0.9704, 1e-4);
    }

    void TestSpheroidStatistics()
    {
        EXIT_IF_PARALLEL;

        // Set up mesh
        HoneycombMeshGenerator generator(5, 5, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Set up cells
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(ApoptoticCellProperty, p_apoptotic_state);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);

        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            SimpleOxygenBasedCellCycleModel* p_model = new SimpleOxygenBasedCellCycleModel();
            p_model->SetDimension(2);

            // Use non-default G1 durations
            p_model->SetStemCellG1Duration(8.0);
            p_model->SetTransitCellG1Duration(8.0);

            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetCellProliferativeType(p_stem_type);
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
        cell_population.AddPopulationWriter<CellPopulationAreaWriter>(); // record the spheroid radius and apoptotic radius

        typedef RadialCellDataDistributionWriter<2,2> RadWriter;
        MAKE_PTR(RadWriter, p_radial_writer);
        p_radial_writer->SetVariableName("oxygen");
        p_radial_writer->SetNumRadialBins(5);
        cell_population.AddPopulationWriter(p_radial_writer);

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestSpheroidStatistics");
        simulator.SetEndTime(1.0/120.0);

        // Create PDE and boundary condition objects
        MAKE_PTR_ARGS(UniformSourceEllipticPde<2>, p_pde, (-0.1));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(EllipticGrowingDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false));
        p_pde_modifier->SetDependentVariableName("oxygen");

        simulator.AddSimulationModifier(p_pde_modifier);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetCutOffLength(1.5);
        simulator.AddForce(p_linear_force);

        // Add an oxygen-dependent cell killer to the cell-based simulation
        MAKE_PTR_ARGS(ApoptoticCellKiller<2>, p_killer, (&cell_population));
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
         * The area of a single hexagonal cell is sqrt(3.0)/2, so the
         * correct spheroid radius is given by sqrt((16*sqrt(3.0)/2)/pi).
         *
         * Since there are 3 apoptotic cells, the correct apoptotic radius
         * is given by sqrt((3*sqrt(3.0)/2)/pi).
         */
        OutputFileHandler handler("TestSpheroidStatistics", false);

        std::string areas_results_file = handler.GetOutputDirectoryFullPath() + "results_from_time_0/cellpopulationareas.dat";
        FileComparison(areas_results_file, "cell_based/test/data/TestSpheroidStatistics/cellpopulationareas.dat").CompareFiles();

        std::string results_dir = handler.GetOutputDirectoryFullPath() + "results_from_time_0/";
        NumericFileComparison comparison(results_dir + "radial_dist.dat", "cell_based/test/data/TestSpheroidStatistics/radial_dist.dat");
        TS_ASSERT(comparison.CompareFiles(5e-3));
    }

    void TestCoarseSourceMesh()
    {
        EXIT_IF_PARALLEL;

        // Set up mesh
        HoneycombMeshGenerator generator(5, 5, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Set up cells
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            SimpleOxygenBasedCellCycleModel* p_model = new SimpleOxygenBasedCellCycleModel();
            p_model->SetDimension(2);

            // Use non-default G1 durations
            p_model->SetStemCellG1Duration(8.0);
            p_model->SetTransitCellG1Duration(8.0);

            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetCellProliferativeType(p_stem_type);
            double birth_time = -RandomNumberGenerator::Instance()->ranf()*18.0;
            p_cell->SetBirthTime(birth_time);

            cells.push_back(p_cell);
        }

        // Set up cell population
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestCoarseSourceMesh");
        simulator.SetEndTime(0.05);

        // Create PDE and boundary condition objects
        MAKE_PTR_ARGS(AveragedSourceEllipticPde<2>, p_pde, (cell_population, -0.1));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));

        MAKE_PTR_ARGS(AveragedSourceEllipticPde<2>, p_pde2, (cell_population, -0.5));

        // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
        c_vector<double,2> centroid = cell_population.GetCentroidOfCellPopulation();
        ChastePoint<2> lower(centroid(0)-25.0, centroid(1)-25.0);
        ChastePoint<2> upper(centroid(0)+25.0, centroid(1)+25.0);
        MAKE_PTR_ARGS(ChasteCuboid<2>, p_cuboid, (lower, upper));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false, p_cuboid, 10.0));
        p_pde_modifier->SetDependentVariableName("oxygen");
        p_pde_modifier->SetBcsOnBoxBoundary(false);

        simulator.AddSimulationModifier(p_pde_modifier);

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier<2>, p_pde_modifier2, (p_pde2, p_bc, false, p_cuboid, 10.0));
        p_pde_modifier2->SetDependentVariableName("dunno");
        p_pde_modifier2->SetBcsOnBoxBoundary(false);

        simulator.AddSimulationModifier(p_pde_modifier2);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetCutOffLength(1.5);
        simulator.AddForce(p_linear_force);

        // Set up cell killer and pass into simulation
        MAKE_PTR_ARGS(ApoptoticCellKiller<2>, p_killer, (&cell_population));
        simulator.AddCellKiller(p_killer);

        // Find centre of cell population
        c_vector<double,2> centre_of_cell_population = cell_population.GetCentroidOfCellPopulation();

        // Find centre of coarse PDE mesh
        c_vector<double,2> centre_of_coarse_pde_mesh = zero_vector<double>(2);

        TetrahedralMesh<2,2>* p_coarse_mesh = p_pde_modifier->GetFeMesh();

        for (unsigned i=0; i<p_coarse_mesh->GetNumNodes(); i++)
        {
            centre_of_coarse_pde_mesh += p_coarse_mesh->GetNode(i)->rGetLocation();
        }
        centre_of_coarse_pde_mesh /= p_coarse_mesh->GetNumNodes();

        // Test that the two centres match
        c_vector<double,2> centre_diff = centre_of_cell_population - centre_of_coarse_pde_mesh;
        TS_ASSERT_DELTA(norm_2(centre_diff), 0.0, 1e-4);

        // Run cell-based simulation
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());

        TS_ASSERT(p_coarse_mesh != NULL);

        ReplicatableVector pde_solution0(p_pde_modifier->GetSolution());
        ReplicatableVector pde_solution1(p_pde_modifier2->GetSolution());

        TS_ASSERT_EQUALS(pde_solution0.GetSize(), pde_solution1.GetSize());

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
            unsigned elem_index = p_pde_modifier->GetFeMesh()->GetContainingElementIndex(cell_population.GetLocationOfCellCentre(*cell_iter));
            Element<2,2>* p_element = p_coarse_mesh->GetElement(elem_index);

            double max0 = std::max(pde_solution0[p_element->GetNodeGlobalIndex(0)], pde_solution0[p_element->GetNodeGlobalIndex(1)]);
            max0 = std::max(max0, pde_solution0[p_element->GetNodeGlobalIndex(2)]);

            double max1 = std::max(pde_solution1[p_element->GetNodeGlobalIndex(0)], pde_solution1[p_element->GetNodeGlobalIndex(1)]);
            max1 = std::max(max1, pde_solution1[p_element->GetNodeGlobalIndex(2)]);

            double min0 = std::min(pde_solution0[p_element->GetNodeGlobalIndex(0)], pde_solution0[p_element->GetNodeGlobalIndex(1)]);
            min0 = std::min(min0, pde_solution0[p_element->GetNodeGlobalIndex(2)]);

            double min1 = std::min(pde_solution1[p_element->GetNodeGlobalIndex(0)], pde_solution1[p_element->GetNodeGlobalIndex(1)]);
            min1 = std::min(min1, pde_solution1[p_element->GetNodeGlobalIndex(2)]);

            double value0_at_cell = cell_iter->GetCellData()->GetItem("oxygen");
            double value1_at_cell = cell_iter->GetCellData()->GetItem("dunno");

            TS_ASSERT_LESS_THAN_EQUALS(value1_at_cell, value0_at_cell);
            TS_ASSERT_LESS_THAN_EQUALS(min0, value0_at_cell + DBL_EPSILON);
            TS_ASSERT_LESS_THAN_EQUALS(value0_at_cell, max0 + DBL_EPSILON);
            TS_ASSERT_LESS_THAN_EQUALS(min1, value1_at_cell + DBL_EPSILON);
            TS_ASSERT_LESS_THAN_EQUALS(value1_at_cell, max1 + DBL_EPSILON);
        }
    }

    void TestCoarseSourceMeshWithGhostNodes()
    {
        EXIT_IF_PARALLEL;

        // Set up mesh
        HoneycombMeshGenerator generator(5, 5, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumNodes()-1);

        std::vector<unsigned> cell_location_indices;
        for (unsigned i=0; i<cells.size(); i++)
        {
            cell_location_indices.push_back(i);
        }

        // Passes as the cell population constructor automatically works out which
        // cells are ghost nodes using the mesh and cell_location_indices
        MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, cell_location_indices);

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestCoarseSourceMeshWithGhostNodes");
        simulator.SetEndTime(0.05);

        // Create PDE and boundary condition objects
        MAKE_PTR_ARGS(AveragedSourceEllipticPde<2>, p_pde, (cell_population, 0.0));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));

        // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
        c_vector<double,2> centroid = cell_population.GetCentroidOfCellPopulation();
        ChastePoint<2> lower(centroid(0)-25.0, centroid(1)-25.0);
        ChastePoint<2> upper(centroid(0)+25.0, centroid(1)+25.0);
        MAKE_PTR_ARGS(ChasteCuboid<2>, p_cuboid, (lower, upper));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false, p_cuboid));
        p_pde_modifier->SetDependentVariableName("nutrient");

        simulator.AddSimulationModifier(p_pde_modifier);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetCutOffLength(1.5);
        simulator.AddForce(p_linear_force);

        // Run cell-based simulation
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());

        // Test solution is constant
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            double analytic_solution = 1.0;
            // Test that PDE solver is working correctly
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("nutrient"), analytic_solution, 1e-2);
        }
    }

    void TestArchivingWithSimplePde()
    {
        EXIT_IF_PARALLEL;

        // Set up mesh
        HoneycombMeshGenerator generator(5, 5, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Set up cells
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            SimpleOxygenBasedCellCycleModel* p_model = new SimpleOxygenBasedCellCycleModel();
            p_model->SetDimension(2);

            // Use non-default G1 durations
            p_model->SetStemCellG1Duration(8.0);
            p_model->SetTransitCellG1Duration(8.0);

            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetCellProliferativeType(p_stem_type);
            double birth_time = -1.0 - ((double) i/p_mesh->GetNumNodes())*18.0;
            p_cell->SetBirthTime(birth_time);

            cells.push_back(p_cell);
        }

        // Set up cell population
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("OffLatticeSimulationWithPdesSaveAndLoad");
        simulator.SetEndTime(0.2);

        // Create PDE and boundary condition objects
        MAKE_PTR_ARGS(UniformSourceEllipticPde<2>, p_pde, (-0.1));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(EllipticGrowingDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false));
        p_pde_modifier->SetDependentVariableName("oxygen");

        simulator.AddSimulationModifier(p_pde_modifier);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetCutOffLength(1.5);
        simulator.AddForce(p_linear_force);

        // Set up cell killer and pass into simulation
        MAKE_PTR_ARGS(ApoptoticCellKiller<2>, p_killer, (&cell_population));
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
        TS_ASSERT_DELTA(node_5_location[0], 0.4987, 1e-4);
        TS_ASSERT_DELTA(node_5_location[1], 0.8694, 1e-4);

        std::vector<double> node_15_location = p_simulator->GetNodeLocation(15);
        TS_ASSERT_DELTA(node_15_location[0], 0.5106, 1e-4);
        TS_ASSERT_DELTA(node_15_location[1], 2.6052, 1e-4);

        // Test the CellData result
        TS_ASSERT_DELTA((p_simulator->rGetCellPopulation().GetCellUsingLocationIndex(5))->GetCellData()->GetItem("oxygen"), 0.9605, 1e-4);
        TS_ASSERT_DELTA((p_simulator->rGetCellPopulation().GetCellUsingLocationIndex(15))->GetCellData()->GetItem("oxygen"), 0.9582, 1e-4);

        // Tidy up
        delete p_simulator;
    }

    /**
     * This test demonstrates how to archive a OffLatticeSimulation
     * in the case where the PDE has the cell population as a member variable.
     */
    void TestArchivingWithCellwisePde()
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
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            SimpleOxygenBasedCellCycleModel* p_model = new SimpleOxygenBasedCellCycleModel();
            p_model->SetDimension(2);

            // Use non-default G1 durations
            p_model->SetStemCellG1Duration(8.0);
            p_model->SetTransitCellG1Duration(8.0);

            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetCellProliferativeType(p_stem_type);
            double birth_time = -1.0 - ((double) i/p_mesh->GetNumNodes())*18.0;
            p_cell->SetBirthTime(birth_time);

            cells.push_back(p_cell);
        }

        // Set up cell population
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory(output_directory);
        simulator.SetEndTime(end_time);

        // Create PDE and boundary condition objects
        MAKE_PTR_ARGS(CellwiseSourceEllipticPde<2>, p_pde, (cell_population, -0.03));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(EllipticGrowingDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false));
        p_pde_modifier->SetDependentVariableName("oxygen");

        simulator.AddSimulationModifier(p_pde_modifier);

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
    }

    void Test3DOffLatticeSimulationWithPdes()
    {
        EXIT_IF_PARALLEL;

        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_136_elements");
        MutableMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        TrianglesMeshWriter<3,3> mesh_writer("TestSolveMethodSpheroidSimulation3DMesh", "StartMesh");
        mesh_writer.WriteFilesUsingMesh(mesh);

        // Set up cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 3> generator;
        generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Set some model parameters for the cell-cycle model
        for (unsigned index=0; index < cells.size(); index++)
        {
            static_cast<FixedG1GenerationalCellCycleModel*>(cells[index]->GetCellCycleModel())->SetTransitCellG1Duration(8.0);
            static_cast<FixedG1GenerationalCellCycleModel*>(cells[index]->GetCellCycleModel())->SetStemCellG1Duration(8.0);
        }

        // Set up cell population
        MeshBasedCellPopulation<3> cell_population(mesh, cells);

        // Set up cell-based simulation
        OffLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory("OffLatticeSimulationWithOxygen3d");
        simulator.SetSamplingTimestepMultiple(12);
        simulator.SetEndTime(0.1);

        // Create PDE and boundary condition objects
        MAKE_PTR_ARGS(UniformSourceEllipticPde<3>, p_pde, (-0.03));
        MAKE_PTR_ARGS(ConstBoundaryCondition<3>, p_bc, (1.0));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(EllipticGrowingDomainPdeModifier<3>, p_pde_modifier, (p_pde, p_bc, false));
        p_pde_modifier->SetDependentVariableName("oxygen");

        simulator.AddSimulationModifier(p_pde_modifier);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<3>, p_linear_force);
        p_linear_force->SetCutOffLength(1.5);
        simulator.AddForce(p_linear_force);

        // Set up cell killer and pass into simulation
        MAKE_PTR_ARGS(ApoptoticCellKiller<3>, p_killer, (&cell_population));
        simulator.AddCellKiller(p_killer);

        // Run cell-based simulation
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());
    }

    /**
     * Test that the simulation runs correctly when a PDE with time-dependent
     * boundary condition is used. Here we test the PDE solver using the problem
     * del squared C = 1 on the unit disc, with boundary condition C=t on r=1,
     * which has analytic solution C = t - 0.25*(1-r^2).
     */
    void TestWithBoundaryConditionVaryingInTime()
    {
        EXIT_IF_PARALLEL;

        // Set up mesh
        MutableMesh<2,2> mesh;
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_522_elements");
        mesh.ConstructFromMeshReader(mesh_reader);

        // Set up cells
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            SimpleOxygenBasedCellCycleModel* p_model = new SimpleOxygenBasedCellCycleModel();
            p_model->SetDimension(2);
            p_model->SetHypoxicConcentration(0.9);
            p_model->SetQuiescentConcentration(0.9);

            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetCellProliferativeType(p_diff_type);
            p_cell->SetBirthTime(-1.0);
            cells.push_back(p_cell);
        }

        // Set up cell population
        MeshBasedCellPopulation<2> cell_population(mesh, cells);

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestUpdateAtEndOfTimeStepMethod");

        // Create PDE and boundary condition objects
        MAKE_PTR(SimplePdeForTesting, p_pde);
        MAKE_PTR_ARGS(FunctionalBoundaryCondition<2>, p_functional_bc, (&bc_func));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(EllipticGrowingDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_functional_bc, false));
        p_pde_modifier->SetDependentVariableName("oxygen");

        simulator.AddSimulationModifier(p_pde_modifier);

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
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("oxygen"), analytic_solution, 0.02);
        }
    }

    void TestOffLatticeSimulationWithPdesParameterOutputMethods()
    {
        EXIT_IF_PARALLEL;

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
        TS_ASSERT_EQUALS(simulator.GetIdentifier(), "OffLatticeSimulation-2-2");

        // Create PDE and boundary condition objects
        MAKE_PTR_ARGS(CellwiseSourceEllipticPde<2>, p_pde, (cell_population, -0.03));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(EllipticGrowingDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false));
        p_pde_modifier->SetDependentVariableName("oxygen");

        simulator.AddSimulationModifier(p_pde_modifier);

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
        FileComparison( results_dir + "cell_based_sim_with_pde_results.parameters", "cell_based/test/data/TestOffLatticeSimulationOutputParameters/cell_based_sim_with_pde_results.parameters").CompareFiles();

        ///\todo check output of simulator.OutputSimulationSetup();
    }

    void TestNodeBasedWithCoarseMesh()
    {
        EXIT_IF_PARALLEL;

        // Create a simple mesh - larger to use coarse graining.
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_522_elements");
        TetrahedralMesh<2,2> temp_mesh;
        temp_mesh.ConstructFromMeshReader(mesh_reader);
        temp_mesh.Scale(20.0,20.0);

        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(temp_mesh, 1.5);

        // Set up cells
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel();
            p_model->SetDimension(2);

            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetCellProliferativeType(p_diff_type);
            double birth_time = -RandomNumberGenerator::Instance()->ranf()*18.0;
            p_cell->SetBirthTime(birth_time);

            cells.push_back(p_cell);
        }

        // Set up cell population
        NodeBasedCellPopulation<2> cell_population(mesh, cells);

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestNodeBasedCellPopulationWithCoarseMeshPDE");
        simulator.SetEndTime(0.01);

        // Create PDE and boundary condition objects
        MAKE_PTR_ARGS(AveragedSourceEllipticPde<2>, p_pde, (cell_population, 0.0));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));

        // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
        c_vector<double,2> centroid = cell_population.GetCentroidOfCellPopulation();
        ChastePoint<2> lower(centroid(0)-25.0, centroid(1)-25.0);
        ChastePoint<2> upper(centroid(0)+25.0, centroid(1)+25.0);
        MAKE_PTR_ARGS(ChasteCuboid<2>, p_cuboid, (lower, upper));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false, p_cuboid));
        p_pde_modifier->SetDependentVariableName("nutrient");

        simulator.AddSimulationModifier(p_pde_modifier);

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
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("nutrient"), analytic_solution, 1e-2);
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
        TetrahedralMesh<2,2>* p_coarse_mesh = p_pde_modifier->GetFeMesh();
        for (unsigned i=0; i<p_coarse_mesh->GetNumNodes(); i++)
        {
            centre_of_coarse_pde_mesh += p_coarse_mesh->GetNode(i)->rGetLocation();
        }
        centre_of_coarse_pde_mesh /= p_coarse_mesh->GetNumNodes();

        // Test that the two centres match
        c_vector<double,2> centre_diff = centre_of_cell_population - centre_of_coarse_pde_mesh;
        TS_ASSERT_DELTA(norm_2(centre_diff), 0.0, 1e-4);
    }

    void TestVertexBasedWithCoarseMesh()
    {
        EXIT_IF_PARALLEL;

        /// Create a simple 2D VertexMesh
        HoneycombVertexMeshGenerator generator(5, 3);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements());

        // Create cell population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestVertexBasedCellPopulationWithCoarseMeshPDE");
        simulator.SetEndTime(0.01);

        // Create PDE and boundary condition objects (zero uptake to check analytic solution)
        MAKE_PTR_ARGS(AveragedSourceEllipticPde<2>, p_pde, (cell_population, 0.0));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));

        // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
        c_vector<double,2> centroid = cell_population.GetCentroidOfCellPopulation();
        ChastePoint<2> lower(centroid(0)-25.0, centroid(1)-25.0);
        ChastePoint<2> upper(centroid(0)+25.0, centroid(1)+25.0);
        MAKE_PTR_ARGS(ChasteCuboid<2>, p_cuboid, (lower, upper));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false, p_cuboid));
        p_pde_modifier->SetDependentVariableName("nutrient");

        simulator.AddSimulationModifier(p_pde_modifier);

        // Create a force law and pass it to the simulation
        MAKE_PTR(NagaiHondaForce<2>, p_nagai_honda_force);
        simulator.AddForce(p_nagai_honda_force);

        // A NagaiHondaForce has to be used together with an AbstractTargetAreaModifier
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        // Find centre of coarse PDE mesh
        c_vector<double,2> centre_of_coarse_pde_mesh = zero_vector<double>(2);
        TetrahedralMesh<2,2>* p_coarse_mesh = p_pde_modifier->GetFeMesh();
        for (unsigned i=0; i<p_coarse_mesh->GetNumNodes(); i++)
        {
            centre_of_coarse_pde_mesh += p_coarse_mesh->GetNode(i)->rGetLocation();
        }
        centre_of_coarse_pde_mesh /= p_coarse_mesh->GetNumNodes();

        // Find centre of cell population
        c_vector<double,2> centre_of_cell_population = cell_population.GetCentroidOfCellPopulation();

        // Test that the two centres match
        c_vector<double,2> centre_diff = centre_of_cell_population - centre_of_coarse_pde_mesh;
        TS_ASSERT_DELTA(norm_2(centre_diff), 0.0, 1e-4);

        // Solve the system
        simulator.Solve();

        // Test solution is constant
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            double analytic_solution = 1.0;
            // Test that PDE solver is working correctly
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("nutrient"), analytic_solution, 1e-2);
        }
    }

    void TestVolumeDependentAveragedPde()
    {
        EXIT_IF_PARALLEL;

        // Create a simple mesh
        TetrahedralMesh<2,2> temp_mesh;
        temp_mesh.ConstructRegularSlabMesh(2.0,2,2);

        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(temp_mesh, 1.5);

        // Set some cell radius
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            mesh.GetNode(i)->SetRadius( 0.9);
        }

        // Set up cells
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel();
            p_model->SetDimension(2);

            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetCellProliferativeType(p_diff_type);
            double birth_time = -RandomNumberGenerator::Instance()->ranf()*18.0;
            p_cell->SetBirthTime(birth_time);

            cells.push_back(p_cell);
        }

        // Set up cell population
        NodeBasedCellPopulation<2> cell_population(mesh, cells);

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestCoarsePdeSolutionOnNodeBased");
        simulator.SetEndTime(0.01);

        // Create PDE and boundary condition objects (uniform uptake at each cell)
        MAKE_PTR_ARGS(VolumeDependentAveragedSourceEllipticPde<2>, p_pde, (cell_population, -0.01));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));

        // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
        c_vector<double,2> centroid = cell_population.GetCentroidOfCellPopulation();
        ChastePoint<2> lower(centroid(0)-25.0, centroid(1)-25.0);
        ChastePoint<2> upper(centroid(0)+25.0, centroid(1)+25.0);
        MAKE_PTR_ARGS(ChasteCuboid<2>, p_cuboid, (lower, upper));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false, p_cuboid));
        p_pde_modifier->SetDependentVariableName("nutrient");

        simulator.AddSimulationModifier(p_pde_modifier);

        // Solve the system
        simulator.Solve();

        TetrahedralMesh<2,2>* p_coarse_mesh = p_pde_modifier->GetFeMesh();

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
            unsigned containing_element_index = p_pde_modifier->mCellPdeElementMap[*cell_iter];
            cells_in_each_coarse_element[containing_element_index] += 1;
        }

        for (unsigned i=0; i<num_elements_in_coarse_mesh; i++)
        {
            c_matrix<double, 2, 2> jacobian;
            double det;
            p_coarse_mesh->GetElement(i)->CalculateJacobian(jacobian, det);
            double element_volume = p_coarse_mesh->GetElement(i)->GetVolume(det);
            double uptake_rate = p_pde->GetUptakeRateForElement(i);
            double expected_uptake = 0.9*0.9*(cells_in_each_coarse_element[i]/element_volume);

            TS_ASSERT_DELTA(uptake_rate, expected_uptake, 1e-4);
        }
    }

    void TestCoarsePdeSolutionOnNodeBased1d()
    {
        EXIT_IF_PARALLEL;

        // Create mesh
        std::vector<Node<1>*> nodes;
        nodes.push_back(new Node<1>(0, true,  0.0));

        NodesOnlyMesh<1> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);

        // Set up differentiated cells
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel();
            p_model->SetDimension(1);

            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetCellProliferativeType(p_diff_type);
            double birth_time = -RandomNumberGenerator::Instance()->ranf()*18.0;
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        // Set up cell population
        NodeBasedCellPopulation<1> cell_population(mesh, cells);

        // Set up cell-based simulation
        OffLatticeSimulation<1> simulator(cell_population);
        simulator.SetEndTime(0.01);

        // Create PDE and boundary condition objects
        MAKE_PTR_ARGS(AveragedSourceEllipticPde<1>, p_pde, (cell_population, -1.0));
        MAKE_PTR_ARGS(ConstBoundaryCondition<1>, p_bc, (0.0));

        // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
        ChastePoint<1> lower(0.0);
        ChastePoint<1> upper(2.0);
        MAKE_PTR_ARGS(ChasteCuboid<1>, p_cuboid, (lower, upper));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier<1>, p_pde_modifier, (p_pde, p_bc, false, p_cuboid));
        p_pde_modifier->SetDependentVariableName("nutrient");

        simulator.AddSimulationModifier(p_pde_modifier);

        simulator.SetOutputDirectory("TestCoarsePdeSolutionOnNodeBased1d");

        // Solve the system
        simulator.Solve();

        // When the mesh goes out of scope, then it's a different set of nodes that get destroyed
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

    void TestCoarsePdeSolutionOnNodeBased2d()
    {
        EXIT_IF_PARALLEL;

        // Create a simple mesh
        TetrahedralMesh<2,2> temp_mesh;
        temp_mesh.ConstructRectangularMesh(10,10);

        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(temp_mesh, 1.5);

        // Set up cells
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel();
            p_model->SetDimension(2);

            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetCellProliferativeType(p_diff_type);
            double birth_time = -RandomNumberGenerator::Instance()->ranf()*18.0;
            p_cell->SetBirthTime(birth_time);

            cells.push_back(p_cell);
        }

        // Set up cell population
        NodeBasedCellPopulation<2> cell_population(mesh, cells);

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestCoarsePdeSolutionOnNodeBased2d");
        simulator.SetEndTime(0.01);

        // Create PDE and boundary condition objects
        MAKE_PTR_ARGS(AveragedSourceEllipticPde<2>, p_pde, (cell_population, -0.01));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));

        // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
        c_vector<double,2> centroid = cell_population.GetCentroidOfCellPopulation();
        ChastePoint<2> lower(centroid(0)-25.0, centroid(1)-25.0);
        ChastePoint<2> upper(centroid(0)+25.0, centroid(1)+25.0);
        MAKE_PTR_ARGS(ChasteCuboid<2>, p_cuboid, (lower, upper));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false, p_cuboid, 10.0));
        p_pde_modifier->SetDependentVariableName("nutrient");
        p_pde_modifier->SetBcsOnBoxBoundary(false);

        simulator.AddSimulationModifier(p_pde_modifier);

        // Solve the system
        simulator.Solve();

        // Test solution is unchanged at first two cells
        TS_ASSERT_DELTA(  (cell_population.Begin())->GetCellData()->GetItem("nutrient"), 0.9543, 1e-4);
        TS_ASSERT_DELTA((++cell_population.Begin())->GetCellData()->GetItem("nutrient"), 0.9589, 1e-4);
    }

    void TestCoarsePdeSolutionOnNodeBased3d()
    {
        EXIT_IF_PARALLEL;

        // Create a simple mesh
        TetrahedralMesh<3,3> temp_mesh;
        temp_mesh.ConstructCuboid(5,5,5);

        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(temp_mesh, 1.5);

        // Set up cells
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel();
            p_model->SetDimension(3);

            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetCellProliferativeType(p_diff_type);
            double birth_time = -RandomNumberGenerator::Instance()->ranf()*18.0;
            p_cell->SetBirthTime(birth_time);

            cells.push_back(p_cell);
        }

        // Set up cell population
        NodeBasedCellPopulation<3> cell_population(mesh, cells);

        // Set up cell-based simulation
        OffLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory("TestCoarsePdeSolutionOnNodeBased3d");
        simulator.SetEndTime(0.01);

        // Create PDE and boundary condition objects
        MAKE_PTR_ARGS(AveragedSourceEllipticPde<3>, p_pde, (cell_population, -0.01));
        MAKE_PTR_ARGS(ConstBoundaryCondition<3>, p_bc, (1.0));

        // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
        c_vector<double,3> centroid = cell_population.GetCentroidOfCellPopulation();
        ChastePoint<3> lower(centroid(0)-25.0, centroid(1)-25.0, centroid(2)-25.0);
        ChastePoint<3> upper(centroid(0)+25.0, centroid(1)+25.0, centroid(2)+25.0);
        MAKE_PTR_ARGS(ChasteCuboid<3>, p_cuboid, (lower, upper));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier<3>, p_pde_modifier, (p_pde, p_bc, false, p_cuboid, 10.0));
        p_pde_modifier->SetDependentVariableName("nutrient");
        p_pde_modifier->SetBcsOnBoxBoundary(false);

        simulator.AddSimulationModifier(p_pde_modifier);

        // Solve the system
        simulator.Solve();

        // Test solution is unchanged at first two cells
        TS_ASSERT_DELTA(  (cell_population.Begin())->GetCellData()->GetItem("nutrient"), 1.0, 1e-4);
        TS_ASSERT_DELTA((++cell_population.Begin())->GetCellData()->GetItem("nutrient"), 1.0, 1e-4);
    }
};

#endif /*TESTOFFLATTICESIMULATIONWITHPDES_HPP_*/
