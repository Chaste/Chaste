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
/*
 *
 *  Chaste tutorial - this page gets automatically changed to a wiki page
 *  DO NOT remove the comments below, and if the code has to be changed in
 *  order to run, please check the comments are still accurate
 *
 *
 */

#ifndef TESTRUNNINGPOTTSBASEDSIMULATIONSTUTORIAL_HPP_
#define TESTRUNNINGPOTTSBASEDSIMULATIONSTUTORIAL_HPP_

/*
 * = Examples showing how to create, run and visualize Potts-based simulations =
 *
 * == Introduction ==
 *
 * In this tutorial we show how Chaste can be used to create, run and visualize Potts-based simulations.
 * Full details of the mathematical model can be found in Graner, F. and Glazier, J. A. (1992). Simulation
 * of biological cell sorting using a two-dimensional extended Potts model. Phys. Rev. Lett., 69(13):2015–2016.
 *
 * == The test ==
 *
 * As in previous cell-based Chaste tutorials, we begin by including the necessary header files.
 */
#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "PetscSetupAndFinalize.hpp"

/* The remaining header files define classes that will be used in the cell population
 * simulation test. We have encountered some of these header files in previous cell-based
 * Chaste tutorials. */
#include "CellsGenerator.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "SmartPointers.hpp"
#include "UniformCellCycleModel.hpp"
/* The next header file defines a helper class for generating a suitable mesh. */
#include "PottsMeshGenerator.hpp"
/* The next header file defines the class that simulates the evolution of an on lattice {{{CellPopulation}}}. */
#include "OnLatticeSimulation.hpp"
/* The next header file defines a{{{CellPopulation}}} class for implementing a cellular Potts model.*/
#include "PottsBasedCellPopulation.hpp"
/* The next header files define some update rules for describing the Hamiltonian used to define the Potts simulations. */
#include "VolumeConstraintPottsUpdateRule.hpp"
#include "AdhesionPottsUpdateRule.hpp"
#include "DifferentialAdhesionPottsUpdateRule.hpp"
#include "TransitCellProliferativeType.hpp"
/* Finally these headers allow us to output cell labels. */
#include "CellLabel.hpp"
#include "CellLabelWriter.hpp"

/*
 * Next, we define the test class, which inherits from {{{AbstractCellBasedTestSuite}}}
 * and defines some test methods.
 */
class TestRunningPottsBasedSimulationsTutorial : public AbstractCellBasedWithTimingsTestSuite
{
public:
    /* EMPTYLINE
     *
     * == Test 1 - a basic Potts-based simulation ==
     *
     * EMPTYLINE
     *
     * In the first test, we run a simple Potts-based simulation, in which we create a monolayer
     * of cells, using a Potts mesh. Each cell is assigned a stochastic cell-cycle model.
     */
    void TestMonolayer()
    {
        /** The next line is needed because we cannot currently run Potts simulations in parallel. */
        EXIT_IF_PARALLEL;

        /* First, we generate a Potts mesh. To create a {{{PottsMesh}}}, we can use
         * the {{{PottsMeshGenerator}}}. This generates a regular square-shaped mesh,
         * in which all elements are the same size.
         * Here the first three arguments specify the domain width; the number of elements across; and the width of
         * elements. The second set of three arguments specify the domain height; the number of elements up; and
         * the height of individual elements.
         * We have chosen a 2 by 2 block of elements, each consisting of 4 by 4  ( = 16) lattice sites.
         */
        PottsMeshGenerator<2> generator(50, 2, 4, 50, 2, 4);  // Parameters are: lattice sites across; num elements across; element width; lattice sites up; num elements up; and element height
        PottsMesh<2>* p_mesh = generator.GetMesh();

        /* Having created a mesh, we now create a {{{std::vector}}} of {{{CellPtr}}}s.
         * To do this, we the `CellsGenerator` helper class, which is templated over the type
         * of cell model required (here {{{UniformCellCycleModel}}})
         * and the dimension. We create an empty vector of cells and pass this into the
         * method along with the mesh. The second argument represents the size of that the vector
         * {{{cells}}} should become - one cell for each element. Third argument makes all cells
         * proliferate.*/
        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_transit_type);

        /* Now we have a mesh and a set of cells to go with it, we can create a {{{CellPopulation}}}.
         * In general, this class associates a collection of cells with a mesh.
         * For this test, because we have a {{{PottsMesh}}}, we use a particular type of
         * cell population called a {{{PottsBasedCellPopulation}}}.
         */
        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);

        /*
         * We can set the "Temperature" to be used in the Potts Simulation using the optional command below.
         * The default value is 0.1.
         */
        cell_population.SetTemperature(0.1);
        /*
         * By default the Potts simulation will make 1 sweep over the whole domain per timestep.  To use a different
         * number of sweeps per timestep use the command.
         */
        cell_population.SetNumSweepsPerTimestep(1);

        /* We then pass in the cell population into an {{{OnLatticeSimulation}}},
         * and set the output directory and end time.*/
        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("PottsBasedMonolayer");
        simulator.SetEndTime(50.0);
        /*
         * The default timestep is 0.1, but can be changed using the below command. The timestep is used in conjunction with the "Temperature" and
         * number of sweeps per timestep to specify the relationship between cell movement and proliferation. We also set the simulation to only output
         * every 10 steps i.e. once per hour.
         */
        simulator.SetDt(0.1);
        simulator.SetSamplingTimestepMultiple(10);

        /* We must now create one or more update rules, which determine the Hamiltonian
         * in the Potts simulation. For this test, we use two update rules based upon
         * a volume constraint ({{{VolumeConstraintPottsUpdateRule}}}) and adhesion between cells ({{{AdhesionPottsUpdateRule}}}) and pass them to the {{{OnLatticeSimulation}}}.
         * For a list of possible update rules see subclasses of {{{AbstractPottsUpdateRule}}}.
         * These can be found in the inheritance diagram, here, [class:AbstractPottsUpdateRule AbstractPottsUpdateRule].
         *
         * Similarly to specifying forces for off lattice simulations we use the {{{MAKE_PTR}}} macro
         * to make a boost shared pointer to our required update rule before specifying parameters and passing to the simulation as follows
         */
        MAKE_PTR(VolumeConstraintPottsUpdateRule<2>, p_volume_constraint_update_rule);
        /*
         * Set an appropriate target volume in number of lattice sites. Here we use the default value of 16 lattice sites.
         */
        p_volume_constraint_update_rule->SetMatureCellTargetVolume(16);
        /*
         * You can also vary the deformation energy parameter. The larger the parameter
         * the more cells will try to maintain target volume. Here we use the default value of 0.2.
         */
        p_volume_constraint_update_rule->SetDeformationEnergyParameter(0.2);
        /*
         * Finally we add the update rule to the simulator.
         */
        simulator.AddUpdateRule(p_volume_constraint_update_rule);
        /*
         * We repeat the process for any other update rules.
         */
        MAKE_PTR(AdhesionPottsUpdateRule<2>, p_adhesion_update_rule);
        simulator.AddUpdateRule(p_adhesion_update_rule);

        /* To run the simulation, we call {{{Solve()}}}. */
        simulator.Solve();

        /* The next two lines are for test purposes only and are not part of this tutorial. If different simulation input parameters are being explored
         * the lines should be removed.*/
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 64u);
        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), 50.0, 1e-10);
    }

    /*
     * EMPTYLINE
     *
     * To visualize the results, open a new terminal, {{{cd}}} to the Chaste directory,
     * then {{{cd}}} to {{{anim}}}. Then do: {{{java Visualize2dVertexCells /tmp/$USER/testoutput/PottsBasedMonolayer/results_from_time_0}}}.
     * We may have to do: {{{javac Visualize2dVertexCells.java}}} beforehand to create the
     * java executable.
     *
     * We could also visualize the results using paraview.
     *
     * See UserTutorials/VisualizingWithParaview for more information.
     *
     * Load the file {{{/tmp/$USER/testoutput/PottsBasedMonolayer/results_from_time_0/results.pvd}}}, and click apply.
     *
     * Add box "Glyphs" to represent lattice sites. You will need to adjust the size so they don't overlap.
     *
     * Note that, for larger simulations, you may need to unclick "Mask Points" (or similar) so as not to limit the number of glyphs
     * displayed by Paraview.
     *
     * Select the "Display" tab and select "color by" cell index to see individual cells.
     *
     * Add a "Threshold" filter, filter by cell type and make the lower threshold 0 or greater (unoccupied lattice sites are labelled with -1). This will allow you to view only the cells.
     *
     * Load the files {{{/tmp/$USER/testoutput/PottsBasedMonolayer/results_from_time_0/outlines_..vtu}}}, and click apply.
     *
     * In order to see the cell outlines you will need to select "Surface With Edges" in the drop down menu.
     *
     * Click play to see the evolution of the simulation.
     *
     * You should see that the cells sort into ones of the same type.
     *
     * EMPTYLINE
     *
     * == Test 2 - Cell Sorting ==
     *
     * EMPTYLINE
     *
     * The next test generates a collection of cells, there are two types of cells, labelled ones
     * and non labelled ones, there is differential adhesion between the cell types. For the
     * parameters specified, the cells sort into separate types.
     *
     * Parameters are taken from Graner, F. and Glazier, J. A. (1992). Simulation of biological
     * cell sorting using a two-dimensional extended Potts model. Phys. Rev. Lett., 69(13):2015–2016.
     *
     */
    void TestPottsMonolayerCellSorting()
    {
        /** The next line is needed because we cannot currently run Potts simulations in parallel. */
        EXIT_IF_PARALLEL;

        /* First, we generate a Potts mesh. To create a {{{PottsMesh}}}, we can use
         * the {{{PottsMeshGenerator}}}. This generates a regular square-shaped mesh,
         * in which all elements are the same size.
         * We have chosen an 8 by 8 block of elements each consisting of 4 by 4  ( = 16) lattice sites.
         */
        PottsMeshGenerator<2> generator(50, 8, 4, 50, 8, 4);  // Parameters are: lattice sites across; num elements across; element width; lattice sites up; num elements up; and element height
        PottsMesh<2>* p_mesh = generator.GetMesh();

        /* Having created a mesh, we now create a {{{std::vector}}} of {{{CellPtr}}}s.
         * To do this, we the `CellsGenerator` helper class, as before but this time
         * the third argument is set to make all cells non-proliferative. */
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_diff_type);

        /* Before we make a {{{CellPopulation}}} we make a boost shared pointer to a cell label and then assign this
         * label to some randomly chosen cells. */
        MAKE_PTR(CellLabel, p_label);
        for (unsigned i = 0; i<cells.size(); i++)
        {
            if (RandomNumberGenerator::Instance()->ranf() < 0.5)
            {
                cells[i]->AddCellProperty(p_label);
            }
        }

        /* Now we have a mesh and a set of cells to go with it, we can create a {{{CellPopulation}}}.
         */
        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);

        /* In order to visualize labelled cells we need to use the following command.*/
        cell_population.AddCellWriter<CellLabelWriter>();

        /* We then pass in the cell population into an {{{OnLatticeSimulation}}},
         * and set the output directory and end time. */
        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("PottsMonolayerCellSorting");
        simulator.SetEndTime(20.0);
        simulator.SetSamplingTimestepMultiple(10);

        /* We must now create one or more update rules, which determine the Hamiltonian
         * in the Potts simulation. For this test, we use two update rules based upon
         * a volume constraint ({{{VolumeConstraintPottsUpdateRule}}}) and differential adhesion between cells ({{{DifferentialAdhesionPottsUpdateRule}}}), set appropriate parameters, and pass them to the {{{OnLatticeSimulation}}}.
         */
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
        /*
         * These parameters cause the cells to sort, for different values you can get different patterns.
         *
         * To run the simulation, we call {{{Solve()}}}. */
        simulator.Solve();

        /* The next two lines are for test purposes only and are not part of this tutorial.
         */
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 64u);
        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), 20.0, 1e-10);
    }

    /*
     * EMPTYLINE
     *
     * To visualize the results, open a new terminal, {{{cd}}} to the Chaste directory,
     * then {{{cd}}} to {{{anim}}}. Then do: {{{java Visualize2dVertexCells /tmp/$USER/testoutput/PottsMonolayerCellSorting/results_from_time_0}}}.
     * We may have to do: {{{javac Visualize2dVertexCells.java}}} beforehand to create the
     * java executable.
     *
     *  You could also visualize in paraview as above.
     *
     * EMPTYLINE
     *
     * == Test 3 - 3D Cell Sorting ==
     *
     * EMPTYLINE
     *
     * The next test extends the previous example to three dimensions.
     *
     */
    void TestPottsSpheroidCellSorting()
    {
        /** The next line is needed because we cannot currently run Potts simulations in parallel. */
        EXIT_IF_PARALLEL;

        /* First, we generate a Potts mesh. To create a {{{PottsMesh}}}, we can use
         * the {{{PottsMeshGenerator}}}. This generates a regular square-shaped mesh,
         * in which all elements are the same size.
         *
         * Here the first three arguments specify the domain width; the number of elements across; and the width of
         * elements. The second set of three arguments specify the domain height; the number of elements up; and
         * the height of individual elements.  The third set of three arguments specify the domain depth; the number of elements deep; and
         * the depth of individual elements.
         * We have chosen an 4 by 4 by 4 ( = 64) block of elements each consisting of 2 by 2 by 2 ( = 8) lattice sites.
         */
        PottsMeshGenerator<3> generator(10, 4, 2, 10, 4, 2, 10, 4, 2);  // Parameters are: lattice sites across; num elements across; element width; lattice sites up; num elements up; and element height; lattice sites deep; num elements deep; and element depth
        PottsMesh<3>* p_mesh = generator.GetMesh();

        /* Having created a mesh, we now create a {{{std::vector}}} of {{{CellPtr}}}s.
         * To do this, we the `CellsGenerator` helper class, as before but this time
         * the third argument is set to make all cells non-proliferative.*/
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<UniformCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_diff_type);

        /* As for the 2D case before we make a {{{CellPopulation}}} we make a pointer to a cell label and then assign this
         * label to some randomly chosen cells. */
        MAKE_PTR(CellLabel, p_label);
        for (unsigned i = 0; i<cells.size(); i++)
        {
            if (RandomNumberGenerator::Instance()->ranf() < 0.5)
            {
                cells[i]->AddCellProperty(p_label);
            }
        }

        /* Now we have a mesh and a set of cells to go with it, we can create a {{{CellPopulation}}}.
         * In general, this class associates a collection of cells with a set of elements or a mesh.
         * For this test, because we have a {{{PottsMesh}}}, we use a particular type of
         * cell population called a {{{PottsBasedCellPopulation}}}.
         */
        PottsBasedCellPopulation<3> cell_population(*p_mesh, cells);

        /* In order to visualize labelled cells we need to use the following command.*/
        cell_population.AddCellWriter<CellLabelWriter>();

        /* We then pass in the cell population into an {{{OnLatticeSimulation}}},
         * and set the output directory and end time. */
        OnLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory("PottsCellSorting3D");
        simulator.SetEndTime(20.0);
        simulator.SetSamplingTimestepMultiple(10);

        /* We must now create one or more update rules, which determine the Hamiltonian
         * in the Potts simulation. For this test, we use two update rules based upon
         * an area constraint and differential adhesion between cells and pass them to the {{{OnLatticeSimulation}}}.
         */
        MAKE_PTR(VolumeConstraintPottsUpdateRule<3>, p_volume_constraint_update_rule);
        /*
         * Now set the target volume to be appropriate for this 3D simulation.
         */
        p_volume_constraint_update_rule->SetMatureCellTargetVolume(8.0);
        p_volume_constraint_update_rule->SetDeformationEnergyParameter(0.2);
        simulator.AddUpdateRule(p_volume_constraint_update_rule);

        /*
         * We use the same differential adhesion parameters as in the 2D case.
         */
        MAKE_PTR(DifferentialAdhesionPottsUpdateRule<3>, p_differential_adhesion_update_rule);
        p_differential_adhesion_update_rule->SetLabelledCellLabelledCellAdhesionEnergyParameter(0.16);
        p_differential_adhesion_update_rule->SetLabelledCellCellAdhesionEnergyParameter(0.11);
        p_differential_adhesion_update_rule->SetCellCellAdhesionEnergyParameter(0.02);
        p_differential_adhesion_update_rule->SetLabelledCellBoundaryAdhesionEnergyParameter(0.16);
        p_differential_adhesion_update_rule->SetCellBoundaryAdhesionEnergyParameter(0.16);
        simulator.AddUpdateRule(p_differential_adhesion_update_rule);

        /* To run the simulation, we call {{{Solve()}}}. */
        simulator.Solve();

        /* The next two lines are for test purposes only and are not part of this tutorial.
         */
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 64u);
        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), 20.0, 1e-10);
    }
    /*
     * EMPTYLINE
     *
     * To visualize the results, we need to use Paraview. Note that we don't output the cell boundaries (outlines) in 3D.
     * See UserTutorials/VisualizingWithParaview for more information.
     *
     * Load the file {{{/tmp/$USER/testoutput/PottsCellSorting3D/results_from_time_0/results.pvd}}}, and click apply.
     *
     * Add box "Glyphs" to represent lattice sites. You will need to adjust the size so they don't overlap.
     *
     * Note that, for larger simulations, you may need to unclick "Mask Points" (or similar) so as not to limit the number of glyphs
     * displayed by Paraview.
     *
     * Select the "Display" tab and select "color by" cell label (you can also "color by" cell index to see individual cells).
     *
     * Add a "Threshold" filter, filter by cell type and make the lower threshold 0 or greater (unoccupied lattice sites are labelled with -1). This will allow you to view only the cells.
     *
     * Click play to see the evolution of the simulation.
     *
     * You should see that the cells sort into ones of the same type.
     *
     * EMPTYLINE
     */
};

#endif /* TESTRUNNINGPOTTSBASEDSIMULATIONSTUTORIAL_HPP_ */
