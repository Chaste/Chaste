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

#ifndef TESTRUNNIGCONTACTINHIBITIONSIMULATIONSTUTORIAL_HPP_
#define TESTRUNNIGCONTACTINHIBITIONSIMULATIONSTUTORIAL_HPP_

/*
 * = An example showing how to use a contact inhibition cell cycle model and volume tracking simulation modifier =
 *
 * == Introduction ==
 *
 * In this tutorial, we will show how to use a simple implementation of the contact inhibition cell-cycle model,
 * which prevents a cell from undergoing division when its volume is smaller than a critical value.
 *
 * We first consider two mesh-based populations in 2D with cells trapped in a square box. In the first population,
 * all the cells are contact inhibited and we study the effect of the critical volume upon the global cell density. In the
 * second population, we consider a mix of 'normal' cells (contact inhibited) and 'tumour' cells (not inhibited) and study the growth of the
 * tumour cells within the box.
 *
 * We then go on to consider the behaviour of a vertex-based population in a box that experiences contact inhibition.
 *
 * == Including header files ==
 *
 * We begin by including the necessary header files. */
#include <cxxtest/TestSuite.h>
/* This header needs to be included here to ensure archiving of {{{CelwiseData}}} works on all Boost versions. */
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"

/* The next header includes the Boost {{{shared_ptr}}} smart pointer, and defines some useful
 * macros to save typing when using it. */
#include "SmartPointers.hpp"
/* The next header include the {{{NEVER_REACHED}}} macro, which is used in one of the methods below. */
#include "Exception.hpp"

/*
 * The next header file defines the contact inhibition cell-cycle model that inherits from {{{AbstractCellCycleModel}}}.
 * The duration of the G1 phase depends on the deviation from a 'target' volume (or area/length in 2D/1D): if a cell's volume is
 * lower than a given fraction of its target volume, the G1 phase continues.
 * This model of cell-cycle progression allows for quiescence imposed by transient periods of high stress, followed by relaxation. Note that
 * in this cell cycle model, quiescence is implemented only by extending the G1 phase. Therefore, if a cell
 * is compressed during G2 or S phases then it will still divide, and thus cells whose volumes are smaller
 * than the given threshold may still divide.
 *
 * The target volume and the critical fraction are specified using the methods {{{SetEquilibriumVolume()}}} and {{{SetQuiescentVolumeFraction()}}} respectively.
 * Within the {{{ContactInhibitionCellCycleModel}}}'s {{{UpdateCellCyclePhase()}}} method these parameters are compared to the actual cell volumes, which are stored
 * using the cell property {{{CellData}}}.
 */
#include "ContactInhibitionCellCycleModel.hpp"

/*
 * The next header defines the simulation class modifier corresponding to the contact inhibition cell-cycle model.
 * This modifier updates the {{{CellData}}} cell property at each timestep with the volume of each cell.
 */
#include "VolumeTrackingModifier.hpp"

/* The remaining header files define classes that will be also be used and are presented in other tutorials. */
#include "OffLatticeSimulation.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "CellMutationStatesCountWriter.hpp"
#include "UniformCellCycleModel.hpp"
#include "WildTypeCellMutationState.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "OutputFileHandler.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "NagaiHondaForce.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "SimulationTime.hpp"
#include "CellLabel.hpp"
#include "MutableMesh.hpp"
#include "MutableVertexMesh.hpp"
#include "PlaneBoundaryCondition.hpp"
#include "FakePetscSetup.hpp"

/* We first define the global test class that inherits from {{{AbstractCellBasedTestSuite}}}. */
class TestRunningContactInhibitionSimulationsTutorial : public AbstractCellBasedTestSuite
{
public:
    /*
     * == Testing healthy cell contact inhibition with mesh-based population ==
     *
     * In this first test we show how to simulate the behaviour of cells healthy cells trapped in a box.
     * Each cell will only divide if there is sufficient room.
     */
    void TestContactInhibitionInBox()
    {
        /* We use the honeycomb mesh generator to create a honeycomb mesh and
         * the associated mutable mesh. */
        HoneycombMeshGenerator generator(3, 3);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        /* We now create a vector of cell pointers. */
        std::vector<CellPtr> cells;

        /* We then define the mutation state of the cells we are working with. We will just consider
         * wild type mutations here. */
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);

        /* We now create a cell-cycle (only contact inhibited) model for these cells and loop over the
         * nodes of the mesh to create as many elements in the vector of cell pointers as there are
         * in the initial mesh. */
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            ContactInhibitionCellCycleModel* p_cycle_model = new ContactInhibitionCellCycleModel();
            p_cycle_model->SetDimension(2);
            p_cycle_model->SetBirthTime(-2.0*(double)i);
            p_cycle_model->SetQuiescentVolumeFraction(0.5);
            p_cycle_model->SetEquilibriumVolume(1.0);

            CellPtr p_cell(new Cell(p_state, p_cycle_model));
            p_cell->SetCellProliferativeType(p_transit_type);
            p_cell->InitialiseCellCycleModel();

            cells.push_back(p_cell);
        }

        /* We now create a cell population, that takes several inputs: the mesh (for the position); and
         * the vector of cell pointers (for cycles and states)*/
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        /* In order to visualize labelled cells (i.e. those that are inhibited from division) you need to use the following command.*/
        cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();

        /* Here we create a simulation as before. We also set up the output directory, the end time and the output multiple.*/
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestContactInhibitionInBox");
        simulator.SetSamplingTimestepMultiple(12);
        simulator.SetEndTime(20.0);

        /* Then, we define the modifier class, which automatically updates the volumes of the cells in {{{CellData}}} and passes it to the simulation.*/
        MAKE_PTR(VolumeTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        /* Next, we create a force law (springs) to be applied between cell centres and set up a
         * cut-off length beyond which cells stop interacting. We then pass this to the {{{VolumeTrackedOffLatticeSimulation}}}. */
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        p_force->SetCutOffLength(1.5);
        simulator.AddForce(p_force);

        /*
         * To study the behaviour of the cells with varying volume, we trap them in the square domain [0,2.5]x[0,2.5].
         * This is implemented using four {{{PlaneBoundaryCondition}}} objects.
         * These planes are indicated by a point and a normal and then passed to the {{{VolumeTrackedOffLatticeSimulation}}}.
         * The domain is chosen to be quite small so as to make the test run in a short amount of time.
         */

        /* First we impose a wall at x=0: */
        c_vector<double,2> point = zero_vector<double>(2);
        c_vector<double,2> normal = zero_vector<double>(2);
        normal(0) = -1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc1, (&cell_population, point, normal));
        simulator.AddCellPopulationBoundaryCondition(p_bc1);
        /* Then we impose a wall at x<=2.5: */
        point(0) = 2.5;
        normal(0) = 1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc2, (&cell_population, point, normal));
        simulator.AddCellPopulationBoundaryCondition(p_bc2);
        /* Then we impose a wall at y>0: */
        point(0) = 0.0;
        point(1) = 0.0;
        normal(0) = 0.0;
        normal(1) = -1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc3, (&cell_population, point, normal));
        simulator.AddCellPopulationBoundaryCondition(p_bc3);
        /* Finally we impose a wall at y<2.5: */
        point(1) = 2.5;
        normal(1) = 1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc4, (&cell_population, point, normal));
        simulator.AddCellPopulationBoundaryCondition(p_bc4);

        /* To run the simulation, we call {{{Solve()}}}. */
        simulator.Solve();
    }
    /*
     * EMPTYLINE
     *
     * To visualize the results, open a new terminal, {{{cd}}} to the Chaste directory,
     * then {{{cd}}} to {{{anim}}}. Then do: {{{java Visualize2dCentreCells /tmp/$USER/testoutput/TestContactInhibitionInBox/results_from_time_0}}}.
     * We may have to do: {{{javac Visualize2dCentreCells.java}}} beforehand to create the
     * java executable.
     *
     * You will notice that once the cells are below a certain size they no longer proliferate and turn dark blue in the visualisation.
     *
     * EMPTYLINE
     *
     * == Testing normal and tumour cells with mesh-based population ==
     *
     * We now test the behaviour of a mixture of healthy and tumour cells in a Box. In this test healthy cells will only
     * divide if there is sufficient room whereas tumour cells will divide regardless.
     */
    void TestContactInhibitionInBoxWithMutants()
    {
        /* Just as before we create a simple mesh. */
        HoneycombMeshGenerator generator(3, 3);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        /* We again create the cells. The difference here is that one of the cells is not contact-inhibited, but rather
         * is defined by a {{{UniformCellCycleModel}}}. */
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        std::vector<CellPtr> cells;
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            if (i==1)
            {
                UniformCellCycleModel* p_cycle_model = new UniformCellCycleModel();
                p_cycle_model->SetBirthTime(-14.0);

                CellPtr p_cell(new Cell(p_state, p_cycle_model));
                p_cell->SetCellProliferativeType(p_stem_type);
                p_cell->SetBirthTime(0.0);
                cells.push_back(p_cell);
            }
            else
            {
                ContactInhibitionCellCycleModel* p_cycle_model = new ContactInhibitionCellCycleModel();
                p_cycle_model->SetDimension(2);
                p_cycle_model->SetBirthTime(-2.0*(double)i);
                p_cycle_model->SetQuiescentVolumeFraction(0.8);
                p_cycle_model->SetEquilibriumVolume(1.0);

                CellPtr p_cell(new Cell(p_state, p_cycle_model));
                p_cell->SetCellProliferativeType(p_transit_type);
                p_cell->InitialiseCellCycleModel();
                cells.push_back(p_cell);
            }
        }

        /* We now create a cell population, that takes several inputs: the mesh; and
         * the vector of cell pointers*/
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        /* In order to visualize labelled cells (i.e those that are inhibited from division) you need to use the following command.*/
        cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();

        /* Here we create a simulation as before. We also set up the output directory, the end time and the output multiple.*/
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestContactInhibitionTumourInBox");
        simulator.SetSamplingTimestepMultiple(12);
        simulator.SetEndTime(20.0);

        /* Then, we define the modifier class, which automatically updates the volumes of the cells in {{{CellData}}} and passes it to the simulation.*/
        MAKE_PTR(VolumeTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        /* Next, we create a force law (springs) to be applied between cell centres and set up a
         * cut-off length beyond which cells stop interacting. We then pass this to the {{{VolumeTrackedOffLatticeSimulation}}} */
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        p_force->SetCutOffLength(1.5);
        simulator.AddForce(p_force);

        /* As in the previous test, we trap the cells in the square domain [0,2.5]x[0,2.5]: */
        c_vector<double,2> point = zero_vector<double>(2);
        c_vector<double,2> normal = zero_vector<double>(2);
        normal(0) = -1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc1, (&cell_population, point, normal));
        simulator.AddCellPopulationBoundaryCondition(p_bc1);

        point(0) = 2.5;
        normal(0) = 1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc2, (&cell_population, point, normal));
        simulator.AddCellPopulationBoundaryCondition(p_bc2);

        point(0) = 0.0;
        point(1) = 0.0;
        normal(0) = 0.0;
        normal(1) = -1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc3, (&cell_population, point, normal));
        simulator.AddCellPopulationBoundaryCondition(p_bc3);

        point(1) = 2.5;
        normal(1) = 1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc4, (&cell_population, point, normal));
        simulator.AddCellPopulationBoundaryCondition(p_bc4);

        /* Finally, to run the simulation, we call {{{Solve()}}}. */
        simulator.Solve();
    }
    /*
     * EMPTYLINE
     *
     * To visualize the results, open a new terminal, {{{cd}}} to the Chaste directory,
     * then {{{cd}}} to {{{anim}}}. Then do: {{{java Visualize2dCentreCells /tmp/$USER/testoutput/TestContactInhibitionTumourInBox/results_from_time_0}}}.
     * We may have to do: {{{javac Visualize2dCentreCells.java}}} beforehand to create the
     * java executable.
     *
     * You will notice that once the healthy cells (yellow) are below a certain size they no longer proliferate and turn dark blue in the visualisation.
     * Whereas Tumour cells (light blue) on the other hand will continue to proliferate. You may want to run the simulation for longer to see this more clearly.
     *
     * EMPTYLINE
     *
     * == Testing contact inhibition in vertex-based monolayer ==
     *
     * We now test the behaviour of normal contact inhibited cells for a vertex-based population.
     * The example we use is a growing monolayer.
     *
     */
    void TestContactInhibitionWithVertex()
    {
        /* First we create a simple 2D MutableVertexMesh.*/
        HoneycombVertexMeshGenerator generator(2, 2);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        /*
         * We then create cells as before, only this time we need one per element. We also create the cell population (a {{{VertexBasedCellPopulation}}}).
         */
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        std::vector<CellPtr> cells;
        for (unsigned i=0; i<p_mesh->GetNumElements(); i++)
        {
            ContactInhibitionCellCycleModel* p_cycle_model = new ContactInhibitionCellCycleModel();
            p_cycle_model->SetDimension(2);
            p_cycle_model->SetBirthTime(-(double)i - 2.0); // So all out of M phase
            p_cycle_model->SetQuiescentVolumeFraction(0.9);
            p_cycle_model->SetEquilibriumVolume(1.0);

            CellPtr p_cell(new Cell(p_state, p_cycle_model));
            p_cell->SetCellProliferativeType(p_transit_type);
            p_cell->InitialiseCellCycleModel();
            cells.push_back(p_cell);
        }

        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();

        /* Here we create a simulation as before. We also set up the output directory, the end time and the output multiple.*/
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestVertexContactInhibition");
        simulator.SetSamplingTimestepMultiple(50);
        simulator.SetEndTime(10.0);

        /* Then, we define the modifier class, which automatically updates the volumes of the cells in {{{CellData}}} and passes it to the simulation.*/
        MAKE_PTR(VolumeTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        /* Next, we create a force law, {{{NagaiHondaForce}}}, to be applied to vertices.
         * We then pass this to the {{{VolumeTrackedOffLatticeSimulation}}}. */
        MAKE_PTR(NagaiHondaForce<2>, p_nagai_honda_force);
        simulator.AddForce(p_nagai_honda_force);

        /* A {{{NagaiHondaForce}}} assumes that each cell has a target area assigned to it. In order to
         * assign target areas to each cell and define how they grow, we add a {{{SimpleTargetAreaModifier}}}
         * to the simulator.
         */
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        /* To run the simulation, we call {{{Solve()}}}. */
        simulator.Solve();
    }
    /*
     * EMPTYLINE
     *
     * To visualize the results, open a new terminal, {{{cd}}} to the Chaste directory,
     * then {{{cd}}} to {{{anim}}}. Then do: {{{java Visualize2dVertexCells /tmp/$USER/testoutput/TestVertexContactInhibition/results_from_time_0}}}.
     * We may have to do: {{{javac Visualize2dVertexCells.java}}} beforehand to create the
     * java executable.
     *
     * You will notice that once the healthy cells (yellow) are below a certain size they no longer proliferate and turn dark blue in the visualisation.
     * If you run the simulation for a long time these cells occur primarily towards the centre of the monolayer.
     *
     * EMPTYLINE
     */
};

#endif /*TESTRUNNIGCONTACTINHIBITIONSIMULATIONSTUTORIAL_HPP_*/
