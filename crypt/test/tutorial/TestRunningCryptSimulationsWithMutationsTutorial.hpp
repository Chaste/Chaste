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

#ifndef TESTRUNNINGCRYPTSIMULATIONSWITHMUTATIONSTUTORIAL_HPP_
#define TESTRUNNINGCRYPTSIMULATIONSWITHMUTATIONSTUTORIAL_HPP_

/*
 * = Examples showing how to run crypt simulations with various mutations =
 *
 * == Introduction ==
 *
 * This tutorial assumes you have already read UserTutorials/RunningMeshBasedCryptSimulations.
 *
 * In this tutorial we show how Chaste can be used to simulate a cylindrical model of an
 * intestinal crypt with mutations using both mesh and vertex-based simulations.
 * Full details of the computational model can be found in the paper by
 * Osborne ''et al.'' (2010) [10.1098/rsta.2010.0173].
 *
 * As in previous cell-based Chaste tutorials, we begin by including the necessary header files.
 */
#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "SmartPointers.hpp"
#include "AbstractCellBasedTestSuite.hpp"

/* The next header file defines a helper class for generating cells for crypt simulations. */
#include "CryptCellsGenerator.hpp"
/*
 * The next header file defines a {{{WntCellCycleModel}}}, where the proliferative behaviour of a cell is
 * dependent on the concentration of Wnt at that point in space. Cells proliferate where there is a plentiful level of Wnt
 * and cease proliferation below a given threshold.
 */
#include "SimpleWntCellCycleModel.hpp"
/* The next header file defines a helper class for generating a suitable triangular mesh
 * for the crypt simulation, such that the cell corresponding to each node is initially
 * in mechanical equilibrium with its neighours and periodic boundary conditions are applied
 * at the left- and right-hand sides of the mesh (hence the "cylindrical"). */
#include "CylindricalHoneycombMeshGenerator.hpp"
/* The next header file defines a {{{CellPopulation}}} class that uses a triangular mesh, and allows
 * for the inclusion of 'ghost nodes'. These are nodes in the mesh that do not correspond
 * to cells; instead they help ensure that a sensible Delaunay triangulation is generated
 * at each timestep. This is because the triangulation algorithm requires a convex hull. */
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
/*
 * The next header file defines a force law, based on a linear spring, for describing
 * the mechanical interactions between neighbouring cells in the crypt.
 */
#include "GeneralisedLinearSpringForce.hpp"
/*
 * The next header file defines the class that simulates the evolution of a {{{CellPopulation}}},
 * specialized to deal with the cylindrical crypt geometry.
 */
#include "CryptSimulation2d.hpp"
/*
 * The next header file defines a Wnt singleton class, which (if used) deals with the
 * imposed Wnt gradient in our crypt model. This affects cell proliferation in the case
 * where we construct each cell with a {{{WntCellCycleModel}}}.
 */
#include "WntConcentration.hpp"
/*
 * The next header file defines a cell killer class, which implements sloughing of cells
 * into the lumen once they reach the top of the crypt.
 */
#include "SloughingCellKiller.hpp"
/* These headers are used for defining and recording mutations. */
#include "ApcTwoHitCellMutationState.hpp"
#include "CellMutationStatesCountWriter.hpp"
/* The final header ensures that this test is only ever run sequentially, not in parallel. */
#include "FakePetscSetup.hpp"

/*
 * Next, we define the test class, which inherits from {{{AbstractCellBasedTestSuite}}}
 * and defines some test methods.
 */
class TestRunningCryptSimulationsWithMutationsTutorial : public AbstractCellBasedTestSuite
{
public:
    /* == Test 1: a mesh-based crypt simulation with mutations ==
     *
     * In the first test, we demonstrate how to introduce mutations into a simulation of a crypt.
     */
    void TestMeshBasedCryptWithMutations()
    {
        /* Note that time is re-initialized to zero and the random number generator is re-seeded to zero in the {{{AbstractCellBasedTestSuite}}}.
         *
         * We first create a cylindrical mesh, and get the cell location indices, exactly as before. */
        CylindricalHoneycombMeshGenerator generator(6, 9, 2);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        /* We create the cells, using the same method as before. Here, though, we use a {{{SimpleWntCellCycleModel}}}.*/
        std::vector<CellPtr> cells;
        CryptCellsGenerator<SimpleWntCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, location_indices, true);

        /* We now create boost shared pointers to any mutations we wish to use.
         * We need to do this using the {{{CellPropertyRegistry}}}, otherwise
         * the numbers of each type of mutation aren't correctly tracked. For
         * a list of possible mutations, see subclasses of {{{AbstractCellMutationState}}}.
         * These can be found in the inheritance diagram, here,
         * [class:AbstractCellMutationState AbstractCellMutationState].
         * Each mutation has a different effect on the cell cycle models; see the class
         * documentation for details.
         */
        boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<ApcTwoHitCellMutationState>());

        /* We create the cell population, as before. */
        MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);

        /* In order to visualize mutant cells and to count how many cells there are of each type we need to use the following command.*/
        cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();

        /*
         * We set the height of the crypt. As well as passing this variable into the {{{SloughingCellKiller}}},
         * we will pass it to the {{{WntConcentration}}} object (see below).
         */
        double crypt_height = 8.0;

        /*
         * When using a {{{SimpleWntCellCycleModel}}}, we need a way of telling each cell what the Wnt concentration
         * is at its location. To do this, we set up a {{{WntConcentration}}} object. Like {{{SimulationTime}}},
         * {{{WntConcentration}}} is a singleton class, so when instantiated it is accessible from anywhere in
         * the code (and in particular, all cells and cell-cycle models can access it). We need to say what
         * the profile of the Wnt concentation should be up the crypt: here, we say it is {{{LINEAR}}} (linear
         * decreasing from 1 to 0 from the bottom of the crypt to the top). We also need to inform the
         * {{{WntConcentration}}} of the cell population and the height of the crypt.
         */
        WntConcentration<2>::Instance()->SetType(LINEAR);
        WntConcentration<2>::Instance()->SetCellPopulation(cell_population);
        WntConcentration<2>::Instance()->SetCryptLength(crypt_height);

        /* Create a simulator as before (except setting a different output directory). */
        CryptSimulation2d simulator(cell_population);
        simulator.SetOutputDirectory("MeshBasedCryptWithMutations");
        simulator.SetSamplingTimestepMultiple(12);
        simulator.SetEndTime(10);

        /* As before, we create a force law and cell killer and pass these objects to the simulator, then call
         * Solve(). */
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        simulator.AddForce(p_linear_force);
        MAKE_PTR_ARGS(SloughingCellKiller<2>, p_killer, (&cell_population, crypt_height));
        simulator.AddCellKiller(p_killer);

        simulator.Solve();

        /*
         * Now we have run the simulation to a steady state (where the initial regular configuration is lost) we select a cell to become mutant.
         * We select one of the cells and set the mutation state to {{{ApcTwoHitCellMutationState}}} (i.e. p_state).
         */
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            unsigned node_index = cell_population.GetLocationIndexUsingCell(*cell_iter);

            if (node_index == 74) // Chosen from looking at the results from steady state
            {
                cell_iter->SetMutationState(p_state);
            }
        }

       /*
        * We also change the value of the damping constant for mutant cells to be
        * ten times the normal value.
        */
       double normal_damping_constant = cell_population.GetDampingConstantNormal();
       cell_population.SetDampingConstantMutant(10*normal_damping_constant);

       /* Next we reset the end time to some later time. */
       simulator.SetEndTime(20);

       /* Run the simulation to the new end time. */
       simulator.Solve();

       /*
        * Finally, we must tidy up by destroying the {{{WntConcentration}}}
        * singleton object. This avoids memory leaks occurring.
        */
       WntConcentration<2>::Destroy();
    }
    /* To visualize the results, open a new terminal, {{{cd}}} to the Chaste directory,
     * then {{{cd}}} to {{{anim}}}. Then do: {{{java Visualize2dCentreCells /tmp/$USER/testoutput/MeshBasedCryptWithMutations/results_from_time_0}}}.
     *
     * These are the results before we add the mutations do: {{{java Visualize2dCentreCells /tmp/$USER/testoutput/MeshBasedCryptWithMutations/results_from_time_10}}}
     * to see the results from after the mutation has been added.
     *
     * We may have to do: {{{javac Visualize2dCentreCells.java}}} beforehand to create the
     * java executable.
     *
     * In the results folder there is also a file {{{cellmutationstates.dat}}} which tracks the numbers of each mutation type in the simulation.
     * These results are just tab separated columns so may be visualized by using gnuplot, Matlab or similar.
     */
};

#endif /*TESTRUNNINGCRYPTSIMULATIONSWITHMUTATIONSTUTORIAL_HPP_*/
