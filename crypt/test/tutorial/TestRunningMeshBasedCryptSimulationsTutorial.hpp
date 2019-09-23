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

#ifndef TESTRUNNINGMESHBASEDCRYPTSIMULATIONSTUTORIAL_HPP_
#define TESTRUNNINGMESHBASEDCRYPTSIMULATIONSTUTORIAL_HPP_

/*
 * = Examples showing how to run crypt simulations on periodic meshes with different cell-cycle models =
 *
 * EMPTYLINE
 *
 * == Introduction ==
 *
 * EMPTYLINE
 *
 * In this tutorial we show how Chaste can be used to simulate a cylindrical model of an
 * intestinal crypt. Full details of the computational model can be found in the paper by
 * van Leeuwen ''et al.'' (2009) [doi:10.1111/j.1365-2184.2009.00627.x].
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
 * The next two header files define two different types of cell-cycle model.
 * In a {{{FixedG1GenerationalCellCycleModel}}}, the duration of each phase
 * of the cell cycle is fixed. In a {{{WntCellCycleModel}}}, the duration of a cell's G1 phase
 * is determined by a system of nonlinear ODEs describing a cell's response to the local
 * concentration of Wnt,
 * a secreted cell–cell signalling molecule that is known to play a key role in cell
 * proliferation in the crypt. In our crypt simulations, we impose a fixed gradient of
 * Wnt up the axis of the crypt.
 */
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "WntCellCycleModel.hpp"
/* The next header file defines a helper class for generating a suitable triangular mesh
 * for the crypt simulation, such that the cell corresponding to each node is initially
 * in mechanical equilibrium with its neighours and periodic boundary conditions are applied
 * at the left- and right-hand sides of the mesh (hence the "cylindrical"). */
#include "CylindricalHoneycombMeshGenerator.hpp"
/* The next two header files were encountered in UserTutorials/RunningMeshBasedSimulations.
 * The first header
 * defines a {{{CellPopulation}}} class that uses a triangular mesh, and allows
 * for the inclusion of 'ghost nodes': these are nodes in the mesh that do not correspond
 * to cells, but help ensure that a sensible Delaunay triangulation is generated
 * at each timestep; this is because the triangulation algorithm requires a convex hull.
 * The next header file defines a force law, based on a linear spring, for describing
 * the mechanical interactions between neighbouring cells in the crypt.*/
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
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
 * The final header file defines a cell killer class, which implements sloughing of cells
 * into the lumen once they reach the top of the crypt.
 */
#include "SloughingCellKiller.hpp"

/* This header ensures that this test is only run on one process, since it doesn't support parallel execution. */
#include "FakePetscSetup.hpp"

/* Next, we define the test class. */
class TestRunningMeshBasedCryptSimulationsTutorial : public AbstractCellBasedTestSuite
{
public:
    /* EMPTYLINE
     *
     * == Test 1: a basic crypt simulation ==
     *
     * EMPTYLINE
     *
     * In the first test, we demonstrate how to create a crypt simulation using a
     * cylindrical mesh, with each cell progressing through a fixed cell-cycle model,
     * and sloughing enforced at the top of the crypt.
     */
    void TestCryptWithFixedCellCycle()
    {
        /* First, we generate a mesh. The basic Chaste mesh is a {{{TetrahedralMesh}}}.
         * To enforce periodicity at the left- and right-hand sides of the mesh, we
         * use a subclass called {{{Cylindrical2dMesh}}}, which has extra methods for
         * maintaining periodicity. To create a {{{Cylindrical2dMesh}}}, we can use a helper class called
         * {{{CylindricalHoneycombMeshGenerator}}}. This generates a periodic honeycomb-shaped mesh,
         * in which all nodes are equidistant to their neighbours. Here the first and second arguments
         * define the size of the mesh - we have chosen a mesh that is 6 nodes (i.e.
         * cells) wide, and 9 nodes high. The third argument indicates that we require
         * a double layer of ghost nodes around the mesh (technically, just above
         * and below the mesh, since it is periodic). We call {{{GetCylindricalMesh()}}} on the {{{CylindricalHoneycombMeshGenerator}}} to
         * return our {{{Cylindrical2dMesh}}}, and call {{{ GetCellLocationIndices()}}}
         * to return a {{{std::vector}}} of indices of nodes in the mesh that correspond to real cells (as opposed
         * to ghost nodes).
         */
        CylindricalHoneycombMeshGenerator generator(6, 9, 2);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        /*
         * Having created a mesh, we now create a {{{std::vector}}} of {{{CellPtr}}}s.
         * To do this, we use the `CryptCellsGenerator` helper class, which is templated over the type
         * of cell-cycle model required (here {{{FixedG1GenerationalCellCycleModel}}})
         * and the dimension. We create an empty vector of cells and pass this into the
         * method {{{Generate()}}} along with the mesh. The fourth argument 'true' indicates that the cells
         * should be assigned random birth times, to avoid synchronous division. The
         * {{{cells}}} vector is populated once the method {{{Generate()}}} is
         * called. Note that we only ever deal with shared pointers to cells, named {{{CellPtr}}}s.
         */
        std::vector<CellPtr> cells;
        CryptCellsGenerator<FixedG1GenerationalCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, location_indices, true);

        /*
         * Now we have a mesh, a set of cells to go with it, and a vector of node indices
         * corresponding to real cells, we can create a {{{CellPopulation}}} object. In general,
         * this class associates a collection of cells with a set of nodes or a mesh.
         * For this test, because we have a mesh and ghost nodes, we use a particular type of
         * cell population called a {{{MeshBasedCellPopulationWithGhostNodes}}}.
         */
        MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);

        /*
         * Next we use the {{{CellPopulation}}} object to construct a {{{CryptSimulation2d}}} object,
         * which will be used to simulate the crypt model. */
        CryptSimulation2d simulator(cell_population);
        /*
         * We must set the output directory on the simulator (relative to
         * "/tmp/<USER_NAME>/testoutput") and the end time (in hours).
         */
        simulator.SetOutputDirectory("CryptTutorialFixedCellCycle");
        simulator.SetEndTime(1);
        /*
         * For longer simulations, we may not want to output the results
         * every time step. In this case we can use the following method,
         * to print results every 12 time steps instead. As the time step
         * used by the simulator, is 30 seconds, this method will cause the
         * simulator to print results every 6 minutes.
         */
        simulator.SetSamplingTimestepMultiple(12);

        /*
         * Before running the simulation, we must add one or more force laws, which determine the mechanical
         * behaviour of the cell population. For this test, we use a {{{GeneralisedLinearSpringForce}}}, which assumes
         * that every cell experiences a force from each of its neighbours that can be represented as a linear overdamped
         * spring.
         */
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        simulator.AddForce(p_linear_force);

        /*
         * We also add a cell killer to the simulator. This object
         * dictates under what conditions cells die. For this test, we use
         * a {{{SloughingCellKiller}}}, which kills cells above a certain
         * height (passed as an argument to the constructor).
         */
        double crypt_height = 8.0;
        MAKE_PTR_ARGS(SloughingCellKiller<2>, p_killer, (&cell_population, crypt_height));
        simulator.AddCellKiller(p_killer);

        /* To run the simulation, we call {{{Solve()}}}. */
        simulator.Solve();
    }

    /*
     * EMPTYLINE
     *
     * Finally, to visualize the results, we open a new terminal, {{{cd}}} to the Chaste directory,
     * then {{{cd}}} to {{{anim}}}. Then we do: {{{java Visualize2dCentreCells /tmp/$USER/testoutput/CryptTutorialFixedCellCycle/results_from_time_0}}}.
     * It may be necessary to do: {{{javac Visualize2dCentreCells.java}}} beforehand to create the
     * java executable. Further details on visualization can be found on the Chaste wiki page
     * For further details on visualization, see ChasteGuides/RunningCellBasedVisualization.
     *
     * EMPTYLINE
     *
     * == Test 2: a Wnt-dependent crypt simulation ==
     *
     * EMPTYLINE
     *
     * The next test is very similar to Test 1, except that instead of
     * using a fixed cell-cycle model, we use a Wnt-dependent cell-cycle model,
     * with the Wnt concentration varying within the crypt in a predefined manner.
     */
    void TestCryptWithWntCellCycle()
    {
        /* First we create a cylindrical mesh, and get the cell location indices, exactly as before.
         * Note that time is re-initialized to zero and random number generator is re-seeded to zero in the {{{AbstractCellBasedTestSuite}}}.*/
        CylindricalHoneycombMeshGenerator generator(6, 9, 2);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        /* We create the cells, using the same method as before. Here, though, we use a {{{WntCellCycleModel}}}.*/
        std::vector<CellPtr> cells;
        CryptCellsGenerator<WntCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, location_indices, true);

        /* We create the cell population, as before. */
        MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);

        /*
         * We set the height of the crypt. As well as passing this variable into the {{{sloughingCellKiller}}},
         * we will pass it to the {{{WntConcentration}}} object (see below).
         */
        double crypt_height = 8.0;

        /*
         * When using a {{{WntCellCycleModel}}}, we need a way of telling each cell what the Wnt concentration
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
        simulator.SetOutputDirectory("CryptTutorialWntCellCycle");
        simulator.SetEndTime(1);

        /* As before, we create a force law and cell killer and pass these objects to the simulator, then call
         * Solve(). */
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        simulator.AddForce(p_linear_force);
        MAKE_PTR_ARGS(SloughingCellKiller<2>, p_killer, (&cell_population, crypt_height));
        simulator.AddCellKiller(p_killer);

        simulator.Solve();

        /*
         * Finally, we must tidy up by destroying the {{{WntConcentration}}}
         * singleton object. This avoids memory leaks occurring. */
        WntConcentration<2>::Destroy();
    }
};
    /*
     * EMPTYLINE
     *
     * The results of this test can be visualized as in Test 1, with the correct output directory.
     */
#endif /*TESTRUNNINGMESHBASEDCRYPTSIMULATIONSTUTORIAL_HPP_*/
