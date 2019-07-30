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

#ifndef TESTRUNNINGMESHBASEDSIMULATIONSTUTORIAL_HPP_
#define TESTRUNNINGMESHBASEDSIMULATIONSTUTORIAL_HPP_

/*
 * = Examples showing how to create, run and visualize mesh-based simulations =
 *
 * == Introduction ==
 *
 * In this tutorial we show how Chaste can be used to create, run and visualize mesh-based simulations.
 * Full details of the mathematical model can be found in van Leeuwen ''et al.'' (2009) [doi:10.1111/j.1365-2184.2009.00627.x].
 *
 * == The test ==
 *
 * We begin by including the necessary header files. The first thing to do is include the
 * following header file, which allows us to use certain methods in our test. This header
 * file must be included in any Chaste test.
 */
#include <cxxtest/TestSuite.h>
/* The following header is usually included in all cell-based test suites.  It enables us to write tests where the
 * {{{SimulationTime}}} is handled automatically and simplifies the tests. It also sets up the random number generator
 * and the {{{CellPropertyRegistry}}}. You will learn about both of them in later tutorials.
 */
#include "AbstractCellBasedTestSuite.hpp"
/* Any test in which the {{{GetIdentifier()}}} method is used, even via the main
 * `cell_based` code (through calls to {{{AbstractCellPopulation}}} output methods),
 * must also include {{{CheckpointArchiveTypes.hpp}}} or {{{CellBasedSimulationArchiver.hpp}}}
 * as the first Chaste header file.
 */
#include "CheckpointArchiveTypes.hpp"
/* The next header includes the Boost shared_ptr smart pointer, and defines some useful
 * macros to save typing when using it.
 */
#include "SmartPointers.hpp"
/* The remaining header files define classes that will be used in the cell population
 * simulation test. The first defines a helper class for generating a suitable collection
 * of cells. */
#include "CellsGenerator.hpp"
#include "TransitCellProliferativeType.hpp"
/* The next header file defines a stochastic cell-cycle model class. */
#include "UniformCellCycleModel.hpp"
/* The next header file defines a helper class for generating a suitable mesh. */
#include "HoneycombMeshGenerator.hpp"
/* The next header file defines the class that simulates the evolution of an off-lattice {{{CellPopulation}}}. */
#include "OffLatticeSimulation.hpp"
/* The next header files define classes for mesh-based {{{CellPopulation}}}s with and without ghost nodes.*/
#include "MeshBasedCellPopulation.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
/* The next header file defines a force law for describing the mechanical interactions
 * between neighbouring cells in the cell population.
 */
#include "GeneralisedLinearSpringForce.hpp"
/* The next header file defines a class for writing output that can be visualized in Paraview. */
#include "VoronoiDataWriter.hpp"
/* Finally the following header ensures that the test never runs in parallel. */
#include "FakePetscSetup.hpp"
/* Next, we define the test class.
 */
class TestRunningMeshBasedSimulationsTutorial : public AbstractCellBasedTestSuite
{
public:
    /*
     * == Test 1 - a basic mesh-based simulation ==
     *
     * In the first test, we run a simple mesh-based simulation, in which we create a monolayer
     * of cells, using a mutable mesh. Each cell is assigned a stochastic cell-cycle model.
     */
    void TestMonolayer()
    {
        /* Next, we generate a mutable mesh. To create a {{{MutableMesh}}}, we can use
         * the {{{HoneycombMeshGenerator}}}. This generates a honeycomb-shaped mesh,
         * in which all nodes are equidistant. Here the first and second arguments
         * define the size of the mesh - we have chosen a mesh that is 2 nodes (i.e.
         * cells) wide, and 2 nodes high.
         */
        HoneycombMeshGenerator generator(2, 2);    // Parameters are: cells across, cells up
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        /* Having created a mesh, we now create a {{{std::vector}}} of {{{CellPtr}}}s.
         * To do this, we use the `CellsGenerator` helper class, which is templated over the type
         * of cell cycle model required (here {{{UniformCellCycleModel}}})
         * and the dimension.
         * For a list of possible cell cycle models see subclasses of {{{AbstractCellCycleModel}}}.
         * These can be found in the inheritance diagram, here, [class:AbstractCellCycleModel AbstractCellCycleModel].
         * Note that some of these models will require information on the surrounding medium such as Oxygen concentration to work,
         * see specific class documentation for details. Some of these will be covered in later tutorials (UserTutorials/RunningContactInhibitionSimulations,
         * UserTutorials/RunningDeltaNotchSimulations, and UserTutorials/RunningTumourSpheroidSimulations).
         * We create an empty vector of cells and pass this into the
         * method along with the mesh. The second argument represents the size of that the vector
         * {{{cells}}} should become - one cell for each node, the third argument specifies
         * the proliferative type of the cell. */
        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes(), p_transit_type);

        /* Now we have a mesh and a set of cells to go with it, we can create a {{{CellPopulation}}}.
         * In general, this class associates a collection of cells with a mesh.
         * For this test, because we have a {{{MutableMesh}}}, we use a particular type of
         * cell population called a {{{MeshBasedCellPopulation}}}.
         */
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        /* To view the results of this and the next test in Paraview it is necessary to explicitly
        * generate the required .vtu files. This is detailed in the [wiki:UserTutorials/VisualizingWithParaview] tutorial.
        * Note that the results in Paraview may appear different to those in the java based visualizer. This is related
        * to the different methods used to generate voronoi tesselations in each and is resolved through the use of
        * 'ghost nodes', as shown in the next test. */
        cell_population.AddPopulationWriter<VoronoiDataWriter>();

        /* We then pass in the cell population into an {{{OffLatticeSimulation}}},
         * and set the output directory and end time. */
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("MeshBasedMonolayer");
        simulator.SetEndTime(10.0);

        /*
         * For longer simulations, we may not want to output the results
         * every time step. In this case we can use the following method,
         * to print results every 12 time steps instead. As the default time step
         * used by the simulator is 30 seconds, this method will cause the
         * simulator to print results every 6 minutes (or 0.1 hours).
         */
        simulator.SetSamplingTimestepMultiple(12);

        /* We must now create one or more force laws, which determine the mechanics of the centres
         * of each cell in a cell population. For this test, we use one force law, based on the
         * spring based model, and pass it to the {{{OffLatticeSimulation}}}.
         * For a list of possible forces see subclasses of {{{AbstractForce}}}.
         * These can be found in the inheritance diagram, here, [class:AbstractForce AbstractForce].
         * Note that some of these forces are not compatible with mesh-based simulations,
         * see the specific class documentation for details.  If you try to use an incompatible class
         * then you will receive a warning.
         */
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        simulator.AddForce(p_force);

        /* To run the simulation, we call {{{Solve()}}}. */
        simulator.Solve();

        /* The next two lines are for test purposes only and are not part of this tutorial. If different simulation input parameters are being explored
         * the lines should be removed.*/
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 8u);
        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), 10.0, 1e-10);
    }

    /*
     * To visualize the results, open a new terminal, {{{cd}}} to the Chaste directory,
     * then {{{cd}}} to {{{anim}}}. Then do: {{{java Visualize2dCentreCells /tmp/$USER/testoutput/MeshBasedMonolayer/results_from_time_0}}}.
     * We may have to do: {{{javac Visualize2dCentreCells.java}}} beforehand to create the
     * java executable.
     *
     * For further details on visualization, see [wiki:ChasteGuides/RunningCellBasedVisualization].
     *
     * You will notice that half of each cell cell around the edge is missing.
     * This is because the Voronoi region for nodes on the edge of the mesh can be
     * infinite, therefore we only visualize the part inside the mesh.
     * This also means there may be "long" edges in the mesh which can cause the cells
     * to move due long range interactions resulting in an artificially rounded shape.
     *
     * There are two solutions to this. The first is to define a cut off length on the force,
     * which can be done by using the command
     * {{{p_force->SetCutOffLength(1.5);}}}
     * on the {{{GeneralisedLinearSpringForce}}}. Here there will be no forces exerted
     * on any "springs" which are longer than 1.5 cell radii.
     *
     * The second solution is to use 'ghost nodes'. Ghost nodes can be added to mesh-based
     * simulations to remove infinite voronoi regions and long edges. To do this, a set of
     * nodes (known as ghost nodes) are added around the original mesh which exert forces
     * on each other but do not exert forces on the nodes of the original mesh (known as
     * real nodes). In addition real nodes exert forces on ghost nodes so the ghost nodes
     * remain surrounding the cell population.
     *
     * == Test 2 - a basic mesh-based simulation with ghost nodes ==
     *
     * In the second test, we run a simple mesh-based simulation with ghost nodes, in which we
     * create a monolayer of cells, using a mutable mesh.
     * Each cell is assigned a stochastic cell-cycle model.
     */
    void TestMonolayerWithGhostNodes()
    {
        /* We start by generating a mutable mesh. To create a {{{MutableMesh}}}, we can use
         * the {{{HoneycombMeshGenerator}}} as before. Here the first and second arguments
         * define the size of the mesh - we have chosen a mesh that is 2 nodes (i.e.
         * cells) wide, and 2 nodes high.  The third argument specifies the number of layers
         * of ghost nodes to make.
         */
        HoneycombMeshGenerator generator(2, 2, 2);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        /* We only want to create cells to attach to real nodes, so we
         * use the method {{{GetCellLocationIndices}}} to get the indices
         * of the real nodes in the mesh. This will be passed in to the
         * cell population later on.
         */
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        /* Having created a mesh, we now create a {{{std::vector}}} of {{{CellPtr}}}s.
         * To do this, we the `CellsGenerator` helper class again. This time the second
         * argument is different and is the number of real nodes in the mesh.
         * As before all cells have {{{TransitCellProliferativeType}}}. */
        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, location_indices.size(), p_transit_type);

        /* Now we have a mesh and a set of cells to go with it, we can create a {{{CellPopulation}}}.
         * In general, this class associates a collection of cells with a set of elements or a mesh.
         * For this test, because we have a {{{MutableMesh}}}, and ghost nodes we use a particular type of
         * cell population called a {{{MeshBasedCellPopulationWithGhostNodes}}}. The third
         * argument of the constructor takes a vector of the indices of the real nodes and should be the
         * same length as the vector of cell pointers.
         */
        MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices); //**Changed**//

        /* Again Paraview output is explicitly requested.*/
        cell_population.AddPopulationWriter<VoronoiDataWriter>();

        /* We then pass in the cell population into an {{{OffLatticeSimulation}}},
         * and set the output directory, output multiple and end time. */
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("MeshBasedMonolayerWithGhostNodes");
        simulator.SetSamplingTimestepMultiple(12);
        simulator.SetEndTime(10.0);

        /* Again we create a force law, and pass it to the {{{OffLatticeSimulation}}}. This
         * force law ensures that ghost nodes don't exert forces on real nodes but real nodes
         * exert forces on ghost nodes.*/
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        simulator.AddForce(p_force);

        /* To run the simulation, we call {{{Solve()}}}. */
        simulator.Solve();

        /* The next two lines are for test purposes only and are not part of this tutorial.
         */
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 8u);
        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), 10.0, 1e-10);
    }
    /*
     * To visualize the results, open a new terminal, {{{cd}}} to the Chaste directory,
     * then {{{cd}}} to {{{anim}}}. Then do: {{{java Visualize2dCentreCells /tmp/$USER/testoutput/MeshBasedMonolayerWithGhostNodes/results_from_time_0}}}.
     */
};

#endif /* TESTRUNNINGMESHBASEDSIMULATIONSTUTORIAL_HPP_ */
