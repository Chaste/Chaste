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

#ifndef TESTRUNNINGNODEBASEDSIMULATIONSTUTORIAL_HPP_
#define TESTRUNNINGNODEBASEDSIMULATIONSTUTORIAL_HPP_

/*
 * = Examples showing how to create, run and visualize node-based simulations =
 *
 * EMPTYLINE
 *
 * == Introduction ==
 *
 * EMPTYLINE
 *
 * In this tutorial we show how Chaste can be used to create, run and visualize node-based simulations.
 * Full details of the mechanical model can be found in Pathamathan et al "A computational study of
 * discrete mechanical tissue models", Physical Biology. Vol. 6. No. 3. 2009.. DOI (10.1088/1478-3975/6/3/036001).
 *
 * EMPTYLINE
 *
 * == The test ==
 *
 * EMPTYLINE
 *
 * As in previous cell-based Chaste tutorials (UserTutorials/RunningMeshBasedSimulations), we begin by including the necessary header files.
 */
#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"

/* The following header is usually included in all cell-based test suites. It enables us to write tests where the {{{SimulationTime}}} is handled automatically and simplifies the tests.*/
#include "AbstractCellBasedTestSuite.hpp"
#include "PetscSetupAndFinalize.hpp"
/* The remaining header files define classes that will be used in the cell population
 * simulation test. We encountered some of these header files in
 * UserTutorials/RunningMeshBasedSimulations. */
#include "CellsGenerator.hpp"
#include "TransitCellProliferativeType.hpp"
#include "UniformCellCycleModel.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "OffLatticeSimulation.hpp"
#include "SmartPointers.hpp"
/* The next header file defines the class for storing the spatial information of cells. */
#include "NodesOnlyMesh.hpp"
/* The next header file defines a node-based {{{CellPopulation}}} class.*/
#include "NodeBasedCellPopulation.hpp"
/* The next header file defines a boundary condition to be used in the third test.*/
#include "SphereGeometryBoundaryCondition.hpp"
/* Next, we define the test class.
 */
class TestRunningNodeBasedSimulationsTutorial : public AbstractCellBasedTestSuite
{
public:
    /* EMPTYLINE
     *
     * == Test 1 - a basic node-based simulation ==
     *
     * EMPTYLINE
     *
     * In the first test, we run a simple node-based simulation, in which we create a monolayer
     * of cells, using a nodes only mesh. Each cell is assigned a stochastic cell-cycle model.
     */
    void TestMonolayer()
    {
        /** The next line is needed because HoneycombMeshGenerator is not designed to be run in parallel */
        EXIT_IF_PARALLEL;

        /* The first thing we do is generate a nodes only mesh. To do this we first create a {{{MutableMesh}}}
         * to use as a generating mesh.
         * To do this we can use the {{{HoneycombMeshGenerator}}}. This generates a honeycomb-shaped mesh,
         * in which all nodes are equidistant. Here the first and second arguments
         * define the size of the mesh - we have chosen a mesh that is 2 nodes (i.e.
         * cells) wide, and 2 nodes high.
         */
        HoneycombMeshGenerator generator(2, 2);
        MutableMesh<2,2>* p_generating_mesh = generator.GetMesh();
        /* Once we have a {{{MutableMesh}}} we can generate a {{{NodesOnlyMesh}}} from it using the
         * following commands. Note you can also generate the {{{NodesOnlyMesh}}} from a collection of
         * nodes, see  [class:NodesOnlyMesh NodesOnlyMesh] for details.
         */
        NodesOnlyMesh<2> mesh;
        /* To run node-based simulations you need to define a cut off length (second argument in
         * {{{ConstructNodesWithoutMesh}}}), which defines the connectivity of the nodes by defining
         * a radius of interaction. */
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh, 1.5);

        /* Having created a mesh, we now create a {{{std::vector}}} of {{{CellPtr}}}s.
         * To do this, we the `CellsGenerator` helper class, which is templated over the type
         * of cell model required (here {{{UniformCellCycleModel}}})
         * and the dimension. We create an empty vector of cells and pass this into the
         * method along with the mesh. The second argument represents the size of that the vector
         * {{{cells}}} should become - one cell for each node, the third argument specifies
         * the proliferative type of the cell. */
        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_transit_type);

        /* Now we have a mesh and a set of cells to go with it, we can create a {{{CellPopulation}}}.
        * In general, this class associates a collection of cells with a mesh.
        * For this test, because we have a {{{NodesOnlyMesh}}}, we use a particular type of
        * cell population called a {{{NodeBasedCellPopulation}}}.
        */
        NodeBasedCellPopulation<2> cell_population(mesh, cells);

        /* We then pass in the cell population into an {{{OffLatticeSimulation}}},
         * and set the output directory, output multiple and end time. */
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("NodeBasedMonolayer");
        simulator.SetSamplingTimestepMultiple(12);
        simulator.SetEndTime(10.0);

        /* We now pass a force law to the simulation. */
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
     * EMPTYLINE
     *
     * To visualize the results, open a new terminal, {{{cd}}} to the Chaste directory,
     * then {{{cd}}} to {{{anim}}}. Then do: {{{java Visualize2dCentreCells /tmp/$USER/testoutput/NodeBasedMonolayer/results_from_time_0}}}.
     * we need to select the 'Cells as circles` option to be able to visualize the cells, as opposed
     * to just the centres.
     * We may have to do: {{{javac Visualize2dCentreCells.java}}} beforehand to create the
     * java executable.
     *
     * Alternatively, to view in Paraview, load the file {{{/tmp/$USER/testoutput/NodeBasedMonolayer/results_from_time_0/results.pvd}}}
     * and add glyphs to represent cells. An option is to use 3D spherical glyphs and then make a planar cut.
     * Note that, for larger simulations, you may need to unclick "Mask Points" (or similar) so as not to limit the number of glyphs
     * displayed by Paraview.
     *
     *
     *
     * EMPTYLINE
     *
     * == Test 2 - a basic node-based simulation in 3D ==
     *
     * EMPTYLINE
     *
     * In the second test we run a simple node-based simulation in 3D. This is very similar
     * to the 2D test with the dimension template (<2,2> and <2>) changed from 2 to 3 and instead of using a mesh
     * generator we generate the nodes directly.
     */
    void TestSpheroid()
    {
        /** The next line is needed because we cannot currently run node based simulations in parallel. */
        EXIT_IF_PARALLEL;

        /*
         * First, we generate a nodes only mesh. This time we specify the nodes manually by first
         * creating a vector of nodes. */
        std::vector<Node<3>*> nodes;
        /* We then create some nodes to add to this vector. */
        nodes.push_back(new Node<3>(0u,  false,  0.5, 0.0, 0.0));
        nodes.push_back(new Node<3>(1u,  false,  -0.5, 0.0, 0.0));
        nodes.push_back(new Node<3>(2u,  false,  0.0, 0.5, 0.0));
        nodes.push_back(new Node<3>(3u,  false,  0.0, -0.5, 0.0));
        /* Finally a {{{NodesOnlyMesh}}} is created and the vector of nodes is passed to
         * the {{{ConstructNodesWithoutMesh}}} method. */
        NodesOnlyMesh<3> mesh;
        /* To run node-based simulations you need to define a cut off length (second argument in
         * {{{ConstructNodesWithoutMesh}}}), which defines the connectivity of the nodes by defining
         * a radius of interaction. */
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);

        /*
         * Having created a mesh, we now create a {{{std::vector}}} of {{{CellPtr}}}s.
         * As before, we do this with the `CellsGenerator` helper class (this time with dimension 3).
         */
        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<UniformCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_transit_type);

        /* We make a {{{NodeBasedCellPopulation}}} (this time with dimension 3) as before.
         */
        NodeBasedCellPopulation<3> cell_population(mesh, cells);

        /* We then pass in the cell population into an {{{OffLatticeSimulation}}},
         * (this time with dimension 3) and set the output directory, output multiple and end time. */
        OffLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory("NodeBasedSpheroid");
        simulator.SetSamplingTimestepMultiple(12);
        simulator.SetEndTime(10.0);

        /* Again we create a force law (this time with dimension 3), and pass it to the {{{OffLatticeSimulation}}}.*/
        MAKE_PTR(GeneralisedLinearSpringForce<3>, p_force);
        simulator.AddForce(p_force);

        /* To run the simulation, we call {{{Solve()}}}. */
        simulator.Solve();

        /* The next two lines are for test purposes only and are not part of this tutorial.
         */
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 8u);
        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), 10.0, 1e-10);

        /* To avoid memory leaks, we conclude by deleting any pointers that we created in the test.*/
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

    /*
     * EMPTYLINE
     *
     * Note that you '''cannot view the results of a 3D simulation using the Java visualiser''' but
     * to visualize the results, use Paraview. See the UserTutorials/VisualizingWithParaview tutorial for more information.
     *
     * Load the file {{{/tmp/$USER/testoutput/NodeBasedSpheroid/results_from_time_0/results.pvd}}},
     * and add spherical glyphs to represent cells.
     *
     * EMPTYLINE
     *
     * == Test 3 - a node-based simulation on a restricted geometry ==
     *
     * EMPTYLINE
     *
     * In the third test we run a node-based simulation restricted to the surface of a sphere.
     */
    void TestOnSurfaceOfSphere()
    {
        /** The next line is needed because we cannot currently run node based simulations in parallel. */
        EXIT_IF_PARALLEL;

        /*
         * We begin with exactly the same code as the previous test: we create a cell population
         * from a mesh and vector of cells, and use this in turn to create
         * a simulation object.
         */

        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0u,  false,  0.5, 0.0, 0.0));
        nodes.push_back(new Node<3>(1u,  false,  -0.5, 0.0, 0.0));
        nodes.push_back(new Node<3>(2u,  false,  0.0, 0.5, 0.0));
        nodes.push_back(new Node<3>(3u,  false,  0.0, -0.5, 0.0));
        NodesOnlyMesh<3> mesh;
        /* To run node-based simulations you need to define a cut off length (second argument in
         * {{{ConstructNodesWithoutMesh}}}), which defines the connectivity of the nodes by defining
         * a radius of interaction. */
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);

        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<UniformCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_transit_type);

        NodeBasedCellPopulation<3> cell_population(mesh, cells);

        OffLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory("NodeBasedOnSphere");
        simulator.SetSamplingTimestepMultiple(12);
        simulator.SetEndTime(10.0);

        /* As before, we create a linear spring force and pass it to the simulation object. */
        MAKE_PTR(GeneralisedLinearSpringForce<3>, p_force);
        simulator.AddForce(p_force);

        /*
         * This time we create a {{{CellPopulationBoundaryCondition}}} and pass this to
         * the {{{OffLatticeSimulation}}}. Here we use a {{{SphereGeometryBoundaryCondition}}}
         * which restricts cells to lie on a sphere (in 3D) or circle (in 2D).
         *
         * For a list of possible boundary conditions see subclasses of {{{AbstractCellPopulationBoundaryCondition}}}.
         * These can be found in the inheritance diagram, here, [class:AbstractCellPopulationBoundaryCondition AbstractCellPopulationBoundaryCondition].
         * Note that some of these boundary conditions are not compatible with node-based simulations see the specific class documentation for details,
         * if you try to use an incompatible class then you will receive a warning.
         *
         * First we set the centre (0,0,1) and radius of the sphere (1).
         */
        c_vector<double,3> centre = zero_vector<double>(3);
        centre(2) = 1.0;
        double radius = 1.0;
        /* We then make a pointer to the boundary condition using the MAKE_PTR_ARGS macro, and pass
         * it to the {{{OffLatticeSimulation}}}. */
        MAKE_PTR_ARGS(SphereGeometryBoundaryCondition<3>, p_boundary_condition, (&cell_population, centre, radius));
        simulator.AddCellPopulationBoundaryCondition(p_boundary_condition);

        /* To run the simulation, we call {{{Solve()}}}. */
        simulator.Solve();

        /* The next two lines are for test purposes only and are not part of this tutorial.
         */
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 8u);
        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), 10.0, 1e-10);

        /* To avoid memory leaks, we conclude by deleting any pointers that we created in the test.*/
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }
    /*
     * EMPTYLINE
     *
     * To visualize the results, use Paraview. See the UserTutorials/VisualizingWithParaview tutorial for more information.
     *
     * Load the file {{{/tmp/$USER/testoutput/NodeBasedOnSphere/results_from_time_0/results.pvd}}},
     * and add spherical glyphs to represent cells.
     *
     * EMPTYLINE
     */
};

#endif /* TESTRUNNINGNODEBASEDSIMULATIONSTUTORIAL_HPP_ */
