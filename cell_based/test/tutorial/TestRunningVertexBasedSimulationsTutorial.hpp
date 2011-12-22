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
/*
 *
 *  Chaste tutorial - this page gets automatically changed to a wiki page
 *  DO NOT remove the comments below, and if the code has to be changed in
 *  order to run, please check the comments are still accurate
 *
 *
 */

#ifndef TESTRUNNINGVERTEXBASEDSIMULATIONSTUTORIAL_HPP_
#define TESTRUNNINGVERTEXBASEDSIMULATIONSTUTORIAL_HPP_

/*
 * = Examples showing how to create, run and visualize vertex-based simulations =
 *
 * EMPTYLINE
 *
 * == Introduction ==
 *
 * EMPTYLINE
 *
 * In this tutorial we show how Chaste can be used to create, run and visualize vertex-based simulations.
 * Full details of the mechanical model proposed by T. Nagai and H. Honda ("A dynamic cell model for
 * the formation of epithelial tissues", Philosophical Magazine Part B 81:699-719).
 *
 * EMPTYLINE
 *
 * == The test ==
 *
 * EMPTYLINE
 *
 * As in previous cell-based Chaste tutorials, we begin by including the necessary header files.
 */
#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"

/* The remaining header files define classes that will be used in the cell-based
 * simulation. We have encountered some of these header files in previous cell-based
 * Chaste tutorials. */
#include "CellsGenerator.hpp"
#include "OffLatticeSimulation.hpp"
#include "SmartPointers.hpp"
/* The next header file defines the cell cycle model. */
#include "StochasticDurationCellCycleModel.hpp"
/* The next two header files define a helper class for generating suitable meshes: one planar and one periodic. */
#include "HoneycombVertexMeshGenerator.hpp"
#include "CylindricalHoneycombVertexMeshGenerator.hpp"
/* The next header file defines a vertex-based {{{CellPopulation}}} class.*/
#include "VertexBasedCellPopulation.hpp"
/* The next header file defines a force law for describing the mechanical interactions
 * between neighbouring cells in the cell population, subject to each vertex.
 */
#include "NagaiHondaForce.hpp"
/* The next header file defines a boundary condition for the cells.*/
#include "PlaneBoundaryCondition.hpp"
/* The next header file defines a cell killer, which specifies how cells are removed from the simulation.*/
#include "PlaneBasedCellKiller.hpp"

/* Next, we define the test class, which inherits from {{{AbstractCellBasedTestSuite}}}
 * and defines some test methods.
 */
class TestRunningVertexBasedSimulationsTutorial : public AbstractCellBasedTestSuite
{
public:
    /* EMPTYLINE
    *
    * == Test 1 - a basic vertex-based simulation ==
    *
    * EMPTYLINE
    *
    * In the first test, we run a simple vertex-based simulation, in which we create a monolayer
    * of cells, using a mutable vertex mesh. Each cell is assigned a stochastic cell-cycle model.
    */
    void TestMonolayer() throw(Exception)
    {
        /* First, we generate a vertex mesh. To create a {{{MutableVertexMesh}}}, we can use
        * the {{{HoneycombVertexMeshGenerator}}}. This generates a honeycomb-shaped mesh,
        * in which all nodes are equidistant. Here the first and second arguments
        * define the size of the mesh - we have chosen a mesh that is 2 elements (i.e.
        * cells) wide, and 2 elements high.
        */
        HoneycombVertexMeshGenerator generator(2, 2);    // Parameters are: cells across, cells up
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        /* Having created a mesh, we now create a {{{std::vector}}} of {{{CellPtr}}}s.
        * To do this, we use the `CellsGenerator` helper class, which is templated over the type
        * of cell model required (here {{{StochasticDurationCellCycleModel}}})
        * and the dimension. We create an empty vector of cells and pass this into the
        * method along with the mesh. The second argument represents the size of that the vector
        * {{{cells}}} should become - one cell for each element, the third argument specifies
        * the proliferative type of the cell STEM, TRANSIT or DIFFERENTIATED. */
        std::vector<CellPtr> cells;
        CellsGenerator<StochasticDurationCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(),TRANSIT);

        /* Now we have a mesh and a set of cells to go with it, we can create a {{{CellPopulation}}}.
        * In general, this class associates a collection of cells with a mesh.
        * For this test, because we have a {{{MutableVertexMesh}}}, we use a particular type of
        * cell population called a {{{VertexBasedCellPopulation}}}.
        */
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        /* We then pass the cell population into an {{{OffLatticeSimulation}}},
         * and set the output directory and end time. */
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("VertexBasedMonolayer");
        simulator.SetEndTime(1.0);

        /*
         * For longer simulations, we may not want to output the results
         * every time step. In this case we can use the following method,
         * to print results every 50 time steps instead. As the default time step
         * used by the simulator (for vertex based simulations), is 0.02 hours, this method will cause the
         * simulator to print results every 6 minutes (i.e. 0.1 hours).
         */
        simulator.SetSamplingTimestepMultiple(50);

        /* We must now create one or more force laws, which determine the mechanics of the vertices
        * of each cell in a cell population. For this test, we use one force law, based on the
        * Nagai-Honda mechanics, and pass it to the {{{OffLatticeSimulation}}}.
        * For a list of possible forces see subclasses of {{{AbstractForce}}}.
        * These can be found in the inheritance diagram, here, [class:AbstractForce AbstractForce].
        * Note that some of these forces are not compatible with vertex-based simulations see the specific class documentation for details,
        * if you try to use an incompatible class then you will receive a warning.
        */
        MAKE_PTR(NagaiHondaForce<2>, p_force);
        simulator.AddForce(p_force);

        /* To run the simulation, we call {{{Solve()}}}. */
        simulator.Solve();
    }

    /*
    * EMPTYLINE
    *
    * To visualize the results, open a new terminal, {{{cd}}} to the Chaste directory,
    * then {{{cd}}} to {{{anim}}}. Then do: {{{java Visualize2dVertexCells /tmp/$USER/testoutput/VertexBasedMonolayer/results_from_time_0}}}.
    * We may have to do: {{{javac Visualize2dVertexCells.java}}} beforehand to create the
    * java executable.
    *
    * EMPTYLINE
    *
    * == Test 2 - introducing periodicity, boundaries and cell killers ==
    *
    * EMPTYLINE
    *
    * In the second test, we run a simple vertex-based simulation, in which we create a monolayer
    * of cells in a periodic geometry, using a cylindrical vertex mesh. We also include a fixed
    * boundary which cells can't pass through and a cell killer which removes cells once they leave
    * a region. As before each cell is assigned a stochastic cell-cycle model.
    */
    void TestPeriodicMonolayer() throw(Exception)
    {
        /* First, we generate a periodic vertex mesh. To create a {{{Cylindrical2dVertexMesh}}}, we can use
         * the {{{CylindricalHoneycombVertexMeshGenerator}}}. This generates a honeycomb-shaped mesh,
         * in which all nodes are equidistant and the right hand side is associated with the left hand side.
         * Here the first and second arguments define the size of the mesh - we have chosen a mesh that
         * is 4 elements (i.e. cells) wide, and 4 elements high.
         */
        CylindricalHoneycombVertexMeshGenerator generator(4, 4);    // Parameters are: cells across, cells up
        Cylindrical2dVertexMesh* p_mesh = generator.GetCylindricalMesh();

        /* Having created a mesh, we now create a {{{std::vector}}} of {{{CellPtr}}}s.
        * This is exactly the same as the above test. */
        std::vector<CellPtr> cells;
        CellsGenerator<StochasticDurationCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(),TRANSIT);

        /* Now we have a mesh and a set of cells to go with it, we can create a {{{CellPopulation}}}.
         * This is also the same as in the above test.
         */
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        /* As always we then pass the cell population into an {{{OffLatticeSimulation}}},
         * and set the output directory, output multiple and end time. */
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("VertexBasedPeriodicMonolayer");
        simulator.SetSamplingTimestepMultiple(50);
        simulator.SetEndTime(1.0);

        /* We now make a pointer to an appropriate force and pass it to the
         * {{{OffLatticeSimulation}}}.
         */
        MAKE_PTR(NagaiHondaForce<2>, p_force);
        simulator.AddForce(p_force);

        /* We now create one or more {{{CellPopulationBoundaryCondition}}}s, which determine
         * any conditions which each cell in a cell population must satisfy. For this test,
         * we use a {{{PlaneBoundaryCondition}}}, and pass it to the {{{OffLatticeSimulation}}}.
         * For a list of possible boundary condition see subclasses of {{{AbstractCellPopulationBoundaryCondition}}}.
         * These can be found in the inheritance diagram, here,
         * [class:AbstractCellPopulationBoundaryCondition AbstractCellPopulationBoundaryCondition].
         * Note that some of these boundary conditions are not compatible with vertex-based
         * simulations see the specific class documentation for details, if you try to use an
         * incompatible class then you will receive a warning.
         *
         * The first step is to define a point on the plane boundary and a normal to the plane.
         */
        c_vector<double,2> point = zero_vector<double>(2);
        c_vector<double,2> normal = zero_vector<double>(2);
        normal(1) = -1.0;
        /* We can now make a pointer to a {{{PlaneBoundaryCondition}}} (passing the point
         * and normal to the plane) and pass it to the {{{OffLatticeSimulation}}}.*/
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc, (&cell_population, point, normal));
        simulator.AddCellPopulationBoundaryCondition(p_bc);

        /* We now create one or more {{{CellKiller}}}s, which determine how cells are removed
         * from the simulation. For this test, we use a {{{PlaneBasedCellKiller}}}, and pass
         * it to the {{{OffLatticeSimulation}}}. For a list of possible cell killers see subclasses
         * of {{{AbstractCellKiller}}}. These can be found in the inheritance diagram, here,
         * [class:AbstractCellKiller AbstractCellKiller].
         *
         * The first step is to define a point on the plane boundary and a normal to the plane.
         * We reuse the point and normal from the {{{PlaneBoundaryCondition}}}.
         */
        point(1) = 3.0;
        normal(1) = 1.0;
        /* Finally we now make a pointer to a {{{PlaneBasedCellKiller}}} (passing the point
         * and normal to the plane) and pass it to the {{{OffLatticeSimulation}}}.*/
        MAKE_PTR_ARGS(PlaneBasedCellKiller<2>, p_killer, (&cell_population, point, normal));
        simulator.AddCellKiller(p_killer);

        /* To run the simulation, we call {{{Solve()}}}. */
        simulator.Solve();
    }
    /*
     * EMPTYLINE
     *
     * To visualize the results, open a new terminal, {{{cd}}} to the Chaste directory,
     * then {{{cd}}} to {{{anim}}}. Then do: {{{java Visualize2dVertexCells /tmp/$USER/testoutput/VertexBasedPeriodicMonolayer/results_from_time_0}}}.
     *
     * You should see that the edges of the mesh are identical on both sides; cells no
     * longer pass through the line y=0; and cells are removed at y=3.
     *
     * EMPTYLINE
     */
};

#endif /* TESTRUNNINGVERTEXBASEDSIMULATIONSTUTORIAL_HPP_ */
