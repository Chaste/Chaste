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

#ifndef TESTVISUALIZINGWITHPARAVIEWTUTORIAL_HPP_
#define TESTVISUALIZINGWITHPARAVIEWTUTORIAL_HPP_

/*
 * = Examples showing how to visualize simulations in Paraview =
 *
 * EMPTYLINE
 *
 * == Introduction ==
 *
 * EMPTYLINE
 *
 * In this tutorial we show how Chaste can be used to generate simulations
 * that can be viewed in Paraview, and how to use Paraview itself. Four examples
 * are provided: the first two use a `MeshBasedCellPopulation` (a cell-centre model
 * based on a Delaunay triangulation description of cell neighbours);
 * the third uses a `NodeBasedCellPopulation` (a cell-centre model based on an
 * 'overlapping spheres' description of cell neighbours); and the fourth uses a
 * `VertexBasedCellPopulation` (in which each cell is represented by a polygon).
 * To be able to view these simulations, we must first have
 * downloaded and installed VTK and Paraview, and updated our hostconfig file
 * to ensure that it knows to use VTK.
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

/* The remaining header files define classes that will be used in the cell population
 * simulation test. We have encountered each of these header files in previous cell-based
 * Chaste tutorials. */
#include "UniformCellCycleModel.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "CylindricalHoneycombMeshGenerator.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "CellsGenerator.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "NagaiHondaForce.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "OffLatticeSimulation.hpp"
#include "TransitCellProliferativeType.hpp"
#include "SmartPointers.hpp"
#include "VoronoiDataWriter.hpp"

#include "FakePetscSetup.hpp"

/* Next, we define the test class, which inherits from {{{AbstractCellBasedTestSuite}}}
 * and defines some test methods.
 */
class TestVisualizingWithParaviewTutorial : public AbstractCellBasedTestSuite
{
public:
    /* EMPTYLINE
     *
     * == Test 1 - a mesh-based cell centre monolayer simulation ==
     *
     * EMPTYLINE
     *
     * In the first test, we run a simple cell-based simulation using a `MeshBasedCellPopulation`,
     * in which we use
     * a honeycomb mesh with ghost nodes, and give each cell a stochastic cell-cycle model.
     */
    void Test2DMeshBasedMonolayerSimulationForVisualizing()
    {
        /* In a similar way to previous cell-based Chaste tutorials,
         * we create a mesh-based cell population in which cells are defined by their centres,
         * and cell proliferation is governed by a stochastic generation-based cell-cycle model
         * with no differentiation.
         */
        HoneycombMeshGenerator generator(10, 10, 2);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, location_indices.size(), p_transit_type);

        MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);

        /*
         * The following line tells the cell population to write data to .vtu files with cells
         * not as points, but as polytopes. This is the default setting: we include the call
         * here to highlight this option. If writing point data, we may choose the shape used
         * to visualize each cell in Paraview using glyphs.
         */
        cell_population.SetWriteVtkAsPoints(false);

        /* In order to output the .vtu files required for Paraview, we explicitly
         * instruct the simulation to output the data we need.
         */
        cell_population.AddPopulationWriter<VoronoiDataWriter>();

        /* We then pass in the cell population into an {{{OffLatticeSimulation}}},
         * and set the output directory and end time. */
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("Test2DMeshBasedMonolayerSimulationForVisualizing");
        simulator.SetEndTime(1.0);

        /* We create a force law and pass it to the {{{OffLatticeSimulation}}}. */
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetCutOffLength(1.5);
        simulator.AddForce(p_linear_force);

        /* To run the simulation, we call {{{Solve()}}}. */
        simulator.Solve();

        /* The next two lines are for test purposes only and are not part of this tutorial.
         * We are checking that we reached the end time of the simulation
         * with the correct number of cells. If different simulation input parameters are being explored
         * the lines should be removed.
         */
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 108u);
        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), 1.0, 1e-10);
    }

    /*
    * EMPTYLINE
    *
    * To visualize the results, we must first open Paraview. We open the folder containing our test output using the 'file' menu at
    * the top. The output will be located in {{{/tmp/$USER/testoutput/Test2DMeshBasedMonolayerSimulationForVisualizing/results_from_time_0}}}.
    * There will be a .vtu file generated for every timestep, which must all be opened at once to view the simulation. To do this,
    * simply select {{{results.pvd}}}. We should now see {{{results.pvd}}}  in the pipeline browser. We click {{{Apply}}} in the properties tab
    * of the object inspector, and we should now see a visualization in the right hand window.  (An alternative to opening the {{{results.pvd}}}
    * file is to open all the time steps en masse where we open {{{results_..vtu}}} and see {{{results_*}}} appear in the pipeline browser.)
    *
    * At this stage, it will be necessary to refine how we wish to view this particular visualisation. The viewing styles can be edited using
    * the display tab of the object inspector. In particular, under {{{Style}}}, the representation drop down menu allows us to view
    * the cells as a surface with edges, or as simply a wireframe. It is advisable at this point to familiarize ourselves with the different
    * viewing options, colour and size settings.
    *
    * At this stage, the viewer is showing all cells in the simulation, including the ghost nodes. In order to view only real cells, we must
    * apply a threshold. This is achieved using the threshold button on the third toolbar (the icon is a cube with a green 'T' inside). Once you
    * click the threshold button, you will see a new threshold appear below your results in the pipeline browser. Go to the properties tab and
    * reset the lower threshold to be less than 0, and the upper threshold to be between 0 and 1, ensuring that the 'Non-ghosts' option is
    * selected in the 'Scalars' drop down menu. Once we have edited this, we click apply (we may need to click it twice), and the visualisation on the
    * right window will have changed to eliminate ghost nodes.
    *
    * To view the simulation, simply use the animation buttons located on the top toolbar. We can also save a screenshot, or an animation, using
    * the appropriate options from the file menu. Next to the threshold button are two other useful options, 'slice' and 'clip', but these will
    * only be applicable for 3D visualisations.
    *
    * EMPTYLINE
    *
    * == Test 2 - a periodic mesh-based cell centre monolayer simulation ==
     *
     * EMPTYLINE
     *
     * In the second test, similar to the first test, we run a simple cell-based simulation using a `MeshBasedCellPopulation`,
     * in which we use
     * a honeycomb mesh with ghost nodes, and give each cell a stochastic cell-cycle model. However here we impose periodic boundaries.
     * The only difference in this test is the generation of the mesh
     */
    void Test2DPeriodicMeshBasedMonolayerSimulationForVisualizing()
    {
        /*
         * We setup the simulation in the same way as above but
         * here we use a cylindrical mesh as we wish to enforce periodicity
         * in the x direction.
         */
        CylindricalHoneycombMeshGenerator generator(10, 10, 2);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, location_indices.size(), p_transit_type);

        MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);


        cell_population.SetWriteVtkAsPoints(false);
        cell_population.AddPopulationWriter<VoronoiDataWriter>();

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("Test2DPeriodicMeshBasedMonolayerSimulationForVisualizing");
        simulator.SetEndTime(1.0);

        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetCutOffLength(1.5);
        simulator.AddForce(p_linear_force);

        simulator.Solve();

        /* The next two lines are for test purposes only and are not part of this tutorial.
         * We are checking that we reached the end time of the simulation
         * with the correct number of cells. If different simulation input parameters are being explored
         * the lines should be removed.
         */
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 108u);
        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), 1.0, 1e-10);
    }

    /*
    * EMPTYLINE
    *
    * To visualize the results, we follow the instructions above for the first simulation, ensuring that we open the
    * test output from the new folder, {{{Test2DPeriodicMeshBasedMonolayerSimulationForVisualizing}}}. You will see that the left an righ sides
    * of the monolayer are the same. As in the first test you can threshold to remove the ghost nodes.
    *
    * EMPTYLINE
    *
    * == Test 3 - a node-based simulation ==
    *
     * EMPTYLINE
     *
    * We next run a similar simulation to the first two examples, but now use a `NodeBasedCellPopulation`,
    * in which cells are represented as 'overlapping spheres'.
    */
    void Test2DNodeBasedMonolayerSimulationForVisualizing()
    {
        /* We set up the simulation in much the same way as above, except now using a `NodesOnlyMesh` and
         * `NodeBasedCellPopulation`. Further details on how to set up a node-based simulation can be found in
         * UserTutorials/RunningNodeBasedSimulations.
         */
        HoneycombMeshGenerator generator(10, 10, 0);
        TetrahedralMesh<2,2>* p_generating_mesh = generator.GetMesh();

        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh, 1.5);

        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes());

        NodeBasedCellPopulation<2> cell_population(mesh, cells);

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("Test2DNodeBasedMonolayerSimulationForVisualizing");
        simulator.SetEndTime(1.0);

        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetCutOffLength(1.5);
        simulator.AddForce(p_linear_force);

        simulator.Solve();

        /* The next two lines are for test purposes only and are not part of this tutorial.
         * We are checking that we reached the end time of the simulation
         * with the correct number of cells. */
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 108u);
        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), 1.0, 1e-10);
    }

    /*
    * EMPTYLINE
    *
    * To visualize the results, we follow the instructions above for the first simulation, ensuring that we open the
    * test output from the new folder, {{{Test2DNodeBasedMonolayerSimulationForVisualizing}}}.
    * After opening Paraview, load the file {{{results.pvd}}}, then click "Apply" in the object inspector panel.
    * As this simulation uses a `NodeBasedCellPopulation`, you must use glyphs to visualize cells: click the button
    * marked "Glyph" in the toolbar of common filters; specify cells to be displayed as spheres; then click "Apply".
    *
     * Note that, for larger simulations, you may need to unclick "Mask Points" (or similar) so as not to limit the number of glyphs
     * displayed by Paraview.
    *
    * EMPTYLINE
    *
    * == Test 4 - a basic vertex-based simulation ==
    *
     * EMPTYLINE
     *
    * Here, we run a simple vertex-based simulation, in which we create a monolayer
    * of cells using a mutable vertex mesh. Each cell is assigned a fixed cell-cycle model.
    */
    void Test2DVertexBasedMonolayerSimulationForVisualizing()
    {
        /* In this test, we create a vertex-based cell population in which cells are defined
         * by their vertices, and cell proliferation is governed by a fixed generation-based
         * cell-cycle model (with differentiation after a default number of generations).
         */
        HoneycombVertexMeshGenerator generator(6, 9);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements());

        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        /* We then pass in the cell population into an {{{OffLatticeSimulation}}},
         * and set the output directory and end time. */
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("Test2DVertexMonolayerSimulationForVisualizing");
        simulator.SetEndTime(0.1);

        /* We create a force law and pass it to the {{{OffLatticeSimulation}}}. */
        MAKE_PTR(NagaiHondaForce<2>, p_nagai_honda_force);
        simulator.AddForce(p_nagai_honda_force);

        /* We also make a pointer to target area modifier and add it to the simulator.
         * The target area modifier assigns target areas to cells throughout the simulation.
         */
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        /* To run the simulation, we call {{{Solve()}}}. */
        simulator.Solve();

        /* The next two lines are for test purposes only and are not part of this tutorial.
         * We are checking that we reached the end time of the simulation
         * with the correct number of cells. */
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 84u);
        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), 0.1, 1e-10);
    }
    /*
    * EMPTYLINE
    *
    * To visualize the results, we follow the instructions above for the first simulation, ensuring that we open the
    * test output from the new folder, {{{Test2DVertexMonolayerSimulationForVisualizing}}}.
    *
    */
};

#endif /* TESTVISUALIZINGWITHPARAVIEWTUTORIAL_HPP_ */
