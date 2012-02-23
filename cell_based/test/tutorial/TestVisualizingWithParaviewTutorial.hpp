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
 * In this tutorial we show how Chaste is used to generate simulations
 * that can be viewed in Paraview, and how to use Paraview itself. Two examples
 * are provided - one using a cell-centre based model, and the second using
 * a vertex model. To be able to view these simulations, we must first have
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
#include "StochasticDurationCellCycleModel.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "CellsGenerator.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "NagaiHondaForce.hpp"
#include "OffLatticeSimulation.hpp"
#include "SmartPointers.hpp"

/* Next, we define the test class, which inherits from {{{AbstractCellBasedTestSuite}}}
 * and defines some test methods.
 */
class TestVisualizingWithParaviewTutorial : public AbstractCellBasedTestSuite
{
public:
    /* EMPTYLINE
     *
     * == Test 1 - a cell-based monolayer simulation ==
     *
     * EMPTYLINE
     *
     * In the first test, we run a simple cell-based simulation, in which we use
     * a honeycomb mesh with ghost nodes, and give each cell a stochastic cell-cycle model.
     */
    void Test2DMeshBasedMonolayerSimulationForVisualizing() throw (Exception)
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
        CellsGenerator<StochasticDurationCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, location_indices.size(), TRANSIT);

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
        cell_population.SetOutputVoronoiData(true);

        /* We then pass in the cell population into an {{{OffLatticeSimulation}}},
         * and set the output directory and end time. */
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("Test2DMonolayerSimulationForVisualizing");
        simulator.SetEndTime(1.0);

        /* We create a force law and pass it to the {{{OffLatticeSimulation}}}. */
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetCutOffLength(1.5);
        simulator.AddForce(p_linear_force);

        /* To run the simulation, we call {{{Solve()}}}. */
        simulator.Solve();
    }

    /*
    * EMPTYLINE
    *
    * To visualize the results, we must first open Paraview. We open the folder containing our test output using the 'file' menu at
    * the top. The output will be located in {{{/tmp/$USER/testoutput/Test2DMonolayerSimulationForVisualizing/results_from_time_0}}}.
    * There will be a .vtu file generated for every timestep, which must all be opened at once to view the simulation. To do this,
    * simply select {{{results_..vtu}}}. We should now see {{{results_*}}} in the pipeline browser. We click {{{Apply}}} in the properties tab
    * of the object inspector, and we should now see a visualization in the right hand window.
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
    * == Test 2 - a basic vertex-based simulation ==
    *
    * EMPTYLINE
    *
    * Here, we run a simple vertex-based simulation, in which we create a monolayer
    * of cells using a mutable vertex mesh. Each cell is assigned a fixed cell-cycle model.
    */
    void Test2DVertexBasedMonolayerSimulationForVisualizing() throw(Exception)
    {
        /* In this test, we create a vertex-based cell population in which cells are defined
         * by their vertices, and cell proliferation is governed by a fixed generation-based
         * cell-cycle model (with differentiation after a default number of generations).
         */
        HoneycombVertexMeshGenerator generator(6, 9);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
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

        /* To run the simulation, we call {{{Solve()}}}. */
        simulator.Solve();
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
