/*

Copyright (c) 2005-2024, University of Oxford.
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

#ifndef TESTCOMMANDLINEVERTEXMESHTUTORIAL_HPP_
#define TESTCOMMANDLINEVERTEXMESHTUTORIAL_HPP_

/*
 * ## Examples showing how to create, run and visualize vertex-based simulations using command line arguments.
 *
 * ### Introduction
 *
 * In this tutorial we show how Chaste can be used to create, run and visualize vertex-based simulations with variables from command line arguments.
 * Full details of the mechanical model proposed by T. Nagai and H. Honda, 2000, "A dynamic cell model for
 * the formation of epithelial tissues", Philosophical Magazine Part B 81:699-719, doi:[10.1103/PhysRevLett.69.2013](https://doi.org/10.1103/PhysRevLett.69.2013).
 *
 * ### The test
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
#include "TransitCellProliferativeType.hpp"
#include "SmartPointers.hpp"
/* The next header file defines the cell cycle model. */
#include "UniformG1GenerationalCellCycleModel.hpp"
/* The next two header files define a helper class for generating suitable meshes: one planar and one periodic. */
#include "HoneycombVertexMeshGenerator.hpp"
#include "CylindricalHoneycombVertexMeshGenerator.hpp"
/* The next header file defines a vertex-based `CellPopulation` class.*/
#include "VertexBasedCellPopulation.hpp"
/* The next header file defines a force law for describing the mechanical interactions
 * between neighbouring cells in the cell population, subject to each vertex.
 */
#include "NagaiHondaForce.hpp"
/* This force law assumes that cells possess a "target area" property which determines the size of each
 * cell in the simulation. In order to assign target areas to cells and update them in each time step, we need
 * the next header file.
 */
#include "SimpleTargetAreaModifier.hpp"
/* The next header file defines a boundary condition for the cells.*/
#include "PlaneBoundaryCondition.hpp"
/* The next header file defines a cell killer, which specifies how cells are removed from the simulation.*/
#include "PlaneBasedCellKiller.hpp"

/* The next header file enforces that we are running this test only on one process. */
#include "FakePetscSetup.hpp"
/* Finally, we include our header file that allows for taking command line arguments into our test.*/
#include "CommandLineArguments.hpp"

/* NOTE: This tutorial utilises command line arguments which are implemented in a BASH script. 
 * If this test is called utilising the standard method previously used such as ctest and not
 * provided the correct command line arguments it will fail. The corresponding BASH script can be 
 * found at Chaste/cell_based/test/runcommandlinemeshtutorials.sh .*/

/* Next, we define the test class, which inherits from `AbstractCellBasedTestSuite`
 * and defines some test methods.
 */
class TestCommandLineVertexMeshTutorial : public AbstractCellBasedTestSuite
{
public:
    /*
    * Test  - a basic vertex-based simulation which takes in command line arguements for the height
    * and width of the tissue from a bash script.  
    * This test for the most part operates the same as TestRunningVertexBasedSimulationsTutorial
    * with slight additions for command line arguements. Additionaly, as our bash script executes a for loop to 
    * run multiple simulations, this tutorial will output a total of 9 meshes of varying height and widths with
    * correspodning paraview files.
    */
    void TestCommandLineMeshTutorial1()
    {
        /* First, we generate a vertex mesh. To create a `MutableVertexMesh`, we can use
        * the `HoneycombVertexMeshGenerator`. This generates a honeycomb-shaped mesh,
        * in which all nodes are equidistant. Here the first and second arguments
        * define the size of the mesh - in this test we utilise our command line arguements
        * as our tissues height and width. Thus, if outp = 2 our mesh will have a width and height 
        * of 2 cells for a total of 4 cells.
        */

        // We take in the command line arguement defined in the accompanying bash script (runcommandlinemeshtutorial.sh).
        // In this tutorial we take in our command lined argument as an unsigned interger. 
        unsigned outp = CommandLineArguments::Instance()->GetUnsignedCorrespondingToOption("-opt");
        // We also need to make sure that our outp variable has taken in a value from the command line.
        // Additionally, as this arguement will be defining our tissue height and width it needs to be a positive number 
        // greater than 1.
        TS_ASSERT_LESS_THAN(1, outp);

        // We define our mesh width and height based on the provided command line argument.
        HoneycombVertexMeshGenerator generator(outp, outp);    // Parameters are: cells across, cells up
        boost::shared_ptr<MutableVertexMesh<2,2> > p_mesh = generator.GetMesh();

        /* Having created a mesh, we now create a `std::vector` of `CellPtr`s.
        * To do this, we use the `CellsGenerator` helper class, which is templated over the type
        * of cell model required (here `UniformG1GenerationalCellCycleModel`)
        * and the dimension. We create an empty vector of cells and pass this into the
        * method along with the mesh. The second argument represents the size of that the vector
        * `cells` should become - one cell for each element, the third argument specifies
        * the proliferative type of the cell. */
        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<UniformG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_transit_type);

        /* Now we have a mesh and a set of cells to go with it, we can create a `CellPopulation`.
        * In general, this class associates a collection of cells with a mesh.
        * For this test, because we have a `MutableVertexMesh`, we use a particular type of
        * cell population called a `VertexBasedCellPopulation`.
        */
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        /* We then pass the cell population into an `OffLatticeSimulation`,
         * and set the output directory and end time. */
        OffLatticeSimulation<2> simulator(cell_population);
        // As we will be running multiple simulations, the results need to be saved to different locations.
        // Here we choose to utilise the height and width of the tissue as our title of each output file.
        simulator.SetOutputDirectory("CommandLineMeshTutorial/tissue_height_and_width_"+std::to_string(outp));
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
        * Nagai-Honda mechanics, and pass it to the `OffLatticeSimulation`.
        * For a list of possible forces see subclasses of `AbstractForce`.
        * These can be found in the inheritance diagram, here,
        * [AbstractForce](/doxygen-latest/classAbstractForce).
        * Note that some of these forces are not compatible with vertex-based simulations see the specific class documentation for details,
        * if you try to use an incompatible class then you will receive a warning.
        */
        MAKE_PTR(NagaiHondaForce<2>, p_force);
        simulator.AddForce(p_force);

        /* A `NagaiHondaForce` assumes that each cell has a target area. The target areas of cells are used to determine pressure
         * forces on each vertex and eventually determine the size of each cell in the simulation. In order to assign target areas to cells
         * and update them in each time step to model growth, we add a `SimpleTargetAreaModifier` to the simulation, which inherits from
         *  `AbstractTargetAreaModifier`.
         */
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        /* To run the simulation, we call `Solve()`. */
        simulator.Solve();

        /* The next two lines are for test purposes only and are not part of this tutorial. If different simulation input parameters are being explored
         * the lines should be removed.*/
        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), 1.0, 1e-10);
    }
};
/*
 *
 * To visualize the results, open paraview and click File->Open, navigate to /build/testoutput/CommandLineMeshTutorial/,
 * you should now see several "tissue_height_and_width_" folders. Open one of these folders and select the vtu file to open.
 *
 * You should now see a mesh of cells of height and width equal to the number after tissue_height_and_width_. To make visualisation
 * easier find the drop down menu with "Surface" currently selected and change this to "Surface and edges".
 */

#endif /* TESTCOMMANDLINEMESHTUTORIAL_HPP_ */
