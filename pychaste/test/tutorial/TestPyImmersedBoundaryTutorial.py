"""Copyright (c) 2005-2024, University of Oxford.
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
"""

#ifndef
#define TRIGGER_WIKI

## ## Introduction
## This tutorial is a demonstration of the immersed boundary method, a technique
## for simulating fluid-structure interactions. We can use the immersed boundary
## method to simulate a cell as a structure with its outer **boundary immersed**
## in a fluid. There is a two-way coupling between the fluid and the structure:
## the flow of the fluid exerts a force on the structure, and the structure
## influences the flow of the fluid.
##
## In this tutorial, we demonstrate:
## 1. Building single-cell immersed boundary capable simulations.
## 2. Building multi-cellular immersed boundary simulations.
## 3. Adding and manipulating immersed boundary fluid sources.
##
## ## Imports

import unittest

import chaste

chaste.init()  # setup MPI

from chaste.cell_based import (
    AbstractCellBasedTestSuite,
    CellsGeneratorUniformCellCycleModel_2,
    DifferentiatedCellProliferativeType,
    ForwardEulerNumericalMethod2_2,
    ImmersedBoundaryCellPopulation2,
    ImmersedBoundaryLinearInteractionForce2,
    ImmersedBoundaryLinearMembraneForce2,
    ImmersedBoundarySimulationModifier2,
    OffLatticeSimulation2_2,
    SetupNotebookTest,
    SimulationTime,
    TearDownNotebookTest,
)

from chaste.mesh import FluidSource2, ImmersedBoundaryPalisadeMeshGenerator

from chaste.visualization import (
    JupyterNotebookManager,
    JupyterSceneModifier2,
    VtkScene2
)

class TestImmersedBoundaryTutorial(AbstractCellBasedTestSuite):

    ## ### 1. Simple Immersed Boundary Simulations
    ## We begin by exploring simulations containing a single cell. This will
    ## familiarise you with how to generate immersed boundary cells, the steps
    ## involved in setting up an immersed boundary simulation, and the options
    ## available for controlling how the cells are generated and behave.
    ##

    ## Immersed boundary simulations operate over a square domain, with `x` and `y`
    ## coordinates lying in the range `0` to `1`. The domain wraps on both axes -
    ## this means that if a cell moves off the right hand edge of the domain,
    ## the segment will appear on the left hand side. This is not purely visual;
    ## forces are also transmitted across these boundaries.
    ##

    ## **Tip** Make sure all your coordinates are between `0` and `1`.

    def test_simple_immersed_boundary_simulation(self):

        ## Setup the simulation environment in the notebook

        # JUPYTER_SETUP

        ## Next, we define the necessary geometry by generating a mesh to 
        ## contain a single cell.

        gen = ImmersedBoundaryPalisadeMeshGenerator(1, 128, 0.1, 2.0, 0.0, False)
        mesh = gen.GetMesh()

        ## The first line of code defines an `ImmersedBoundaryPalisadeMeshGenerator`
        ## called `gen`. The 3rd parameter controls the exponent of the superellipse(`0.1`) 
        ## and the 4th parameter controls the aspect ratio of the cell(`2.0`). You can
        ## experiment with modifying these to change the initial shape of the cell.
        ##

        ## The second line of code instructs the mesh generator to generate a mesh.
        ## Checking the type of mesh with `type(mesh)` will show it as
        ## `ImmersedBoundaryMesh2_2`. The `2_2` suffix denotes that we are using
        ## a 2-dimensional space, and 2-dimensional elements to define the mesh.
        ##

        ## We now set the fluid grid resolution. The following code specifies
        ## that we are using a 64x64 grid to simulate our fluid over.

        mesh.SetNumGridPtsXAndY(64)

        ## Next, we generate the cells. We specify a cell type and cell cycle model. 
        ## These can be changed to modify the life cycle of the cells. The
        ## cell generator then constructs the necessary information for each
        ## of the elements in the mesh.

        cell_type = DifferentiatedCellProliferativeType()
        cell_generator = CellsGeneratorUniformCellCycleModel_2()
        cells = cell_generator.GenerateBasicRandom(mesh.GetNumElements(), cell_type)

        ## Finally, we construct the cell population. We then specify whether the
        ## population has active fluid sources or not. For now, we are not
        ## using any fluid sources, so we set this to `False`

        cell_population = ImmersedBoundaryCellPopulation2(mesh, cells)
        cell_population.SetIfPopulationHasActiveSources(False)

        ## We can make a quick visualization of the cell population

        scene = VtkScene2()
        scene.SetCellPopulation(cell_population)
        nb_manager = JupyterNotebookManager()
        nb_manager.vtk_show(scene, height=300)

        ## Next, we create an `OffLatticeSimulation` simulator to control the
        ## simulation. Although the fluid is simulated on a lattice (grid),
        ## the nodes/cells are not bound to a lattice.

        simulator = OffLatticeSimulation2_2(cell_population)
        simulator.SetNumericalMethod(ForwardEulerNumericalMethod2_2())
        simulator.GetNumericalMethod().SetUseUpdateNodeLocation(True)

        ## As we have an off-lattice simulation, we need a way to model the
        ## fluid. This is handled by the `ImmersedBoundarySimulationModifier`.
        ## Modifiers in Chaste are classes that can be attached to simulations
        ## to perform some additional custom functionality each timestep.
        ## In this case, the modifier is responsible for solving the
        ## Navier-Stokes equations and propagating forces between the nodes and
        ## the fluid.

        ib_modifier = ImmersedBoundarySimulationModifier2()
        simulator.AddSimulationModifier(ib_modifier)

        ## We must also provide the modifier with a force model to govern
        ## interactions between the nodes forming the cell membrane.
        ## Note that these forces only act between nodes in the same cell; 
        ## they do not control interactions between cells.

        membrane_force = ImmersedBoundaryLinearMembraneForce2()
        membrane_force.SetElementSpringConst(1.0 * 1e7)
        ib_modifier.AddImmersedBoundaryForce(membrane_force)

        ## The `ImmersedBoundaryLinearMembraneForce` models forces between
        ## membrane nodes using linear springs i.e, the force applied is
        ## proportional to the deviation of the distance between nodes
        ## from a rest length. The spring constant(`1.0 * 1e7`) defines how 
        ## stiff the cell boundary is.
        ##

        ## **Practice** Experiment with adjusting the spring constant to 
        ## change the force behaviour between nodes of the cell boundary.
        ##

        ## Next, we set the simulation properties

        dt = 0.05
        simulator.SetOutputDirectory("Python/TestImmersedBoundary_1")
        simulator.SetDt(dt)
        simulator.SetSamplingTimestepMultiple(4)
        simulator.SetEndTime(1000 * dt)
        
        ## We can add a modifier to visualize the cell population while the 
        ## simulation is in progress

        scene_modifier = JupyterSceneModifier2(nb_manager)
        scene_modifier.SetVtkScene(scene)
        scene_modifier.SetUpdateFrequency(1000)
        simulator.AddSimulationModifier(scene_modifier)

        ## Finally, to run the simulation we call the `Solve()` method.

        simulator.Solve()

        ## Reset the simulation environment in the notebook

        # JUPYTER_TEARDOWN

    ## ### 2. Adding More Cells

    def test_multicell_immersed_boundary_simulation(self):

        ## #### Multiple Cells

        ## Setup the simulation environment in the notebook

        # JUPYTER_SETUP

        ## We can use the mesh generator to generate multiple cells. The first
        ## parameter of the mesh generator constructor controls the number of
        ## cells.
        ##

        ## **Practice** Try increasing the number of cells by adjusting the 
        ## parameter value. A sensible range for this tutorial is 4-10 cells.

        gen = ImmersedBoundaryPalisadeMeshGenerator(5, 128, 0.1, 2.0, 0.0, False)

        ## #### Laminas
        ## In addition to the cells we have seen so far, we can introduce
        ## laminas to the simulation. Laminas are surfaces with reduced
        ## dimensionality. For 3D elements, a lamina is a 2D surface. For the
        ## 2D elements we are currently working with, laminas are lines. 
        ## Changing the last parameter of the mesh generator constructor from `False`
        ## to `True` will generate a basal lamina spanning the palisade cells.
        ## Laminas can also interact with the fluid field, and can be made
        ## "leaky" to allow some flow across their boundary. This can be used
        ## to model a permeable boundary.
        ##

        ## **Practice** Try changing the 6th constructor parameter to create a lamina.
        ##


        ## #### Cell Variations
        ## Apart from using the 3rd and 4th constructor parameters to modify
        ## the cell shapes, we can also introduce variation between cells by
        ## modifying the 5th parameter.
        ##

        ## **Practice** Try adjusting the 3rd and 4th constructor parameters to
        ## introduce cell variations.
        ##

        ## Next, we generate the mesh and set the fluid grid resolution

        mesh = gen.GetMesh()
        mesh.SetNumGridPtsXAndY(64)

        ## Below, we generate the cells

        cell_type = DifferentiatedCellProliferativeType()
        cell_generator = CellsGeneratorUniformCellCycleModel_2()
        cells = cell_generator.GenerateBasicRandom(mesh.GetNumElements(), cell_type)

        ## Then we set up the cell population with no active fluid sources

        cell_population = ImmersedBoundaryCellPopulation2(mesh, cells)
        cell_population.SetIfPopulationHasActiveSources(False)

        ## We can visualize the cell population below

        scene = VtkScene2()
        scene.SetCellPopulation(cell_population)
        nb_manager = JupyterNotebookManager()
        nb_manager.vtk_show(scene, height=300)

        ## Now we create a simulator to manage the simulation

        simulator = OffLatticeSimulation2_2(cell_population)
        simulator.SetNumericalMethod(ForwardEulerNumericalMethod2_2())
        simulator.GetNumericalMethod().SetUseUpdateNodeLocation(True)

        ## We add an immersed boundary simulation modifier to the simulator

        ib_modifier = ImmersedBoundarySimulationModifier2()
        simulator.AddSimulationModifier(ib_modifier)

        ## We then add a force law to the simulation modifier to model the 
        ## behaviour of the cell membrane

        membrane_force = ImmersedBoundaryLinearMembraneForce2()
        membrane_force.SetElementSpringConst(1.0 * 1e7)
        ib_modifier.AddImmersedBoundaryForce(membrane_force)

        ## #### Inter-cellular Interactions
        ## So far, we have encountered forces that act to maintain the shape
        ## of the cell membrane. We can also introduce an inter-cellular
        ## force law using `ImmersedBoundaryLinearInteractionForce`.
        ## This has a `SetSpringConst` method instead of a `SetElementSpringConst` 
        ## method. It also has a `SetRestLength` method that we can use to 
        ## modify the rest length.

        interaction_force = ImmersedBoundaryLinearInteractionForce2()
        interaction_force.SetSpringConst(1.0 * 1e6)
        ib_modifier.AddImmersedBoundaryForce(interaction_force)

        ## Next, we set the simulation properties

        dt = 0.05
        simulator.SetOutputDirectory("Python/TestImmersedBoundary_2")
        simulator.SetDt(dt)
        simulator.SetSamplingTimestepMultiple(4)
        simulator.SetEndTime(1000 * dt)

        ## Finally, we run the simulation

        simulator.Solve()

        ## We can visualize the end state of the cell population

        nb_manager.vtk_show(scene, height=300)

        ## Reset the simulation environment in the notebook

        # JUPYTER_TEARDOWN

    ## ### 3. Adding Fluid Sources
    ## Now that we are familiar with how to generate the cells, we will
    ## introduce fluid sources.

    def test_fluid_source_immersed_boundary_simulation(self):

        ## #### Adding a Fluid Source

        ## Setup the simulation environment in the notebook

        # JUPYTER_SETUP

        ## We begin by constructing a fluid source object:

        source = FluidSource2(0, 0.5, 0.7)

        ## This constructs a `FluidSource` object in 2 dimensions. The first
        ## parameter supplies the index of the fluid source. Each source we
        ## create must have a unique index. The next two parameters are the
        ## `x` and `y` coordinates of the source. Fluid sources in Chaste are
        ## point-like, that is to say they do not have any area/volume.
        ##
        ## Having created the fluid source, we set its strength:

        source.SetStrength(0.012)

        ## Next, we create the mesh

        gen = ImmersedBoundaryPalisadeMeshGenerator(5, 128, 0.1, 2.0, 0.0, False)
        mesh = gen.GetMesh()
        mesh.SetNumGridPtsXAndY(64)

        ## We must associate the source with an element in the simulation
        ## so that the simulation is aware of the source.

        mesh.GetElement(0).SetFluidSource(source)

        ## We now generate the cells

        cell_type = DifferentiatedCellProliferativeType()
        cell_generator = CellsGeneratorUniformCellCycleModel_2()
        cells = cell_generator.GenerateBasicRandom(mesh.GetNumElements(), cell_type)

        ## Then we set up the cell population

        cell_population = ImmersedBoundaryCellPopulation2(mesh, cells)

        ## Finally, we must tell the cell population that fluid sources are present.

        cell_population.SetIfPopulationHasActiveSources(True)

        ## #### Varying the Source Location and Strength
        ## **Practice** You can experiment with the source location. Try moving it 
        ## closer to and further away from the cells.
        ##

        ## **Practice** Try modifying the source strength to see what impact this 
        ## has on the cell shapes.
        ##

        ## Below, we visualize the cell population

        scene = VtkScene2()
        scene.SetCellPopulation(cell_population)
        nb_manager = JupyterNotebookManager()
        nb_manager.vtk_show(scene, height=300)

        ## Create a simulator to manage the simulation

        simulator = OffLatticeSimulation2_2(cell_population)
        simulator.SetNumericalMethod(ForwardEulerNumericalMethod2_2())
        simulator.GetNumericalMethod().SetUseUpdateNodeLocation(True)

        ## Add an immersed boundary simulation modifier

        ib_modifier = ImmersedBoundarySimulationModifier2()
        simulator.AddSimulationModifier(ib_modifier)

        ## #### Fluid-Cell Interaction
        ## **Practice** Try modifying the spring constant of the
        ## `ImmersedBoundaryLinearMembraneForce` to see how this changes the
        ## effect of the fluid source on the cells.

        membrane_force = ImmersedBoundaryLinearMembraneForce2()
        membrane_force.SetElementSpringConst(1.0 * 1e7)
        ib_modifier.AddImmersedBoundaryForce(membrane_force)

        ## Add an inter-cellular force law

        interaction_force = ImmersedBoundaryLinearInteractionForce2()
        interaction_force.SetSpringConst(1.0 * 1e6)
        ib_modifier.AddImmersedBoundaryForce(interaction_force)

        ## #### Adding More Sources
        ## **Practice** Try adding a second fluid source. You will need to 
        ## use a unique index, and attach it to a different element as 
        ## each element can only manage a single fluid source.
        ##

        ## Next, we set the simulation properties

        dt = 0.05
        simulator.SetOutputDirectory("Python/TestImmersedBoundary_3")
        simulator.SetDt(dt)
        simulator.SetSamplingTimestepMultiple(4)
        simulator.SetEndTime(300 * dt)

        ## Finally, we run the simulation

        simulator.Solve()

        ## Then we visualize the end state

        nb_manager.vtk_show(scene, height=300)
        
        ## Reset the simulation environment in the notebook
        # JUPYTER_TEARDOWN

        ## #### Further Exercises
        ## * Try integrating a different cell cycle model to introduce cell
        ##   division. See how the presence of a fluid source impacts the
        ##   structure that is formed.
        ## * Use one of the cell writers to collect some statistics


if __name__ == "__main__":
    unittest.main(verbosity=2)

#endif END_WIKI
