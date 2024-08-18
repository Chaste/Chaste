
---
title : "Test Scratch Assay Tutorial"
summary: ""
draft: false
images: []
toc: true
layout: "single"
---

This tutorial is automatically generated from [TestScratchAssayTutorial](https://github.com/Chaste/PyChaste/blob/develop/test/python/cell_based/tutorials/TestScratchAssayTutorial.py) at revision [f810861a](https://github.com/Chaste/PyChaste/commit/f810861afe376ba19bd791e14e85f29583993205).
Note that the code is given in full at the bottom of the page.


## Introduction
This tutorial is an example of modelling a scratch assay using a simple cellular automaton
representation of cells. It will cover the following techniques:

 * Setting up a regular mesh (or lattice)
 * Visualizing the mesh
 * Working with file-based output
 * Generating cells and adding them to the mesh
 * Simulating cell migration on the mesh
 * Real-time visualization of the cell population and plotting of population statistics
 
## The Test

```python
import unittest # Python testing framework
import matplotlib.pyplot as plt # Plotting
import numpy as np # Matrix tools
import chaste # The PyChaste module
chaste.init() # Set up MPI
import chaste.cell_based # Contains cell populations
import chaste.mesh # Contains meshes
import chaste.visualization # Visualization tools

class TestScratchAssayTutorial(chaste.cell_based.AbstractCellBasedTestSuite):

```
### Test 1 - Scratch Assay
In this test we will create a scratch along the middle of a domain and quantify the migration
of cells into the region. Cells will migrate by random walk on the their regular mesh  (lattice).

```python
    def test_single_scratch(self):

        # JUPYTER_SETUP

```
Chaste is based on the concept of `Cells` and `Meshes`. 'Cells' do not store their position in space,
or connectivity, these are managed by a `Mesh`. The first step in most Chaste simulations is to
set up a mesh, on which we can locate cells. A collection of `Cells` and a `Mesh` are a `CellPopulation`
in Chaste terminology. The most simple `CellPopulation` is the `CaBasedCellPopulation` which corresponds
to cells occupying discrete locations on a regular mesh (lattice). Our first step is to set up the mesh.
Here we set up a 2D lattice.

```python
        num_points_in_x = 100
        num_points_in_y = 12
        generator = chaste.mesh.PottsMeshGenerator2(num_points_in_x, 0, 0, num_points_in_y, 0, 0)
        mesh = generator.GetMesh()

```
Note that we are using a `PottsMeshGenerator2` to set up the grid and we are setting some terms to 0. Chaste
design is based on re-use of components, the `PottsMeshGenerator` can be used to set up other types of
cell population which require these extra terms. Note also the '2' at the end of the class name. This
tells us that we are working in 2D. Most Chaste classes are specialized (templated) for spatial dimension,
so we need to make sure we are consistent in the dimensionality of the classes we are using.

Next we set up some cells. We create and empty container `VectorSharedPtrCell` (which will behave like a Python list)
and will fill it with cells of our chosen type. In Chaste cells can be assinged a number of proliferative types
(Default, Differentiated, Stem, Transit or User Defined). These types will define how cells behave in certain
simulations, for example whether they will proliferate. We just want our cells to migrate in this example, so
we set a DifferentiatedCellProliferativeType.

```python
        cells = []
        differentiated_type = chaste.cell_based.DifferentiatedCellProliferativeType()

```
We are not interested in cell cycling so we specialize the generator to NoCellCycleModel.

```python
        cell_generator = chaste.cell_based.CellsGeneratorNoCellCycleModel_2()

```
We want two sets of cells, starting on opposite sides of the mesh. We use `location_indices` to map cells onto
locations (or Nodes) on the mesh. For our regular mesh the Node indices increase fastest in x, then y. We will
add four layers of cells to each side of the mesh.

```python
        num_cell_layers = 4
        bottom_location_indices = list(range(num_cell_layers*num_points_in_x))
        num_grid_points = num_points_in_x*num_points_in_y
        top_location_indices = list(range(num_grid_points-1, num_grid_points -
                                     num_cell_layers*num_points_in_x-1, -1))
        cells = cell_generator.GenerateGivenLocationIndices(
                                                    bottom_location_indices + top_location_indices,
                                                    differentiated_type)

```
Now we have a mesh and a set of cells to go with it, we can create a CellPopulation.

```python
        cell_population = chaste.cell_based.CaBasedCellPopulation2(mesh, cells,
                                                                   bottom_location_indices +
                                                                   top_location_indices)

```
Next, we set up an `OffLatticeSimulation` which will manage the solver. We need to add some custom rules to
this solver to specify how we want the cells to migrate.

```python
        simulator = chaste.cell_based.OnLatticeSimulation2(cell_population)
        simulator.SetOutputDirectory("Python/TestScratchAssayTutorial")
        simulator.SetEndTime(10.0)
        simulator.SetDt(0.1)
        simulator.SetSamplingTimestepMultiple(1)

```
We must now create a rule for cell migration. We will use an existing diffusion type rule.

```python
        diffusion_update_rule = chaste.cell_based.DiffusionCaUpdateRule2()
        simulator.AddUpdateRule(diffusion_update_rule)

```
PyChaste can do simple 3D rendering with VTK. We set up a `VtkScene` so that we can see the population
evovle in real time.

```python
        scene= chaste.visualization.VtkScene2()
        scene.SetCellPopulation(cell_population)
        scene.GetCellPopulationActorGenerator().SetShowCellCentres(True)
        # JUPYTER_SHOW_FIRST
        scene.Start()  # JUPYTER_SHOW

```
We add the scene to the simulation for real-time updating using a `VtkSceneModifier`. Such
modifiers are called by the simulator at regular periods during the main time loop and
have access to the cell population. We will use a similar idea in a moment to record cell
positions for real time plotting.

```python
        scene_modifier = chaste.cell_based.VtkSceneModifier2()
        scene_modifier.SetVtkScene(scene)
        scene_modifier.SetUpdateFrequency(10)
        simulator.AddSimulationModifier(scene_modifier)

```
Chaste and PyChaste use object oriented programming. This may require some background reading,
but allows for great flexibility in terms of modifying existing functionality. In
order to pull the data we want out of the simulation as it runs we will create our own
simulation modifier class and use it for real time plotting. This Python class over-rides
one of the built-in classes, giving us access to the quantities we want during the simulation.
Usually we would define such a class in a different module and import it, it is placed
here for the purposes of the tutorial.

```python
        class PlottingModifier(chaste.cell_based.PythonSimulationModifier2):

            """ Class for real time plotting of cell numbers using Matplotlib
            """

            def __init__(self, num_points_in_x, num_points_in_y):
                super(PlottingModifier, self).__init__()

                # Set up a figure for plotting
                plt.ioff()
                self.fig = plt.figure()
                self.fig.ax = self.fig.add_subplot(111)
                self.fig.ax.set_xlabel("y - Position (Cell Lengths)")
                self.fig.ax.set_ylabel("Number Of Cells")
                self.plot_frequency = 10 # only plot every 10 steps
                self.num_points_in_x = num_points_in_x
                self.num_points_in_y = num_points_in_y

            def UpdateAtEndOfTimeStep(self, cell_population):

                """ Plot the number of cells at each lattice point and time-point

                Use the SimulationTime singleton to determine when to plot.
                """

                num_increments = chaste.cell_based.SimulationTime.Instance().GetTimeStepsElapsed()
                if num_increments % self.plot_frequency == 0:
                    y_locations = np.linspace(0, num_points_in_y, num_points_in_y)
                    num_cells = []
                    for idx in range(num_points_in_y):
                        counter = 0
                        for jdx in range(num_points_in_x):
                            if cell_population.IsCellAttachedToLocationIndex(jdx +
                                                                             idx*num_points_in_x):
                                counter +=1
                        num_cells.append(counter)

                    self.fig.ax.plot(y_locations, num_cells, color='black')
                    self.fig.canvas.draw()
                    #display.display(self.fig)
                    #display.clear_output(wait=True)

            def SetupSolve(self, cell_population, output_directory):

                """ Ensure the cell population is in the correct state at the start of the simulation
                """

                cell_population.Update()

        plotting_modifier = PlottingModifier(num_points_in_x, num_points_in_y)
        simulator.AddSimulationModifier(plotting_modifier)

```
To run the simulation, we call `Solve()` and optionally set up interactive plotting. We will see the cells
migrate and the population distribution gradually become more uniform.

```python
        scene.Start()
        plt.ion()
        plt.show()
        simulator.Solve()

        # JUPYTER_TEARDOWN

if __name__ == '__main__':
    unittest.main(verbosity=2)

```


## Full code 


**File name:** `TestScratchAssayTutorial.py` 

```python
import unittest # Python testing framework
import matplotlib.pyplot as plt # Plotting
import numpy as np # Matrix tools
import chaste # The PyChaste module
chaste.init() # Set up MPI
import chaste.cell_based # Contains cell populations
import chaste.mesh # Contains meshes
import chaste.visualization # Visualization tools

class TestScratchAssayTutorial(chaste.cell_based.AbstractCellBasedTestSuite):

    def test_single_scratch(self):

        # JUPYTER_SETUP

        num_points_in_x = 100
        num_points_in_y = 12
        generator = chaste.mesh.PottsMeshGenerator2(num_points_in_x, 0, 0, num_points_in_y, 0, 0)
        mesh = generator.GetMesh()

        cells = []
        differentiated_type = chaste.cell_based.DifferentiatedCellProliferativeType()

        cell_generator = chaste.cell_based.CellsGeneratorNoCellCycleModel_2()

        num_cell_layers = 4
        bottom_location_indices = list(range(num_cell_layers*num_points_in_x))
        num_grid_points = num_points_in_x*num_points_in_y
        top_location_indices = list(range(num_grid_points-1, num_grid_points -
                                     num_cell_layers*num_points_in_x-1, -1))
        cells = cell_generator.GenerateGivenLocationIndices(
                                                    bottom_location_indices + top_location_indices,
                                                    differentiated_type)

        cell_population = chaste.cell_based.CaBasedCellPopulation2(mesh, cells,
                                                                   bottom_location_indices +
                                                                   top_location_indices)

        simulator = chaste.cell_based.OnLatticeSimulation2(cell_population)
        simulator.SetOutputDirectory("Python/TestScratchAssayTutorial")
        simulator.SetEndTime(10.0)
        simulator.SetDt(0.1)
        simulator.SetSamplingTimestepMultiple(1)

        diffusion_update_rule = chaste.cell_based.DiffusionCaUpdateRule2()
        simulator.AddUpdateRule(diffusion_update_rule)

        scene= chaste.visualization.VtkScene2()
        scene.SetCellPopulation(cell_population)
        scene.GetCellPopulationActorGenerator().SetShowCellCentres(True)
        # JUPYTER_SHOW_FIRST
        scene.Start()  # JUPYTER_SHOW

        scene_modifier = chaste.cell_based.VtkSceneModifier2()
        scene_modifier.SetVtkScene(scene)
        scene_modifier.SetUpdateFrequency(10)
        simulator.AddSimulationModifier(scene_modifier)

        class PlottingModifier(chaste.cell_based.PythonSimulationModifier2):

            """ Class for real time plotting of cell numbers using Matplotlib
            """

            def __init__(self, num_points_in_x, num_points_in_y):
                super(PlottingModifier, self).__init__()

                # Set up a figure for plotting
                plt.ioff()
                self.fig = plt.figure()
                self.fig.ax = self.fig.add_subplot(111)
                self.fig.ax.set_xlabel("y - Position (Cell Lengths)")
                self.fig.ax.set_ylabel("Number Of Cells")
                self.plot_frequency = 10 # only plot every 10 steps
                self.num_points_in_x = num_points_in_x
                self.num_points_in_y = num_points_in_y

            def UpdateAtEndOfTimeStep(self, cell_population):

                """ Plot the number of cells at each lattice point and time-point

                Use the SimulationTime singleton to determine when to plot.
                """

                num_increments = chaste.cell_based.SimulationTime.Instance().GetTimeStepsElapsed()
                if num_increments % self.plot_frequency == 0:
                    y_locations = np.linspace(0, num_points_in_y, num_points_in_y)
                    num_cells = []
                    for idx in range(num_points_in_y):
                        counter = 0
                        for jdx in range(num_points_in_x):
                            if cell_population.IsCellAttachedToLocationIndex(jdx +
                                                                             idx*num_points_in_x):
                                counter +=1
                        num_cells.append(counter)

                    self.fig.ax.plot(y_locations, num_cells, color='black')
                    self.fig.canvas.draw()
                    #display.display(self.fig)
                    #display.clear_output(wait=True)

            def SetupSolve(self, cell_population, output_directory):

                """ Ensure the cell population is in the correct state at the start of the simulation
                """

                cell_population.Update()

        plotting_modifier = PlottingModifier(num_points_in_x, num_points_in_y)
        simulator.AddSimulationModifier(plotting_modifier)

        scene.Start()
        plt.ion()
        plt.show()
        simulator.Solve()

        # JUPYTER_TEARDOWN

if __name__ == '__main__':
    unittest.main(verbosity=2)

```

