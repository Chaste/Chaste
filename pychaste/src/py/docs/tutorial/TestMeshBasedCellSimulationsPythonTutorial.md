
---
title : "Test Mesh Based Cell Simulations Python Tutorial"
summary: ""
draft: false
images: []
toc: true
layout: "single"
---

This tutorial is automatically generated from [TestMeshBasedCellSimulationsPythonTutorial](https://github.com/Chaste/PyChaste/blob/develop/test/python/cell_based/tutorials/TestMeshBasedCellSimulationsPythonTutorial.py) at revision [f810861a](https://github.com/Chaste/PyChaste/commit/f810861afe376ba19bd791e14e85f29583993205).
Note that the code is given in full at the bottom of the page.


## Introduction
In this tutorial we show how Chaste can be used to create, run and visualize mesh-based simulations.
Full details of the mathematical model can be found in van Leeuwen et al. (2009) [doi:10.1111/j.1365-2184.2009.00627.x].

## The Test

```python
import unittest  # Python testing framework
import matplotlib.pyplot as plt  # Plotting
import numpy as np  # Matrix tools
import chaste  # The PyChaste module
chaste.init()  # Set up MPI
import chaste.cell_based  # Contains cell populations
import chaste.mesh  # Contains meshes
import chaste.visualization  # Visualization tools

class TestRunningMeshBasedSimulationsTutorial(chaste.cell_based.AbstractCellBasedTestSuite):

```
### Test 1 - a basic mesh-based simulation
In the first test, we run a simple mesh-based simulation,
in which we create a monolayer of cells, using a mutable mesh. Each cell is assigned a stochastic cell-cycle model.

```python
    def test_monolayer(self):

        # JUPYTER_SETUP

```
Next, we generate a mutable mesh. To create a `MutableMesh`, we can use the `HoneycombMeshGenerator`.
This generates a honeycomb-shaped mesh, in which all nodes are equidistant. Here the first and second arguments define the size of the mesh -
we have chosen a mesh that is 4 nodes (i.e. cells) wide, and 4 nodes high.

```python
        chaste.core.OutputFileHandler("Python/TestMeshBasedCellSimulationsTutorial")
        generator = chaste.mesh.HoneycombMeshGenerator(4, 4)
        mesh = generator.GetMesh()

```
Having created a mesh, we now create some cells. To do this, we use the `CellsGenerator` helper class,
which is specialized by the type of cell cycle model required (here `UniformCellCycleModel`) and the dimension.
For a list of possible cell cycle models see subclasses of `AbstractCellCycleModel`.
Note that some of these models will require information on the surrounding medium such as Oxygen concentration to work,
see specific class documentation for details. We create an empty vector of cells and pass this into the method along with the mesh.
The second argument represents the size of that the list of cells should become - one cell for each node,
the third argument specifies the proliferative type of the cell.

```python
        transit_type = chaste.cell_based.TransitCellProliferativeType()
        cell_generator = chaste.cell_based.CellsGeneratorUniformCellCycleModel_2()
        cells = cell_generator.GenerateBasicRandom(mesh.GetNumNodes(),
                                                   transit_type)

```
Now we have a mesh and a set of cells to go with it, we can create a `CellPopulation`.
In general, this class associates a collection of cells with a mesh. For this test, because we have a `MutableMesh`,
we use a particular type of cell population called a `MeshBasedCellPopulation`.

```python
        cell_population = chaste.cell_based.MeshBasedCellPopulation2_2(mesh,
                                                                       cells)

```
To view the results of this and the next test in Paraview it is necessary to explicitly
generate the required .vtu files.

```python
        cell_population.AddPopulationWriterVoronoiDataWriter()

```
We can set up a `VtkScene` to do a quick visualization of the population before running the analysis.

```python
        scene = chaste.visualization.VtkScene2()
        scene.SetCellPopulation(cell_population)
        # JUPYTER_SHOW_FIRST
        scene.Start()  # JUPYTER_SHOW

```
We then pass in the cell population into an `OffLatticeSimulation`, and set the output directory and end time.

```python
        simulator = chaste.cell_based.OffLatticeSimulation2_2(cell_population)
        simulator.SetOutputDirectory("Python/TestMeshBasedCellSimulationsTutorial")
        simulator.SetEndTime(10.0)

```
For longer simulations, we may not want to output the results every time step. In this case we can use the following method,
to print results every 12 time steps instead. As the default time step used by the simulator is 30 seconds,
this method will cause the simulator to print results every 6 minutes (or 0.1 hours).

```python
        simulator.SetSamplingTimestepMultiple(12)

```
We must now create one or more force laws, which determine the mechanics of the centres of each cell in a cell population.
For this test, we use one force law, based on the spring based model, and pass it to the `OffLatticeSimulation`.
For a list of possible forces see subclasses of `AbstractForce`. Note that some of these forces are not compatible with mesh-based simulations,
see the specific class documentation for details. If you try to use an incompatible class then you will receive a warning.

```python
        force = chaste.cell_based.GeneralisedLinearSpringForce2_2()
        simulator.AddForce(force)

```
Save snapshot images of the population during the simulation

```python
        scene_modifier = chaste.cell_based.VtkSceneModifier2()
        scene_modifier.SetVtkScene(scene)
        scene_modifier.SetUpdateFrequency(100)
        simulator.AddSimulationModifier(scene_modifier)

```
To run the simulation, we call `Solve()`. We can again do a quick rendering of the population at the end of the simulation

```python
        scene.Start()
        simulator.Solve()
        scene.End()

        # JUPYTER_TEARDOWN

```
Full results can be visualized in Paraview from the `file_handler.GetOutputDirectoryFullPath()` directory.

### Test 2 -  a basic mesh-based simulation with ghost nodes
In the second test, we run a simple mesh-based simulation with ghost nodes, in which we create a monolayer of cells, using a mutable mesh.
Each cell is assigned a stochastic cell-cycle model.

```python
    def test_monolayer_with_ghost_nodes(self):

        # JUPYTER_SETUP

```
We start by generating a mutable mesh. To create a `MutableMesh`, we can use the `HoneycombMeshGenerator` as before.
Here the first and second arguments define the size of the mesh - we have chosen a mesh that is 2 nodes (i.e. cells) wide,
and 2 nodes high. The third argument specifies the number of layers of ghost nodes to make.

```python
        chaste.core.OutputFileHandler("Python/TestMeshBasedCellPopulationWithGhostNodes")
        generator = chaste.mesh.HoneycombMeshGenerator(5, 5, 2)
        mesh = generator.GetMesh()

```
We only want to create cells to attach to real nodes, so we use the method `GetCellLocationIndices` to get the
indices of the real nodes in the mesh. This will be passed in to the cell population later on.

```python
        locs = generator.GetCellLocationIndices()

```
Having created a mesh, we now create some cells. To do this, we use the `CellsGenerator` helper class again.
This time the second argument is different and is the number of real nodes in the mesh.
As before all cells have `TransitCellProliferativeType`.

```python
        transit_type = chaste.cell_based.TransitCellProliferativeType()
        cell_generator = chaste.cell_based.CellsGeneratorUniformCellCycleModel_2()
        cells = cell_generator.GenerateBasicRandom(len(locs),
                                                   transit_type)

```
Now we have a mesh and a set of cells to go with it, we can create a `CellPopulation`.
In general, this class associates a collection of cells with a set of elements or a mesh.
For this test, because we have a `MutableMesh`, and ghost nodes we use a particular type of cell population called
a `MeshBasedCellPopulationWithGhostNodes`. The third argument of the constructor takes a vector of the indices of the real nodes
and should be the same length as the vector of cell pointers.

```python
        cell_population = chaste.cell_based.MeshBasedCellPopulationWithGhostNodes2(mesh,
                                                                                   cells,
                                                                                   locs)

```
Again Paraview output is explicitly requested.

```python
        cell_population.AddPopulationWriterVoronoiDataWriter()

```
We can set up a `VtkScene` to do a quick visualization of the population before running the analysis.

```python
        scene = chaste.visualization.VtkScene2()
        scene.SetCellPopulation(cell_population)
        scene.GetCellPopulationActorGenerator().SetShowVoronoiMeshEdges(True)
        # JUPYTER_SHOW

```
We then pass in the cell population into an `OffLatticeSimulation`, and set the output directory, output multiple and end time.

```python
        simulator = chaste.cell_based.OffLatticeSimulation2_2(cell_population)
        simulator.SetOutputDirectory("Python/TestMeshBasedCellPopulationWithGhostNodes")
        simulator.SetEndTime(10.0)
        simulator.SetSamplingTimestepMultiple(12)

```
Save snapshot images of the population during the simulation

```python
        scene_modifier = chaste.cell_based.VtkSceneModifier2()
        scene_modifier.SetVtkScene(scene)
        scene_modifier.SetUpdateFrequency(300)
        simulator.AddSimulationModifier(scene_modifier)

```
Again we create a force law, and pass it to the `OffLatticeSimulation`.
This force law ensures that ghost nodes don't exert forces on real nodes but real nodes exert forces on ghost nodes.

```python
        force = chaste.cell_based.GeneralisedLinearSpringForce2_2()
        simulator.AddForce(force)

```
To run the simulation, we call `Solve()`.

```python
        scene.Start()
        simulator.Solve()

```
The next two lines are for test purposes only and are not part of this tutorial.
If different simulation input parameters are being explored the lines should be removed.

```python
        self.assertEqual(cell_population.GetNumRealCells(), 48)
        self.assertAlmostEqual(chaste.cell_based.SimulationTime.Instance().GetTime(), 10.0, 6)

        # JUPYTER_TEARDOWN

```
Full results can be visualized in Paraview from the `file_handler.GetOutputDirectoryFullPath()` directory.

```python
if __name__ == '__main__':
    unittest.main(verbosity=2)

```


## Full code 


**File name:** `TestMeshBasedCellSimulationsPythonTutorial.py` 

```python
import unittest  # Python testing framework
import matplotlib.pyplot as plt  # Plotting
import numpy as np  # Matrix tools
import chaste  # The PyChaste module
chaste.init()  # Set up MPI
import chaste.cell_based  # Contains cell populations
import chaste.mesh  # Contains meshes
import chaste.visualization  # Visualization tools

class TestRunningMeshBasedSimulationsTutorial(chaste.cell_based.AbstractCellBasedTestSuite):

    def test_monolayer(self):

        # JUPYTER_SETUP

        chaste.core.OutputFileHandler("Python/TestMeshBasedCellSimulationsTutorial")
        generator = chaste.mesh.HoneycombMeshGenerator(4, 4)
        mesh = generator.GetMesh()

        transit_type = chaste.cell_based.TransitCellProliferativeType()
        cell_generator = chaste.cell_based.CellsGeneratorUniformCellCycleModel_2()
        cells = cell_generator.GenerateBasicRandom(mesh.GetNumNodes(),
                                                   transit_type)

        cell_population = chaste.cell_based.MeshBasedCellPopulation2_2(mesh,
                                                                       cells)

        cell_population.AddPopulationWriterVoronoiDataWriter()

        scene = chaste.visualization.VtkScene2()
        scene.SetCellPopulation(cell_population)
        # JUPYTER_SHOW_FIRST
        scene.Start()  # JUPYTER_SHOW

        simulator = chaste.cell_based.OffLatticeSimulation2_2(cell_population)
        simulator.SetOutputDirectory("Python/TestMeshBasedCellSimulationsTutorial")
        simulator.SetEndTime(10.0)

        simulator.SetSamplingTimestepMultiple(12)

        force = chaste.cell_based.GeneralisedLinearSpringForce2_2()
        simulator.AddForce(force)

        scene_modifier = chaste.cell_based.VtkSceneModifier2()
        scene_modifier.SetVtkScene(scene)
        scene_modifier.SetUpdateFrequency(100)
        simulator.AddSimulationModifier(scene_modifier)

        scene.Start()
        simulator.Solve()
        scene.End()

        # JUPYTER_TEARDOWN

    def test_monolayer_with_ghost_nodes(self):

        # JUPYTER_SETUP

        chaste.core.OutputFileHandler("Python/TestMeshBasedCellPopulationWithGhostNodes")
        generator = chaste.mesh.HoneycombMeshGenerator(5, 5, 2)
        mesh = generator.GetMesh()

        locs = generator.GetCellLocationIndices()

        transit_type = chaste.cell_based.TransitCellProliferativeType()
        cell_generator = chaste.cell_based.CellsGeneratorUniformCellCycleModel_2()
        cells = cell_generator.GenerateBasicRandom(len(locs),
                                                   transit_type)

        cell_population = chaste.cell_based.MeshBasedCellPopulationWithGhostNodes2(mesh,
                                                                                   cells,
                                                                                   locs)

        cell_population.AddPopulationWriterVoronoiDataWriter()

        scene = chaste.visualization.VtkScene2()
        scene.SetCellPopulation(cell_population)
        scene.GetCellPopulationActorGenerator().SetShowVoronoiMeshEdges(True)
        # JUPYTER_SHOW

        simulator = chaste.cell_based.OffLatticeSimulation2_2(cell_population)
        simulator.SetOutputDirectory("Python/TestMeshBasedCellPopulationWithGhostNodes")
        simulator.SetEndTime(10.0)
        simulator.SetSamplingTimestepMultiple(12)

        scene_modifier = chaste.cell_based.VtkSceneModifier2()
        scene_modifier.SetVtkScene(scene)
        scene_modifier.SetUpdateFrequency(300)
        simulator.AddSimulationModifier(scene_modifier)

        force = chaste.cell_based.GeneralisedLinearSpringForce2_2()
        simulator.AddForce(force)

        scene.Start()
        simulator.Solve()

        self.assertEqual(cell_population.GetNumRealCells(), 48)
        self.assertAlmostEqual(chaste.cell_based.SimulationTime.Instance().GetTime(), 10.0, 6)

        # JUPYTER_TEARDOWN

if __name__ == '__main__':
    unittest.main(verbosity=2)

```

