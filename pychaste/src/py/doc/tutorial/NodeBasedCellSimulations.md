---
title : "Node Based Cell Simulations"
summary: ""
draft: false
images: []
toc: true
layout: "single"
---
This tutorial is automatically generated from [TestPyNodeBasedCellSimulationsTutorial.py](https://github.com/Chaste/Chaste/blob/develop/pychaste/test/tutorial/TestPyNodeBasedCellSimulationsTutorial.py) at revision [4045f91a83f5](https://github.com/Chaste/Chaste/commit/4045f91a83f55dc4a97f2ca4f97b0c32f4e43a4a).

Note that the code is given in full at the bottom of the page.

## Introduction
In this tutorial we show how Chaste can be used to create, run and visualize node-based simulations. Full details of the mechanical model can be found in Pathamathan et
al "A computational study of discrete mechanical tissue models", Physical Biology. Vol. 6. No. 3. 2009.. DOI (10.1088/1478-3975/6/3/036001).

## The Test

```python
import unittest  # Python testing framework

import numpy as np  # Matrix tools

import chaste  # The PyChaste module
import chaste.cell_based  # Contains cell populations
import chaste.mesh  # Contains meshes
import chaste.visualization  # Visualization tools

from chaste.cell_based import AbstractCellBasedTestSuite

class TestPyNodeBasedCellSimulationsTutorial(AbstractCellBasedTestSuite):
```
### Test 1 - A basic node-based simulation
In the first test, we run a simple node-based simulation, in which we create a monolayer of cells,
using a nodes only mesh. Each cell is assigned a uniform cell-cycle model.

```python
    def test_monolayer(self):

        # JUPYTER_SETUP
```
The first thing we do is generate a nodes only mesh. To do this we first create a `MutableMesh` to use as a generating mesh.
To do this we can use the `HoneycombMeshGenerator`. This generates a honeycomb-shaped mesh, in which all nodes are equidistant.
Here the first and second arguments define the size of the mesh - we have chosen a mesh that is 2 nodes (i.e. cells) wide, and 2 nodes high.

```python
        chaste.core.OutputFileHandler("Python/TestNodeBasedCellSimulationsTutorial")
        generator = chaste.mesh.HoneycombMeshGenerator(2, 2)
        generating_mesh = generator.GetMesh()
```
Once we have a MutableMesh we can generate a NodesOnlyMesh from it using the following commands.
Note you can also generate the NodesOnlyMesh from a collection of nodes.

```python
        mesh = chaste.mesh.NodesOnlyMesh[2]()
```
To run node-based simulations you need to define a cut off length (second argument in `ConstructNodesWithoutMesh`),
which defines the connectivity of the nodes by defining a radius of interaction.

```python
        mesh.ConstructNodesWithoutMesh(generating_mesh, 1.5)
```
Having created a mesh, we now create a (wrapped) vector of CellPtrs. To do this, we use the `CellsGenerator` helper class,
which is specialized for the type of cell model required (here `UniformCellCycleModel`) and the dimension.
We create an empty vector of cells and pass this into the method along with the mesh.
The second argument represents the size of that the vector cells should become - one cell for each node,
the third argument specifies the proliferative type of the cell.

```python
        transit_type = chaste.cell_based.TransitCellProliferativeType()
        cell_generator = chaste.cell_based.CellsGenerator["UniformCellCycleModel", 2]()
        cells = cell_generator.GenerateBasicRandom(mesh.GetNumNodes(), transit_type)
```
Now we have a mesh and a set of cells to go with it, we can create a `CellPopulation`.
In general, this class associates a collection of cells with a mesh. For this test,
because we have a `NodesOnlyMesh`, we use a particular type of cell population called a `NodeBasedCellPopulation`.

```python
        cell_population = chaste.cell_based.NodeBasedCellPopulation[2](mesh, cells)
```
We can set up a `VtkScene` to do a quick visualization of the population before running the analysis.

```python
        scene = chaste.visualization.VtkScene[2]()
        scene.SetCellPopulation(cell_population)
        # JUPYTER_SHOW_FIRST
        scene.Start()  # JUPYTER_SHOW
```
We then pass in the cell population into an `OffLatticeSimulation`, and set the output directory, output multiple and end time

```python
        simulator = chaste.cell_based.OffLatticeSimulation[2, 2](cell_population)
        simulator.SetOutputDirectory("Python/TestNodeBasedCellSimulationsTutorial")
        simulator.SetSamplingTimestepMultiple(100)
        simulator.SetEndTime(10.0)
```
We now pass a force law to the simulation.

```python
        force = chaste.cell_based.GeneralisedLinearSpringForce[2, 2]()
        simulator.AddForce(force)
```
Save snapshot images of the population during the simulation

```python
        scene_modifier = chaste.cell_based.VtkSceneModifier[2]()
        scene_modifier.SetVtkScene(scene)
        scene_modifier.SetUpdateFrequency(100)
        simulator.AddSimulationModifier(scene_modifier)
```
To run the simulation, we call `Solve()`. We can again do a quick rendering of the population at the end of the simulation

```python
        scene.Start()
        simulator.Solve()
        scene.End()
```
The next two lines are for test purposes only and are not part of this tutorial.
If different simulation input parameters are being explored the lines should be removed.

```python
        self.assertEqual(cell_population.GetNumRealCells(), 8)
        self.assertAlmostEqual(chaste.cell_based.SimulationTime.Instance().GetTime(), 10.0, 6)

        # JUPYTER_TEARDOWN
```
### Test 2 - a basic node-based simulation in 3D
In the second test we run a simple node-based simulation in 3D. This is very similar to the 2D test with the dimension changed from 2 to 3 and
instead of using a mesh generator we generate the nodes directly.

```python
    def test_spheroid(self):

        # JUPYTER_SETUP
```
First, we generate a nodes only mesh. This time we specify the nodes manually by first creating a vector of nodes

```python
        chaste.core.OutputFileHandler(
            "Python/TestNodeBasedCellSimulationsSpheroidTutorial"
        )
        nodes = []
        nodes.append(chaste.mesh.Node[3](0, False, 0.5, 0.0, 0.0))
        nodes.append(chaste.mesh.Node[3](1, False, -0.5, 0.0, 0.0))
        nodes.append(chaste.mesh.Node[3](2, False, 0.0, 0.5, 0.0))
        nodes.append(chaste.mesh.Node[3](3, False, 0.0, -0.5, 0.0))
```
Finally a NodesOnlyMesh is created and the vector of nodes is passed to the ConstructNodesWithoutMesh method.

```python
        mesh = chaste.mesh.NodesOnlyMesh[3]()
```
To run node-based simulations you need to define a cut off length (second argument in ConstructNodesWithoutMesh),
which defines the connectivity of the nodes by defining a radius of interaction.

```python
        mesh.ConstructNodesWithoutMesh(nodes, 1.5)
```
Having created a mesh, we now create a std::vector of CellPtrs.
As before, we do this with the CellsGenerator helper class (this time with dimension 3).

```python
        transit_type = chaste.cell_based.TransitCellProliferativeType()
        cell_generator = chaste.cell_based.CellsGenerator["UniformCellCycleModel", 3]()
        cells = cell_generator.GenerateBasicRandom(mesh.GetNumNodes(), transit_type)
```
Now we have a mesh and a set of cells to go with it, we can create a `CellPopulation`.
In general, this class associates a collection of cells with a mesh. For this test,
because we have a `NodesOnlyMesh`, we use a particular type of cell population called a `NodeBasedCellPopulation`.

```python
        cell_population = chaste.cell_based.NodeBasedCellPopulation[3](mesh, cells)
```
We can set up a `VtkScene` to do a quick visualization of the population before running the analysis.

```python
        scene = chaste.visualization.VtkScene[3]()
        scene.SetCellPopulation(cell_population)
        scene.Start()  # JUPYTER_SHOW
```
We then pass in the cell population into an `OffLatticeSimulation`, and set the output directory, output multiple and end time

```python
        simulator = chaste.cell_based.OffLatticeSimulation[3, 3](cell_population)
        simulator.SetOutputDirectory("Python/TestNodeBasedCellSimulationsSpheroidTutorial")
        simulator.SetSamplingTimestepMultiple(12)
        simulator.SetEndTime(10.0)
```
We now pass a force law to the simulation.

```python
        force = chaste.cell_based.GeneralisedLinearSpringForce[3, 3]()
        simulator.AddForce(force)
```
Save snapshot images of the population during the simulation

```python
        scene_modifier = chaste.cell_based.VtkSceneModifier[3]()
        scene_modifier.SetVtkScene(scene)
        scene_modifier.SetUpdateFrequency(100)
        simulator.AddSimulationModifier(scene_modifier)
```
To run the simulation, we call `Solve()`. We can again do a quick rendering of the population at the end of the simulation

```python
        scene.Start()
        simulator.Solve()
        scene.End()
```
The next two lines are for test purposes only and are not part of this tutorial.
If different simulation input parameters are being explored the lines should be removed.

```python
        self.assertEqual(cell_population.GetNumRealCells(), 8)
        self.assertAlmostEqual(chaste.cell_based.SimulationTime.Instance().GetTime(), 10.0, 6)

        # JUPYTER_TEARDOWN
```
### Test 3 - a node-based simulation on a restricted geometry
In the second test we run a simple node-based simulation in 3D. This is very similar to the 2D test with the dimension changed from 2 to 3 and
instead of using a mesh generator we generate the nodes directly.

```python
    def test_spheroid_on_sphere(self):

        # JUPYTER_SETUP
```
In the third test we run a node-based simulation restricted to the surface of a sphere.

```python
        chaste.core.OutputFileHandler("Python/TestNodeBasedCellSimulationsRestrictedSpheroidTutorial")
        nodes = []
        nodes.append(chaste.mesh.Node[3](0, False, 0.5, 0.0, 0.0))
        nodes.append(chaste.mesh.Node[3](1, False, -0.5, 0.0, 0.0))
        nodes.append(chaste.mesh.Node[3](2, False, 0.0, 0.5, 0.0))
        nodes.append(chaste.mesh.Node[3](3, False, 0.0, -0.5, 0.0))
        mesh = chaste.mesh.NodesOnlyMesh[3]()
```
To run node-based simulations you need to define a cut off length (second argument in ConstructNodesWithoutMesh),
which defines the connectivity of the nodes by defining a radius of interaction.

```python
        mesh.ConstructNodesWithoutMesh(nodes, 1.5)

        transit_type = chaste.cell_based.TransitCellProliferativeType()
        cell_generator = chaste.cell_based.CellsGenerator["UniformCellCycleModel", 3]()
        cells = cell_generator.GenerateBasicRandom(mesh.GetNumNodes(), transit_type)
        cell_population = chaste.cell_based.NodeBasedCellPopulation[3](mesh, cells)
```
We can set up a `VtkScene` to do a quick visualization of the population before running the analysis.

```python
        scene = chaste.visualization.VtkScene[3]()
        scene.SetCellPopulation(cell_population)
        scene.Start()  # JUPYTER_SHOW

        simulator = chaste.cell_based.OffLatticeSimulation[3, 3](cell_population)
        simulator.SetOutputDirectory("Python/TestNodeBasedCellSimulationsRestrictedSpheroidTutorial")
        simulator.SetSamplingTimestepMultiple(12)
        simulator.SetEndTime(10.0)
```
We now pass a force law to the simulation.

```python
        force = chaste.cell_based.GeneralisedLinearSpringForce[3, 3]()
        simulator.AddForce(force)
```
This time we create a CellPopulationBoundaryCondition and pass this to the OffLatticeSimulation.
Here we use a SphereGeometryBoundaryCondition which restricts cells to lie on a sphere (in 3D) or circle (in 2D).
For a list of possible boundary conditions see subclasses of AbstractCellPopulationBoundaryCondition.
Note that some of these boundary conditions are not compatible with node-based simulations see the specific class documentation
for details, if you try to use an incompatible class then you will receive a warning.
First we set the centre (0,0,1) and radius of the sphere (1).

```python
        centre = np.array([0.0, 0.0, 1.0])
        radius = 5.0
        point2 = chaste.mesh.ChastePoint[3](centre)
        boundary_condition = chaste.cell_based.SphereGeometryBoundaryCondition[3](cell_population, point2.rGetLocation(), radius)
        simulator.AddCellPopulationBoundaryCondition(boundary_condition)
```
Save snapshot images of the population during the simulation
scene_modifier = chaste.cell_based.VtkSceneModifier[3]()

```python
        scene_modifier.SetVtkScene(scene)
        scene_modifier.SetUpdateFrequency(100)
        simulator.AddSimulationModifier(scene_modifier)
```
To run the simulation, we call `Solve()`. We can again do a quick rendering of the population at the end of the simulation

```python
        scene.Start()
        simulator.Solve()
        scene.End()
```
The next two lines are for test purposes only and are not part of this tutorial.
If different simulation input parameters are being explored the lines should be removed.

```python
        self.assertEqual(cell_population.GetNumRealCells(), 8)
        self.assertAlmostEqual(chaste.cell_based.SimulationTime.Instance().GetTime(), 10.0, 6)

        # JUPYTER_TEARDOWN

if __name__ == "__main__":
    unittest.main(verbosity=2)
```

## Full code

```python
import unittest  # Python testing framework

import numpy as np  # Matrix tools

import chaste  # The PyChaste module
import chaste.cell_based  # Contains cell populations
import chaste.mesh  # Contains meshes
import chaste.visualization  # Visualization tools

from chaste.cell_based import AbstractCellBasedTestSuite

class TestPyNodeBasedCellSimulationsTutorial(AbstractCellBasedTestSuite):
    def test_monolayer(self):

        # JUPYTER_SETUP

        chaste.core.OutputFileHandler("Python/TestNodeBasedCellSimulationsTutorial")
        generator = chaste.mesh.HoneycombMeshGenerator(2, 2)
        generating_mesh = generator.GetMesh()

        mesh = chaste.mesh.NodesOnlyMesh[2]()

        mesh.ConstructNodesWithoutMesh(generating_mesh, 1.5)

        transit_type = chaste.cell_based.TransitCellProliferativeType()
        cell_generator = chaste.cell_based.CellsGenerator["UniformCellCycleModel", 2]()
        cells = cell_generator.GenerateBasicRandom(mesh.GetNumNodes(), transit_type)

        cell_population = chaste.cell_based.NodeBasedCellPopulation[2](mesh, cells)

        scene = chaste.visualization.VtkScene[2]()
        scene.SetCellPopulation(cell_population)
        # JUPYTER_SHOW_FIRST
        scene.Start()  # JUPYTER_SHOW

        simulator = chaste.cell_based.OffLatticeSimulation[2, 2](cell_population)
        simulator.SetOutputDirectory("Python/TestNodeBasedCellSimulationsTutorial")
        simulator.SetSamplingTimestepMultiple(100)
        simulator.SetEndTime(10.0)

        force = chaste.cell_based.GeneralisedLinearSpringForce[2, 2]()
        simulator.AddForce(force)

        scene_modifier = chaste.cell_based.VtkSceneModifier[2]()
        scene_modifier.SetVtkScene(scene)
        scene_modifier.SetUpdateFrequency(100)
        simulator.AddSimulationModifier(scene_modifier)

        scene.Start()
        simulator.Solve()
        scene.End()

        self.assertEqual(cell_population.GetNumRealCells(), 8)
        self.assertAlmostEqual(chaste.cell_based.SimulationTime.Instance().GetTime(), 10.0, 6)

        # JUPYTER_TEARDOWN

    def test_spheroid(self):

        # JUPYTER_SETUP

        chaste.core.OutputFileHandler(
            "Python/TestNodeBasedCellSimulationsSpheroidTutorial"
        )
        nodes = []
        nodes.append(chaste.mesh.Node[3](0, False, 0.5, 0.0, 0.0))
        nodes.append(chaste.mesh.Node[3](1, False, -0.5, 0.0, 0.0))
        nodes.append(chaste.mesh.Node[3](2, False, 0.0, 0.5, 0.0))
        nodes.append(chaste.mesh.Node[3](3, False, 0.0, -0.5, 0.0))

        mesh = chaste.mesh.NodesOnlyMesh[3]()

        mesh.ConstructNodesWithoutMesh(nodes, 1.5)

        transit_type = chaste.cell_based.TransitCellProliferativeType()
        cell_generator = chaste.cell_based.CellsGenerator["UniformCellCycleModel", 3]()
        cells = cell_generator.GenerateBasicRandom(mesh.GetNumNodes(), transit_type)

        cell_population = chaste.cell_based.NodeBasedCellPopulation[3](mesh, cells)

        scene = chaste.visualization.VtkScene[3]()
        scene.SetCellPopulation(cell_population)
        scene.Start()  # JUPYTER_SHOW

        simulator = chaste.cell_based.OffLatticeSimulation[3, 3](cell_population)
        simulator.SetOutputDirectory("Python/TestNodeBasedCellSimulationsSpheroidTutorial")
        simulator.SetSamplingTimestepMultiple(12)
        simulator.SetEndTime(10.0)

        force = chaste.cell_based.GeneralisedLinearSpringForce[3, 3]()
        simulator.AddForce(force)

        scene_modifier = chaste.cell_based.VtkSceneModifier[3]()
        scene_modifier.SetVtkScene(scene)
        scene_modifier.SetUpdateFrequency(100)
        simulator.AddSimulationModifier(scene_modifier)

        scene.Start()
        simulator.Solve()
        scene.End()

        self.assertEqual(cell_population.GetNumRealCells(), 8)
        self.assertAlmostEqual(chaste.cell_based.SimulationTime.Instance().GetTime(), 10.0, 6)

        # JUPYTER_TEARDOWN

    def test_spheroid_on_sphere(self):

        # JUPYTER_SETUP

        chaste.core.OutputFileHandler("Python/TestNodeBasedCellSimulationsRestrictedSpheroidTutorial")
        nodes = []
        nodes.append(chaste.mesh.Node[3](0, False, 0.5, 0.0, 0.0))
        nodes.append(chaste.mesh.Node[3](1, False, -0.5, 0.0, 0.0))
        nodes.append(chaste.mesh.Node[3](2, False, 0.0, 0.5, 0.0))
        nodes.append(chaste.mesh.Node[3](3, False, 0.0, -0.5, 0.0))
        mesh = chaste.mesh.NodesOnlyMesh[3]()

        mesh.ConstructNodesWithoutMesh(nodes, 1.5)

        transit_type = chaste.cell_based.TransitCellProliferativeType()
        cell_generator = chaste.cell_based.CellsGenerator["UniformCellCycleModel", 3]()
        cells = cell_generator.GenerateBasicRandom(mesh.GetNumNodes(), transit_type)
        cell_population = chaste.cell_based.NodeBasedCellPopulation[3](mesh, cells)

        scene = chaste.visualization.VtkScene[3]()
        scene.SetCellPopulation(cell_population)
        scene.Start()  # JUPYTER_SHOW

        simulator = chaste.cell_based.OffLatticeSimulation[3, 3](cell_population)
        simulator.SetOutputDirectory("Python/TestNodeBasedCellSimulationsRestrictedSpheroidTutorial")
        simulator.SetSamplingTimestepMultiple(12)
        simulator.SetEndTime(10.0)

        force = chaste.cell_based.GeneralisedLinearSpringForce[3, 3]()
        simulator.AddForce(force)

        centre = np.array([0.0, 0.0, 1.0])
        radius = 5.0
        point2 = chaste.mesh.ChastePoint[3](centre)
        boundary_condition = chaste.cell_based.SphereGeometryBoundaryCondition[3](cell_population, point2.rGetLocation(), radius)
        simulator.AddCellPopulationBoundaryCondition(boundary_condition)

        scene_modifier.SetVtkScene(scene)
        scene_modifier.SetUpdateFrequency(100)
        simulator.AddSimulationModifier(scene_modifier)

        scene.Start()
        simulator.Solve()
        scene.End()

        self.assertEqual(cell_population.GetNumRealCells(), 8)
        self.assertAlmostEqual(chaste.cell_based.SimulationTime.Instance().GetTime(), 10.0, 6)

        # JUPYTER_TEARDOWN

if __name__ == "__main__":
    unittest.main(verbosity=2)
```
