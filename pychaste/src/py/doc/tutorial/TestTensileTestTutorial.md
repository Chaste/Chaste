
---
title : "Test Tensile Test Tutorial"
summary: ""
draft: false
images: []
toc: true
layout: "single"
---

This tutorial is automatically generated from [TestTensileTestTutorial](https://github.com/Chaste/PyChaste/blob/develop/test/python/cell_based/tutorials/TestTensileTestTutorial.py) at revision [f810861a](https://github.com/Chaste/PyChaste/commit/f810861afe376ba19bd791e14e85f29583993205).
Note that the code is given in full at the bottom of the page.


## Introduction
In this tutorial we will demonstrate a simulated tensile test on an epithelial sheet. This test
demonstrates:
 * Working with vertex based off lattice populations
 * Applying boundary conditions
 * Working with forces
 
## The Test

```python
import unittest  # Python testing framework
import numpy as np  # Matrix tools
import chaste  # The PyChaste module
import chaste.mesh  # Contains meshes
import chaste.cell_based  # Contains cell populations
import chaste.visualization  # Visualization tools
chaste.init()  # Set up MPI

class TestTensileTestTutorial(chaste.cell_based.AbstractCellBasedTestSuite):

```
### Test 1 - A 2D test

```python
    def test_monolayer(self):

        # JUPYTER_SETUP

```
First, we generate a vertex mesh using a HoneycombVertexMeshGenerator.

```python
        generator = chaste.mesh.HoneycombVertexMeshGenerator(5, 15)
        mesh = generator.GetMesh()

```
Now set up the cells, again we want to avoid proliferation.

```python
        differentiated_type = chaste.cell_based.DifferentiatedCellProliferativeType()
        cell_generator = chaste.cell_based.CellsGeneratorUniformG1GenerationalCellCycleModel_2()
        cells = cell_generator.GenerateBasicRandom(mesh.GetNumElements(), differentiated_type)

```
Next, create the cell population

```python
        cell_population = chaste.cell_based.VertexBasedCellPopulation2(mesh,
                                                                       cells)

```
Pass the cell population into an `OffLatticeSimulation`, and set the output directory, output multiple and end time

```python
        simulator = chaste.cell_based.OffLatticeSimulation2_2(cell_population)
        simulator.SetOutputDirectory("Python/TestTensileTest")
        simulator.SetEndTime(1.0)
        simulator.SetSamplingTimestepMultiple(1000)

```
Now create a force law

```python
        force = chaste.cell_based.NagaiHondaForce2()
        simulator.AddForce(force)

```
A `NagaiHondaForce` assumes that each cell has a target area. The target areas of cells are used to determine
pressure forces on each vertex and eventually determine the size of each cell in the simulation.
In order to assign target areas to cells and update them in each time step we add a `SimpleTargetAreaModifier`
to the simulation, which inherits from `AbstractTargetAreaModifier`.

```python
        growth_modifier = chaste.cell_based.SimpleTargetAreaModifier2()
        simulator.AddSimulationModifier(growth_modifier)

```
For our tensile test we will fix the bottom of the sheet and subject the top to an applied displacement. We neglect
fixing lateral degress of freedom for simplicity, since we are using an over-damped mechanical model.

```python
        my_point = np.array([0.0, 0.0])
        normal = np.array([0.0, -1.0])
        bc = chaste.cell_based.AttractingPlaneBoundaryCondition2_2(cell_population, my_point, normal)
        simulator.AddCellPopulationBoundaryCondition(bc)
        point = np.array([0.0, 15.5])
        normal = np.array([0.0, -1.0])
        bc2 = chaste.cell_based.AttractingPlaneBoundaryCondition2_2(cell_population, point, normal)
        simulator.AddCellPopulationBoundaryCondition(bc2)

```
We want to displace our top boundary over time. We could write a custom boundary condition class to do this.
A more simple alternative is to modify the the position of the point describing our boundary plane in `bc2`
as the simulation progresses. As per earlier tutorials we make a new `SimulationModifier` class to do this.

```python
        class BoundaryConditionModifier(chaste.cell_based.PythonSimulationModifier2):

            """ Class for time varying boundary conditions
            """

            def __init__(self, boundary_condition):
                self.boundary_condition = boundary_condition
                self.original_location = boundary_condition.rGetPointOnPlane()
                self.velocity = 0.5  # cell lengths per time
                super(BoundaryConditionModifier, self).__init__()

            def UpdateAtEndOfTimeStep(self, cell_population):

                """ Move the boundary upwards at the specified velocity
                """

                total_time = chaste.cell_based.SimulationTime.Instance().GetTime()
                new_location = [self.original_location[0],
                                self.original_location[1] + self.velocity*total_time]
                self.boundary_condition.SetPointOnPlane(np.array(new_location))

            def SetupSolve(self, cell_population, output_directory):

                """ Make sure the cell population is in the correct state at the start of the simulation
                """

                cell_population.Update()

        bc_modifier = BoundaryConditionModifier(bc2)
        simulator.AddSimulationModifier(bc_modifier)

```
PyChaste can do simple 3D rendering with VTK. We set up a `VtkScene` so that we can see the population
evovle in real time.

```python
        scene = chaste.visualization.VtkScene2()
        scene.SetCellPopulation(cell_population)
        # JUPYTER_SHOW_FIRST

        scene_modifier = chaste.cell_based.VtkSceneModifier2()
        scene_modifier.SetVtkScene(scene)
        scene_modifier.SetUpdateFrequency(1000)
        simulator.AddSimulationModifier(scene_modifier)

```
To run the simulation, we call `Solve()`.

```python
        scene.Start()
        simulator.Solve()

        # JUPYTER_TEARDOWN

if __name__ == '__main__':
    unittest.main(verbosity=2)

```


## Full code 


**File name:** `TestTensileTestTutorial.py` 

```python
import unittest  # Python testing framework
import numpy as np  # Matrix tools
import chaste  # The PyChaste module
import chaste.mesh  # Contains meshes
import chaste.cell_based  # Contains cell populations
import chaste.visualization  # Visualization tools
chaste.init()  # Set up MPI

class TestTensileTestTutorial(chaste.cell_based.AbstractCellBasedTestSuite):

    def test_monolayer(self):

        # JUPYTER_SETUP

        generator = chaste.mesh.HoneycombVertexMeshGenerator(5, 15)
        mesh = generator.GetMesh()

        differentiated_type = chaste.cell_based.DifferentiatedCellProliferativeType()
        cell_generator = chaste.cell_based.CellsGeneratorUniformG1GenerationalCellCycleModel_2()
        cells = cell_generator.GenerateBasicRandom(mesh.GetNumElements(), differentiated_type)

        cell_population = chaste.cell_based.VertexBasedCellPopulation2(mesh,
                                                                       cells)

        simulator = chaste.cell_based.OffLatticeSimulation2_2(cell_population)
        simulator.SetOutputDirectory("Python/TestTensileTest")
        simulator.SetEndTime(1.0)
        simulator.SetSamplingTimestepMultiple(1000)

        force = chaste.cell_based.NagaiHondaForce2()
        simulator.AddForce(force)

        growth_modifier = chaste.cell_based.SimpleTargetAreaModifier2()
        simulator.AddSimulationModifier(growth_modifier)

        my_point = np.array([0.0, 0.0])
        normal = np.array([0.0, -1.0])
        bc = chaste.cell_based.AttractingPlaneBoundaryCondition2_2(cell_population, my_point, normal)
        simulator.AddCellPopulationBoundaryCondition(bc)
        point = np.array([0.0, 15.5])
        normal = np.array([0.0, -1.0])
        bc2 = chaste.cell_based.AttractingPlaneBoundaryCondition2_2(cell_population, point, normal)
        simulator.AddCellPopulationBoundaryCondition(bc2)

        class BoundaryConditionModifier(chaste.cell_based.PythonSimulationModifier2):

            """ Class for time varying boundary conditions
            """

            def __init__(self, boundary_condition):
                self.boundary_condition = boundary_condition
                self.original_location = boundary_condition.rGetPointOnPlane()
                self.velocity = 0.5  # cell lengths per time
                super(BoundaryConditionModifier, self).__init__()

            def UpdateAtEndOfTimeStep(self, cell_population):

                """ Move the boundary upwards at the specified velocity
                """

                total_time = chaste.cell_based.SimulationTime.Instance().GetTime()
                new_location = [self.original_location[0],
                                self.original_location[1] + self.velocity*total_time]
                self.boundary_condition.SetPointOnPlane(np.array(new_location))

            def SetupSolve(self, cell_population, output_directory):

                """ Make sure the cell population is in the correct state at the start of the simulation
                """

                cell_population.Update()

        bc_modifier = BoundaryConditionModifier(bc2)
        simulator.AddSimulationModifier(bc_modifier)

        scene = chaste.visualization.VtkScene2()
        scene.SetCellPopulation(cell_population)
        # JUPYTER_SHOW_FIRST

        scene_modifier = chaste.cell_based.VtkSceneModifier2()
        scene_modifier.SetVtkScene(scene)
        scene_modifier.SetUpdateFrequency(1000)
        simulator.AddSimulationModifier(scene_modifier)

        scene.Start()
        simulator.Solve()

        # JUPYTER_TEARDOWN

if __name__ == '__main__':
    unittest.main(verbosity=2)

```

