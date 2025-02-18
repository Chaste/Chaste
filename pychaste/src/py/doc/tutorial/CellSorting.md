---
title : "Cell Sorting"
summary: ""
draft: false
images: []
toc: true
layout: "single"
---
This tutorial is automatically generated from [TestPyCellSortingTutorial.py](https://github.com/Chaste/Chaste/blob/develop/pychaste/test/tutorial/TestPyCellSortingTutorial.py) at revision [4045f91a83f5](https://github.com/Chaste/Chaste/commit/4045f91a83f55dc4a97f2ca4f97b0c32f4e43a4a).

Note that the code is given in full at the bottom of the page.

## Introduction
This test is a demonstration of cell sorting using a Cellular Potts based framework.
It shows:
 * How to set up a Potts simulation
 * Working with labels
 
## The Test

```python
import unittest  # Python testing framework

import chaste  # The PyChaste module
import chaste.cell_based  # Contains cell populations
import chaste.mesh  # Contains meshes
import chaste.visualization  # Visualization tools
import matplotlib.pyplot as plt  # Plotting
import numpy as np  # Matrix tools

class TestPyCellSortingTutorial(chaste.cell_based.AbstractCellBasedTestSuite):
```
### Test 1 - Cell sorting
The next test generates a collection of cells, there are two types of cells, labelled ones and non labelled ones,
there is differential adhesion between the cell types. For the parameters specified, the cells sort into separate types.

```python
    def test_potts_monolayer_cell_sorting(self):

        # JUPYTER_SETUP
```
First, we generate a `Potts` mesh. To create a `PottsMesh`, we can use the `PottsMeshGenerator`.
This generates a regular square-shaped mesh, in which all elements are the same size.
We have chosen an 8 by 8 block of elements each consisting of 4 by 4 ( = 16) lattice sites.

```python
        generator = chaste.mesh.PottsMeshGenerator[2](50, 8, 4, 50, 8, 4)
        mesh = generator.GetMesh()
```
Having created a mesh, we now create some cells. To do this, we the `CellsGenerator` helper class,
as before but this time the third argument is set to make all cells non-proliferative.

```python
        differentiated_type = chaste.cell_based.DifferentiatedCellProliferativeType()
        cell_generator = chaste.cell_based.CellsGenerator["UniformCellCycleModel", 2]()
        cells = cell_generator.GenerateBasicRandom(mesh.GetNumElements(), differentiated_type)
```
Before we make a CellPopulation we make a cell label and then assign this label to some randomly chosen cells.

```python
        label = chaste.cell_based.CellLabel()
        for eachCell in cells:
            if chaste.core.RandomNumberGenerator.Instance().ranf() < 0.5:
                eachCell.AddCellProperty(label)
```
Now we have a mesh and a set of cells to go with it, we can create a `CellPopulation`.

```python
        cell_population = chaste.cell_based.PottsBasedCellPopulation[2](mesh, cells)
```
In order to visualize labelled cells we need to use the following command.

```python
        cell_population.AddCellWriterCellLabelWriter()
```
PyChaste can do simple 3D rendering with VTK. We set up a VtkScene so that we can
see the population evovle in real time.

```python
        scene = chaste.visualization.VtkScene[2]()
        scene.SetCellPopulation(cell_population)
        scene.GetCellPopulationActorGenerator().SetShowPottsMeshEdges(True)
        # JUPYTER_SHOW_FIRST
        scene.Start()  # JUPYTER_SHOW
```
We then pass in the cell population into an `OffLatticeSimulation`, and set the output directory and end time

```python
        simulator = chaste.cell_based.OnLatticeSimulation[2](cell_population)
        simulator.SetOutputDirectory("Python/TestCellSorting")
        simulator.SetEndTime(20.0)
        simulator.SetSamplingTimestepMultiple(10)
```
We must now create one or more update rules, which determine the Hamiltonian in the Potts simulation.
For this test, we use two update rules based upon a volume constraint (`VolumeConstraintPottsUpdateRule`) and
differential adhesion between cells (`DifferentialAdhesionPottsUpdateRule`), set appropriate parameters, and
pass them to the `OnLatticeSimulation`.

```python
        volume_constraint_update_rule = chaste.cell_based.VolumeConstraintPottsUpdateRule[2]()
        volume_constraint_update_rule.SetMatureCellTargetVolume(16)
        volume_constraint_update_rule.SetDeformationEnergyParameter(0.2)
        simulator.AddUpdateRule(volume_constraint_update_rule)
```
We repeat the process for any other update rules.

```python
        differential_adhesion_update_rule = chaste.cell_based.DifferentialAdhesionPottsUpdateRule[2]()
        differential_adhesion_update_rule.SetLabelledCellLabelledCellAdhesionEnergyParameter(0.16)
        differential_adhesion_update_rule.SetLabelledCellCellAdhesionEnergyParameter(0.11)
        differential_adhesion_update_rule.SetCellCellAdhesionEnergyParameter(0.02)
        differential_adhesion_update_rule.SetLabelledCellBoundaryAdhesionEnergyParameter(0.16)
        differential_adhesion_update_rule.SetCellBoundaryAdhesionEnergyParameter(0.16)
        simulator.AddUpdateRule(differential_adhesion_update_rule)
```
Set up plotting

```python
        scene_modifier = chaste.cell_based.VtkSceneModifier[2]()
        scene_modifier.SetVtkScene(scene)
        scene_modifier.SetUpdateFrequency(1000)
        simulator.AddSimulationModifier(scene_modifier)
```
To run the simulation, we call `Solve()`.

```python
        scene.Start()
        simulator.Solve()

        # JUPYTER_TEARDOWN

if __name__ == "__main__":
    unittest.main(verbosity=2)
```

## Full code

```python
import unittest  # Python testing framework

import chaste  # The PyChaste module
import chaste.cell_based  # Contains cell populations
import chaste.mesh  # Contains meshes
import chaste.visualization  # Visualization tools
import matplotlib.pyplot as plt  # Plotting
import numpy as np  # Matrix tools

class TestPyCellSortingTutorial(chaste.cell_based.AbstractCellBasedTestSuite):

    def test_potts_monolayer_cell_sorting(self):

        # JUPYTER_SETUP

        generator = chaste.mesh.PottsMeshGenerator[2](50, 8, 4, 50, 8, 4)
        mesh = generator.GetMesh()

        differentiated_type = chaste.cell_based.DifferentiatedCellProliferativeType()
        cell_generator = chaste.cell_based.CellsGenerator["UniformCellCycleModel", 2]()
        cells = cell_generator.GenerateBasicRandom(mesh.GetNumElements(), differentiated_type)

        label = chaste.cell_based.CellLabel()
        for eachCell in cells:
            if chaste.core.RandomNumberGenerator.Instance().ranf() < 0.5:
                eachCell.AddCellProperty(label)

        cell_population = chaste.cell_based.PottsBasedCellPopulation[2](mesh, cells)

        cell_population.AddCellWriterCellLabelWriter()

        scene = chaste.visualization.VtkScene[2]()
        scene.SetCellPopulation(cell_population)
        scene.GetCellPopulationActorGenerator().SetShowPottsMeshEdges(True)
        # JUPYTER_SHOW_FIRST
        scene.Start()  # JUPYTER_SHOW

        simulator = chaste.cell_based.OnLatticeSimulation[2](cell_population)
        simulator.SetOutputDirectory("Python/TestCellSorting")
        simulator.SetEndTime(20.0)
        simulator.SetSamplingTimestepMultiple(10)

        volume_constraint_update_rule = chaste.cell_based.VolumeConstraintPottsUpdateRule[2]()
        volume_constraint_update_rule.SetMatureCellTargetVolume(16)
        volume_constraint_update_rule.SetDeformationEnergyParameter(0.2)
        simulator.AddUpdateRule(volume_constraint_update_rule)

        differential_adhesion_update_rule = chaste.cell_based.DifferentialAdhesionPottsUpdateRule[2]()
        differential_adhesion_update_rule.SetLabelledCellLabelledCellAdhesionEnergyParameter(0.16)
        differential_adhesion_update_rule.SetLabelledCellCellAdhesionEnergyParameter(0.11)
        differential_adhesion_update_rule.SetCellCellAdhesionEnergyParameter(0.02)
        differential_adhesion_update_rule.SetLabelledCellBoundaryAdhesionEnergyParameter(0.16)
        differential_adhesion_update_rule.SetCellBoundaryAdhesionEnergyParameter(0.16)
        simulator.AddUpdateRule(differential_adhesion_update_rule)

        scene_modifier = chaste.cell_based.VtkSceneModifier[2]()
        scene_modifier.SetVtkScene(scene)
        scene_modifier.SetUpdateFrequency(1000)
        simulator.AddSimulationModifier(scene_modifier)

        scene.Start()
        simulator.Solve()

        # JUPYTER_TEARDOWN

if __name__ == "__main__":
    unittest.main(verbosity=2)
```
