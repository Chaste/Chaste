
---
title : "Test Potts Based Cell Simulations Python Tutorial"
summary: ""
draft: false
images: []
toc: true
layout: "single"
---

This tutorial is automatically generated from [TestPottsBasedCellSimulationsPythonTutorial](https://github.com/Chaste/PyChaste/blob/develop/test/python/cell_based/tutorials/TestPottsBasedCellSimulationsPythonTutorial.py) at revision [f810861a](https://github.com/Chaste/PyChaste/commit/f810861afe376ba19bd791e14e85f29583993205).
Note that the code is given in full at the bottom of the page.


## Introduction
In this tutorial we show how Chaste can be used to create, run and visualize Potts-based simulations.
Full details of the mathematical model can be found in Graner, F. and Glazier, J. A. (1992).

## The Test

```python
import unittest  # Python testing framework
import chaste  # The PyChaste module
import chaste.cell_based  # Contains cell populations
import chaste.mesh  # Contains meshes
import chaste.visualization  # Visualization tools
chaste.init()  # Set up MPI

class TestRunningPottsBasedSimulationsTutorial(chaste.cell_based.AbstractCellBasedTestSuite):
```
### Test 1 - A basic node-based simulation
In the first test, we run a simple Potts-based simulation, in which we create a monolayer of cells, using a Potts mesh.
Each cell is assigned a stochastic cell-cycle model.

```python
    def test_monolayer(self):

        # JUPYTER_SETUP

```
First, we generate a Potts mesh. To create a PottsMesh, we can use the PottsMeshGenerator.
This generates a regular square-shaped mesh, in which all elements are the same size.
Here the first three arguments specify the domain width; the number of elements across; and the width of elements.
The second set of three arguments specify the domain height; the number of elements up; and the height of individual elements.
We have chosen a 2 by 2 block of elements, each consisting of 4 by 4 ( = 16) lattice sites.

```python
        chaste.core.OutputFileHandler("Python/TestPottsBasedCellSimulationsTutorial")
        generator = chaste.mesh.PottsMeshGenerator2(50, 2, 4,
                                                    50, 2, 4)
        mesh = generator.GetMesh()

```
Having created a mesh, we now create a vector of CellPtrs. To do this, we the CellsGenerator helper class,
which is templated over the type of cell model required and the dimension.
We create an empty vector of cells and pass this into the method along with the mesh.
The second argument represents the size of that the vector cells should become - one cell for each element.
Third argument makes all cells proliferate.

```python
        cell_generator = chaste.cell_based.CellsGeneratorUniformCellCycleModel_2()
        cells = cell_generator.GenerateBasic(mesh.GetNumElements())

```
Now we have a mesh and a set of cells to go with it, we can create a CellPopulation.
In general, this class associates a collection of cells with a mesh. For this test, because we have a PottsMesh,
we use a particular type of cell population called a PottsBasedCellPopulation.

```python
        cell_population = chaste.cell_based.PottsBasedCellPopulation2(mesh,
                                                                      cells)

```
We can set the "Temperature" to be used in the Potts Simulation using the optional command below. The default value is 0.1.

```python
        cell_population.SetTemperature(0.1)

```
By default the Potts simulation will make 1 sweep over the whole domain per timestep.
To use a different number of sweeps per timestep use the command.

```python
        cell_population.SetNumSweepsPerTimestep(1)

```
We can set up a `VtkScene` to do a quick visualization of the population before running the analysis.

```python
        scene = chaste.visualization.VtkScene2()
        scene.SetCellPopulation(cell_population)
        scene.GetCellPopulationActorGenerator().SetShowPottsMeshEdges(True)
        # JUPYTER_SHOW_FIRST
        scene.Start()  # JUPYTER_SHOW

```
We then pass in the cell population into an `OffLatticeSimulation`, and set the output directory and end time

```python
        simulator = chaste.cell_based.OnLatticeSimulation2(cell_population)
        simulator.SetOutputDirectory("Python/TestPottsBasedCellSimulationsTutorial")
        simulator.SetEndTime(50.0)

```
The default timestep is 0.1, but can be changed using the below command. The timestep is used in conjunction with the "Temperature"
and number of sweeps per timestep to specify the relationship between cell movement and proliferation.
We also set the simulation to only output every 10 steps i.e. once per hour.

```python
        simulator.SetDt(0.1)
        simulator.SetSamplingTimestepMultiple(10)

```
We must now create one or more update rules, which determine the Hamiltonian in the Potts simulation.
For this test, we use two update rules based upon a volume constraint (VolumeConstraintPottsUpdateRule)
and adhesion between cells (AdhesionPottsUpdateRule) and pass them to the OnLatticeSimulation.
For a list of possible update rules see subclasses of AbstractPottsUpdateRule.

```python
        volume_constraint_update_rule = chaste.cell_based.VolumeConstraintPottsUpdateRule2()

```
Set an appropriate target volume in number of lattice sites. Here we use the default value of 16 lattice sites.

```python
        volume_constraint_update_rule.SetMatureCellTargetVolume(16)

```
You can also vary the deformation energy parameter.
The larger the parameter the more cells will try to maintain target volume. Here we use the default value of 0.2.

```python
        volume_constraint_update_rule.SetDeformationEnergyParameter(0.2)

```
Finally we add the update rule to the simulator.

```python
        simulator.AddUpdateRule(volume_constraint_update_rule)

```
We repeat the process for any other update rules.

```python
        adhesion_update_rule = chaste.cell_based.AdhesionPottsUpdateRule2()
        simulator.AddUpdateRule(adhesion_update_rule)

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

```
The next two lines are for test purposes only and are not part of this tutorial.
If different simulation input parameters are being explored the lines should be removed.

```python
        self.assertEqual(cell_population.GetNumRealCells(), 41)
        self.assertAlmostEqual(chaste.cell_based.SimulationTime.Instance().GetTime(),
                               50.0, 6)

        # JUPYTER_TEARDOWN

```
### Test 2 - Cell sorting
The next test generates a collection of cells, there are two types of cells, labelled ones and non labelled ones,
there is differential adhesion between the cell types. For the parameters specified, the cells sort into separate types.

```python
    def test_potts_monolayer_cell_sorting(self):

        # JUPYTER_SETUP

```
First, we generate a Potts mesh. To create a PottsMesh, we can use the PottsMeshGenerator.
This generates a regular square-shaped mesh, in which all elements are the same size.
We have chosen an 8 by 8 block of elements each consisting of 4 by 4 ( = 16) lattice sites.

```python
        generator = chaste.mesh.PottsMeshGenerator2(50, 8, 4,
                                                    50, 8, 4)
        mesh = generator.GetMesh()

```
Having created a mesh, we now create a VectorSharedPtrCells. To do this, we the CellsGenerator helper class,
as before but this time the third argument is set to make all cells non-proliferative.

```python
        differentiated_type = chaste.cell_based.DifferentiatedCellProliferativeType()
        cell_generator = chaste.cell_based.CellsGeneratorUniformCellCycleModel_2()
        cells = cell_generator.GenerateBasicRandom(mesh.GetNumElements(),
                                              differentiated_type)

```
Before we make a CellPopulation we make a cell label and then assign this label to some randomly chosen cells.

```python
        label = chaste.cell_based.CellLabel()
        for eachCell in cells:
            if(chaste.core.RandomNumberGenerator.Instance().ranf() < 0.5):
                eachCell.AddCellProperty(label)

```
Now we have a mesh and a set of cells to go with it, we can create a CellPopulation.

```python
        cell_population = chaste.cell_based.PottsBasedCellPopulation2(mesh,
                                                                      cells)

```
In order to visualize labelled cells we need to use the following command.

```python
        cell_population.AddCellWriterCellLabelWriter()

```
We then pass in the cell population into an `OffLatticeSimulation`, and set the output directory and end time

```python
        simulator = chaste.cell_based.OnLatticeSimulation2(cell_population)
        simulator.SetOutputDirectory("Python/TestPottsBasedCellSorting")
        simulator.SetEndTime(20.0)
        simulator.SetSamplingTimestepMultiple(10)

```
We must now create one or more update rules, which determine the Hamiltonian in the Potts simulation.
For this test, we use two update rules based upon a volume constraint (VolumeConstraintPottsUpdateRule) and
differential adhesion between cells (DifferentialAdhesionPottsUpdateRule), set appropriate parameters, and
pass them to the OnLatticeSimulation.

```python
        volume_constraint_update_rule = chaste.cell_based.VolumeConstraintPottsUpdateRule2()
        volume_constraint_update_rule.SetMatureCellTargetVolume(16)
        volume_constraint_update_rule.SetDeformationEnergyParameter(0.2)
        simulator.AddUpdateRule(volume_constraint_update_rule)

```
We repeat the process for any other update rules.

```python
        differential_adhesion_update_rule = chaste.cell_based.DifferentialAdhesionPottsUpdateRule2()
        differential_adhesion_update_rule.SetLabelledCellLabelledCellAdhesionEnergyParameter(0.16)
        differential_adhesion_update_rule.SetLabelledCellCellAdhesionEnergyParameter(0.11)
        differential_adhesion_update_rule.SetCellCellAdhesionEnergyParameter(0.02)
        differential_adhesion_update_rule.SetLabelledCellBoundaryAdhesionEnergyParameter(0.16)
        differential_adhesion_update_rule.SetCellBoundaryAdhesionEnergyParameter(0.16)
        simulator.AddUpdateRule(differential_adhesion_update_rule)

```
To run the simulation, we call `Solve()`.

```python
        simulator.Solve()

```
The next two lines are for test purposes only and are not part of this tutorial.
If different simulation input parameters are being explored the lines should be removed.

```python
        self.assertEqual(cell_population.GetNumRealCells(), 64)
        self.assertAlmostEqual(chaste.cell_based.SimulationTime.Instance().GetTime(),
                               20.0, 6)

        # JUPYTER_TEARDOWN

```
### Test 3 - 3D Cell Sorting
The next test extends the previous example to three dimensions.

```python
    def test_potts_spheroid_cell_sorting(self):

        # JUPYTER_SETUP

```
First, we generate a Potts mesh. To create a PottsMesh, we can use the PottsMeshGenerator.
This generates a regular square-shaped mesh, in which all elements are the same size.
Here the first three arguments specify the domain width; the number of elements across; and the width of elements.
The second set of three arguments specify the domain height; the number of elements up; and the height of individual elements.
The third set of three arguments specify the domain depth; the number of elements deep; and the depth of individual elements.
We have chosen an 4 by 4 by 4 ( = 64) block of elements each consisting of 2 by 2 by 2 ( = 8) lattice sites.

```python
        generator = chaste.mesh.PottsMeshGenerator3(10, 4, 2,
                                                    10, 4, 2,
                                                    10, 4, 2)
        mesh = generator.GetMesh()

```
Having created a mesh, we now create a VectorSharedPtrCells. To do this, we the CellsGenerator helper class,
as before but this time the third argument is set to make all cells non-proliferative.

```python
        differentiated_type = chaste.cell_based.DifferentiatedCellProliferativeType()
        cell_generator = chaste.cell_based.CellsGeneratorUniformCellCycleModel_3()
        cells = cell_generator.GenerateBasicRandom(mesh.GetNumElements(),
                                                   differentiated_type)

```
As for the 2D case before we make a CellPopulation we make a pointer to a cell label and then assign this label to some randomly chosen cells.

```python
        label = chaste.cell_based.CellLabel()
        for eachCell in cells:
            if(chaste.core.RandomNumberGenerator.Instance().ranf() < 0.5):
                eachCell.AddCellProperty(label)

```
Now we have a mesh and a set of cells to go with it, we can create a CellPopulation.

```python
        cell_population = chaste.cell_based.PottsBasedCellPopulation3(mesh,
                                                                      cells)

```
In order to visualize labelled cells we need to use the following command.

```python
        cell_population.AddCellWriterCellLabelWriter()

```
We then pass in the cell population into an `OffLatticeSimulation`, and set the output directory and end time

```python
        simulator = chaste.cell_based.OnLatticeSimulation3(cell_population)
        simulator.SetOutputDirectory("Python/TestPottsBasedCellSorting3D")
        simulator.SetEndTime(20.0)
        simulator.SetSamplingTimestepMultiple(10)

```
We must now create one or more update rules, which determine the Hamiltonian in the Potts simulation.
Now set the target volume to be appropriate for this 3D simulation.

```python
        volume_constraint_update_rule = chaste.cell_based.VolumeConstraintPottsUpdateRule3()
        volume_constraint_update_rule.SetMatureCellTargetVolume(8)
        volume_constraint_update_rule.SetDeformationEnergyParameter(0.2)
        simulator.AddUpdateRule(volume_constraint_update_rule)

```
We use the same differential adhesion parameters as in the 2D case.

```python
        differential_adhesion_update_rule = chaste.cell_based.DifferentialAdhesionPottsUpdateRule3()
        differential_adhesion_update_rule.SetLabelledCellLabelledCellAdhesionEnergyParameter(0.16)
        differential_adhesion_update_rule.SetLabelledCellCellAdhesionEnergyParameter(0.11)
        differential_adhesion_update_rule.SetCellCellAdhesionEnergyParameter(0.02)
        differential_adhesion_update_rule.SetLabelledCellBoundaryAdhesionEnergyParameter(0.16)
        differential_adhesion_update_rule.SetCellBoundaryAdhesionEnergyParameter(0.16)
        simulator.AddUpdateRule(differential_adhesion_update_rule)

```
To run the simulation, we call `Solve()`.

```python
        simulator.Solve()

```
The next two lines are for test purposes only and are not part of this tutorial.
If different simulation input parameters are being explored the lines should be removed.

```python
        self.assertEqual(cell_population.GetNumRealCells(), 64)
        self.assertAlmostEqual(chaste.cell_based.SimulationTime.Instance().GetTime(),
                               20.0, 6)

        # JUPYTER_TEARDOWN

if __name__ == '__main__':
    unittest.main(verbosity=2)

```


## Full code 


**File name:** `TestPottsBasedCellSimulationsPythonTutorial.py` 

```python
import unittest  # Python testing framework
import chaste  # The PyChaste module
import chaste.cell_based  # Contains cell populations
import chaste.mesh  # Contains meshes
import chaste.visualization  # Visualization tools
chaste.init()  # Set up MPI

class TestRunningPottsBasedSimulationsTutorial(chaste.cell_based.AbstractCellBasedTestSuite):
    def test_monolayer(self):

        # JUPYTER_SETUP

        chaste.core.OutputFileHandler("Python/TestPottsBasedCellSimulationsTutorial")
        generator = chaste.mesh.PottsMeshGenerator2(50, 2, 4,
                                                    50, 2, 4)
        mesh = generator.GetMesh()

        cell_generator = chaste.cell_based.CellsGeneratorUniformCellCycleModel_2()
        cells = cell_generator.GenerateBasic(mesh.GetNumElements())

        cell_population = chaste.cell_based.PottsBasedCellPopulation2(mesh,
                                                                      cells)

        cell_population.SetTemperature(0.1)

        cell_population.SetNumSweepsPerTimestep(1)

        scene = chaste.visualization.VtkScene2()
        scene.SetCellPopulation(cell_population)
        scene.GetCellPopulationActorGenerator().SetShowPottsMeshEdges(True)
        # JUPYTER_SHOW_FIRST
        scene.Start()  # JUPYTER_SHOW

        simulator = chaste.cell_based.OnLatticeSimulation2(cell_population)
        simulator.SetOutputDirectory("Python/TestPottsBasedCellSimulationsTutorial")
        simulator.SetEndTime(50.0)

        simulator.SetDt(0.1)
        simulator.SetSamplingTimestepMultiple(10)

        volume_constraint_update_rule = chaste.cell_based.VolumeConstraintPottsUpdateRule2()

        volume_constraint_update_rule.SetMatureCellTargetVolume(16)

        volume_constraint_update_rule.SetDeformationEnergyParameter(0.2)

        simulator.AddUpdateRule(volume_constraint_update_rule)

        adhesion_update_rule = chaste.cell_based.AdhesionPottsUpdateRule2()
        simulator.AddUpdateRule(adhesion_update_rule)

        scene_modifier = chaste.cell_based.VtkSceneModifier2()
        scene_modifier.SetVtkScene(scene)
        scene_modifier.SetUpdateFrequency(100)
        simulator.AddSimulationModifier(scene_modifier)

        scene.Start()
        simulator.Solve()

        self.assertEqual(cell_population.GetNumRealCells(), 41)
        self.assertAlmostEqual(chaste.cell_based.SimulationTime.Instance().GetTime(),
                               50.0, 6)

        # JUPYTER_TEARDOWN

    def test_potts_monolayer_cell_sorting(self):

        # JUPYTER_SETUP

        generator = chaste.mesh.PottsMeshGenerator2(50, 8, 4,
                                                    50, 8, 4)
        mesh = generator.GetMesh()

        differentiated_type = chaste.cell_based.DifferentiatedCellProliferativeType()
        cell_generator = chaste.cell_based.CellsGeneratorUniformCellCycleModel_2()
        cells = cell_generator.GenerateBasicRandom(mesh.GetNumElements(),
                                              differentiated_type)

        label = chaste.cell_based.CellLabel()
        for eachCell in cells:
            if(chaste.core.RandomNumberGenerator.Instance().ranf() < 0.5):
                eachCell.AddCellProperty(label)

        cell_population = chaste.cell_based.PottsBasedCellPopulation2(mesh,
                                                                      cells)

        cell_population.AddCellWriterCellLabelWriter()

        simulator = chaste.cell_based.OnLatticeSimulation2(cell_population)
        simulator.SetOutputDirectory("Python/TestPottsBasedCellSorting")
        simulator.SetEndTime(20.0)
        simulator.SetSamplingTimestepMultiple(10)

        volume_constraint_update_rule = chaste.cell_based.VolumeConstraintPottsUpdateRule2()
        volume_constraint_update_rule.SetMatureCellTargetVolume(16)
        volume_constraint_update_rule.SetDeformationEnergyParameter(0.2)
        simulator.AddUpdateRule(volume_constraint_update_rule)

        differential_adhesion_update_rule = chaste.cell_based.DifferentialAdhesionPottsUpdateRule2()
        differential_adhesion_update_rule.SetLabelledCellLabelledCellAdhesionEnergyParameter(0.16)
        differential_adhesion_update_rule.SetLabelledCellCellAdhesionEnergyParameter(0.11)
        differential_adhesion_update_rule.SetCellCellAdhesionEnergyParameter(0.02)
        differential_adhesion_update_rule.SetLabelledCellBoundaryAdhesionEnergyParameter(0.16)
        differential_adhesion_update_rule.SetCellBoundaryAdhesionEnergyParameter(0.16)
        simulator.AddUpdateRule(differential_adhesion_update_rule)

        simulator.Solve()

        self.assertEqual(cell_population.GetNumRealCells(), 64)
        self.assertAlmostEqual(chaste.cell_based.SimulationTime.Instance().GetTime(),
                               20.0, 6)

        # JUPYTER_TEARDOWN

    def test_potts_spheroid_cell_sorting(self):

        # JUPYTER_SETUP

        generator = chaste.mesh.PottsMeshGenerator3(10, 4, 2,
                                                    10, 4, 2,
                                                    10, 4, 2)
        mesh = generator.GetMesh()

        differentiated_type = chaste.cell_based.DifferentiatedCellProliferativeType()
        cell_generator = chaste.cell_based.CellsGeneratorUniformCellCycleModel_3()
        cells = cell_generator.GenerateBasicRandom(mesh.GetNumElements(),
                                                   differentiated_type)

        label = chaste.cell_based.CellLabel()
        for eachCell in cells:
            if(chaste.core.RandomNumberGenerator.Instance().ranf() < 0.5):
                eachCell.AddCellProperty(label)

        cell_population = chaste.cell_based.PottsBasedCellPopulation3(mesh,
                                                                      cells)

        cell_population.AddCellWriterCellLabelWriter()

        simulator = chaste.cell_based.OnLatticeSimulation3(cell_population)
        simulator.SetOutputDirectory("Python/TestPottsBasedCellSorting3D")
        simulator.SetEndTime(20.0)
        simulator.SetSamplingTimestepMultiple(10)

        volume_constraint_update_rule = chaste.cell_based.VolumeConstraintPottsUpdateRule3()
        volume_constraint_update_rule.SetMatureCellTargetVolume(8)
        volume_constraint_update_rule.SetDeformationEnergyParameter(0.2)
        simulator.AddUpdateRule(volume_constraint_update_rule)

        differential_adhesion_update_rule = chaste.cell_based.DifferentialAdhesionPottsUpdateRule3()
        differential_adhesion_update_rule.SetLabelledCellLabelledCellAdhesionEnergyParameter(0.16)
        differential_adhesion_update_rule.SetLabelledCellCellAdhesionEnergyParameter(0.11)
        differential_adhesion_update_rule.SetCellCellAdhesionEnergyParameter(0.02)
        differential_adhesion_update_rule.SetLabelledCellBoundaryAdhesionEnergyParameter(0.16)
        differential_adhesion_update_rule.SetCellBoundaryAdhesionEnergyParameter(0.16)
        simulator.AddUpdateRule(differential_adhesion_update_rule)

        simulator.Solve()

        self.assertEqual(cell_population.GetNumRealCells(), 64)
        self.assertAlmostEqual(chaste.cell_based.SimulationTime.Instance().GetTime(),
                               20.0, 6)

        # JUPYTER_TEARDOWN

if __name__ == '__main__':
    unittest.main(verbosity=2)

```

