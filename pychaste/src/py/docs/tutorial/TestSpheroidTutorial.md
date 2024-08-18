
---
title : "Test Spheroid Tutorial"
summary: ""
draft: false
images: []
toc: true
layout: "single"
---

This tutorial is automatically generated from [TestSpheroidTutorial](https://github.com/Chaste/PyChaste/blob/develop/test/python/cell_based/tutorials/TestSpheroidTutorial.py) at revision [f810861a](https://github.com/Chaste/PyChaste/commit/f810861afe376ba19bd791e14e85f29583993205).
Note that the code is given in full at the bottom of the page.


## Introduction
This tutorial is an example of modelling spheroid growth with a nutrient.
It covers:
 * Setting up an off-lattice cell population
 * Setting up a cell cycle model with oxygen dependence
 * Setting up and solving an oxygen transport PDE
 * Setting up a cell killer
 
## The Test

```python
import unittest # Python testing framework
import matplotlib.pyplot as plt # Plotting
import numpy as np # Matrix tools
import chaste # The PyChaste module
import chaste.mesh # Contains meshes
import chaste.pde # PDEs
import chaste.cell_based # Contains cell populations
import chaste.visualization # Visualization tools
chaste.init() # Set up MPI

class TestSpheroidTutorial(chaste.cell_based.AbstractCellBasedTestSuite):

```
### Test 1 - a 2D mesh-based spheroid
In this test we set up a spheroid with a plentiful supply of oxygen on the boundary and watch it grow
over time. Cells can gradually become apoptotic if the oxygen tension is too low.

```python
    def test_spheroid(self):

        # JUPYTER_SETUP

```
This time we will use on off-lattice `MeshBased` cell population. Cell centres are joined with
springs with a Delauney Triangulation used to identify neighbours. Cell area is given by the dual
(Voronoi Tesselation). We start off with a small number of cells. We use a `MutableMesh` which
can change connectivity over time and a `HoneycombMeshGenerator` to set it up with a simple
honeycomb pattern. Here the first and second arguments define the size of the mesh -
we have chosen a mesh that is 5 nodes (i.e. cells) wide, and 5 nodes high. The extra '2' argument puts
two layers of non-cell elements around the mesh, which help to form a nicer voronoi tesselation
for area calculations.

```python
        chaste.core.OutputFileHandler("Python/TestSpheroidTutorial")
        generator = chaste.mesh.HoneycombMeshGenerator(5, 5)
        mesh = generator.GetMesh()

```
We create some cells next, with a stem-like proliferative type. This means they will continually
proliferate if there is enough oxygen, similar to how a tumour spheroid may behave.

```python
        stem_type = chaste.cell_based.StemCellProliferativeType()
        cell_generator = chaste.cell_based.CellsGeneratorSimpleOxygenBasedCellCycleModel_2()
        cells = cell_generator.GenerateBasicRandom(mesh.GetNumNodes(), stem_type)

```
Define when cells become apoptotic

```python
        for eachCell in cells:
            cell_cycle_model = eachCell.GetCellCycleModel()
            eachCell.GetCellData().SetItem("oxygen", 30.0)
            cell_cycle_model.SetDimension(2)
            cell_cycle_model.SetStemCellG1Duration(4.0)
            cell_cycle_model.SetHypoxicConcentration(0.1)
            cell_cycle_model.SetQuiescentConcentration(0.3)
            cell_cycle_model.SetCriticalHypoxicDuration(8)
            g1_duration = cell_cycle_model.GetStemCellG1Duration()
            sg2m_duration = cell_cycle_model.GetSG2MDuration()
            rnum = chaste.core.RandomNumberGenerator.Instance().ranf()
            birth_time = -rnum * (g1_duration + sg2m_duration)
            eachCell.SetBirthTime(birth_time)

```
Now we have a mesh and a set of cells to go with it, we can create a `CellPopulation` as before.

```python
        cell_population = chaste.cell_based.MeshBasedCellPopulation2_2(mesh, cells)

```
To view the results of this and the next test in Paraview it is necessary to explicitly generate the required .vtu files.

```python
        cell_population.AddPopulationWriterVoronoiDataWriter()

```
We then pass in the cell population into an `OffLatticeSimulation`, and set the output directory and end time.

```python
        simulator = chaste.cell_based.OffLatticeSimulation2_2(cell_population)
        simulator.SetOutputDirectory("Python/TestSpheroidTutorial")
        simulator.SetEndTime(5.0)

```
We ask for output every 12 increments

```python
        simulator.SetSamplingTimestepMultiple(100)

```
We define how the springs between cells behave using a force law.

```python
        force = chaste.cell_based.GeneralisedLinearSpringForce2_2()
        simulator.AddForce(force)

```
We set up a PDE for oxygen diffusion and consumption by cells, setting the rate of consumption to 0.1

```python
        pde = chaste.cell_based.CellwiseSourceEllipticPde2(cell_population, -0.5)

```
We set a constant amount of oxygen on the edge of the spheroid

```python
        bc = chaste.pde.ConstBoundaryCondition2(1.0)
        is_neumann_bc = False

```
Set up a pde modifier to solve the PDE at each simulation time step

```python
        #pde_modifier = chaste.cell_based.EllipticGrowingDomainPdeModifier2(pde, bc, is_neumann_bc)
        #pde_modifier.SetDependentVariableName("oxygen")
        #simulator.AddSimulationModifier(pde_modifier)

```
As before, we set up a scene modifier for real-time visualization

```python
        scene = chaste.visualization.VtkScene2()
        scene.SetCellPopulation(cell_population)
        scene.GetCellPopulationActorGenerator().SetColorByCellData(True)
        scene.GetCellPopulationActorGenerator().SetDataLabel("oxygen")
        scene.GetCellPopulationActorGenerator().SetShowCellCentres(True)
        scene.GetCellPopulationActorGenerator().SetShowVoronoiMeshEdges(False)
        # JUPYTER_SHOW_FIRST

        scene_modifier = chaste.cell_based.VtkSceneModifier2()
        scene_modifier.SetVtkScene(scene)
        scene_modifier.SetUpdateFrequency(100)
        simulator.AddSimulationModifier(scene_modifier)

```
Eventually remove apoptotic cells

```python
        killer = chaste.cell_based.ApoptoticCellKiller2(cell_population)
        simulator.AddCellKiller(killer)

```
To run the simulation, we call `Solve()`. We can again do a quick rendering of the population at the end of the simulation

```python
        scene.Start()
        simulator.Solve()

        # JUPYTER_TEARDOWN

```
Full results can be visualized in Paraview from the `file_handler.GetOutputDirectoryFullPath()` directory.

```python
if __name__ == '__main__':
    unittest.main(verbosity=2)

```


## Full code 


**File name:** `TestSpheroidTutorial.py` 

```python
import unittest # Python testing framework
import matplotlib.pyplot as plt # Plotting
import numpy as np # Matrix tools
import chaste # The PyChaste module
import chaste.mesh # Contains meshes
import chaste.pde # PDEs
import chaste.cell_based # Contains cell populations
import chaste.visualization # Visualization tools
chaste.init() # Set up MPI

class TestSpheroidTutorial(chaste.cell_based.AbstractCellBasedTestSuite):

    def test_spheroid(self):

        # JUPYTER_SETUP

        chaste.core.OutputFileHandler("Python/TestSpheroidTutorial")
        generator = chaste.mesh.HoneycombMeshGenerator(5, 5)
        mesh = generator.GetMesh()

        stem_type = chaste.cell_based.StemCellProliferativeType()
        cell_generator = chaste.cell_based.CellsGeneratorSimpleOxygenBasedCellCycleModel_2()
        cells = cell_generator.GenerateBasicRandom(mesh.GetNumNodes(), stem_type)

        for eachCell in cells:
            cell_cycle_model = eachCell.GetCellCycleModel()
            eachCell.GetCellData().SetItem("oxygen", 30.0)
            cell_cycle_model.SetDimension(2)
            cell_cycle_model.SetStemCellG1Duration(4.0)
            cell_cycle_model.SetHypoxicConcentration(0.1)
            cell_cycle_model.SetQuiescentConcentration(0.3)
            cell_cycle_model.SetCriticalHypoxicDuration(8)
            g1_duration = cell_cycle_model.GetStemCellG1Duration()
            sg2m_duration = cell_cycle_model.GetSG2MDuration()
            rnum = chaste.core.RandomNumberGenerator.Instance().ranf()
            birth_time = -rnum * (g1_duration + sg2m_duration)
            eachCell.SetBirthTime(birth_time)

        cell_population = chaste.cell_based.MeshBasedCellPopulation2_2(mesh, cells)

        cell_population.AddPopulationWriterVoronoiDataWriter()

        simulator = chaste.cell_based.OffLatticeSimulation2_2(cell_population)
        simulator.SetOutputDirectory("Python/TestSpheroidTutorial")
        simulator.SetEndTime(5.0)

        simulator.SetSamplingTimestepMultiple(100)

        force = chaste.cell_based.GeneralisedLinearSpringForce2_2()
        simulator.AddForce(force)

        pde = chaste.cell_based.CellwiseSourceEllipticPde2(cell_population, -0.5)

        bc = chaste.pde.ConstBoundaryCondition2(1.0)
        is_neumann_bc = False

        #pde_modifier = chaste.cell_based.EllipticGrowingDomainPdeModifier2(pde, bc, is_neumann_bc)
        #pde_modifier.SetDependentVariableName("oxygen")
        #simulator.AddSimulationModifier(pde_modifier)

        scene = chaste.visualization.VtkScene2()
        scene.SetCellPopulation(cell_population)
        scene.GetCellPopulationActorGenerator().SetColorByCellData(True)
        scene.GetCellPopulationActorGenerator().SetDataLabel("oxygen")
        scene.GetCellPopulationActorGenerator().SetShowCellCentres(True)
        scene.GetCellPopulationActorGenerator().SetShowVoronoiMeshEdges(False)
        # JUPYTER_SHOW_FIRST

        scene_modifier = chaste.cell_based.VtkSceneModifier2()
        scene_modifier.SetVtkScene(scene)
        scene_modifier.SetUpdateFrequency(100)
        simulator.AddSimulationModifier(scene_modifier)

        killer = chaste.cell_based.ApoptoticCellKiller2(cell_population)
        simulator.AddCellKiller(killer)

        scene.Start()
        simulator.Solve()

        # JUPYTER_TEARDOWN

if __name__ == '__main__':
    unittest.main(verbosity=2)

```

