
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
## In this tutorial we show how Chaste can be used to create, run and visualize vertex-based simulations. 
## Full details of the mechanical model proposed by T. Nagai and H. Honda ("A dynamic cell model for the formation of epithelial tissues", 
## Philosophical Magazine Part B 81:699-719).
##
## ## The Test

import unittest  # Python testing framework
import matplotlib.pyplot as plt  # Plotting
import numpy as np  # Matrix tools
import chaste  # The PyChaste module
import chaste.cell_based  # Contains cell populations
import chaste.mesh  # Contains meshes
import chaste.visualization  # Visualization tools
chaste.init()  # Set up MPI

class TestRunningVertexBasedSimulationsTutorial(chaste.cell_based.AbstractCellBasedTestSuite):
    ## ### Test 1 - A basic vertex-based simulation
    ## In the first test, we run a simple vertex-based simulation, in which we create a monolayer of cells, 
    ## using a mutable vertex mesh. Each cell is assigned a stochastic cell-cycle model.

    def test_monolayer(self):

        # JUPYTER_SETUP

        ## First, we generate a vertex mesh. To create a MutableVertexMesh, we can use the HoneycombVertexMeshGenerator. 
        ## This generates a honeycomb-shaped mesh, in which all nodes are equidistant. Here the first and second arguments 
        ## define the size of the mesh - we have chosen a mesh that is 2 elements (i.e. cells) wide, and 2 elements high.

        chaste.core.OutputFileHandler("Python/TestVertexBasedCellSimulationsTutorial")
        generator = chaste.mesh.HoneycombVertexMeshGenerator(2, 2)
        mesh = generator.GetMesh()

        ## Having created a mesh, we now create a std::vector of CellPtrs. To do this, we use the CellsGenerator helper class, 
        ## which is templated over the type of cell model required 
        ## and the dimension. We create an empty vector of cells and pass this into the method along with the mesh. 
        ## The second argument represents the size of that the vector cells should become - one cell for each element, 
        ## the third argument specifies the proliferative type of the cell.

        transit_type = chaste.cell_based.TransitCellProliferativeType()
        cell_generator = chaste.cell_based.CellsGeneratorUniformG1GenerationalCellCycleModel_2()
        cells = cell_generator.GenerateBasicRandom(mesh.GetNumElements(), transit_type)

        ## Now we have a mesh and a set of cells to go with it, we can create a CellPopulation. 
        ## In general, this class associates a collection of cells with a mesh. For this test, because we have a MutableVertexMesh, 
        ## we use a particular type of cell population called a VertexBasedCellPopulation.

        cell_population = chaste.cell_based.VertexBasedCellPopulation2(mesh, cells)

        ## We can set up a `VtkScene` to do a quick visualization of the population before running the analysis.

        scene = chaste.visualization.VtkScene2()
        scene.SetCellPopulation(cell_population)
        # JUPYTER_SHOW_FIRST
        scene.Start()  # JUPYTER_SHOW

        ## We then pass in the cell population into an `OffLatticeSimulation`, and set the output directory, output multiple and end time

        simulator = chaste.cell_based.OffLatticeSimulation2_2(cell_population)
        simulator.SetOutputDirectory("Python/TestVertexBasedCellSimulationsTutorial")
        simulator.SetEndTime(5.0)

        ## For longer simulations, we may not want to output the results every time step. 
        ## In this case we can use the following method, to print results every 50 time steps instead. 
        ## As the default time step used by the simulator (for vertex based simulations), is 0.02 hours, this method 
        ## will cause the simulator to print results every 6 minutes (i.e. 0.1 hours).

        simulator.SetSamplingTimestepMultiple(50)

        ## We must now create one or more force laws, which determine the mechanics of the vertices of each cell in a cell population. 
        ## For this test, we use one force law, based on the Nagai-Honda mechanics, and pass it to the OffLatticeSimulation. 
        ## For a list of possible forces see subclasses of AbstractForce. 
        ## Note that some of these forces are not compatible with vertex-based simulations see the specific class documentation for details, 
        ## if you try to use an incompatible class then you will receive a warning.

        force = chaste.cell_based.NagaiHondaForce2()
        simulator.AddForce(force)

        ## A NagaiHondaForce assumes that each cell has a target area. The target areas of cells are used to determine 
        ## pressure forces on each vertex and eventually determine the size of each cell in the simulation. 
        ## In order to assign target areas to cells and update them in each time step we add a SimpleTargetAreaModifier 
        ## to the simulation, which inherits from AbstractTargetAreaModifier.

        growth_modifier = chaste.cell_based.SimpleTargetAreaModifier2()
        simulator.AddSimulationModifier(growth_modifier)

        ## Save snapshot images of the population during the simulation

        scene_modifier = chaste.cell_based.VtkSceneModifier2()
        scene_modifier.SetVtkScene(scene)
        scene_modifier.SetUpdateFrequency(100)
        simulator.AddSimulationModifier(scene_modifier)

        ## To run the simulation, we call `Solve()`. We can again do a quick rendering of the population at the end of the simulation

        scene.Start()
        simulator.Solve()
        scene.End()

        ## The next two lines are for test purposes only and are not part of this tutorial. 
        ## If different simulation input parameters are being explored the lines should be removed.

        self.assertEqual(cell_population.GetNumRealCells(), 7)
        self.assertAlmostEqual(chaste.cell_based.SimulationTime.Instance().GetTime(), 5.0, 6)

        # JUPYTER_TEARDOWN

    ## ### Test 2 - introducing periodicity, boundaries and cell killers
    ## In the second test, we run a simple vertex-based simulation, in which we create a monolayer of cells in a periodic geometry, 
    ## using a cylindrical vertex mesh. We also include a fixed boundary which cells can't pass through and a cell killer which removes
    ## cells once they leave a region. As before each cell is assigned a stochastic cell-cycle model.

    def test_periodic_monolayer(self):

        # JUPYTER_SETUP

        ## First, we generate a periodic vertex mesh. To create a Cylindrical2dVertexMesh, we can use the CylindricalHoneycombVertexMeshGenerator. 
        ## This generates a honeycomb-shaped mesh, in which all nodes are equidistant and the right hand side is associated with the left hand side. 
        ## Here the first and second arguments define the size of the mesh - we have chosen a mesh that is
        ## 4 elements (i.e. cells) wide, and 4 elements high.

        generator = chaste.mesh.CylindricalHoneycombVertexMeshGenerator(4, 4)
        mesh = generator.GetCylindricalMesh()

        ## Having created a mesh, we now create a VectorSharedPtrCells. This is exactly the same as the above test.

        transit_type = chaste.cell_based.TransitCellProliferativeType()
        cell_generator = chaste.cell_based.CellsGeneratorUniformG1GenerationalCellCycleModel_2()
        cells = cell_generator.GenerateBasicRandom(mesh.GetNumElements(), transit_type)

        ## Now we have a mesh and a set of cells to go with it, we can create a CellPopulation. This is also the same as in the above test.

        cell_population = chaste.cell_based.VertexBasedCellPopulation2(mesh, cells)

        ## We then pass in the cell population into an `OffLatticeSimulation`, and set the output directory, output multiple and end time

        simulator = chaste.cell_based.OffLatticeSimulation2_2(cell_population)
        simulator.SetOutputDirectory("Python/TestPeriodicVertexBasedCellPopulation")
        simulator.SetEndTime(1.0)
        simulator.SetSamplingTimestepMultiple(50)

        ## We now make a pointer to an appropriate force and pass it to the OffLatticeSimulation.

        force = chaste.cell_based.NagaiHondaForce2()
        simulator.AddForce(force)

        ## We also make a pointer to the target area modifier and add it to the simulator.

        growth_modifier = chaste.cell_based.SimpleTargetAreaModifier2()
        simulator.AddSimulationModifier(growth_modifier)

        ## We now create one or more CellPopulationBoundaryConditions, which determine any conditions which each cell in a cell population must satisfy. 
        ## For this test, we use a PlaneBoundaryCondition, and pass it to the OffLatticeSimulation. For a list of possible boundary condition
        ## see subclasses of AbstractCellPopulationBoundaryCondition. 
        ## Note that some of these boundary conditions are not compatible with vertex-based simulations see the specific class documentation 
        ## for details, if you try to use an incompatible class then you will receive a warning.
        ## The first step is to define a point on the plane boundary and a normal to the plane.

        point = np.array([0.0, 0.0])
        normal = np.array([0.0, -1.0])

        ## We can now make a PlaneBoundaryCondition (passing the point and normal to the plane) and pass it to the OffLatticeSimulation.
        bc = chaste.cell_based.PlaneBoundaryCondition2_2(cell_population, point, normal)
        simulator.AddCellPopulationBoundaryCondition(bc)

        ## We now create one or more CellKillers, which determine how cells are removed from the simulation. 
        ## For this test, we use a PlaneBasedCellKiller, and pass it to the OffLatticeSimulation. 
        ## For a list of possible cell killers see subclasses of AbstractCellKiller.
        ## The first step is to define a point on the plane boundary and a normal to the plane.

        point = np.array([0.0, 3.0])
        normal = np.array([0.0, 1.0])

        ## Finally we now make a PlaneBasedCellKiller (passing the point and normal to the plane) and pass it to the OffLatticeSimulation.
        killer = chaste.cell_based.PlaneBasedCellKiller2(cell_population, point, normal)
        simulator.AddCellKiller(killer)

        ## To run the simulation, we call `Solve()`.
        simulator.Solve()

        ## The next two lines are for test purposes only and are not part of this tutorial. 
        ## If different simulation input parameters are being explored the lines should be removed.

        self.assertEqual(cell_population.GetNumRealCells(), 12)
        self.assertAlmostEqual(chaste.cell_based.SimulationTime.Instance().GetTime(), 1.0, 6)

        # JUPYTER_TEARDOWN 


if __name__ == '__main__':
    unittest.main(verbosity=2)

#endif END_WIKI