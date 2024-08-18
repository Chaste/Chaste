
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
## This tutorial is an example of modelling spheroid growth with a nutrient.
## It covers:
## * Setting up an off-lattice cell population
## * Setting up a cell cycle model with oxygen dependence
## * Setting up and solving an oxygen transport PDE
## * Setting up a cell killer
##
## ## The Test

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

    ## ### Test 1 - a 2D mesh-based spheroid
    ## In this test we set up a spheroid with a plentiful supply of oxygen on the boundary and watch it grow
    ## over time. Cells can gradually become apoptotic if the oxygen tension is too low.

    def test_spheroid(self):

        # JUPYTER_SETUP

        ## This time we will use on off-lattice `MeshBased` cell population. Cell centres are joined with
        ## springs with a Delauney Triangulation used to identify neighbours. Cell area is given by the dual
        ## (Voronoi Tesselation). We start off with a small number of cells. We use a `MutableMesh` which
        ## can change connectivity over time and a `HoneycombMeshGenerator` to set it up with a simple
        ## honeycomb pattern. Here the first and second arguments define the size of the mesh - 
        ## we have chosen a mesh that is 5 nodes (i.e. cells) wide, and 5 nodes high. The extra '2' argument puts
        ## two layers of non-cell elements around the mesh, which help to form a nicer voronoi tesselation
        ## for area calculations.

        chaste.core.OutputFileHandler("Python/TestSpheroidTutorial")
        generator = chaste.mesh.HoneycombMeshGenerator(5, 5)
        mesh = generator.GetMesh()

        ## We create some cells next, with a stem-like proliferative type. This means they will continually
        ## proliferate if there is enough oxygen, similar to how a tumour spheroid may behave.

        stem_type = chaste.cell_based.StemCellProliferativeType()
        cell_generator = chaste.cell_based.CellsGeneratorSimpleOxygenBasedCellCycleModel_2()
        cells = cell_generator.GenerateBasicRandom(mesh.GetNumNodes(), stem_type)

        ## Define when cells become apoptotic

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

        ## Now we have a mesh and a set of cells to go with it, we can create a `CellPopulation` as before. 

        cell_population = chaste.cell_based.MeshBasedCellPopulation2_2(mesh, cells)

        ## To view the results of this and the next test in Paraview it is necessary to explicitly generate the required .vtu files. 

        cell_population.AddPopulationWriterVoronoiDataWriter()

        ## We then pass in the cell population into an `OffLatticeSimulation`, and set the output directory and end time.

        simulator = chaste.cell_based.OffLatticeSimulation2_2(cell_population)
        simulator.SetOutputDirectory("Python/TestSpheroidTutorial")
        simulator.SetEndTime(5.0)

        ## We ask for output every 12 increments

        simulator.SetSamplingTimestepMultiple(100)

        ## We define how the springs between cells behave using a force law.

        force = chaste.cell_based.GeneralisedLinearSpringForce2_2()
        simulator.AddForce(force)

        ## We set up a PDE for oxygen diffusion and consumption by cells, setting the rate of consumption to 0.1

        pde = chaste.cell_based.CellwiseSourceEllipticPde2(cell_population, -0.5)

        ## We set a constant amount of oxygen on the edge of the spheroid

        bc = chaste.pde.ConstBoundaryCondition2(1.0)
        is_neumann_bc = False

        ## Set up a pde modifier to solve the PDE at each simulation time step

        #pde_modifier = chaste.cell_based.EllipticGrowingDomainPdeModifier2(pde, bc, is_neumann_bc)
        #pde_modifier.SetDependentVariableName("oxygen")
        #simulator.AddSimulationModifier(pde_modifier)

        ## As before, we set up a scene modifier for real-time visualization

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

        ## Eventually remove apoptotic cells

        killer = chaste.cell_based.ApoptoticCellKiller2(cell_population)
        simulator.AddCellKiller(killer)

        ## To run the simulation, we call `Solve()`. We can again do a quick rendering of the population at the end of the simulation

        scene.Start()
        simulator.Solve()

        # JUPYTER_TEARDOWN 

        ## Full results can be visualized in Paraview from the `file_handler.GetOutputDirectoryFullPath()` directory.


if __name__ == '__main__':
    unittest.main(verbosity=2)

#endif END_WIKI