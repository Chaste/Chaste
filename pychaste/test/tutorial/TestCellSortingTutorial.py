
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
## This test is a demonstration of cell sorting using a Cellular Potts based framework.
## It shows:
## * How to set up a Potts simulation
## * Working with labels
##
## ## The Test

import unittest # Python testing framework
import matplotlib.pyplot as plt # Plotting
import numpy as np # Matrix tools
import chaste # The PyChaste module
chaste.init() # Set up MPI
import chaste.cell_based # Contains cell populations
import chaste.mesh # Contains meshes
import chaste.visualization # Visualization tools

class TestCellSortingTutorial(chaste.cell_based.AbstractCellBasedTestSuite):
    
    ## ### Test 1 - Cell sorting
    ## The next test generates a collection of cells, there are two types of cells, labelled ones and non labelled ones, 
    ## there is differential adhesion between the cell types. For the parameters specified, the cells sort into separate types.
    
    def test_potts_monolayer_cell_sorting(self):
        
        # JUPYTER_SETUP 
        
        ## First, we generate a `Potts` mesh. To create a `PottsMesh`, we can use the `PottsMeshGenerator`. 
        ## This generates a regular square-shaped mesh, in which all elements are the same size. 
        ## We have chosen an 8 by 8 block of elements each consisting of 4 by 4 ( = 16) lattice sites.
        
        generator = chaste.mesh.PottsMeshGenerator2(50, 8, 4, 50, 8, 4)
        mesh = generator.GetMesh()
          
        ## Having created a mesh, we now create some cells. To do this, we the `CellsGenerator` helper class, 
        ## as before but this time the third argument is set to make all cells non-proliferative.
        
        differentiated_type = chaste.cell_based.DifferentiatedCellProliferativeType()
        cell_generator = chaste.cell_based.CellsGeneratorUniformCellCycleModel_2()
        cells = cell_generator.GenerateBasicRandom(mesh.GetNumElements(), differentiated_type)
        
        ## Before we make a CellPopulation we make a cell label and then assign this label to some randomly chosen cells.
        
        label = chaste.cell_based.CellLabel()
        for eachCell in cells:
            if(chaste.core.RandomNumberGenerator.Instance().ranf()<0.5):
                eachCell.AddCellProperty(label)
          
        ## Now we have a mesh and a set of cells to go with it, we can create a `CellPopulation`. 
        
        cell_population = chaste.cell_based.PottsBasedCellPopulation2(mesh, cells)
        
        ## In order to visualize labelled cells we need to use the following command.
        
        cell_population.AddCellWriterCellLabelWriter()
        
        ## PyChaste can do simple 3D rendering with VTK. We set up a VtkScene so that we can 
        ## see the population evovle in real time.
        
        scene= chaste.visualization.VtkScene2()
        scene.SetCellPopulation(cell_population)
        scene.GetCellPopulationActorGenerator().SetShowPottsMeshEdges(True)
        # JUPYTER_SHOW_FIRST
        scene.Start()  # JUPYTER_SHOW
        
        ## We then pass in the cell population into an `OffLatticeSimulation`, and set the output directory and end time

        simulator = chaste.cell_based.OnLatticeSimulation2(cell_population)
        simulator.SetOutputDirectory("Python/TestCellSorting")
        simulator.SetEndTime(20.0)
        simulator.SetSamplingTimestepMultiple(10)

        ## We must now create one or more update rules, which determine the Hamiltonian in the Potts simulation. 
        ## For this test, we use two update rules based upon a volume constraint (`VolumeConstraintPottsUpdateRule`) and 
        ## differential adhesion between cells (`DifferentialAdhesionPottsUpdateRule`), set appropriate parameters, and 
        ## pass them to the `OnLatticeSimulation`.
        
        volume_constraint_update_rule = chaste.cell_based.VolumeConstraintPottsUpdateRule2()
        volume_constraint_update_rule.SetMatureCellTargetVolume(16)
        volume_constraint_update_rule.SetDeformationEnergyParameter(0.2)
        simulator.AddUpdateRule(volume_constraint_update_rule)
        
        ## We repeat the process for any other update rules.

        differential_adhesion_update_rule = chaste.cell_based.DifferentialAdhesionPottsUpdateRule2() 
        differential_adhesion_update_rule.SetLabelledCellLabelledCellAdhesionEnergyParameter(0.16)
        differential_adhesion_update_rule.SetLabelledCellCellAdhesionEnergyParameter(0.11)
        differential_adhesion_update_rule.SetCellCellAdhesionEnergyParameter(0.02)
        differential_adhesion_update_rule.SetLabelledCellBoundaryAdhesionEnergyParameter(0.16)
        differential_adhesion_update_rule.SetCellBoundaryAdhesionEnergyParameter(0.16)
        simulator.AddUpdateRule(differential_adhesion_update_rule)
        
        ## Set up plotting
        
        scene_modifier = chaste.cell_based.VtkSceneModifier2()
        scene_modifier.SetVtkScene(scene)
        scene_modifier.SetUpdateFrequency(1000)
        simulator.AddSimulationModifier(scene_modifier)
        
        ## To run the simulation, we call `Solve()`.
    
        scene.Start()
        simulator.Solve()
        
        # JUPYTER_TEARDOWN 
        
if __name__ == '__main__':
    unittest.main(verbosity=2)
    
#endif END_WIKI