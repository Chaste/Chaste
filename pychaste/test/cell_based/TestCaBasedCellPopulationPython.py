
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

import unittest
import chaste.mesh
import chaste.cell_based
import chaste.visualization
chaste.init()


class TestPottsMesh(chaste.cell_based.AbstractCellBasedTestSuite):

    def test_construct(self):
        work_dir = "Python/TestCaBasedPopulationPython"
        file_handler = chaste.core.OutputFileHandler(work_dir)
        generator = chaste.mesh.PottsMeshGenerator3(10, 0, 0,
                                                    10, 0, 0,
                                                    10, 0, 0)
        mesh = generator.GetMesh()

        locs = range(5)
        cell_generator = chaste.cell_based.CellsGeneratorUniformCellCycleModel_3()
        cells = cell_generator.GenerateBasic(len(locs))
        cell_population = chaste.cell_based.CaBasedCellPopulation3(mesh,
                                                                   cells,
                                                                   locs)

        scene = chaste.visualization.VtkScene3()
        scene.SetCellPopulation(cell_population)
        scene.SetSaveAsImages(True)
        scene.SetOutputFilePath(file_handler.GetOutputDirectoryFullPath() +
                                "/cell_population")
        scene.GetCellPopulationActorGenerator().SetShowPottsMeshEdges(True)
        scene.GetCellPopulationActorGenerator().SetVolumeOpacity(1.0)
        scene_modifier = chaste.cell_based.VtkSceneModifier3()
        scene_modifier.SetVtkScene(scene)
        scene.Start()

        simulator = chaste.cell_based.OnLatticeSimulation3(cell_population)
        simulator.SetOutputDirectory(work_dir)
        simulator.SetDt(10.0)
        simulator.SetEndTime(300.0)
        simulator.AddSimulationModifier(scene_modifier)
        simulator.Solve()


if __name__ == '__main__':
    unittest.main()
