/*

Copyright (c) 2005-2024, University of Oxford.
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

*/

#ifndef TESTVTKSCENEWITHPOTTSBASEDPOPAULTION_HPP_
#define TESTVTKSCENEWITHPOTTSBASEDPOPAULTION_HPP_

#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "PottsMeshGenerator.hpp"
#include "Cell.hpp"
#include "CellsGenerator.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "UniformCellCycleModel.hpp"
#include "CaBasedCellPopulation.hpp"
#include "CellMutationStatesCountWriter.hpp"
#include "CellProliferativeTypesCountWriter.hpp"
#include "CellProliferativePhasesCountWriter.hpp"
#include "CellProliferativePhasesWriter.hpp"
#include "CellAncestorWriter.hpp"
#include "CellAgesWriter.hpp"
#include "CellVolumesWriter.hpp"
#include "CellMutationStatesWriter.hpp"
#include "OnLatticeSimulation.hpp"
#include "SmartPointers.hpp"
#include "VtkScene.hpp"
#include "VtkSceneModifier.hpp"
#include "FileFinder.hpp"
#include "OutputFileHandler.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "OffLatticeSimulation.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "VoronoiDataWriter.hpp"
#include "AdhesionPottsUpdateRule.hpp"
#include "VolumeConstraintPottsUpdateRule.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestVtkScene : public AbstractCellBasedTestSuite
{
public:

    void TestRenderingPottsBasedPopulation()
    {
        OutputFileHandler file_handler1 = OutputFileHandler("TestVtkSceneWithMeshBasedPopulation/2d/");

        PottsMeshGenerator<2> generator(50, 2, 4, 50, 2, 4);
        boost::shared_ptr<PottsMesh<2> > p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_transit_type);

        auto  p_cell_population = boost::make_shared<PottsBasedCellPopulation<2> >(*p_mesh, cells);
        p_cell_population->SetTemperature(0.1);
        p_cell_population->SetNumSweepsPerTimestep(1);

        OnLatticeSimulation<2> simulator(*p_cell_population);
        simulator.SetOutputDirectory("TestVtkSceneWithMeshBasedPopulation/2d/");
        simulator.SetEndTime(10.0);
        simulator.SetDt(0.1);
        simulator.SetSamplingTimestepMultiple(10);

        auto p_scene = boost::make_shared<VtkScene<2> >();
        p_scene->SetCellPopulation(p_cell_population);
        p_scene->SetIsInteractive(true);
        p_scene->SetSaveAsImages(false);
        p_scene->GetCellPopulationActorGenerator()->SetShowPottsMeshEdges(true);
        p_scene->SetOutputFilePath(file_handler1.GetOutputDirectoryFullPath()+"/cell_population");

        auto p_scene_modifier = boost::make_shared<VtkSceneModifier<2> >();
        p_scene_modifier->SetVtkScene(p_scene);
        p_scene_modifier->SetUpdateFrequency(10);
        p_scene->Start();

        MAKE_PTR(VolumeConstraintPottsUpdateRule<2>, p_volume_constraint_update_rule);
        p_volume_constraint_update_rule->SetMatureCellTargetVolume(16);
        p_volume_constraint_update_rule->SetDeformationEnergyParameter(0.2);
        simulator.AddUpdateRule(p_volume_constraint_update_rule);
        MAKE_PTR(AdhesionPottsUpdateRule<2>, p_adhesion_update_rule);
        simulator.AddUpdateRule(p_adhesion_update_rule);
        simulator.AddSimulationModifier(p_scene_modifier);
        simulator.Solve();
    }

    void TestRenderingPottsBasedPopulation3d()
    {
        OutputFileHandler file_handler1 = OutputFileHandler("TestVtkSceneWithMeshBasedPopulation/3d/");

        PottsMeshGenerator<3> generator(50, 2, 4, 50, 2, 4, 50, 2, 4);
        boost::shared_ptr<PottsMesh<3> > p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<UniformCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_transit_type);

        auto p_cell_population = boost::make_shared<PottsBasedCellPopulation<3> >(*p_mesh, cells);
        p_cell_population->SetTemperature(0.1);
        p_cell_population->SetNumSweepsPerTimestep(1);

        OnLatticeSimulation<3> simulator(*p_cell_population);
        simulator.SetOutputDirectory("TestVtkSceneWithMeshBasedPopulation/3d/");
        simulator.SetEndTime(1.0);
        simulator.SetDt(0.1);
        simulator.SetSamplingTimestepMultiple(100);

        auto p_scene = boost::make_shared<VtkScene<3> >();
        p_scene->SetCellPopulation(p_cell_population);
        p_scene->SetIsInteractive(true);
        p_scene->SetSaveAsImages(false);
        p_scene->GetCellPopulationActorGenerator()->SetShowCellCentres(true);
        p_scene->GetCellPopulationActorGenerator()->SetShowPottsMeshOutlines(true);
        p_scene->SetOutputFilePath(file_handler1.GetOutputDirectoryFullPath()+"/cell_population");

        auto p_scene_modifier = boost::make_shared<VtkSceneModifier<3> >();
        p_scene_modifier->SetVtkScene(p_scene);
        p_scene->Start();

        MAKE_PTR(VolumeConstraintPottsUpdateRule<3>, p_volume_constraint_update_rule);
        p_volume_constraint_update_rule->SetMatureCellTargetVolume(16);
        p_volume_constraint_update_rule->SetDeformationEnergyParameter(0.2);
        simulator.AddUpdateRule(p_volume_constraint_update_rule);
        MAKE_PTR(AdhesionPottsUpdateRule<3>, p_adhesion_update_rule);
        simulator.AddUpdateRule(p_adhesion_update_rule);
        simulator.AddSimulationModifier(p_scene_modifier);
        simulator.Solve();
    }
};
#endif
