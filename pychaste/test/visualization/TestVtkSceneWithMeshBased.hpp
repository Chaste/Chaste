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

#ifndef TESTVTKSCENEWITHMESHBASEDPOPULATION_HPP_
#define TESTVTKSCENEWITHMESHBASEDPOPULATION_HPP_

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

class TestVtkSceneWithMeshBasedPopulation : public AbstractCellBasedTestSuite
{
public:

    void TestRenderingMeshBasedPopulation()
    {
        OutputFileHandler file_handler1 = OutputFileHandler("TestVtkSceneWithMeshBasedPopulation/2d/");

        HoneycombMeshGenerator generator(2, 2, 2);    // Parameters are: cells across, cells up
        boost::shared_ptr<MutableMesh<2,2> > p_mesh = generator.GetMesh();
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, location_indices.size(), p_transit_type);
        auto p_cell_population = boost::make_shared<MeshBasedCellPopulationWithGhostNodes<2> >(*p_mesh, cells, location_indices);
        p_cell_population->AddPopulationWriter<VoronoiDataWriter>();

        auto p_scene = boost::make_shared<VtkScene<2> >();
        p_scene->SetCellPopulation(p_cell_population);
        p_scene->SetSaveAsImages(true);
        p_scene->SetOutputFilePath(file_handler1.GetOutputDirectoryFullPath()+"/cell_population");

        auto p_scene_modifier = boost::make_shared<VtkSceneModifier<2> >();
        p_scene_modifier->SetVtkScene(p_scene);
        p_scene->Start();

        OffLatticeSimulation<2> simulator(*p_cell_population);
        simulator.SetOutputDirectory("TestVtkSceneWithMeshBasedPopulation/2d/");
        simulator.SetEndTime(4.0);
        simulator.SetSamplingTimestepMultiple(12);

        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        simulator.AddForce(p_force);
        simulator.AddSimulationModifier(p_scene_modifier);
        simulator.Solve();
    }

    void TestRendering3dMeshBasedPopulation()
    {
        OutputFileHandler file_handler1 = OutputFileHandler("TestVtkSceneWithMeshBasedPopulation/3d/");

        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, true,  0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, true,  1.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(2, true,  1.0, 0.0, 1.0));
        nodes.push_back(new Node<3>(3, true,  0.0, 1.0, 1.0));
        nodes.push_back(new Node<3>(4, false, 0.5, 0.5, 0.5));
        MutableMesh<3,3> mesh(nodes);

        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<UniformCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_transit_type);

        auto p_cell_population = boost::make_shared<MeshBasedCellPopulation<3> >(mesh, cells);
        p_cell_population->SetAbsoluteMovementThreshold(DBL_MAX);
        p_cell_population->AddPopulationWriter<VoronoiDataWriter>();

        auto p_scene = boost::make_shared<VtkScene<3> >();
        p_scene->SetCellPopulation(p_cell_population);
        p_scene->SetSaveAsImages(true);
        p_scene->SetOutputFilePath(file_handler1.GetOutputDirectoryFullPath()+"/cell_population");

        auto p_scene_modifier = boost::make_shared<VtkSceneModifier<3> >();
        p_scene_modifier->SetVtkScene(p_scene);
        p_scene->Start();

        OffLatticeSimulation<3> simulator(*p_cell_population);
        simulator.SetOutputDirectory("TestVtkSceneWithMeshBasedPopulation/3d/");
        simulator.SetEndTime(2.0);
        simulator.SetSamplingTimestepMultiple(5);

        MAKE_PTR(GeneralisedLinearSpringForce<3>, p_force);
        p_force->SetMeinekeSpringStiffness(30.0); // default is 15.0;
        p_force->SetCutOffLength(1.5);
        simulator.AddForce(p_force);
        simulator.AddSimulationModifier(p_scene_modifier);
        simulator.Solve();
    }
};
#endif
