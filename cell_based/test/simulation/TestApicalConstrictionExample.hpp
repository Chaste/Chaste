/*

Copyright (c) 2005-2016, University of Oxford.
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

#ifndef TESTAPICALCONSTRICTIONEXAMPLE_HPP_
#define TESTAPICALCONSTRICTIONEXAMPLE_HPP_

#include <cxxtest/TestSuite.h>
#include "AbstractCellBasedTestSuite.hpp"

#include "CheckpointArchiveTypes.hpp"

#include "PatternedApicalConstrictionForce.hpp"
#include "MonolayerVertexMeshGenerator.hpp"
#include "VoronoiPrism3dVertexMeshGenerator.hpp"
#include "VoronoiVertexMeshGenerator.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "CellsGenerator.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "NoCellCycleModel.hpp"
#include "SmartPointers.hpp"
#include "OffLatticeSimulation.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "CellLabel.hpp"
#include "CellLabelWriter.hpp"
#include "FakePetscSetup.hpp"

#define OUTPUT_NAME "TestApicalConstrictionExample/InitialMesh"


class TestApicalConstrictionExample : public AbstractCellBasedTestSuite
{
public:

    void TestApicalConstriction() throw (Exception)
    {
        // Make a mesh of 10x10
//        const double z_height = 1;
        const double target_area = 1;
        const unsigned num_cells_x = 10;
        const unsigned num_cells_y = 10;
        // There seems to be a bug somewhere in voronoiprism3dVertexMeshGenerator....
//        VoronoiPrism3dVertexMeshGenerator generator(num_cells_x, num_cells_y, z_height, 5, target_area);
//        MutableVertexMesh<3,3>* p_mesh = generator.GetMeshAfterReMesh();
        HoneycombVertexMeshGenerator generator(num_cells_x, num_cells_y, false, 0.1, 0.01, target_area);
//        VoronoiVertexMeshGenerator generator(num_cells_x, num_cells_y, 5, target_area);
        MutableVertexMesh<2, 2>& vertex_2mesh = *(generator.GetMesh());
        MonolayerVertexMeshGenerator builder("ApicalConstriction");
        MutableVertexMesh<3, 3>* p_mesh = builder.MakeMeshUsing2dMesh(vertex_2mesh);
        builder.WriteVtk(OUTPUT_NAME,"Before");

        char tmp_name[50];
        sprintf(tmp_name, "TestApicalConstrictionExample/HoneyTest%dx%d", num_cells_x, num_cells_y);
        VertexMeshWriter<3, 3> vertex_mesh_writer(tmp_name, "InitialMesh", false);
        vertex_mesh_writer.WriteVtkUsingMeshWithCellId(*p_mesh);

        std::vector<CellPtr> cells;
        CellsGenerator<NoCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
        VertexBasedCellPopulation<3> cell_population(*p_mesh, cells);
        cell_population.AddCellWriter<CellLabelWriter>();


        // Label a few cells near the centre of the population
        c_vector<double, 3> population_centre;
        population_centre(0) = 0.5*p_mesh->GetWidth(0);
        population_centre(1) = 0.5*p_mesh->GetWidth(1);
        population_centre(2) = 0.5*p_mesh->GetWidth(2);

        for (AbstractCellPopulation<3>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            c_vector<double, 3> cell_location = cell_population.GetLocationOfCellCentre(*cell_iter);
            c_vector<double, 3> dist = cell_location - population_centre;
            if (dist(0)*dist(0) + dist(1)*dist(1) < 2.0*2.0)
            {
                boost::shared_ptr<AbstractCellProperty> p_label = cell_population.GetCellPropertyRegistry()->Get<CellLabel>();
                cell_iter->AddCellProperty(p_label);
            }
        }

        OffLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory(tmp_name);
        simulator.SetSamplingTimestepMultiple(10);
        const double end_time = 4;
        simulator.SetEndTime(end_time);

        MAKE_PTR(PatternedApicalConstrictionForce, p_force3);
        p_force3->SetApicalParameters(10, 10, 1);
        p_force3->SetBasalParameters(10, 10, 1);
        p_force3->SetLateralParameter(4);
        p_force3->SetVolumeParameters(200, 1);
        p_force3->SetPatternedApicalParameter(20, 20, 0.5);
        simulator.AddForce(p_force3);

        simulator.Solve();

        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 100u);
        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), end_time, 1e-10);
    }

    void TestApicalConstrictionVoronoi() throw (Exception)
    {
        // Make a mesh of 10x10
        const double z_height = 1;
        const double target_area = 1;
        const unsigned num_cells_x = 10;
        const unsigned num_cells_y = 10;
        // There seems to be a bug somewhere in voronoiprism3dVertexMeshGenerator....
        VoronoiPrism3dVertexMeshGenerator generator(num_cells_x, num_cells_y, z_height, 5, target_area);
        MutableVertexMesh<3,3>* p_mesh = generator.GetMeshAfterReMesh();
//        HoneycombVertexMeshGenerator generator(num_cells_x, num_cells_y, false, 0.1, 0.01, target_area);
//        VoronoiVertexMeshGenerator generator(num_cells_x, num_cells_y, 5, target_area);
//        MutableVertexMesh<2, 2>& vertex_2mesh = *(generator.GetMesh());
//        Helper3dVertexMeshBuilder builder("ApicalConstriction");
//        MutableVertexMesh<3, 3>* p_mesh = builder.MakeMeshUsing2dMesh(vertex_2mesh);
//        builder.WriteVtk(OUTPUT_NAME,"Before");

        char tmp_name[50];
        sprintf(tmp_name, "TestApicalConstrictionExample/VoronoiTest%dx%d", num_cells_x, num_cells_y);
        VertexMeshWriter<3, 3> vertex_mesh_writer(tmp_name, "InitialMesh", false);
        vertex_mesh_writer.WriteVtkUsingMeshWithCellId(*p_mesh);

        std::vector<CellPtr> cells;
        CellsGenerator<NoCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
        VertexBasedCellPopulation<3> cell_population(*p_mesh, cells);
        cell_population.AddCellWriter<CellLabelWriter>();


        // Label a few cells near the centre of the population
        c_vector<double, 3> population_centre;
        population_centre(0) = 0.5*p_mesh->GetWidth(0);
        population_centre(1) = 0.5*p_mesh->GetWidth(1);
        population_centre(2) = 0.5*p_mesh->GetWidth(2);

        for (AbstractCellPopulation<3>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            c_vector<double, 3> cell_location = cell_population.GetLocationOfCellCentre(*cell_iter);
            c_vector<double, 3> dist = cell_location - population_centre;
            if (dist(0)*dist(0) + dist(1)*dist(1) < 2.0*2.0)
            {
                boost::shared_ptr<AbstractCellProperty> p_label = cell_population.GetCellPropertyRegistry()->Get<CellLabel>();
                cell_iter->AddCellProperty(p_label);
            }
        }

        OffLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory(tmp_name);
        simulator.SetSamplingTimestepMultiple(10);
        const double end_time = 4;
        simulator.SetEndTime(end_time);

        MAKE_PTR(PatternedApicalConstrictionForce, p_force3);
        p_force3->SetApicalParameters(10, 10, 1);
        p_force3->SetBasalParameters(10, 10, 1);
        p_force3->SetLateralParameter(4);
        p_force3->SetVolumeParameters(200, 1);
        p_force3->SetPatternedApicalParameter(20, 20, 0.5);
        simulator.AddForce(p_force3);

        simulator.Solve();

        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 100u);
        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), end_time, 1e-10);
    }
};

#endif /*TESTAPICALCONSTRICTIONEXAMPLE_HPP_*/
