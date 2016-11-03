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

#ifndef TESTFINDINGGENERALMONOLAYERVERTEXMESHFORCEPARAMETER_HPP_
#define TESTFINDINGGENERALMONOLAYERVERTEXMESHFORCEPARAMETER_HPP_

#include "Debug.hpp"
#include <cxxtest/TestSuite.h>
#include "AbstractCellBasedTestSuite.hpp"

#include "CheckpointArchiveTypes.hpp"
#include "AbstractForce.hpp"

#include "VoronoiPrism3dVertexMeshGenerator.hpp"
#include "VoronoiVertexMeshGenerator.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "CellsGenerator.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "NoCellCycleModel.hpp"
#include "SmartPointers.hpp"
#include "OffLatticeSimulation.hpp"

#include "FakePetscSetup.hpp"

#include "GeneralMonolayerVertexMeshForce.hpp"
#include "MeshBuilderHelper.hpp"

#define OUTPUT_NAME "TestUniaxialLoad"
#define ADDFORCEPARAMETER p_force3->SetVolumeParameter(0.8, 10);
#define CURRENT_TEST std::string("Volume2")
#define END_TIME 0.5

#include "HoneycombVertexMeshGenerator.hpp"

class TestFindingGeneralMonolayerVertexMeshForceParameter : public AbstractCellBasedTestSuite
{
public:

    void TestAll()
    {
        const double z_height = 1;
        const double target_area = 1;
        const unsigned num_cells_x = 5;
        const unsigned num_cells_y = 5;
        HoneycombVertexMeshGenerator generator(num_cells_x, num_cells_y, false, 0.1, 0.01, target_area);
        MutableVertexMesh<2, 2>& vertex_2mesh = *(generator.GetMesh());
        MeshBuilderHelper builder("TestForceParameter");

        std::vector<CellPtr> cells;
        CellsGenerator<NoCellCycleModel, 3> cells_generator;

MARK; TRACE("Volume vs Lateral")
        // Testing for lateral to volume parameter ratio
        char name_pattern_lateral_volume[] ("TestUniaxialLoad/TestForceParameter/VolumeVsLateral/param=%.2f");
        for (unsigned i=0; i<30; i+=3)
        {
            char tmp_name[100];
            sprintf(tmp_name, name_pattern_lateral_volume, double(i));

            MutableVertexMesh<3, 3>* p_mesh = builder.MakeNewMeshUsing2dMesh(vertex_2mesh, z_height);
            cells_generator.GenerateBasicRandom(cells, num_cells_x*num_cells_y);
            VertexBasedCellPopulation<3> cell_population(*p_mesh, cells);

            OffLatticeSimulation<3> simulator(cell_population);
            simulator.SetOutputDirectory(tmp_name);
            simulator.SetSamplingTimestepMultiple(10);
            const double end_time = 0.1;
            simulator.SetEndTime(end_time);

            MAKE_PTR(GeneralMonolayerVertexMeshForce<3>, p_force3);
//            p_force3->SetApicalParameter(10, 10, 1);
//            p_force3->SetBasalParameter(10, 10, 1);
            p_force3->SetLateralParameter(i);
            p_force3->SetVolumeParameter(100, z_height*target_area);
            simulator.AddForce(p_force3);

            simulator.Solve();

            double sum_elem_volume(0);
            std::vector<double> elem_volumes;
            for (unsigned elem_id=0; elem_id<cell_population.GetNumElements(); ++elem_id)
            {
                elem_volumes.push_back(cell_population.rGetMesh().GetVolumeOfElement(elem_id));
                sum_elem_volume += cell_population.rGetMesh().GetVolumeOfElement(elem_id);
            }
            PRINT_VARIABLE(sum_elem_volume/elem_volumes.size())
            PRINT_VECTOR(elem_volumes)

            TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), num_cells_x*num_cells_y);
            TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), end_time, 1e-10);

            SimulationTime::Destroy();
            SimulationTime::Instance()->SetStartTime(0.0);
        }

MARK; TRACE("Apical line vs apical area")
        // Testing for apical line to apical area
        char name_pattern_apical[] ("TestUniaxialLoad/TestForceParameter/ApicalLineVsArea/param=%.2f");
        for (unsigned i=1; i<30; i+=3)
        {
            char tmp_name[100];
            sprintf(tmp_name, name_pattern_apical, double(i));

            MutableVertexMesh<3, 3>* p_mesh = builder.MakeNewMeshUsing2dMesh(vertex_2mesh, z_height);
            cells_generator.GenerateBasicRandom(cells, num_cells_x*num_cells_y);
            VertexBasedCellPopulation<3> cell_population(*p_mesh, cells);

            OffLatticeSimulation<3> simulator(cell_population);
            simulator.SetOutputDirectory(tmp_name);
            simulator.SetSamplingTimestepMultiple(10);
            const double end_time = 1;
            simulator.SetEndTime(end_time);

            MAKE_PTR(GeneralMonolayerVertexMeshForce<3>, p_force3);
            p_force3->SetApicalParameter(15, i, 1);
            simulator.AddForce(p_force3);

            simulator.Solve();

            double sum_elem_volume(0);
            std::vector<double> elem_volumes;
            for (unsigned elem_id=0; elem_id<cell_population.GetNumElements(); ++elem_id)
            {
                elem_volumes.push_back(cell_population.rGetMesh().GetVolumeOfElement(elem_id));
                sum_elem_volume += cell_population.rGetMesh().GetVolumeOfElement(elem_id);
            }
            PRINT_VARIABLE(sum_elem_volume/elem_volumes.size())
            PRINT_VECTOR(elem_volumes)

            TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), num_cells_x*num_cells_y);
            TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), end_time, 1e-10);

            SimulationTime::Destroy();
            SimulationTime::Instance()->SetStartTime(0.0);
        }

        MARK; TRACE("Apical plus volume")
        // Testing for apical line to apical area
        char name_pattern_apical_volume[] ("TestUniaxialLoad/TestForceParameter/ApicalVsVolume/param=%.2f");
        for (unsigned i=1; i<30; i+=3)
        {
            char tmp_name[100];
            sprintf(tmp_name, name_pattern_apical_volume, double(i));

            MutableVertexMesh<3, 3>* p_mesh = builder.MakeNewMeshUsing2dMesh(vertex_2mesh, z_height);
            cells_generator.GenerateBasicRandom(cells, num_cells_x*num_cells_y);
            VertexBasedCellPopulation<3> cell_population(*p_mesh, cells);

            OffLatticeSimulation<3> simulator(cell_population);
            simulator.SetOutputDirectory(tmp_name);
            simulator.SetSamplingTimestepMultiple(10);
            const double end_time = 0.1;
            simulator.SetEndTime(end_time);

            MAKE_PTR(GeneralMonolayerVertexMeshForce<3>, p_force3);
            p_force3->SetApicalParameter(10, i, 1);
            p_force3->SetVolumeParameter(100, z_height*target_area);
            simulator.AddForce(p_force3);
            simulator.Solve();

            double sum_elem_volume(0);
            std::vector<double> elem_volumes;
            for (unsigned elem_id=0; elem_id<cell_population.GetNumElements(); ++elem_id)
            {
                elem_volumes.push_back(cell_population.rGetMesh().GetVolumeOfElement(elem_id));
                sum_elem_volume += cell_population.rGetMesh().GetVolumeOfElement(elem_id);
            }
            PRINT_VARIABLE(sum_elem_volume/elem_volumes.size())
            PRINT_VECTOR(elem_volumes)

            TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), num_cells_x*num_cells_y);
            TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), end_time, 1e-10);

            SimulationTime::Destroy();
            SimulationTime::Instance()->SetStartTime(0.0);
        }
    }
};

#endif /*TESTFINDINGGENERALMONOLAYERVERTEXMESHFORCEPARAMETER_HPP_*/
