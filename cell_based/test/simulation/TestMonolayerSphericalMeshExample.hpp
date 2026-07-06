/*

Copyright (c) 2005-2026, University of Oxford.
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

#ifndef TESTMONOLAYERSPHERICALMESHEXAMPLE_HPP_
#define TESTMONOLAYERSPHERICALMESHEXAMPLE_HPP_

#include "AbstractCellBasedTestSuite.hpp"

#include "GeodesicSphere23Generator.hpp"
#include "MonolayerVertexMeshGenerator.hpp"

#include "CellsGenerator.hpp"
#include "TransitCellProliferativeType.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "VertexBasedCellPopulation.hpp"

#include "SmartPointers.hpp"
#include "GeneralMonolayerVertexMeshForce.hpp"
#include "OffLatticeSimulation.hpp"

#include <string>

#include "FakePetscSetup.hpp"

class TestMonolayerSphericalMeshExample : public AbstractCellBasedTestSuite
{
private:
    static constexpr double target_area = 1;
    static constexpr double end_time = 10;

public:
    void TestOnSphere1()
    {
        /*
        const std::string output_filename = "TestMonolayerSphericalMeshExample/Sphere1_12cells";
        GeodesicSphere23Generator builder;

        MutableVertexMesh<2, 3>* p_dual_mesh = builder.GetDual();
        VertexMeshWriter<2, 3> Writer(output_filename, "Geodesic_Dual", false);
        Writer.WriteVtkUsingMesh(*p_dual_mesh);

        [[maybe_unused]] const unsigned radius = sqrt(p_dual_mesh->GetNumElements() * target_area / 4 / M_PI);
        MonolayerVertexMeshGenerator sBuilder;
        MutableVertexMesh<3, 3>* p_mesh = sBuilder.MakeSphericalMesh33(p_dual_mesh, 5, 0.5);
        sBuilder.WriteVtk(output_filename, "InitialMesh");

        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<UniformG1GenerationalCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_transit_type);
        VertexBasedCellPopulation<3> cell_population(*p_mesh, cells);

        OffLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory(output_filename);
        simulator.SetSamplingTimestepMultiple(10);
        simulator.SetEndTime(end_time);

        MAKE_PTR(GeneralMonolayerVertexMeshForce, p_force3);
        p_force3->SetApicalParameters(20, 20, 0.7);
        p_force3->SetBasalParameters(20, 20, 0.7);
        p_force3->SetLateralParameter(8);
        p_force3->SetVolumeParameters(350, 2);
        simulator.AddForce(p_force3);
        // MAKE_PTR(HorizontalStretchForce<3>, p_force2);
        // p_force2->SetForceMagnitude(1.0);
        // p_force2->SetRelativeWidth(0.15);
        // simulator.AddForce(p_force2);

        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 176u);
        simulator.Solve();

        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 176u);
        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), end_time, 1e-10);
        */
    }

    void TestOnSphere2()
    {
        /*
        const std::string output_filename = "TestMonolayerSphericalMeshExample/Sphere2_42cells";
        GeodesicSphere23Generator builder;
        builder.SubDivide();

        MutableVertexMesh<2, 3>* p_dual_mesh = builder.GetDual();
        VertexMeshWriter<2, 3> Writer(output_filename, "Geodesic_Dual", false);
        Writer.WriteVtkUsingMesh(*p_dual_mesh);

        [[maybe_unused]] const unsigned radius = sqrt(p_dual_mesh->GetNumElements() * target_area / 4 / M_PI);
        MonolayerVertexMeshGenerator sBuilder;
        MutableVertexMesh<3, 3>* p_mesh = sBuilder.MakeSphericalMesh33(p_dual_mesh, 5, 0.5);
        sBuilder.WriteVtk(output_filename, "InitialMesh");

        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<UniformG1GenerationalCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_transit_type);
        VertexBasedCellPopulation<3> cell_population(*p_mesh, cells);

        OffLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory(output_filename);
        simulator.SetSamplingTimestepMultiple(10);
        simulator.SetEndTime(end_time);

        MAKE_PTR(GeneralMonolayerVertexMeshForce, p_force3);
        p_force3->SetApicalParameters(20, 20, 0.7);
        p_force3->SetBasalParameters(20, 20, 0.7);
        p_force3->SetLateralParameter(8);
        p_force3->SetVolumeParameters(350, 2);
        simulator.AddForce(p_force3);
        // MAKE_PTR(HorizontalStretchForce<3>, p_force2);
        // p_force2->SetForceMagnitude(1.0);
        // p_force2->SetRelativeWidth(0.15);
        // simulator.AddForce(p_force2);

        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 176u);
        simulator.Solve();

        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 176u);
        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), end_time, 1e-10);
        */
    }

    void TestOnSphere3()
    {
        const std::string output_filename = "TestMonolayerSphericalMeshExample/Sphere3_162cells";
        GeodesicSphere23Generator builder;
        builder.SubDivide();
        builder.SubDivide();

        MutableVertexMesh<2, 3>* p_dual_mesh = builder.GetDual();
        VertexMeshWriter<2, 3> Writer(output_filename, "Geodesic_Dual", false);
        Writer.WriteVtkUsingMesh(*p_dual_mesh);

        [[maybe_unused]] const unsigned radius = sqrt(p_dual_mesh->GetNumElements() * target_area / 4 / M_PI);
        MonolayerVertexMeshGenerator sBuilder;
        MutableVertexMesh<3, 3>* p_mesh = sBuilder.MakeSphericalMesh33(p_dual_mesh, 5, 0.5);
        sBuilder.WriteVtk(output_filename, "InitialMesh");

        // \todo #2850 Disable asynchronous T1 swaps defensively: this test uses a proliferative
        // cell type, so should a division ever coincide with a scutoid interface node left by an
        // async swap, DivideElementAlongGivenAxis() would fail in GetOppositeNode(). Changes
        // nothing while no division completes. See TestMonolayerCylindricalMeshExample for the
        // full explanation and the "approach B" that would let scutoids and division coexist.
        p_mesh->SetAllowAsynchronousT1Swaps(false);

        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<UniformG1GenerationalCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_transit_type);
        VertexBasedCellPopulation<3> cell_population(*p_mesh, cells);

        OffLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory(output_filename);
        simulator.SetSamplingTimestepMultiple(10);
        simulator.SetEndTime(end_time);

        MAKE_PTR(GeneralMonolayerVertexMeshForce, p_force3);
        p_force3->SetApicalParameters(20, 20, 0.7);
        p_force3->SetBasalParameters(20, 20, 0.7);
        p_force3->SetLateralParameter(8);
        p_force3->SetVolumeParameters(350, 2);
        simulator.AddForce(p_force3);
        // MAKE_PTR(HorizontalStretchForce<3>, p_force2);
        // p_force2->SetForceMagnitude(1.0);
        // p_force2->SetRelativeWidth(0.15);
        // simulator.AddForce(p_force2);

        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 162u);
        simulator.Solve();

        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 162u);
        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), end_time, 1e-10);
    }
};

#endif /*TESTMONOLAYERSPHERICALMESHEXAMPLE_HPP_*/
