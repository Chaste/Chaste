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

#ifndef TESTVERTEXMESH33UNIAXIALLOAD_HPP_
#define TESTVERTEXMESH33UNIAXIALLOAD_HPP_

#include "AbstractCellBasedTestSuite.hpp"

#include "VoronoiVertexMeshGenerator.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "MonolayerVertexMeshGenerator.hpp"

#include "GeodesicSphere23Generator.hpp"

#include "CellsGenerator.hpp"
#include "NoCellCycleModel.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "SmartPointers.hpp"
#include "GeneralMonolayerVertexMeshForce.hpp"
#include "HorizontalStretchForce.hpp"
#include "OffLatticeSimulation.hpp"
#include <boost/lexical_cast.hpp>

#include "FakePetscSetup.hpp"

#include "Debug.hpp"

#include "TransitCellProliferativeType.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"

class TestVertexMesh33UniaxialLoad : public AbstractCellBasedTestSuite
{
private:
    static const double z_height = 1;
    static const double target_area = 1;
    const unsigned num_cells_x = 9;
    const unsigned num_cells_y = 5;
    static const double end_time = 1;

public:
    void TestOnHexagonalMesh() throw (Exception)
    {
        std::string output_filename = "TestUniaxialLoad/HoneyTest" + boost::lexical_cast<std::string>(num_cells_x)
                                    + "x" + boost::lexical_cast<std::string>(num_cells_y);
        HoneycombVertexMeshGenerator generator(num_cells_x, num_cells_y, false, 0.1, 0.01, target_area);
        MutableVertexMesh<2, 2>& vertex_2mesh = *(generator.GetMesh());
        MonolayerVertexMeshGenerator builder;
        MutableVertexMesh<3, 3>* p_mesh = builder.MakeMeshUsing2dMesh(vertex_2mesh, z_height);
        builder.WriteVtk(output_filename, "InitialMesh");

        std::vector<CellPtr> cells;
        CellsGenerator<NoCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
        VertexBasedCellPopulation<3> cell_population(*p_mesh, cells);

        OffLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory(output_filename);
        simulator.SetSamplingTimestepMultiple(10);
        simulator.SetEndTime(end_time);

        MAKE_PTR(GeneralMonolayerVertexMeshForce, p_force3);
//        p_force3->SetApicalParameters(20, 20, 0.7);
//        p_force3->SetBasalParameters(20, 20, 0.7);
//        p_force3->SetLateralParameter(8);
        p_force3->SetVolumeParameters(350, 0.1);
        simulator.AddForce(p_force3);
        MAKE_PTR(HorizontalStretchForce<3>, p_force2);
        p_force2->SetForceMagnitude(1.0);
        p_force2->SetRelativeWidth(0.15);
//        simulator.AddForce(p_force2);

        simulator.Solve();

        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), num_cells_x*num_cells_y);
        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), end_time, 1e-10);
    }

    void TestOnVoronoiMesh() throw (Exception)
    {
        std::string output_filename = "TestUniaxialLoad/VoronoiTest" + boost::lexical_cast<std::string>(num_cells_x)
                                    + "x" + boost::lexical_cast<std::string>(num_cells_y);
        VoronoiVertexMeshGenerator generator(num_cells_x, num_cells_y, 5, target_area);
        MutableVertexMesh<2, 2>& vertex_2mesh = *(generator.GetMesh());
        MonolayerVertexMeshGenerator builder;
        MutableVertexMesh<3, 3>* p_mesh = builder.MakeMeshUsing2dMesh(vertex_2mesh, z_height);
        builder.WriteVtk(output_filename, "InitialMesh");

        std::vector<CellPtr> cells;
        CellsGenerator<NoCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
        VertexBasedCellPopulation<3> cell_population(*p_mesh, cells);

        OffLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory(output_filename);
        simulator.SetSamplingTimestepMultiple(10);
        simulator.SetEndTime(end_time);

        MAKE_PTR(GeneralMonolayerVertexMeshForce, p_force3);
        p_force3->SetApicalParameters(20, 20, 0.7);
        p_force3->SetBasalParameters(20, 20, 0.7);
        p_force3->SetLateralParameter(8);
        p_force3->SetVolumeParameters(350, 1);
        simulator.AddForce(p_force3);
        MAKE_PTR(HorizontalStretchForce<3>, p_force2);
        p_force2->SetForceMagnitude(1.0);
        p_force2->SetRelativeWidth(0.15);
        simulator.AddForce(p_force2);

        simulator.Solve();

        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), num_cells_x*num_cells_y);
        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), end_time, 1e-10);
    }

//    void TestCylindricalMesh()
//    {
//        const unsigned x = 4;
//        const unsigned y = 4;
//        const double a = 2;
//        const double length = 3*sqrt(3)*y+sqrt(3);
//        const double radius = a/M_PI/2*x;
//        HoneycombVertexMeshGenerator generator(x, y, false, 0.1, 0.01, 2*sqrt(3));
//        MutableVertexMesh<2, 2>& vertex_2mesh = *(generator.GetMesh());
//        MonolayerVertexMeshGenerator builder;
//        MutableVertexMesh<3, 3>* p_mesh = builder.MakeMeshUsing2dMesh(vertex_2mesh);
//        std::string output_filename = "TestUniaxialLoad/CylinderTest_SmallVolume" + boost::lexical_cast<std::string>(x)
//                                            + "x" + boost::lexical_cast<std::string>(y);
//        builder.WriteVtk(output_filename, "InitialMesh");
//
//        builder.ConvertMeshToCylinder(2*x, 1, radius, 1, 1);
//        builder.WriteVtk(output_filename, "After");
//
//        std::vector<CellPtr> cells;
//        CellsGenerator<NoCellCycleModel, 3> cells_generator;
//        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
//        VertexBasedCellPopulation<3> cell_population(*p_mesh, cells);
//
//        OffLatticeSimulation<3> simulator(cell_population);
//        simulator.SetOutputDirectory(output_filename);
//        simulator.SetSamplingTimestepMultiple(10);
//        simulator.SetEndTime(end_time);
//
//        MAKE_PTR(GeneralMonolayerVertexMeshForce, p_force3);
//        p_force3->SetApicalParameters(20, 20, 0.7);
//        p_force3->SetBasalParameters(20, 20, 0.7);
//        p_force3->SetLateralParameter(8);
//        p_force3->SetVolumeParameters(100, 1);
//        simulator.AddForce(p_force3);
//        MAKE_PTR(HorizontalStretchForce<3>, p_force2);
//        p_force2->SetForceMagnitude(1.0);
//        p_force2->SetRelativeWidth(0.15);
////        simulator.AddForce(p_force2);
//
//        simulator.Solve();
//
//        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), x*y);
//        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), end_time, 1e-10);
//    }
//
//    void TestCylindricalMesh2()
//    {
//        const unsigned x = 4;
//        const unsigned y = 4;
//        const double a = 2;
//        const double length = 3*sqrt(3)*y+sqrt(3);
//        const double radius = a/M_PI/2*x;
//        HoneycombVertexMeshGenerator generator(x, y, false, 0.1, 0.01, 2*sqrt(3));
//        MutableVertexMesh<2, 2>& vertex_2mesh = *(generator.GetMesh());
//        MonolayerVertexMeshGenerator builder;
//        MutableVertexMesh<3, 3>* p_mesh = builder.MakeMeshUsing2dMesh(vertex_2mesh);
//        std::string output_filename = "TestUniaxialLoad/CylinderTest_Volume" + boost::lexical_cast<std::string>(x)
//                                                        + "x" + boost::lexical_cast<std::string>(y);
//        builder.WriteVtk(output_filename, "InitialMesh");
//
//        builder.ConvertMeshToCylinder(2*x, 1, radius, 1, 1);
//        builder.WriteVtk(output_filename, "After");
//
//        std::vector<CellPtr> cells;
//        CellsGenerator<NoCellCycleModel, 3> cells_generator;
//        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
//        VertexBasedCellPopulation<3> cell_population(*p_mesh, cells);
//
//        OffLatticeSimulation<3> simulator(cell_population);
//        simulator.SetOutputDirectory(output_filename);
//        simulator.SetSamplingTimestepMultiple(10);
//        simulator.SetEndTime(end_time);
//
//        MAKE_PTR(GeneralMonolayerVertexMeshForce, p_force3);
////        p_force3->SetApicalParameters(20, 20, 0.7);
////        p_force3->SetBasalParameters(20, 20, 0.7);
////        p_force3->SetLateralParameter(8);
//        p_force3->SetVolumeParameters(10, 0.5);
//        simulator.AddForce(p_force3);
//        MAKE_PTR(HorizontalStretchForce<3>, p_force2);
//        p_force2->SetForceMagnitude(1.0);
//        p_force2->SetRelativeWidth(0.15);
//        //        simulator.AddForce(p_force2);
//
//        simulator.Solve();
//
//        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), x*y);
//        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), end_time, 1e-10);
//    }
//
//    void TestCylindricalMesh3()
//    {
//        const unsigned x = 4;
//        const unsigned y = 4;
//        const double a = 2;
//        const double length = 3*sqrt(3)*y+sqrt(3);
//        const double radius = a/M_PI/2*x;
//        HoneycombVertexMeshGenerator generator(x, y, false, 0.1, 0.01, 2*sqrt(3));
//        MutableVertexMesh<2, 2>& vertex_2mesh = *(generator.GetMesh());
//        MonolayerVertexMeshGenerator builder;
//        MutableVertexMesh<3, 3>* p_mesh = builder.MakeMeshUsing2dMesh(vertex_2mesh);
//        std::string output_filename = "TestUniaxialLoad/CylinderTest_SmallVolume3" + boost::lexical_cast<std::string>(x)
//                                                        + "x" + boost::lexical_cast<std::string>(y);
//        builder.WriteVtk(output_filename, "InitialMesh");
//
//        builder.ConvertMeshToCylinder(2*x, 1, radius, 1, 1);
//        builder.WriteVtk(output_filename, "After");
//
//        std::vector<CellPtr> cells;
//        CellsGenerator<NoCellCycleModel, 3> cells_generator;
//        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
//        VertexBasedCellPopulation<3> cell_population(*p_mesh, cells);
//
//        OffLatticeSimulation<3> simulator(cell_population);
//        simulator.SetOutputDirectory(output_filename);
//        simulator.SetSamplingTimestepMultiple(10);
//        simulator.SetEndTime(end_time);
//
//        MAKE_PTR(GeneralMonolayerVertexMeshForce, p_force3);
//        p_force3->SetApicalParameters(20, 20, 0.7);
//        p_force3->SetBasalParameters(20, 20, 0.7);
//        p_force3->SetLateralParameter(8);
//        p_force3->SetVolumeParameters(100, 0.1);
//        simulator.AddForce(p_force3);
//        MAKE_PTR(HorizontalStretchForce<3>, p_force2);
//        p_force2->SetForceMagnitude(1.0);
//        p_force2->SetRelativeWidth(0.15);
//        //        simulator.AddForce(p_force2);
//
//        simulator.Solve();
//
//        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), x*y);
//        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), end_time, 1e-10);
//    }

    void TestCheckForceCalculation()
    {
        const unsigned x = 4;
        const unsigned y = 5;
        const double a = 2;
        const double length = 3*sqrt(3)*y+sqrt(3);
        const double radius = a/M_PI/2*x;
        HoneycombVertexMeshGenerator generator(x, y, false, 0.1, 0.01, 2*sqrt(3));
        MutableVertexMesh<2, 2>& vertex_2mesh = *(generator.GetMesh());
        MonolayerVertexMeshGenerator builder;
        MutableVertexMesh<3, 3>* p_mesh = builder.MakeMeshUsing2dMesh(vertex_2mesh);
        std::string output_filename = "TestUniaxialLoad/CylinderTest_CheckForce" + boost::lexical_cast<std::string>(x)
                                                                    + "x" + boost::lexical_cast<std::string>(y);
        builder.WriteVtk(output_filename, "InitialMesh");

        builder.ConvertMeshToCylinder(2*x, 1, radius, 0.5, 1);
        builder.WriteVtk(output_filename, "After");
PRINT_VECTOR(p_mesh->GetNode(0)->rGetLocation())
PRINT_VECTOR(p_mesh->GetNode(48)->rGetLocation())

        std::vector<CellPtr> cells;
        CellsGenerator<NoCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
        VertexBasedCellPopulation<3> cell_population(*p_mesh, cells);

        OffLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory(output_filename);
        simulator.SetSamplingTimestepMultiple(1);
        simulator.SetEndTime(1);

        MAKE_PTR(GeneralMonolayerVertexMeshForce, p_force3);
        p_force3->SetApicalParameters(20, 20, 0.7);
        p_force3->SetBasalParameters(20, 20, 0.7);
        p_force3->SetLateralParameter(20);
        p_force3->SetVolumeParameters(350, 3);
        simulator.AddForce(p_force3);
        MAKE_PTR(HorizontalStretchForce<3>, p_force2);
        p_force2->SetForceMagnitude(1.0);
        p_force2->SetRelativeWidth(0.15);
        //        simulator.AddForce(p_force2);

        simulator.Solve();

PRINT_VECTOR(p_mesh->GetNode(0)->rGetLocation())
PRINT_VECTOR(p_mesh->GetNode(48)->rGetLocation())

        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), x*y);
        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), end_time, 1e-10);
    }

    void TestCellGrowth() throw (Exception)
    {
        // Make a mesh of 10x10
        //        const double z_height = 1;
        const double target_area = 1;
        const unsigned num_cells_x = 6;
        const unsigned num_cells_y = 3;
        std::string output_filename = "TestCellDivision/HoneyTest" + boost::lexical_cast<std::string>(num_cells_x)
                                            + "x" + boost::lexical_cast<std::string>(num_cells_y);
        // There seems to be a bug somewhere in voronoiprism3dVertexMeshGenerator....
        //        VoronoiPrism3dVertexMeshGenerator generator(num_cells_x, num_cells_y, z_height, 5, target_area);
        //        MutableVertexMesh<3,3>* p_mesh = generator.GetMeshAfterReMesh();
        HoneycombVertexMeshGenerator generator(num_cells_x, num_cells_y, false, 0.1, 0.01, target_area);
        //        VoronoiVertexMeshGenerator generator(num_cells_x, num_cells_y, 5, target_area);
        MutableVertexMesh<2, 2>& vertex_2mesh = *(generator.GetMesh());
        MonolayerVertexMeshGenerator builder("CellGrowth");
        MutableVertexMesh<3, 3>* p_mesh = builder.MakeMeshUsing2dMesh(vertex_2mesh);
        builder.WriteVtk(output_filename,"Before");

        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<UniformG1GenerationalCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_transit_type);

        VertexBasedCellPopulation<3> cell_population(*p_mesh, cells);

        OffLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory(output_filename);
        simulator.SetSamplingTimestepMultiple(200);
        simulator.SetEndTime(20.0);

        MAKE_PTR(GeneralMonolayerVertexMeshForce, p_force3);
        p_force3->SetApicalParameters(20, 20, 0.7);
        p_force3->SetBasalParameters(20, 20, 0.7);
        p_force3->SetLateralParameter(8);
        p_force3->SetVolumeParameters(350, 1);
        simulator.AddForce(p_force3);
        MAKE_PTR(HorizontalStretchForce<3>, p_force2);
        p_force2->SetForceMagnitude(0.2);
        p_force2->SetRelativeWidth(0.15);
        simulator.AddForce(p_force2);

        /* Finally we call the {{{Solve}}} method on the simulation to run the simulation.*/
        simulator.Solve();
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 55u);
        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), 20.0, 1e-10);
    }


    void TestOnSphere() throw (Exception)
   {
        const double s_end_time = 1;
        std::string output_filename = "TestUniaxialLoad/SphereTest" + boost::lexical_cast<std::string>(num_cells_x)
                                    + "x" + boost::lexical_cast<std::string>(num_cells_y);
        GeodesicSphere23Generator builder;
        builder.SubDivide();    // n=42
        builder.SubDivide();    // n=162


        MutableVertexMesh<2,3>* p_dual_mesh = builder.GetDual();
        VertexMeshWriter<2, 3> Writer(output_filename, "Geodesic_Dual", false);
        Writer.WriteVtkUsingMeshWithCellId(*p_dual_mesh);

        const unsigned radius = sqrt(p_dual_mesh->GetNumElements()*target_area/4/M_PI);
        MonolayerVertexMeshGenerator sBuilder;
        MutableVertexMesh<3,3>* p_mesh = sBuilder.MakeSphericalMesh33(p_dual_mesh, 5, 0.5);
        sBuilder.WriteVtk(output_filename, "InitialMesh");


        std::vector<CellPtr> cells;
        CellsGenerator<NoCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
        VertexBasedCellPopulation<3> cell_population(*p_mesh, cells);

        OffLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory(output_filename);
        simulator.SetSamplingTimestepMultiple(10);
        simulator.SetEndTime(s_end_time);

        MAKE_PTR(GeneralMonolayerVertexMeshForce, p_force3);
        p_force3->SetApicalParameters(20, 20, 0.7);
        p_force3->SetBasalParameters(20, 20, 0.7);
        p_force3->SetLateralParameter(8);
        p_force3->SetVolumeParameters(350, 1);
        simulator.AddForce(p_force3);
//        MAKE_PTR(HorizontalStretchForce<3>, p_force2);
//        p_force2->SetForceMagnitude(1.0);
//        p_force2->SetRelativeWidth(0.15);
//        simulator.AddForce(p_force2);

        simulator.Solve();

        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 162u);
        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), s_end_time, 1e-10);
   }

};

#endif /*TESTVERTEXMESH33UNIAXIALLOAD_HPP_*/
