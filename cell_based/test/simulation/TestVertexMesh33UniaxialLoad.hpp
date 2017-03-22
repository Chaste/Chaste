/*

Copyright (c) 2005-2017, University of Oxford.
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

#include "HexagonalPrism3dVertexMeshGenerator.hpp"
#include "MonolayerVertexMeshCustomFunctions.hpp"
#include "MonolayerVertexMeshGenerator.hpp"
#include "VoronoiPrism3dVertexMeshGenerator.hpp"

#include "GeodesicSphere23Generator.hpp"
#include "HoneycombVertexMeshGenerator.hpp"

#include <boost/lexical_cast.hpp>
#include "CellsGenerator.hpp"
#include "GeneralMonolayerVertexMeshForce.hpp"
#include "HorizontalStretchForce.hpp"
#include "LateralNodeModifier.hpp"
#include "NoCellCycleModel.hpp"
#include "OffLatticeSimulation.hpp"
#include "SmartPointers.hpp"
#include "VertexBasedCellPopulation.hpp"

#include "Debug.hpp"

#include "TransitCellProliferativeType.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"

#include "FakePetscSetup.hpp"

class BielmeierExternalForce : public AbstractForce<3>
{
    double mEcmSpringConstant;

    double mT_ext;

    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive& archive, const unsigned int version)
    {
        archive& boost::serialization::base_object<AbstractForce<3> >(*this);
    }

public:
    BielmeierExternalForce(const double springConstant = 0.0, const double T_ext = 0.0)
            : AbstractForce<3>(),
              mEcmSpringConstant(springConstant),
              mT_ext(T_ext)
    {
    }

    void AddForceContribution(AbstractCellPopulation<3>& rCellPopulation)
    {
        if (dynamic_cast<VertexBasedCellPopulation<3>*>(&rCellPopulation) == NULL)
        {
            EXCEPTION("GeneralMonolayerVertexMeshForce is to be used with a VertexBasedCellPopulation only");
        }

        // Define some helper variables
        VertexBasedCellPopulation<3>* p_cell_population = static_cast<VertexBasedCellPopulation<3>*>(&rCellPopulation);
        MutableVertexMesh<3, 3>& rMesh = p_cell_population->rGetMesh();

        for (unsigned i = 0; i < rMesh.GetNumNodes(); ++i)
        {
            Node<3>* p_node = rMesh.GetNode(i);

            if (IsBasalNode(p_node))
            {
                c_vector<double, 3> result = -1 * mEcmSpringConstant * Create_c_vector(0, 0, p_node->rGetLocation()[2]);
                p_node->AddAppliedForceContribution(result);
            }
        }

        for (unsigned i = 0; i < rMesh.GetNumFaces(); ++i)
        {
            VertexElement<2, 3>* p_face = rMesh.GetFace(i);

            if (IsApicalFace(p_face))
            {
                for (unsigned j = 0; j < p_face->GetNumNodes(); ++j)
                {
                    c_vector<double, 3> result = rMesh.GetAreaGradientOfFaceAtNode(p_face, j);
                    result *= mT_ext;
                    p_face->GetNode(j)->AddAppliedForceContribution(result);
                }
            }
        }
    }

    void OutputForceParameters(out_stream& rParamsFile)
    {
        *rParamsFile << "\t\t\t<EcmSpringConstant>" << mEcmSpringConstant << "</EcmSpringConstant>\n";
        *rParamsFile << "\t\t\t<T_ext>" << mT_ext << "</T_ext>\n";
        AbstractForce<3>::OutputForceParameters(rParamsFile);
    }
};
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(BielmeierExternalForce)

#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(BielmeierExternalForce)

class TestVertexMesh33UniaxialLoad : public AbstractCellBasedTestSuite
{
private:
    static const double z_height = 1;
    static const double target_area = 1;
    const unsigned num_cells_x = 9;
    const unsigned num_cells_y = 5;
    static const double end_time = 15;

public:
    void TestOnHexagonalMesh() throw(Exception)
    {
        std::string output_filename = "TestUniaxialLoad/HoneyTest" + boost::lexical_cast<std::string>(num_cells_x)
            + "x" + boost::lexical_cast<std::string>(num_cells_y);
        HexagonalPrism3dVertexMeshGenerator generator(num_cells_x, num_cells_y, target_area, z_height);
        MutableVertexMesh<3, 3>* p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        CellsGenerator<NoCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
        VertexBasedCellPopulation<3> cell_population(*p_mesh, cells);

        OffLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory(output_filename);
        simulator.SetSamplingTimestepMultiple(10);
        simulator.SetEndTime(end_time);

        MAKE_PTR(GeneralMonolayerVertexMeshForce, p_force3);
        p_force3->SetApicalParameters(2.5, 2.5, 0.7);
        p_force3->SetBasalParameters(2.5, 2.5, 0.7);
        p_force3->SetLateralParameter(1);
        p_force3->SetVolumeParameters(400, target_area * 1.2);
        simulator.AddForce(p_force3);
        MAKE_PTR(HorizontalStretchForce<3>, p_force2);
        p_force2->SetForceMagnitude(0.5);
        p_force2->SetRelativeWidth(0.15);
        p_force2->SetUpPinnedElements(cell_population);
        simulator.AddForce(p_force2);
        // MAKE_PTR(GeneralMonolayerVertexMeshForce, p_force3);
        // p_force3->SetApicalParameters(0.18, 3.1);
        // p_force3->SetBasalParameters(0.18, 6.95);
        // p_force3->SetLateralParameter(0, 1);
        // p_force3->SetVolumeParameters(1000, target_area);
        // simulator.AddForce(p_force3);

        // MAKE_PTR_ARGS(BielmeierExternalForce, p_force4, (5, -4.2));
        // simulator.AddForce(p_force4);

        // MAKE_PTR(HorizontalStretchForce<3>, p_force2);
        // p_force2->SetForceMagnitude(0.5);
        // p_force2->SetRelativeWidth(0.15);
        // simulator.AddForce(p_force2);

        MAKE_PTR(LateralNodeModifier, p_node_modifier);
        simulator.AddSimulationModifier(p_node_modifier);

        simulator.Solve();

        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), num_cells_x * num_cells_y);
        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), end_time, 1e-10);

        // loweredALLApical
        {
            // std::string output_filename = "TestUniaxialLoad/HoneyTest_loweredALLApical";
            // HoneycombVertexMeshGenerator generator(num_cells_x, num_cells_y, false, 0.1, 0.01, target_area);
            // MutableVertexMesh<2, 2> &vertex_2mesh = *(generator.GetMesh());
            // MonolayerVertexMeshGenerator builder;
            // MutableVertexMesh<3, 3> *p_mesh = builder.MakeMeshUsing2dMesh(vertex_2mesh, z_height);
            // builder.WriteVtk(output_filename, "InitialMesh");

            // std::vector<CellPtr> cells;
            // CellsGenerator<NoCellCycleModel, 3> cells_generator;
            // cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
            // VertexBasedCellPopulation<3> cell_population(*p_mesh, cells);
            
            // OffLatticeSimulation<3> simulator(cell_population);
            // simulator.SetOutputDirectory(output_filename);
            // simulator.SetSamplingTimestepMultiple(10);
            // simulator.SetEndTime(end_time);

            // MAKE_PTR(GeneralMonolayerVertexMeshForce, p_force3);
            // p_force3->SetApicalParameters(5, 5, 0.2);
            // p_force3->SetBasalParameters(20, 20, 0.7);
            // p_force3->SetLateralParameter(8);
            // p_force3->SetVolumeParameters(400, target_area * 1.2);
            // simulator.AddForce(p_force3);
            // MAKE_PTR(HorizontalStretchForce<3>, p_force2);
            // p_force2->SetForceMagnitude(0.5);
            // p_force2->SetRelativeWidth(0.15);
            // simulator.AddForce(p_force2);

            // simulator.Solve();

            // TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), num_cells_x * num_cells_y);
            // TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), end_time, 1e-10);
        }

        // LoweredApicalLineTension
        {
            // std::string output_filename = "TestUniaxialLoad/HoneyTest_LoweredApicalLineTension";
            // HoneycombVertexMeshGenerator generator(num_cells_x, num_cells_y, false, 0.1, 0.01, target_area);
            // MutableVertexMesh<2, 2> &vertex_2mesh = *(generator.GetMesh());
            // MonolayerVertexMeshGenerator builder;
            // MutableVertexMesh<3, 3> *p_mesh = builder.MakeMeshUsing2dMesh(vertex_2mesh, z_height);
            // builder.WriteVtk(output_filename, "InitialMesh");

            // std::vector<CellPtr> cells;
            // CellsGenerator<NoCellCycleModel, 3> cells_generator;
            // cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
            // VertexBasedCellPopulation<3> cell_population(*p_mesh, cells);
            
            // OffLatticeSimulation<3> simulator(cell_population);
            // simulator.SetOutputDirectory(output_filename);
            // simulator.SetSamplingTimestepMultiple(10);
            // simulator.SetEndTime(end_time);

            // MAKE_PTR(GeneralMonolayerVertexMeshForce, p_force3);
            // p_force3->SetApicalParameters(5, 20, 0.7);
            // p_force3->SetBasalParameters(20, 20, 0.7);
            // p_force3->SetLateralParameter(8);
            // p_force3->SetVolumeParameters(400, target_area * 1.2);
            // simulator.AddForce(p_force3);
            // MAKE_PTR(HorizontalStretchForce<3>, p_force2);
            // p_force2->SetForceMagnitude(0.5);
            // p_force2->SetRelativeWidth(0.15);
            // simulator.AddForce(p_force2);

            // simulator.Solve();

            // TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), num_cells_x * num_cells_y);
            // TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), end_time, 1e-10);
        }

        // LoweredApicalFaceParam
        {
            // std::string output_filename = "TestUniaxialLoad/HoneyTest_LoweredApicalFaceParam";
            // HoneycombVertexMeshGenerator generator(num_cells_x, num_cells_y, false, 0.1, 0.01, target_area);
            // MutableVertexMesh<2, 2> &vertex_2mesh = *(generator.GetMesh());
            // MonolayerVertexMeshGenerator builder;
            // MutableVertexMesh<3, 3> *p_mesh = builder.MakeMeshUsing2dMesh(vertex_2mesh, z_height);
            // builder.WriteVtk(output_filename, "InitialMesh");

            // std::vector<CellPtr> cells;
            // CellsGenerator<NoCellCycleModel, 3> cells_generator;
            // cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
            // VertexBasedCellPopulation<3> cell_population(*p_mesh, cells);
            
            // OffLatticeSimulation<3> simulator(cell_population);
            // simulator.SetOutputDirectory(output_filename);
            // simulator.SetSamplingTimestepMultiple(10);
            // simulator.SetEndTime(end_time);

            // MAKE_PTR(GeneralMonolayerVertexMeshForce, p_force3);
            // p_force3->SetApicalParameters(20, 5, 0.2);
            // p_force3->SetBasalParameters(20, 20, 0.7);
            // p_force3->SetLateralParameter(8);
            // p_force3->SetVolumeParameters(400, target_area * 1.2);
            // simulator.AddForce(p_force3);
            // MAKE_PTR(HorizontalStretchForce<3>, p_force2);
            // p_force2->SetForceMagnitude(0.5);
            // p_force2->SetRelativeWidth(0.15);
            // simulator.AddForce(p_force2);

            // simulator.Solve();

            // TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), num_cells_x * num_cells_y);
            // TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), end_time, 1e-10);
        }

        // AdjustedBasal (not yet implemented)
        {
            // std::string output_filename = "TestUniaxialLoad/HoneyTest_AdjustedBasal";
            // HoneycombVertexMeshGenerator generator(num_cells_x, num_cells_y, false, 0.1, 0.01, target_area);
            // MutableVertexMesh<2, 2> &vertex_2mesh = *(generator.GetMesh());
            // MonolayerVertexMeshGenerator builder;
            // MutableVertexMesh<3, 3> *p_mesh = builder.MakeMeshUsing2dMesh(vertex_2mesh, z_height);
            // builder.WriteVtk(output_filename, "InitialMesh");

            // std::vector<CellPtr> cells;
            // CellsGenerator<NoCellCycleModel, 3> cells_generator;
            // cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
            // VertexBasedCellPopulation<3> cell_population(*p_mesh, cells);

            // OffLatticeSimulation<3> simulator(cell_population);
            // simulator.SetOutputDirectory(output_filename);
            // simulator.SetSamplingTimestepMultiple(10);
            // simulator.SetEndTime(end_time);

            // MAKE_PTR(GeneralMonolayerVertexMeshForce, p_force3);
            // p_force3->SetApicalParameters(20, 20, 0.7);
            // p_force3->SetBasalParameters(20, 20, 0.7);
            // p_force3->SetLateralParameter(8);
            // p_force3->SetVolumeParameters(400, target_area * 1.2);
            // simulator.AddForce(p_force3);
            // MAKE_PTR(HorizontalStretchForce<3>, p_force2);
            // p_force2->SetForceMagnitude(0.5);
            // p_force2->SetRelativeWidth(0.15);
            // simulator.AddForce(p_force2);

            // simulator.Solve();

            // TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), num_cells_x * num_cells_y);
            // TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), end_time, 1e-10);
        }
    }

    void TestOnVoronoiMesh() throw(Exception)
    {
        std::string output_filename = "TestUniaxialLoad/VoronoiTest" + boost::lexical_cast<std::string>(num_cells_x)
            + "x" + boost::lexical_cast<std::string>(num_cells_y);
        VoronoiPrism3dVertexMeshGenerator generator(num_cells_x, num_cells_y, z_height, 5, target_area);
        MutableVertexMesh<3, 3>* p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        CellsGenerator<NoCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
        VertexBasedCellPopulation<3> cell_population(*p_mesh, cells);

        OffLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory(output_filename);
        simulator.SetSamplingTimestepMultiple(10);
        simulator.SetEndTime(end_time);

        MAKE_PTR(GeneralMonolayerVertexMeshForce, p_force3);
        p_force3->SetApicalParameters(2.5, 2.5, 0.7);
        p_force3->SetBasalParameters(2.5, 2.5, 0.7);
        p_force3->SetLateralParameter(1);
        p_force3->SetVolumeParameters(350, 1);
        simulator.AddForce(p_force3);
        MAKE_PTR(HorizontalStretchForce<3>, p_force2);
        p_force2->SetForceMagnitude(1.0);
        p_force2->SetRelativeWidth(0.15);
        simulator.AddForce(p_force2);

        simulator.Solve();

        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), num_cells_x * num_cells_y);
        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), end_time, 1e-10);
    }

    void TestCellGrowth() throw(Exception)
    {
        // Make a mesh of 10x10
        //        const double z_height = 1;
        const double target_area = 1;
        const unsigned num_cells_x = 6;
        const unsigned num_cells_y = 3;
        std::string output_filename = "TestCellDivision/HoneyTest" + boost::lexical_cast<std::string>(num_cells_x)
            + "x" + boost::lexical_cast<std::string>(num_cells_y);
        HexagonalPrism3dVertexMeshGenerator generator(num_cells_x, num_cells_y, target_area, z_height);
        MutableVertexMesh<3, 3>* p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<UniformG1GenerationalCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_transit_type);
        VertexBasedCellPopulation<3> cell_population(*p_mesh, cells);

        OffLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory(output_filename);
        simulator.SetSamplingTimestepMultiple(10);
        simulator.SetEndTime(end_time);
        p_mesh->SetCellRearrangementThreshold(0.1);

        MAKE_PTR(GeneralMonolayerVertexMeshForce, p_force3);
        p_force3->SetApicalParameters(2.5, 2.5, 0.7);
        p_force3->SetBasalParameters(2.5, 2.5, 0.7);
        p_force3->SetLateralParameter(1);
        p_force3->SetVolumeParameters(350, 1);
        simulator.AddForce(p_force3);
        MAKE_PTR(HorizontalStretchForce<3>, p_force2);
        p_force2->SetForceMagnitude(0.5);
        p_force2->SetRelativeWidth(0.15);
        simulator.AddForce(p_force2);

        simulator.Solve();
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 55u);
        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), end_time, 1e-10);
    }

    void TestSingleCube()
    {
        /*
         * Create a mesh with a square element, as shown below.
         *  _____
         * |     |
         * |     |
         * |_____|
         */

        std::string output_filename = "TestUniaxialLoad/SingleCubeTest";
        const double target_volume = 3;

        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, true, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, true, 1.0, 0.0));
        nodes.push_back(new Node<3>(2, true, 0.0, 1.0));
        nodes.push_back(new Node<3>(3, true, 1.0, 1.0));
        const unsigned node_indices_elem_0[4] = { 0, 1, 3, 2 };

        MonolayerVertexMeshGenerator builder(nodes, "TestGradientDeviation");
        builder.BuildElementWith(4, node_indices_elem_0);
        // A reference variable as mesh is noncopyable
        MutableVertexMesh<3, 3>* p_mesh = builder.GenerateMesh();

        std::vector<CellPtr> cells;
        CellsGenerator<NoCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
        VertexBasedCellPopulation<3> cell_population(*p_mesh, cells);

        OffLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory(output_filename);
        simulator.SetSamplingTimestepMultiple(10);
        simulator.SetEndTime(end_time);

        MAKE_PTR(GeneralMonolayerVertexMeshForce, p_force3);
        // p_force3->SetApicalParameters(20, 20, 0.7);
        // p_force3->SetBasalParameters(20, 20, 0.7);
        // p_force3->SetLateralParameter(8);
        p_force3->SetVolumeParameters(350, target_volume);
        simulator.AddForce(p_force3);

        simulator.Solve();

        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 1);
        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), end_time, 1e-10);
        TS_ASSERT_DELTA(p_mesh->GetVolumeOfElement(0), target_volume, 0.05);

        simulator.SetEndTime(end_time * 2);
        p_force3->SetVolumeParameters(350, target_volume / 2);
        simulator.Solve();
        TS_ASSERT_DELTA(p_mesh->GetVolumeOfElement(0), target_volume / 2, 0.05);
    }

    void TestAsynchronousT1()
    {
        /*
         * Using 3 rows of hexagonal cell to test. See output files.
         */
        std::string output_filename = "TestAsynchronousT1/SecondTest";
        const double end_time = 1;

        HexagonalPrism3dVertexMeshGenerator generator(4, 3, target_area, z_height);
        MutableVertexMesh<3, 3>* p_mesh = generator.GetMesh();

        p_mesh->GetNode(16)->rGetModifiableLocation()[1] = 1.4;
        p_mesh->GetNode(21)->rGetModifiableLocation()[1] = 1.7;
        p_mesh->SetCellRearrangementThreshold(0.40);

        p_mesh->ReMesh();
        p_mesh->SetCellRearrangementThreshold(0.01);
        p_mesh->SetCellRearrangementRatio(1.5);

        std::vector<CellPtr> cells;
        CellsGenerator<NoCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
        VertexBasedCellPopulation<3> cell_population(*p_mesh, cells);

        OffLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory(output_filename);
        simulator.SetSamplingTimestepMultiple(1);
        simulator.SetEndTime(end_time);

        MAKE_PTR(GeneralMonolayerVertexMeshForce, p_force3);
        p_force3->SetApicalParameters(2.5, 2.5, 0.7);
        p_force3->SetBasalParameters(2.5, 2.5, 0.7);
        p_force3->SetLateralParameter(1);
        p_force3->SetVolumeParameters(350, target_area * z_height);
        simulator.AddForce(p_force3);

        MAKE_PTR(LateralNodeModifier, p_node_modifier);
        simulator.AddSimulationModifier(p_node_modifier);

        MAKE_PTR(HorizontalStretchForce<3>, p_force2);
        p_force2->SetForceMagnitude(1);
        p_force2->SetRelativeWidth(0.15);
        simulator.AddForce(p_force2);

        simulator.Solve();

        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 12u);
        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), end_time, 1e-10);
    }
};

#endif /*TESTVERTEXMESH33UNIAXIALLOAD_HPP_*/
