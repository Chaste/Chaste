/*
 * TestLateralNodeModifier.hpp
 *
 *  Created on: 24 Mar 2017
 *      Author: Weijie
 */

#ifndef TESTLATERALNODEMODIFIER_HPP_
#define TESTLATERALNODEMODIFIER_HPP_

#include "AbstractCellBasedTestSuite.hpp"

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include <string>
#include "CellsGenerator.hpp"
#include "GeneralMonolayerVertexMeshForce.hpp"
#include "HexagonalPrism3dVertexMeshGenerator.hpp"
#include "HorizontalStretchForce.hpp"
#include "LateralNodeModifier.hpp"
#include "NoCellCycleModel.hpp"
#include "OffLatticeSimulation.hpp"
#include "SmartPointers.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "VertexMeshWriter.hpp"
#include "FakePetscSetup.hpp"

class TestLateralNodeModifier : public AbstractCellBasedTestSuite
{
private:
    void WriteMesh(const MutableVertexMesh<3, 3>& rMesh, const std::string& additionalTag)
    {
        VertexMeshWriter<3, 3> mesh_writer("TestAsynchronousT1", "mesh", false);
        mesh_writer.WriteVtkUsingMesh(rMesh, additionalTag);
    }

public:
    void TestUpdateCellData()
    {
        /**
         * Use 2 hexagonal cells as model.
         */

        HexagonalPrism3dVertexMeshGenerator generator(2, 2, 3 * sqrt(3) / 2, 1);
        MutableVertexMesh<3, 3>* p_mesh = generator.GetMesh();
        WriteMesh(*p_mesh, "ori0");

        p_mesh->GetNode(3)->rGetModifiableLocation()[1] = 0.95;
        p_mesh->GetNode(6)->rGetModifiableLocation()[1] = 1.05;
        p_mesh->SetCellRearrangementThreshold(0.2);
        p_mesh->SetCellRearrangementRatio(1.5);
        WriteMesh(*p_mesh, "ori1");

        p_mesh->ReMesh();
        TS_ASSERT_EQUALS(p_mesh->GetNumFaces(), 28u);
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 33u);
        WriteMesh(*p_mesh, "ori2");

        std::vector<CellPtr> cells;
        CellsGenerator<NoCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
        VertexBasedCellPopulation<3> cell_population(*p_mesh, cells);

        LateralNodeModifier node_modifier;

        node_modifier.UpdateCellData(cell_population);
        WriteMesh(*p_mesh, "ori3");

        p_mesh->GetNode(3 + 16)->rGetModifiableLocation()[1] = 0.9;
        p_mesh->GetNode(6 + 16)->rGetModifiableLocation()[1] = 1.1;
        p_mesh->SetCellRearrangementThreshold(0.3);
        WriteMesh(*p_mesh, "ori4");

        node_modifier.UpdateCellData(cell_population);
        WriteMesh(*p_mesh, "ori5");

        TS_ASSERT_EQUALS(p_mesh->GetNumFaces(), 27u);
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 32u);
    }

    void TestAsynchronousT1()
    {
        /*
         * Using 3 rows of hexagonal cell to test. See output files.
         */
        std::string output_filename = "TestAsynchronousT1/SecondTest";
        const double end_time = 1;
        const double target_volume = 3 * sqrt(3) / 2;

        HexagonalPrism3dVertexMeshGenerator generator(5, 3, target_volume, 1);
        MutableVertexMesh<3, 3>* p_mesh = generator.GetMesh();

        p_mesh->GetNode(19)->rGetModifiableLocation()[1] = 2.45;
        p_mesh->GetNode(25)->rGetModifiableLocation()[1] = 2.55;
        p_mesh->SetCellRearrangementThreshold(0.30);
        p_mesh->DeleteElementPriorToReMesh(9u);

        p_mesh->ReMesh();
        TS_ASSERT_EQUALS(p_mesh->GetNumFaces(), 86u);
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 89u);
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
        p_force3->SetLateralParameter(0, 2.5);
        p_force3->SetVolumeParameters(400, target_volume);
        simulator.AddForce(p_force3);

        MAKE_PTR(LateralNodeModifier, p_node_modifier);
        simulator.AddSimulationModifier(p_node_modifier);

        MAKE_PTR_ARGS(HorizontalStretchForce<3>, p_force2, (1, 0.25));
        p_force2->SetUpPinnedElements(cell_population);
        simulator.AddForce(p_force2);

        simulator.Solve();

        CellBasedSimulationArchiver<3, OffLatticeSimulation<3> >::Save(&simulator);

        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 14u);
        TS_ASSERT_EQUALS(p_mesh->GetNumFaces(), 85u);
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 88u);
        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), end_time, 1e-10);
    }
};

#endif /*TESTLATERALNODEMODIFIER_HPP_*/
