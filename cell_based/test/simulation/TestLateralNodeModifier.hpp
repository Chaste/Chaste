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

#include "Debug.hpp"
class TestLateralNodeModifier : public AbstractCellBasedTestSuite
{
private:
    void WriteMesh(const MutableVertexMesh<3, 3>& rMesh, const std::string& additionalTag)
    {
        VertexMeshWriter<3, 3> mesh_writer("TestAsynchronousT1", "mesh", false);
        mesh_writer.WriteVtkUsingMesh(rMesh, additionalTag);
    }
    static const double target_area = 1;
    static const double z_height = 1;

public:
    void TestUpdateCellData()
    {
        /**
         * Use 2 hexagonal cells as model.
         */

        HexagonalPrism3dVertexMeshGenerator generator(2, 2, 3 * sqrt(3) / 2, 1);
        MutableVertexMesh<3, 3>* p_mesh = generator.GetMesh();
        WriteMesh(*p_mesh, "ori0");

        PRINT_CONTAINER(p_mesh->GetNode(4)->rGetModifiableLocation());
        PRINT_CONTAINER(p_mesh->GetNode(7)->rGetModifiableLocation());

        p_mesh->GetNode(3)->rGetModifiableLocation()[1] = 0.95;
        p_mesh->GetNode(6)->rGetModifiableLocation()[1] = 1.05;
        p_mesh->SetCellRearrangementThreshold(0.2);
        p_mesh->SetCellRearrangementRatio(1.5);
        WriteMesh(*p_mesh, "ori1");

        p_mesh->ReMesh();
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 33u);
        TS_ASSERT_EQUALS(p_mesh->GetNumFaces(), 28u);
        WriteMesh(*p_mesh, "ori2");

        std::vector<CellPtr> cells;
        CellsGenerator<NoCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
        VertexBasedCellPopulation<3> cell_population(*p_mesh, cells);

        LateralNodeModifier node_modifier;

        node_modifier.UpdateCellData(cell_population);
        WriteMesh(*p_mesh, "ori3");

        p_mesh->GetNode(3 + 32)->rGetModifiableLocation()[1] = 0.9;
        p_mesh->GetNode(6 + 32)->rGetModifiableLocation()[1] = 1.1;
        p_mesh->SetCellRearrangementThreshold(0.3);
        WriteMesh(*p_mesh, "ori4");

        node_modifier.UpdateCellData(cell_population);
        WriteMesh(*p_mesh, "ori5");

        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 33u);
        TS_ASSERT_EQUALS(p_mesh->GetNumFaces(), 28u);
    }

    void TestAsynchronousT1()
    {
        /*
         * Using 3 rows of hexagonal cell to test. See output files.
         */
        std::string output_filename = "TestAsynchronousT1/SecondTest";
        const double end_time = 1;

        HexagonalPrism3dVertexMeshGenerator generator(4, 3, 3 * sqrt(3) / 2, 1);
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
        p_force3->SetLateralParameter(0, 2.5);
        p_force3->SetVolumeParameters(400, target_area * z_height);
        simulator.AddForce(p_force3);

        MAKE_PTR(LateralNodeModifier, p_node_modifier);
        simulator.AddSimulationModifier(p_node_modifier);

        MAKE_PTR(HorizontalStretchForce<3>, p_force2);
        p_force2->SetForceMagnitude(1);
        p_force2->SetRelativeWidth(0.15);
        simulator.AddForce(p_force2);

        simulator.Solve();

        CellBasedSimulationArchiver<3, OffLatticeSimulation<3> >::Save(&simulator);

        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 12u);
        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), end_time, 1e-10);
    }
};

#endif /*TESTLATERALNODEMODIFIER_HPP_*/
