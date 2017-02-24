/*
 * TestMonolayerCylindricalMeshExample.hpp
 *
 *  Created on: 24 Feb 2017
 *      Author: Weijie
 */

#ifndef TESTMONOLAYERCYLINDRICALMESHEXAMPLE_HPP_
#define TESTMONOLAYERCYLINDRICALMESHEXAMPLE_HPP_

#include "AbstractCellBasedTestSuite.hpp"

#include "HoneycombVertexMeshGenerator.hpp"
#include "MonolayerVertexMeshGenerator.hpp"

#include "CellsGenerator.hpp"
#include "TransitCellProliferativeType.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "VertexBasedCellPopulation.hpp"

#include "SmartPointers.hpp"
#include "GeneralMonolayerVertexMeshForce.hpp"
#include "HorizontalStretchForce.hpp"
#include "OffLatticeSimulation.hpp"

#include <string>
#include <sstream>

#include "FakePetscSetup.hpp"

class TestMonolayerCylindricalMeshExample : public AbstractCellBasedTestSuite
{
private:
    static const unsigned x = 10;
    static const unsigned y = 11;
    static const double target_area = 1;
    static const double end_time = 10;

public:
    void TestCylindricalMesh()
    {
        std::ostringstream oss;
        oss << "TestMonolayerCylindricalMeshExample/" << x  << "x" << y;
        const std::string output_filename(oss.str());
        const double a = 2;
        const double length = 3 * sqrt(3) * y + sqrt(3);
        const double radius = a / M_PI / 2 * x;
        HoneycombVertexMeshGenerator generator(x, y, false, 0.1, 0.01, 2 * sqrt(3));
        MutableVertexMesh<2, 2>& vertex_2mesh = *(generator.GetMesh());
        MonolayerVertexMeshGenerator builder;
        MutableVertexMesh<3, 3>* p_mesh = builder.MakeMeshUsing2dMesh(vertex_2mesh);
        builder.WriteVtk(output_filename, "InitialMesh");

        builder.ConvertMeshToCylinder(2 * x, 1, radius*0.8, 1.5, 1);

        for (unsigned i = 0; i < p_mesh->GetNumNodes(); ++i)
        {
            c_vector<double, 3>& tmp_loc = p_mesh->GetNode(i)->rGetModifiableLocation();
            double xx = tmp_loc[0];
            tmp_loc[0] = tmp_loc[1];
            tmp_loc[1] = -xx;
        }
        for (unsigned i = 0; i < p_mesh->GetNumElements(); ++i)
        {
            p_mesh->GetElement(i)->MonolayerElementRearrangeFacesNodes();
        }
        builder.WriteVtk(output_filename, "After");

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
        p_force3->SetApicalParameters(5, 5, 0.7);
        p_force3->SetBasalParameters(5, 5, 0.7);
        p_force3->SetLateralParameter(7);
        p_force3->SetVolumeParameters(100, 6);
        simulator.AddForce(p_force3);
        MAKE_PTR(HorizontalStretchForce<3>, p_force2);
        p_force2->SetForceMagnitude(1);
        p_force2->SetRelativeWidth(0.15);
        simulator.AddForce(p_force2);

        simulator.Solve();

        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), x * y);
        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), end_time, 1e-10);
    }
};

#endif /*TESTMONOLAYERCYLINDRICALMESHEXAMPLE_HPP_*/
