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

#ifndef TESTVERTEXMESH33HORIZONTALSTRETCH_HPP_
#define TESTVERTEXMESH33HORIZONTALSTRETCH_HPP_

#include "Debug.hpp"
#include <cxxtest/TestSuite.h>
#include "AbstractCellBasedTestSuite.hpp"

#include "CheckpointArchiveTypes.hpp"
#include "AbstractForce.hpp"

#include "VoronoiPrism3dVertexMeshGenerator.hpp"
#include "CellsGenerator.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "NoCellCycleModel.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "SmartPointers.hpp"
#include "MisraForce.hpp"
#include "OffLatticeSimulation.hpp"

#include "FakePetscSetup.hpp"

template<unsigned DIM>
class HorizontalStretchForce : public AbstractForce<DIM>
{
private:

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<DIM> >(*this);
        archive & mForceMagnitude;
        archive & mRelativeWidth;
    }

    /**
     * Internal variable to save the magnitude
     */
    double mForceMagnitude;

    /**
     * The relative width from the boundary upon with the force is acting
     */
    double mRelativeWidth;

public:
    /**
     * default constructor
     */
    HorizontalStretchForce(const double ForceMagnitude=1, const double RelativeWidth=0.1)
    : mForceMagnitude(ForceMagnitude),
      mRelativeWidth(RelativeWidth)
    {
MARK
    }

    virtual ~HorizontalStretchForce()
    {
    }

    /**
     * Adding force contributions to the Populations
     * @param rCellPopulation
     */
    virtual void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
    {
MARK
        // Throw an exception message if not using a VertexBasedCellPopulation
        ///\todo check whether this line influences profiling tests - if so, we should remove it.
        if (dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation) == NULL)
        {
            EXCEPTION("MisraForce is to be used with a VertexBasedCellPopulation only");
        }

        // Define some helper variables
        VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);
        MutableVertexMesh<DIM, DIM>& rMesh = p_cell_population->rGetMesh();

        const unsigned num_nodes = p_cell_population->GetNumNodes();
        c_vector<double, DIM> force_in_positive_y = zero_vector<double>(DIM);
        force_in_positive_y[0] = mForceMagnitude;
        c_vector<double, DIM> force_in_negative_y = -force_in_positive_y;

        ChasteCuboid<DIM> bounding_box_tmp = rMesh.CalculateBoundingBox();
        const double min_x = bounding_box_tmp.rGetLowerCorner()[0];
        const double max_x = bounding_box_tmp.rGetUpperCorner()[0];
        assert(max_x > min_x);
        const double width_x = max_x - min_x;
        const double width = width_x*mRelativeWidth;
        const double left_bound = min_x + width;
        const double right_bound = max_x - width;
        // could have built two more Cuboids to check if the point is inside,
        // but will have a lot of wrapping which reduce efficiency.
TRACE("Going into the loop")
PRINT_VECTOR(force_in_positive_y)
PRINT_VECTOR(force_in_negative_y)
        // Iterate over nodes in the cell population
        for (unsigned node_global_index=0; node_global_index<num_nodes; ++node_global_index)
        {
            Node<3>* p_this_node = p_cell_population->GetNode(node_global_index);
//            if ( ! p_this_node->IsBoundaryNode())
//                continue;
            const double loc_x = p_this_node->rGetLocation()[0];
            if (loc_x < left_bound)
            {
TRACE("On the Left")
PRINT_VARIABLE(p_this_node->GetIndex())
                p_this_node->AddAppliedForceContribution(force_in_negative_y);
            }
            else if (loc_x > right_bound)
            {
TRACE("On the Right")
PRINT_VARIABLE(p_this_node->GetIndex())
                p_this_node->AddAppliedForceContribution(force_in_positive_y);
            }
        }
MARK
    }

    /**
     * For the compiler now
     * @param rParamsFile
     */
    virtual void OutputForceParameters(out_stream& rParamsFile)
    {
MARK
        *rParamsFile << "\t\t\t<ForceMagnitude>" << mForceMagnitude << "</ForceMagnitude>\n";
        *rParamsFile << "\t\t\t<RelativeWidth>" << mRelativeWidth << "</RelativeWidth>\n";

        // Call method on direct parent class
        AbstractForce<DIM>::OutputForceParameters(rParamsFile);
    }
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS1(HorizontalStretchForce,3)
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS1(HorizontalStretchForce,3)


class TestVertexMesh33HorizontalStretch : public AbstractCellBasedTestSuite
{
public:

    void TestVertexMesh33HorizontalStretch() throw (Exception)
    {
        // Make a mesh of 10x5
        const unsigned z_height = 1;
        const unsigned target_area = 1;
        VoronoiPrism3dVertexMeshGenerator generator(5, 5, z_height, 5, target_area);
        MutableVertexMesh<3,3>* p_mesh = generator.GetMeshAfterReMesh();
        VertexMeshWriter<3, 3> vertex_mesh_writer("TestVertexMesh33HorizontalStretch/Test1", "InitialMesh", false);
        vertex_mesh_writer.WriteVtkUsingMeshWithCellId(*p_mesh);

        std::vector<CellPtr> cells;
//        MAKE_PTR(DifferentiatedCellProliferativeType, p_transit_type);

        CellsGenerator<NoCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
        VertexBasedCellPopulation<3> cell_population(*p_mesh, cells);

//        VertexBasedCellPopulation<3> cell_population(*p_mesh, cells, false, false);
//        cell_population.InitialiseCells();


        OffLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory("TestUniaxialLoad/Test1");
        simulator.SetSamplingTimestepMultiple(10);
        const double end_time = 0.5;
        simulator.SetEndTime(end_time);

        MAKE_PTR(MisraForce<3>, p_force);
        p_force->SetLateralSurfaceEnergyParameter(0);
        p_force->SetApicalLineTensionParameter(0);
        p_force->SetBasalSurfaceEnergyParameter(0);
        p_force->SetTargetVolume(z_height*target_area*2);
        p_force->SetVolumeCompressibilityParameter(2);
        simulator.AddForce(p_force);
//        MAKE_PTR(HorizontalStretchForce<3>, p_force2);
//        simulator.AddForce(p_force2);
MARK
        /* Finally we call the {{{Solve}}} method on the simulation to run the simulation.*/
        simulator.Solve();
MARK
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 25u);
        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), end_time, 1e-10);
    }

    void TestVertexMesh33UiaxialLoad2() throw (Exception)
    {
        // Make a mesh of 10x5
        const unsigned z_height = 5;
        const unsigned target_area = 1;
        VoronoiPrism3dVertexMeshGenerator generator(7, 3, z_height, 5, target_area);
        MutableVertexMesh<3,3>* p_mesh = generator.GetMeshAfterReMesh();
        VertexMeshWriter<3, 3> vertex_mesh_writer("TestVertexMesh33HorizontalStretch/Test2", "InitialMesh", false);
        vertex_mesh_writer.WriteVtkUsingMeshWithCellId(*p_mesh);

        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_transit_type);
        CellsGenerator<FixedG1GenerationalCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_transit_type);
        VertexBasedCellPopulation<3> cell_population(*p_mesh, cells);

        OffLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory("TestUniaxialLoad/Test2");
        simulator.SetSamplingTimestepMultiple(50);
        const double end_time = 0.5;
        simulator.SetEndTime(end_time);

        MAKE_PTR(MisraForce<3>, p_force);
        MAKE_PTR(HorizontalStretchForce<3>, p_force2);
        p_force->SetLateralSurfaceEnergyParameter(1);
        p_force->SetApicalLineTensionParameter(0);
        p_force->SetBasalSurfaceEnergyParameter(0);
        p_force->SetTargetVolume(z_height*target_area*2);
        p_force->SetVolumeCompressibilityParameter(0);
        simulator.AddForce(p_force);
        simulator.AddForce(p_force2);

        /* Finally we call the {{{Solve}}} method on the simulation to run the simulation.*/
        simulator.Solve();

        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 21u);
        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), end_time, 1e-10);
    }
};

#endif /*TESTVERTEXMESH33HORIZONTALSTRETCH_HPP_*/
