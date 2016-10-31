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

//template<unsigned DIM>
//class HorizontalStretchForce : public AbstractForce<DIM>
//{
//private:
//
//    friend class boost::serialization::access;
//    template<class Archive>
//    void serialize(Archive & archive, const unsigned int version)
//    {
//        archive & boost::serialization::base_object<AbstractForce<DIM> >(*this);
//        archive & mForceMagnitude;
//        archive & mRelativeWidth;
//    }
//
//    /**
//     * Internal variable to save the magnitude
//     */
//    double mForceMagnitude;
//
//    /**
//     * The relative width from the boundary upon with the force is acting
//     */
//    double mRelativeWidth;
//
//public:
//    /**
//     * default constructor
//     */
//    HorizontalStretchForce(const double ForceMagnitude=1, const double RelativeWidth=0.1)
//    : mForceMagnitude(ForceMagnitude),
//      mRelativeWidth(RelativeWidth)
//    {
//MARK
//    }
//
//    virtual ~HorizontalStretchForce()
//    {
//    }
//
//    /**
//     * Adding force contributions to the Populations
//     * @param rCellPopulation
//     */
//    virtual void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
//    {
//        // Throw an exception message if not using a VertexBasedCellPopulation
//        ///\todo check whether this line influences profiling tests - if so, we should remove it.
//        if (dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation) == NULL)
//        {
//            EXCEPTION("MisraForce is to be used with a VertexBasedCellPopulation only");
//        }
//
//        // Define some helper variables
//        VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);
//        MutableVertexMesh<DIM, DIM>& rMesh = p_cell_population->rGetMesh();
//
//        const unsigned num_nodes = p_cell_population->GetNumNodes();
//        c_vector<double, DIM> force_in_positive_y = zero_vector<double>(DIM);
//        force_in_positive_y[0] = mForceMagnitude;
//        c_vector<double, DIM> force_in_negative_y = -force_in_positive_y;
//
//        ChasteCuboid<DIM> bounding_box_tmp = rMesh.CalculateBoundingBox();
//        const double min_x = bounding_box_tmp.rGetLowerCorner()[0];
//        const double max_x = bounding_box_tmp.rGetUpperCorner()[0];
//        assert(max_x > min_x);
//        const double width_x = max_x - min_x;
//        const double width = width_x*mRelativeWidth;
//        const double left_bound = min_x + width;
//        const double right_bound = max_x - width;
//        // could have built two more Cuboids to check if the point is inside,
//        // but will have a lot of wrapping which reduce efficiency.
//        // Iterate over nodes in the cell population
//        for (unsigned node_global_index=0; node_global_index<num_nodes; ++node_global_index)
//        {
//            Node<3>* p_this_node = p_cell_population->GetNode(node_global_index);
// //            if ( ! p_this_node->IsBoundaryNode())
// //                continue;
//            const double loc_x = p_this_node->rGetLocation()[0];
//            if (loc_x < left_bound)
//            {
//                p_this_node->AddAppliedForceContribution(force_in_negative_y);
//            }
//            else if (loc_x > right_bound)
//            {
//                p_this_node->AddAppliedForceContribution(force_in_positive_y);
//            }
//        }
//    }
//
//    /**
//     * For the compiler now
//     * @param rParamsFile
//     */
//    virtual void OutputForceParameters(out_stream& rParamsFile)
//    {
//MARK
//        *rParamsFile << "\t\t\t<ForceMagnitude>" << mForceMagnitude << "</ForceMagnitude>\n";
//        *rParamsFile << "\t\t\t<RelativeWidth>" << mRelativeWidth << "</RelativeWidth>\n";
//
//        // Call method on direct parent class
//        AbstractForce<DIM>::OutputForceParameters(rParamsFile);
//    }
//};
//#include "SerializationExportWrapper.hpp"
//EXPORT_TEMPLATE_CLASS1(HorizontalStretchForce,3)
//#include "SerializationExportWrapperForCpp.hpp"
//EXPORT_TEMPLATE_CLASS1(HorizontalStretchForce,3)

template<unsigned DIM>
class From2dGeneral3dVertexModelForce : public AbstractForce<DIM>
{
private:

    double mTargetApicalArea;
    double mApicalAreaParameter;
    double mApicalEdgeParameter;
    double mTargetBasalArea;
    double mBasalAreaParameter;
    double mBasalEdgeParameter;
    double mLateralEdgeParameter;
    double mTargetVolume;
    double mVolumeParameter;

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<DIM> >(*this);
        archive & mTargetApicalArea;
    }

public:
    From2dGeneral3dVertexModelForce()
        : AbstractForce<DIM>(),
          mTargetApicalArea(0),
          mApicalAreaParameter(0),
          mApicalEdgeParameter(0),
          mTargetBasalArea(0),
          mBasalAreaParameter(0),
          mBasalEdgeParameter(0),
          mLateralEdgeParameter(0),
          mTargetVolume(0),
          mVolumeParameter(0)
    {
    }

    void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
    {
        c_vector<double, DIM> force = zero_vector<double>(DIM);
        if (dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation) == NULL)
        {
            EXCEPTION("3dVertexModelForce is to be used with a VertexBasedCellPopulation only");
        }

        // Define some helper variables
        VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);
        MutableVertexMesh<DIM, DIM>& rMesh = p_cell_population->rGetMesh();
        const unsigned num_nodes = p_cell_population->GetNumNodes();
        const unsigned num_elements = p_cell_population->GetNumElements();

        // Begin by computing the volumes of each element in the mesh, to avoid having to do this multiple times
        std::vector<double> element_volumes(num_elements);
        std::vector<double> apical_areas(num_elements);
        std::vector<double> basal_areas(num_elements);
        for (unsigned elem_index = 0 ; elem_index<num_elements ; ++elem_index)
        {
            const VertexElement<DIM, DIM>* p_elem = p_cell_population->GetElement(elem_index);
            assert(elem_index == p_elem->GetIndex());
            element_volumes[elem_index] = rMesh.GetVolumeOfElement(elem_index);
            apical_areas[elem_index] = rMesh.CalculateAreaOfFace(p_elem->GetFace(1));
            basal_areas[elem_index] = rMesh.CalculateAreaOfFace(p_elem->GetFace(0));
        }

        // Iterate over nodes in the cell population
        for (unsigned node_index=0 ; node_index<num_nodes ; node_index++)
        {
            Node<DIM>* p_this_node = p_cell_population->GetNode(node_index);
            assert(node_index == p_this_node->GetIndex());
            // Get the type of node. 1=basal; 2=apical
            const unsigned node_type = unsigned(p_this_node->rGetNodeAttributes()[0]);

            c_vector<double, DIM> basal_face_contribution = zero_vector<double>(DIM);
            c_vector<double, DIM> ab_edge_contribution = zero_vector<double>(DIM);
            c_vector<double, DIM> apical_face_contribution = zero_vector<double>(DIM);
            c_vector<double, DIM> lateral_edge_contribution = zero_vector<double>(DIM);
            c_vector<double, DIM> volume_contribution = zero_vector<double>(DIM);

            const unsigned opposite_node_index = node_index + num_nodes/2*(node_type==1u?1:-1);
            const Node<DIM>* p_opposite_node = p_cell_population->GetNode(opposite_node_index);
            const c_vector<double, DIM> edge_gradient = (p_this_node->rGetLocation() - p_opposite_node->rGetLocation())/norm_2(p_this_node->rGetLocation() - p_opposite_node->rGetLocation());
            lateral_edge_contribution -= edge_gradient*mLateralEdgeParameter;

            // a variable to store such that the apical/basal edge forces are not counted twice for non-boundary edges.
            std::set<unsigned> neighbour_node_indices;

            // Find the indices of the elements owned by this node
            const std::set<unsigned> containing_elem_indices = p_this_node->rGetContainingElementIndices();

            // Iterate over these elements
            for (std::set<unsigned>::iterator iter = containing_elem_indices.begin();
                    iter != containing_elem_indices.end();
                    ++iter)
            {
                // Get this element, its index and its number of nodes
                VertexElement<DIM, DIM>* p_element = p_cell_population->GetElement(*iter);
                std::vector<VertexElement<DIM-1, DIM>*> lateral_faces;

                // Populate the pointers/vector to different types of face
                for (unsigned face_index=0; face_index<p_element->GetNumFaces(); ++face_index)
                {
                    VertexElement<DIM-1,DIM>* p_tmp_face = p_element->GetFace(face_index);
                    switch (unsigned(p_tmp_face->rGetElementAttributes()[0]))
                    {
                    case 1:
                        assert(face_index == 0);
                        break;
                    case 2:
                        assert(face_index == 1);
                        break;
                    case 3:
                        lateral_faces.push_back(p_tmp_face);
                        break;
                    default:
                        NEVER_REACHED;
                    }
                }

                const unsigned elem_index = p_element->GetIndex();

                // Calculating volume contribution
                c_vector<double, DIM> element_volume_gradient = rMesh.GetVolumeGradientofElementAtNode(p_element, node_index);
                // Add the force contribution from this cell's volume compressibility (note the minus sign)
                volume_contribution -= mVolumeParameter*element_volume_gradient*(element_volumes[elem_index] - mTargetVolume);

                // Pointer to the face having the same type as the node
                const VertexElement<DIM-1, DIM>* p_ab_face = p_element->GetFace(node_type - 1);
                const unsigned local_node_index_in_ab_face = p_ab_face->GetNodeLocalIndex(node_index);
                const c_vector<double, DIM> ab_face_gradient = rMesh.GetAreaGradientOfFaceAtNode(p_ab_face, local_node_index_in_ab_face);
                // Calculating apical face contribution
                if (node_type == 2)
                {
                    apical_face_contribution -= mApicalAreaParameter*ab_face_gradient*(apical_areas[elem_index] - mTargetApicalArea);
                }
                // Computing basal face contribution
                if (node_type == 1)
                {
                    basal_face_contribution -= mBasalAreaParameter*ab_face_gradient*(basal_areas[elem_index] - mTargetBasalArea);
                }
                const unsigned num_nodes_in_ab_face = p_ab_face->GetNumNodes();
                neighbour_node_indices.insert(p_ab_face->GetNodeGlobalIndex((local_node_index_in_ab_face+1)%num_nodes_in_ab_face));
                neighbour_node_indices.insert(p_ab_face->GetNodeGlobalIndex((local_node_index_in_ab_face-1+num_nodes_in_ab_face)%num_nodes_in_ab_face));
            }

            for(std::set<unsigned>::iterator it = neighbour_node_indices.begin(); it!=neighbour_node_indices.end() ; ++it)
            {
                Node<DIM>* p_neighbour_node = p_cell_population->GetNode(*it);
                const c_vector<double, DIM> edge_gradient = (p_this_node->rGetLocation() - p_neighbour_node->rGetLocation())/norm_2(p_this_node->rGetLocation() - p_neighbour_node->rGetLocation());
                ab_edge_contribution = - edge_gradient*(node_type==1u ? mBasalEdgeParameter : mApicalEdgeParameter);
            }

            c_vector<double, DIM> force_on_node = basal_face_contribution + ab_edge_contribution + apical_face_contribution
                                                  + lateral_edge_contribution + volume_contribution;
            p_this_node->AddAppliedForceContribution(force_on_node);
        }
    }

    void SetApicalParameter(const double TargetArea, const double AreaParameter, const double LineParameter)
    {
        mTargetApicalArea = TargetArea;
        mApicalAreaParameter = AreaParameter;
        mApicalEdgeParameter = LineParameter;
    }
    void SetBasalParameter(const double TargetArea, const double AreaParameter, const double LineParameter)
    {
        mTargetBasalArea = TargetArea;
        mBasalAreaParameter = AreaParameter;
        mBasalEdgeParameter = LineParameter;
    }
    void SetLateralParameter(const double Parameter)
    {
        mLateralEdgeParameter = Parameter;
    }
    void SetVolumeParameter(const double TargetVolume, const double VolumeParameter)
    {
        mTargetVolume = TargetVolume;
        mVolumeParameter = VolumeParameter;
    }

    void OutputForceParameters(out_stream& rParamsFile)
    {
        *rParamsFile << "\t\t\t<TargetApicalArea>" << mTargetApicalArea << "</TargetApicalArea>\n";
        *rParamsFile << "\t\t\t<TargetApicalArea>" << mTargetApicalArea << "</TargetApicalArea>\n";
        *rParamsFile << "\t\t\t<TargetBasalArea>" << mTargetBasalArea << "</TargetBasalArea>\n";
        *rParamsFile << "\t\t\t<TargetVolume>" << mTargetVolume << "</TargetVolume>\n";
        AbstractForce<3>::OutputForceParameters(rParamsFile);
    }
};
//template class From2dGeneral3dVertexModelForce<3>;
#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS1(From2dGeneral3dVertexModelForce,3)
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS1(From2dGeneral3dVertexModelForce,3)


class TestVertexMesh33UniaxialLoad : public AbstractCellBasedTestSuite
{
public:

    void TestVertexMesh33UiaxialLoad() throw (Exception)
    {
        // Make a mesh of 10x5
        const double z_height = 1;
        const double target_area = 1;
        VoronoiPrism3dVertexMeshGenerator generator(5, 5, z_height, 5, target_area);
        MutableVertexMesh<3,3>* p_mesh = generator.GetMeshAfterReMesh();
        VertexMeshWriter<3, 3> vertex_mesh_writer("TestUniaxialLoad/Test1", "InitialMesh", false);
        vertex_mesh_writer.WriteVtkUsingMeshWithCellId(*p_mesh);

        std::vector<CellPtr> cells;
        CellsGenerator<NoCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
        VertexBasedCellPopulation<3> cell_population(*p_mesh, cells);

        OffLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory("TestUniaxialLoad/Test1");
        simulator.SetSamplingTimestepMultiple(10);
        const double end_time = 0.5;
        simulator.SetEndTime(end_time);

//        MAKE_PTR(MisraForce<3>, p_force);
//        p_force->SetLateralSurfaceEnergyParameter(0);
//        p_force->SetApicalLineTensionParameter(0);
//        p_force->SetBasalSurfaceEnergyParameter(0);
//        p_force->SetTargetVolume(z_height*target_area*2);
//        p_force->SetVolumeCompressibilityParameter(2);
//        simulator.AddForce(p_force);
//        MAKE_PTR(HorizontalStretchForce<3>, p_force2);
//        simulator.AddForce(p_force2);
        MAKE_PTR(From2dGeneral3dVertexModelForce<3>, p_force3);
        p_force3->SetApicalParameter(0, 0, 1);
        simulator.AddForce(p_force3);

        simulator.Solve();

        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 25u);
        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), end_time, 1e-10);
    }

    // Testing and therefore highly un
    void TestApical1() throw (Exception)
    {
        const std::string output_file = "TestUniaxialLoad/Apical1";
        const double z_height = 0.5, target_area = 1;
        VoronoiPrism3dVertexMeshGenerator generator(5, 5, z_height, 5, target_area);
        MutableVertexMesh<3,3>* p_mesh = generator.GetMeshAfterReMesh();
        VertexMeshWriter<3, 3> vertex_mesh_writer(output_file, "InitialMesh", false);
        vertex_mesh_writer.WriteVtkUsingMeshWithCellId(*p_mesh);
        std::vector<CellPtr> cells; CellsGenerator<NoCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
        VertexBasedCellPopulation<3> cell_population(*p_mesh, cells);
        OffLatticeSimulation<3> simulator(cell_population); simulator.SetOutputDirectory(output_file);
        simulator.SetSamplingTimestepMultiple(10); const double end_time = 0.5; simulator.SetEndTime(end_time);
        MAKE_PTR(From2dGeneral3dVertexModelForce<3>, p_force3);

        p_force3->SetApicalParameter(0, 0, 1);
        simulator.AddForce(p_force3); simulator.Solve();
    }
    void TestApical2() throw (Exception)
    {
        const std::string output_file = "TestUniaxialLoad/Apical2";
        const double z_height = 0.5, target_area = 1;
        VoronoiPrism3dVertexMeshGenerator generator(5, 5, z_height, 5, target_area);
        MutableVertexMesh<3,3>* p_mesh = generator.GetMeshAfterReMesh();
        VertexMeshWriter<3, 3> vertex_mesh_writer(output_file, "InitialMesh", false);
        vertex_mesh_writer.WriteVtkUsingMeshWithCellId(*p_mesh);
        std::vector<CellPtr> cells; CellsGenerator<NoCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
        VertexBasedCellPopulation<3> cell_population(*p_mesh, cells);
        OffLatticeSimulation<3> simulator(cell_population); simulator.SetOutputDirectory(output_file);
        simulator.SetSamplingTimestepMultiple(10); const double end_time = 0.5; simulator.SetEndTime(end_time);
        MAKE_PTR(From2dGeneral3dVertexModelForce<3>, p_force3);

        p_force3->SetApicalParameter(1, 1, 0);
        simulator.AddForce(p_force3); simulator.Solve();
    }
    void TestApical3() throw (Exception)
    {
        const std::string output_file = "TestUniaxialLoad/Apical3";
        const double z_height = 0.5, target_area = 1;
        VoronoiPrism3dVertexMeshGenerator generator(5, 5, z_height, 5, target_area);
        MutableVertexMesh<3,3>* p_mesh = generator.GetMeshAfterReMesh();
        VertexMeshWriter<3, 3> vertex_mesh_writer(output_file, "InitialMesh", false);
        vertex_mesh_writer.WriteVtkUsingMeshWithCellId(*p_mesh);
        std::vector<CellPtr> cells; CellsGenerator<NoCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
        VertexBasedCellPopulation<3> cell_population(*p_mesh, cells);
        OffLatticeSimulation<3> simulator(cell_population); simulator.SetOutputDirectory(output_file);
        simulator.SetSamplingTimestepMultiple(10); const double end_time = 0.5; simulator.SetEndTime(end_time);
        MAKE_PTR(From2dGeneral3dVertexModelForce<3>, p_force3);

        p_force3->SetApicalParameter(1, 1, 1);
        simulator.AddForce(p_force3); simulator.Solve();
    }

    void TestBasal1() throw (Exception)
    {
        const std::string output_file = "TestUniaxialLoad/Basal1";
        const double z_height = 0.5, target_area = 1;
        VoronoiPrism3dVertexMeshGenerator generator(5, 5, z_height, 5, target_area);
        MutableVertexMesh<3,3>* p_mesh = generator.GetMeshAfterReMesh();
        VertexMeshWriter<3, 3> vertex_mesh_writer(output_file, "InitialMesh", false);
        vertex_mesh_writer.WriteVtkUsingMeshWithCellId(*p_mesh);
        std::vector<CellPtr> cells; CellsGenerator<NoCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
        VertexBasedCellPopulation<3> cell_population(*p_mesh, cells);
        OffLatticeSimulation<3> simulator(cell_population); simulator.SetOutputDirectory(output_file);
        simulator.SetSamplingTimestepMultiple(10); const double end_time = 0.5; simulator.SetEndTime(end_time);
        MAKE_PTR(From2dGeneral3dVertexModelForce<3>, p_force3);

        p_force3->SetBasalParameter(0, 0, 1);
        simulator.AddForce(p_force3); simulator.Solve();
    }
    void TestBasal2() throw (Exception)
    {
        const std::string output_file = "TestUniaxialLoad/Basal2";
        const double z_height = 0.5, target_area = 1;
        VoronoiPrism3dVertexMeshGenerator generator(5, 5, z_height, 5, target_area);
        MutableVertexMesh<3,3>* p_mesh = generator.GetMeshAfterReMesh();
        VertexMeshWriter<3, 3> vertex_mesh_writer(output_file, "InitialMesh", false);
        vertex_mesh_writer.WriteVtkUsingMeshWithCellId(*p_mesh);
        std::vector<CellPtr> cells; CellsGenerator<NoCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
        VertexBasedCellPopulation<3> cell_population(*p_mesh, cells);
        OffLatticeSimulation<3> simulator(cell_population); simulator.SetOutputDirectory(output_file);
        simulator.SetSamplingTimestepMultiple(10); const double end_time = 0.5; simulator.SetEndTime(end_time);
        MAKE_PTR(From2dGeneral3dVertexModelForce<3>, p_force3);

        p_force3->SetBasalParameter(1, 1, 0);
        simulator.AddForce(p_force3); simulator.Solve();
    }
    void TestBasal3() throw (Exception)
    {
        const std::string output_file = "TestUniaxialLoad/Basal3";
        const double z_height = 0.5, target_area = 1;
        VoronoiPrism3dVertexMeshGenerator generator(5, 5, z_height, 5, target_area);
        MutableVertexMesh<3,3>* p_mesh = generator.GetMeshAfterReMesh();
        VertexMeshWriter<3, 3> vertex_mesh_writer(output_file, "InitialMesh", false);
        vertex_mesh_writer.WriteVtkUsingMeshWithCellId(*p_mesh);
        std::vector<CellPtr> cells; CellsGenerator<NoCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
        VertexBasedCellPopulation<3> cell_population(*p_mesh, cells);
        OffLatticeSimulation<3> simulator(cell_population); simulator.SetOutputDirectory(output_file);
        simulator.SetSamplingTimestepMultiple(10); const double end_time = 0.5; simulator.SetEndTime(end_time);
        MAKE_PTR(From2dGeneral3dVertexModelForce<3>, p_force3);

        p_force3->SetBasalParameter(1, 1, 1);
        simulator.AddForce(p_force3); simulator.Solve();
    }

    void TestApicalAndBasal1() throw (Exception)
    {
        const std::string output_file = "TestUniaxialLoad/ApicalAndBasal1";
        const double z_height = 0.5, target_area = 1;
        VoronoiPrism3dVertexMeshGenerator generator(5, 5, z_height, 5, target_area);
        MutableVertexMesh<3,3>* p_mesh = generator.GetMeshAfterReMesh();
        VertexMeshWriter<3, 3> vertex_mesh_writer(output_file, "InitialMesh", false);
        vertex_mesh_writer.WriteVtkUsingMeshWithCellId(*p_mesh);
        std::vector<CellPtr> cells; CellsGenerator<NoCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
        VertexBasedCellPopulation<3> cell_population(*p_mesh, cells);
        OffLatticeSimulation<3> simulator(cell_population); simulator.SetOutputDirectory(output_file);
        simulator.SetSamplingTimestepMultiple(10); const double end_time = 0.5; simulator.SetEndTime(end_time);
        MAKE_PTR(From2dGeneral3dVertexModelForce<3>, p_force3);

        p_force3->SetApicalParameter(1, 1, 1);
        p_force3->SetBasalParameter(1, 1, 1);
        simulator.AddForce(p_force3); simulator.Solve();
    }
    void TestApicalAndBasal2() throw (Exception)
    {
        const std::string output_file = "TestUniaxialLoad/ApicalAndBasal2";
        const double z_height = 0.5, target_area = 1;
        VoronoiPrism3dVertexMeshGenerator generator(5, 5, z_height, 5, target_area);
        MutableVertexMesh<3,3>* p_mesh = generator.GetMeshAfterReMesh();
        VertexMeshWriter<3, 3> vertex_mesh_writer(output_file, "InitialMesh", false);
        vertex_mesh_writer.WriteVtkUsingMeshWithCellId(*p_mesh);
        std::vector<CellPtr> cells; CellsGenerator<NoCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
        VertexBasedCellPopulation<3> cell_population(*p_mesh, cells);
        OffLatticeSimulation<3> simulator(cell_population); simulator.SetOutputDirectory(output_file);
        simulator.SetSamplingTimestepMultiple(10); const double end_time = 0.5; simulator.SetEndTime(end_time);
        MAKE_PTR(From2dGeneral3dVertexModelForce<3>, p_force3);

        p_force3->SetApicalParameter(1, 1, 1);
        p_force3->SetBasalParameter(2, 2, 2);
        simulator.AddForce(p_force3); simulator.Solve();
    }

    void TestLateral() throw (Exception)
    {
        const std::string output_file = "TestUniaxialLoad/Lateral";
        const double z_height = 0.5, target_area = 1;
        VoronoiPrism3dVertexMeshGenerator generator(5, 5, z_height, 5, target_area);
        MutableVertexMesh<3,3>* p_mesh = generator.GetMeshAfterReMesh();
        VertexMeshWriter<3, 3> vertex_mesh_writer(output_file, "InitialMesh", false);
        vertex_mesh_writer.WriteVtkUsingMeshWithCellId(*p_mesh);
        std::vector<CellPtr> cells; CellsGenerator<NoCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
        VertexBasedCellPopulation<3> cell_population(*p_mesh, cells);
        OffLatticeSimulation<3> simulator(cell_population); simulator.SetOutputDirectory(output_file);
        simulator.SetSamplingTimestepMultiple(10); const double end_time = 0.5; simulator.SetEndTime(end_time);
        MAKE_PTR(From2dGeneral3dVertexModelForce<3>, p_force3);

        p_force3->SetLateralParameter(1);
        simulator.AddForce(p_force3); simulator.Solve();
    }

    void TestVolume1() throw (Exception)
    {
        const std::string output_file = "TestUniaxialLoad/Volume1";
        const double z_height = 0.5, target_area = 1;
        VoronoiPrism3dVertexMeshGenerator generator(5, 5, z_height, 5, target_area);
        MutableVertexMesh<3,3>* p_mesh = generator.GetMeshAfterReMesh();
        VertexMeshWriter<3, 3> vertex_mesh_writer(output_file, "InitialMesh", false);
        vertex_mesh_writer.WriteVtkUsingMeshWithCellId(*p_mesh);
        std::vector<CellPtr> cells; CellsGenerator<NoCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
        VertexBasedCellPopulation<3> cell_population(*p_mesh, cells);
        OffLatticeSimulation<3> simulator(cell_population); simulator.SetOutputDirectory(output_file);
        simulator.SetSamplingTimestepMultiple(10); const double end_time = 0.5; simulator.SetEndTime(end_time);
        MAKE_PTR(From2dGeneral3dVertexModelForce<3>, p_force3);

        p_force3->SetVolumeParameter(0, 1);
        simulator.AddForce(p_force3); simulator.Solve();
    }
    void TestVolume2() throw (Exception)
    {
        const std::string output_file = "TestUniaxialLoad/Volume2";
        const double z_height = 0.5, target_area = 1;
        VoronoiPrism3dVertexMeshGenerator generator(5, 5, z_height, 5, target_area);
        MutableVertexMesh<3,3>* p_mesh = generator.GetMeshAfterReMesh();
        VertexMeshWriter<3, 3> vertex_mesh_writer(output_file, "InitialMesh", false);
        vertex_mesh_writer.WriteVtkUsingMeshWithCellId(*p_mesh);
        std::vector<CellPtr> cells; CellsGenerator<NoCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
        VertexBasedCellPopulation<3> cell_population(*p_mesh, cells);
        OffLatticeSimulation<3> simulator(cell_population); simulator.SetOutputDirectory(output_file);
        simulator.SetSamplingTimestepMultiple(10); const double end_time = 0.5; simulator.SetEndTime(end_time);
        MAKE_PTR(From2dGeneral3dVertexModelForce<3>, p_force3);

        p_force3->SetVolumeParameter(z_height*target_area, 1);
        simulator.AddForce(p_force3); simulator.Solve();
    }
    void TestVolume3() throw (Exception)
    {
        const std::string output_file = "TestUniaxialLoad/Volume3";
        const double z_height = 0.5, target_area = 1;
        VoronoiPrism3dVertexMeshGenerator generator(5, 5, z_height, 5, target_area);
        MutableVertexMesh<3,3>* p_mesh = generator.GetMeshAfterReMesh();
        VertexMeshWriter<3, 3> vertex_mesh_writer(output_file, "InitialMesh", false);
        vertex_mesh_writer.WriteVtkUsingMeshWithCellId(*p_mesh);
        std::vector<CellPtr> cells; CellsGenerator<NoCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
        VertexBasedCellPopulation<3> cell_population(*p_mesh, cells);
        OffLatticeSimulation<3> simulator(cell_population); simulator.SetOutputDirectory(output_file);
        simulator.SetSamplingTimestepMultiple(10); const double end_time = 0.5; simulator.SetEndTime(end_time);
        MAKE_PTR(From2dGeneral3dVertexModelForce<3>, p_force3);

        p_force3->SetVolumeParameter(z_height*target_area*3, 1);
        simulator.AddForce(p_force3); simulator.Solve();
    }
};

#endif /*TESTVERTEXMESH33UNIAXIALLOAD_HPP_*/
