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

#define OUTPUT_NAME "TestApicalConstrictionExample/InitialMesh"
#define ADDFORCEPARAMETER p_force3->SetVolumeParameter(0.8, 10);
#define CURRENT_TEST std::string("Volume2")
#define END_TIME 0.5

#include "HoneycombVertexMeshGenerator.hpp"
#include "MeshBuilderHelper.hpp"
#include "CellLabel.hpp"
#include "CellLabelWriter.hpp"

template<unsigned DIM>
class PatternedApicalConstrictionForce : public GeneralMonolayerVertexMeshForce<DIM>
{
private:

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<GeneralMonolayerVertexMeshForce<DIM> >(*this);
    }

    double mPatternedTargetApicalArea;
    double mPatternedApicalareaParameter;
    double mPatternedApicalEdgeParameter;

public:
    /**
     * Default constructor.
     */
    PatternedApicalConstrictionForce()
        : GeneralMonolayerVertexMeshForce<DIM>(),
          mPatternedTargetApicalArea(0.0),
          mPatternedApicalareaParameter(0.0),
          mPatternedApicalEdgeParameter(0.0)
    {
    }

    virtual ~PatternedApicalConstrictionForce()
    {
    }

    void SetPatternedApicalParameter(const double patternedTargetApicalArea, const double patternedApicalareaParameter, const double patternedApicalEdgeParameter)
    {
        mPatternedTargetApicalArea = patternedTargetApicalArea;
        mPatternedApicalareaParameter = patternedApicalareaParameter;
        mPatternedApicalEdgeParameter = patternedApicalEdgeParameter;
    }

    virtual void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
    {
        c_vector<double, DIM> force = zero_vector<double>(DIM);
        if (dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation) == NULL)
        {
            EXCEPTION("PatternedApicalConstrictionForce is to be used with a VertexBasedCellPopulation only");
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

            // A variable to store such that the apical/basal edge forces are not counted twice for non-boundary edges.
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
                volume_contribution -= this->mVolumeParameter*element_volume_gradient*(element_volumes[elem_index] - this->mTargetVolume);
                // Pointer to the face having the same type as the node
                const VertexElement<DIM-1, DIM>* p_ab_face = p_element->GetFace(node_type - 1);
                const unsigned local_node_index_in_ab_face = p_ab_face->GetNodeLocalIndex(node_index);
                const c_vector<double, DIM> ab_face_gradient = rMesh.GetAreaGradientOfFaceAtNode(p_ab_face, local_node_index_in_ab_face);
                // Calculating apical face contribution
                if (node_type == 2)
                {
                    bool cell_is_labelled = p_cell_population->GetCellUsingLocationIndex(elem_index)->template HasCellProperty<CellLabel>();
                    double apical_target_area = cell_is_labelled ? mPatternedTargetApicalArea : this->mTargetApicalArea;
                    double apical_area_parameter = cell_is_labelled ? mPatternedApicalareaParameter : this->mApicalareaParameter;

                    apical_face_contribution -= apical_area_parameter*ab_face_gradient*(apical_areas[elem_index] - apical_target_area);
                }
                // Computing basal face contribution
                if (node_type == 1)
                {
                    basal_face_contribution -= this->mBasalareaParameter*ab_face_gradient*(basal_areas[elem_index] - this->mTargetBasalArea);
                }
                const unsigned num_nodes_in_ab_face = p_ab_face->GetNumNodes();
                neighbour_node_indices.insert(p_ab_face->GetNodeGlobalIndex((local_node_index_in_ab_face+1)%num_nodes_in_ab_face));
                neighbour_node_indices.insert(p_ab_face->GetNodeGlobalIndex((local_node_index_in_ab_face-1+num_nodes_in_ab_face)%num_nodes_in_ab_face));
            }

            for (std::set<unsigned>::iterator it = neighbour_node_indices.begin();
                 it != neighbour_node_indices.end();
                 ++it)
            {
                Node<DIM>* p_neighbour_node = p_cell_population->GetNode(*it);
                const c_vector<double, DIM> edge_gradient = (p_this_node->rGetLocation() - p_neighbour_node->rGetLocation())/norm_2(p_this_node->rGetLocation() - p_neighbour_node->rGetLocation());
                ab_edge_contribution -= edge_gradient*(node_type==1u ? this->mBasalEdgeParameter : this->mApicalEdgeParameter);
            }

            const unsigned opposite_node_index = node_index + num_nodes/2*(node_type==1u?1:-1);
            const Node<DIM>* p_opposite_node = p_cell_population->GetNode(opposite_node_index);
            const c_vector<double, DIM> edge_gradient = (p_this_node->rGetLocation() - p_opposite_node->rGetLocation())/norm_2(p_this_node->rGetLocation() - p_opposite_node->rGetLocation());
            lateral_edge_contribution -= edge_gradient*this->mLateralEdgeParameter*(containing_elem_indices.size());

            c_vector<double, DIM> force_on_node = basal_face_contribution + ab_edge_contribution + apical_face_contribution
                    + lateral_edge_contribution + volume_contribution;
            p_this_node->AddAppliedForceContribution(force_on_node);
        }
    }

    /**
     * For the compiler now
     * @param rParamsFile
     */
    virtual void OutputForceParameters(out_stream& rParamsFile)
    {
        *rParamsFile << "\t\t\t<PatternedTargetApicalArea>" << mPatternedTargetApicalArea << "</PatternedTargetApicalArea>\n";
        *rParamsFile << "\t\t\t<PatternedApicalareaParameter>" << mPatternedApicalareaParameter << "</PatternedApicalareaParameter>\n";
        *rParamsFile << "\t\t\t<PatternedApicalEdgeParameter>" << mPatternedApicalEdgeParameter << "</PatternedApicalEdgeParameter>\n";

        // Call method on direct parent class
        GeneralMonolayerVertexMeshForce<DIM>::OutputForceParameters(rParamsFile);
    }
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS1(PatternedApicalConstrictionForce,3)
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS1(PatternedApicalConstrictionForce,3)

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
        MeshBuilderHelper builder("ApicalConstriction");
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

        MAKE_PTR(PatternedApicalConstrictionForce<3>, p_force3);
        p_force3->SetApicalParameter(10, 10, 1);
        p_force3->SetBasalParameter(10, 10, 1);
        p_force3->SetLateralParameter(4);
        p_force3->SetVolumeParameter(200, 1);
        p_force3->SetPatternedApicalParameter(8, 10, 1);
        simulator.AddForce(p_force3);

        simulator.Solve();

        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 25u);
        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), end_time, 1e-10);
    }
};

#endif /*TESTAPICALCONSTRICTIONEXAMPLE_HPP_*/
