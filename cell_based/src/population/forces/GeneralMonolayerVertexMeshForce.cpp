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

#include "GeneralMonolayerVertexMeshForce.hpp"
#include "VertexBasedCellPopulation.hpp"

GeneralMonolayerVertexMeshForce::GeneralMonolayerVertexMeshForce()
    : AbstractForce<3>(),
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

GeneralMonolayerVertexMeshForce::~GeneralMonolayerVertexMeshForce()
{
}

void GeneralMonolayerVertexMeshForce::AddForceContribution(AbstractCellPopulation<3>& rCellPopulation)
{
    c_vector<double, 3> force = zero_vector<double>(3);
    if (dynamic_cast<VertexBasedCellPopulation<3>*>(&rCellPopulation) == NULL)
    {
        EXCEPTION("GeneralMonolayerVertexMeshForce is to be used with a VertexBasedCellPopulation only");
    }

    // Define some helper variables
    VertexBasedCellPopulation<3>* p_cell_population = static_cast<VertexBasedCellPopulation<3>*>(&rCellPopulation);
    MutableVertexMesh<3, 3>& rMesh = p_cell_population->rGetMesh();
    const unsigned num_nodes = p_cell_population->GetNumNodes();
    const unsigned num_elements = p_cell_population->GetNumElements();

    // Begin by computing the volumes of each element in the mesh, to avoid having to do this multiple times
    std::vector<double> element_volumes(num_elements);
    std::vector<double> apical_areas(num_elements);
    std::vector<double> basal_areas(num_elements);
    for (unsigned elem_index = 0; elem_index<num_elements; ++elem_index)
    {
        const VertexElement<3, 3>* p_elem = p_cell_population->GetElement(elem_index);
        assert(elem_index == p_elem->GetIndex());
        element_volumes[elem_index] = rMesh.GetVolumeOfElement(elem_index);
        apical_areas[elem_index] = rMesh.CalculateAreaOfFace(p_elem->GetFace(1));
        basal_areas[elem_index] = rMesh.CalculateAreaOfFace(p_elem->GetFace(0));
    }

    // Iterate over nodes in the cell population
    for (unsigned node_index = 0; node_index < num_nodes; node_index++)
    {
        Node<3>* p_this_node = p_cell_population->GetNode(node_index);

        assert(node_index == p_this_node->GetIndex());

        // Get the type of node. 1=basal; 2=apical
        const unsigned node_type = unsigned(p_this_node->rGetNodeAttributes()[0]);
        c_vector<double, 3> basal_face_contribution = zero_vector<double>(3);
        c_vector<double, 3> ab_edge_contribution = zero_vector<double>(3);
        c_vector<double, 3> apical_face_contribution = zero_vector<double>(3);
        c_vector<double, 3> lateral_edge_contribution = zero_vector<double>(3);
        c_vector<double, 3> volume_contribution = zero_vector<double>(3);

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
            const VertexElement<3, 3>* p_element = p_cell_population->GetElement(*iter);
            const unsigned elem_index = p_element->GetIndex();

            // Calculating volume contribution
            if (fabs(mVolumeParameter)>1e-5)
            {
                c_vector<double, 3> element_volume_gradient = rMesh.GetVolumeGradientofElementAtNode(p_element, node_index);
                // Add the force contribution from this cell's volume compressibility (note the minus sign)
                volume_contribution -= mVolumeParameter*element_volume_gradient*(element_volumes[elem_index] - mTargetVolume);
            }

            // Pointer to the face having the same type as the node
            const VertexElement<2, 3>* p_ab_face = p_element->GetFace(node_type - 1);
            const unsigned local_node_index_in_ab_face = p_ab_face->GetNodeLocalIndex(node_index);
            const c_vector<double, 3> ab_face_gradient = rMesh.GetAreaGradientOfFaceAtNode(p_ab_face, local_node_index_in_ab_face);

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

        const c_vector<double, 3>& r_node_location = p_this_node->rGetLocation();

        for (std::set<unsigned>::iterator it = neighbour_node_indices.begin();
             it != neighbour_node_indices.end();
             ++it)
        {
            Node<3>* p_neighbour_node = p_cell_population->GetNode(*it);
            const c_vector<double, 3> edge_gradient = (r_node_location - p_neighbour_node->rGetLocation())/norm_2(p_this_node->rGetLocation() - p_neighbour_node->rGetLocation());
            ab_edge_contribution -= edge_gradient*(node_type==1u ? mBasalEdgeParameter : mApicalEdgeParameter);
        }

        const unsigned opposite_node_index = node_index + num_nodes/2*(node_type==1u?1:-1);
        const c_vector<double, 3>& r_opposite_node_location = p_cell_population->GetNode(opposite_node_index)->rGetLocation();
        const c_vector<double, 3> edge_gradient = (r_node_location- r_opposite_node_location)/norm_2(r_node_location - r_opposite_node_location);
        lateral_edge_contribution -= edge_gradient*mLateralEdgeParameter*(containing_elem_indices.size());

        c_vector<double, 3> force_on_node = basal_face_contribution + ab_edge_contribution + apical_face_contribution
                + lateral_edge_contribution + volume_contribution;
        p_this_node->AddAppliedForceContribution(force_on_node);
    }
}

void GeneralMonolayerVertexMeshForce::SetApicalParameters(const double lineParameter, const double areaParameter, const double targetArea)
{
    mTargetApicalArea = targetArea;
    mApicalAreaParameter = areaParameter;
    mApicalEdgeParameter = lineParameter;
}

void GeneralMonolayerVertexMeshForce::SetBasalParameters(const double lineParameter, const double areaParameter, const double targetArea)
{
    mTargetBasalArea = targetArea;
    mBasalAreaParameter = areaParameter;
    mBasalEdgeParameter = lineParameter;
}

void GeneralMonolayerVertexMeshForce::SetLateralParameter(const double parameter)
{
    mLateralEdgeParameter = parameter;
}

void GeneralMonolayerVertexMeshForce::SetVolumeParameters(const double volumeParameter, const double targetVolume)
{
    mTargetVolume = targetVolume;
    mVolumeParameter = volumeParameter;
}

void GeneralMonolayerVertexMeshForce::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<TargetApicalArea>" << mTargetApicalArea << "</TargetApicalArea>\n";
    *rParamsFile << "\t\t\t<ApicalAreaParameter>" << mApicalAreaParameter << "</ApicalAreaParameter>\n";
    *rParamsFile << "\t\t\t<ApicalEdgeParameter>" << mApicalEdgeParameter << "</ApicalEdgeParameter>\n";

    *rParamsFile << "\t\t\t<TargetBasalArea>" << mTargetBasalArea << "</TargetBasalArea>\n";
    *rParamsFile << "\t\t\t<BasalAreaParameter>" << mBasalAreaParameter << "</BasalAreaParameter>\n";
    *rParamsFile << "\t\t\t<BasalEdgeParameter>" << mBasalEdgeParameter << "</BasalEdgeParameter>\n";

    *rParamsFile << "\t\t\t<LateralEdgeParameter>" << mLateralEdgeParameter << "</LateralEdgeParameter>\n";
    *rParamsFile << "\t\t\t<TargetVolume>" << mTargetVolume << "</TargetVolume>\n";
    *rParamsFile << "\t\t\t<VolumeParameter>" << mVolumeParameter << "</VolumeParameter>\n";

    AbstractForce<3>::OutputForceParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(GeneralMonolayerVertexMeshForce)