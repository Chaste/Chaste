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
#include "MonolayerVertexMeshCustomFunctions.hpp"

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

c_vector<double, 3> CalculateEdgeGradient(const Node<3>* pNode1, const Node<3>* pNode2)
{
    const c_vector<double, 3>& loc1 = pNode1->rGetLocation();
    const c_vector<double, 3>& loc2 = pNode2->rGetLocation();
    return (loc1 - loc2) / norm_2(loc1 - loc2);
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
    for (unsigned elem_index = 0; elem_index < num_elements; ++elem_index)
    {
        const VertexElement<3, 3>* p_elem = p_cell_population->GetElement(elem_index);
        assert(elem_index == p_elem->GetIndex());
        element_volumes[elem_index] = rMesh.GetVolumeOfElement(elem_index);
        apical_areas[elem_index] = rMesh.CalculateAreaOfFace(GetApicalFace(p_elem));
        basal_areas[elem_index] = rMesh.CalculateAreaOfFace(GetBasalFace(p_elem));
    }

    // Iterate over nodes in the cell population
    for (unsigned node_index = 0; node_index < num_nodes; ++node_index)
    {
        Node<3>* p_this_node = p_cell_population->GetNode(node_index);

        assert(node_index == p_this_node->GetIndex());

        // Get the type of node. 1=basal; 2=apical
        const Monolayer::v_type node_type = GetNodeType(p_this_node);
        c_vector<double, 3> volume_contribution = zero_vector<double>(3);

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
            if (fabs(mVolumeParameter) > 1e-5)
            {
                c_vector<double, 3> element_volume_gradient = rMesh.GetVolumeGradientofElementAtNode(p_element, node_index);
                // Add the force contribution from this cell's volume compressibility (note the minus sign)
                volume_contribution -= mVolumeParameter * element_volume_gradient * (element_volumes[elem_index] - mTargetVolume);
            }
        }
        
        p_this_node->AddAppliedForceContribution(volume_contribution);
    }

    // Do lateral face contributions and edge contributions.
    for (unsigned face_id = 0; face_id < rMesh.GetNumFaces(); ++face_id)
    {
        const VertexElement<2, 3>* p_face = rMesh.GetFace(face_id);

        switch (GetFaceType(p_face))
        {
        case Monolayer::ApicalValue:
        {
            // Calculate apical area contribution
            if (fabs(mApicalAreaParameter) > 1e-5)
            {
                const double apical_area = rMesh.CalculateAreaOfFace(p_face);

                for (unsigned node_id = 0; node_id < p_face->GetNumNodes(); ++node_id)
                {
                    c_vector<double, 3> result = rMesh.GetAreaGradientOfFaceAtNode(p_face, node_id);
                    result *= -1 * mApicalAreaParameter * (apical_area - mTargetApicalArea);
                    p_face->GetNode(node_id)->AddAppliedForceContribution(result);
                }
            }
            break;
        }
        case Monolayer::BasalValue:
        {
            // Calculate basal area contribution
            if (fabs(mBasalAreaParameter) > 1e-5)
            {
                const double basal_area = rMesh.CalculateAreaOfFace(p_face);

                for (unsigned node_id = 0; node_id < p_face->GetNumNodes(); ++node_id)
                {
                    c_vector<double, 3> result = rMesh.GetAreaGradientOfFaceAtNode(p_face, node_id);
                    result *= -1 * mBasalAreaParameter * (basal_area - mTargetBasalArea);
                    p_face->GetNode(node_id)->AddAppliedForceContribution(result);
                }
            }
            break;
        }
        case Monolayer::LateralValue:
        {
            Node<3>* p_node1 = p_face->GetNode(p_face->GetNumNodes() - 1);
            // Calculate apical and basal edge contribution here so that it will be counted once.
            for (unsigned i=0; i<p_face->GetNumNodes(); ++i)
            {
                Node<3>* p_node2 = p_face->GetNode(i);
                c_vector<double, 3> result = CalculateEdgeGradient(p_node1, p_node2);
                
                if (GetNodeType(p_node1) == GetNodeType(p_node2))
                {
                    if (IsApicalNode(p_node1))
                    {
                        result *= -1 * mApicalEdgeParameter;
                    }
                    else
                    {
                        assert(IsBasalNode(p_node1));
                        result *= -1 * mBasalEdgeParameter;
                    }
                    
                }
                else
                {
                    result *= -1 * mLateralEdgeParameter;
                }
                
                p_node1->AddAppliedForceContribution(result);
                p_node2->AddAppliedForceContribution(-result);

                p_node1 = p_node2;
            }

            break;
        }
        default:
            NEVER_REACHED;
        }
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
