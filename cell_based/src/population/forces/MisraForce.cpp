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

#include "MisraForce.hpp"

MisraForce::MisraForce()
   : AbstractForce<3>(),
     mApicalLineTensionParameter(1.0), // These parameters are average values in Misra's paper
     mLateralSurfaceEnergyParameter(2.0),
     mBasalSurfaceEnergyParameter(0.98),
     mVolumeCompressibilityParameter(100),
     mTargetVolume(1.0)
{
}

MisraForce::~MisraForce()
{
}

void MisraForce::AddForceContribution(AbstractCellPopulation<3>& rCellPopulation)
{
    // Throw an exception message if not using a VertexBasedCellPopulation
    ///\todo check whether this line influences profiling tests - if so, we should remove it.
    if (dynamic_cast<VertexBasedCellPopulation<3>*>(&rCellPopulation) == NULL)
    {
        EXCEPTION("MisraForce is to be used with a VertexBasedCellPopulation only");
    }

    // Define some helper variables
    VertexBasedCellPopulation<3>* p_cell_population = static_cast<VertexBasedCellPopulation<3>*>(&rCellPopulation);
    MutableVertexMesh<3,3>& rMesh = p_cell_population->rGetMesh();
    const unsigned num_nodes = p_cell_population->GetNumNodes();
    const unsigned num_elements = p_cell_population->GetNumElements();

    // Begin by computing the volumes of each element in the mesh, to avoid having to do this multiple times
    std::vector<double> element_volumes(num_elements);
    for (unsigned elem_index = 0; elem_index<num_elements; ++elem_index)
    {
        assert(elem_index == p_cell_population->GetElement(elem_index)->GetIndex());
        element_volumes[elem_index] = rMesh.GetVolumeOfElement(elem_index);
    }

    // Iterate over nodes in the cell population
    for (unsigned node_global_index=0; node_global_index<num_nodes; node_global_index++)
    {
        Node<3>* p_this_node = p_cell_population->GetNode(node_global_index);
        assert(node_global_index == p_this_node->GetIndex());

        // Get the type of node. 1=basal; 2=apical
        const unsigned node_type = unsigned(p_this_node->rGetNodeAttributes()[0]);

        /*
         * The force on this Node is given by the gradient of the total free
         * energy of the CellPopulation, evaluated at the position of the vertex. This
         * free energy is the sum of the free energies of all CellPtrs in
         * the cell population. The free energy of each CellPtr is comprised of four
         * terms - ///\todo bla bla bla bla.
         *
         * Note that since the movement of this Node only affects the free energy
         * of the CellPtrs containing it, we can just consider the contributions
         * to the free energy gradient from each of these CellPtrs.
         */
        c_vector<double, 3> volume_elasticity_contribution = zero_vector<double>(3);
        c_vector<double, 3> basal_face_contribution = zero_vector<double>(3);
        c_vector<double, 3> apical_line_tension_contribution = zero_vector<double>(3);
        c_vector<double, 3> lateral_face_contribution = zero_vector<double>(3);

        // A helper variable so that "non-boundary" lateral face will not be counted twice
        std::set<unsigned> counted_lateral_face_indices;

        // Find the indices of the elements owned by this node
        const std::set<unsigned> containing_elem_indices = p_this_node->rGetContainingElementIndices();

        // Iterate over these elements
        for (std::set<unsigned>::iterator iter = containing_elem_indices.begin();
             iter != containing_elem_indices.end();
             ++iter)
        {
            // Get this element, its index and its number of nodes
            VertexElement<3,3>* p_element = p_cell_population->GetElement(*iter);
            std::vector<VertexElement<2, 3>*> lateral_faces;
            VertexElement<2,3>* p_basal_face = NULL;
            VertexElement<2,3>* p_apical_face = NULL;
            bool apical_face_orientation = false;

            // Populate the pointers/vector to different types of face
            for (unsigned face_index=0; face_index<p_element->GetNumFaces(); ++face_index)
            {
                VertexElement<2,3>* p_tmp_face = p_element->GetFace(face_index);
                switch (unsigned(p_tmp_face->rGetElementAttributes()[0]))
                {
                    case 1:
                        p_basal_face = p_tmp_face;
                        break;
                    case 2:
                        p_apical_face = p_tmp_face;
                        apical_face_orientation = p_element->FaceIsOrientatedAntiClockwise(face_index);
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
            c_vector<double, 3> element_volume_gradient = rMesh.GetVolumeGradientofElementAtNode(p_element, node_global_index);

            // Add the force contribution from this cell's volume compressibility (note the minus sign)
            volume_elasticity_contribution -= GetVolumeCompressibilityParameter()*(element_volumes[elem_index] - mTargetVolume)*element_volume_gradient;

            // Calculating apical line tension contribution
            ///\todo Maybe a refactoring similar to lateral will eliminate the trouble of getting previous and next
            if (node_type == 2) ///\todo maybe with macro #define BASAL 1; #define APICAL 2; #define LATERAL 3; to increase code's readability?
            {
                // As we create apical faces on top, CCW from the top is CW from the centre of the element
                assert(apical_face_orientation == false);
                const unsigned num_apical_nodes = p_apical_face->GetNumNodes();
                const unsigned local_apical_node_index = p_apical_face->GetNodeLocalIndex(node_global_index);

                // Get the previous and next nodes in this element
                const unsigned previous_apical_node_index = (num_apical_nodes+local_apical_node_index-1)%num_apical_nodes;
                Node<3>* p_previous_apical_node = p_apical_face->GetNode(previous_apical_node_index);

                const unsigned next_apical_node_index = (local_apical_node_index+1)%num_apical_nodes;
                Node<3>* p_next_apical_node = p_apical_face->GetNode(next_apical_node_index);

                // Compute the apical line tension parameter for each of these edges - be aware that this is half of the actual
                // value for internal edges since we are looping over each of the internal edges twice
                const double previous_edge_line_tension_parameter = GetApicalLineTensionParameter(p_previous_apical_node, p_this_node, *p_cell_population);
                const double next_edge_line_tension_parameter = GetApicalLineTensionParameter(p_this_node, p_next_apical_node, *p_cell_population);

                // Compute the gradient of each these edges, computed at the present node
                const c_vector<double, 3> previous_edge_gradient = rMesh.GetPreviousEdgeGradientOfElementAtNode(p_apical_face, local_apical_node_index, apical_face_orientation);
                const c_vector<double, 3> next_edge_gradient = rMesh.GetNextEdgeGradientOfElementAtNode(p_apical_face, local_apical_node_index, apical_face_orientation);

                // Add the force contribution from cell-cell and cell-boundary line tension (note the minus sign)
                apical_line_tension_contribution -= previous_edge_line_tension_parameter*previous_edge_gradient
                                                    + next_edge_line_tension_parameter*next_edge_gradient;
            }

            // Computing basal face contribution
            if (node_type == 1)
            {
                const unsigned localBasalNodeIndex = p_basal_face->GetNodeLocalIndex(node_global_index);
                const c_vector<double, 3> basal_face_gradient = rMesh.GetAreaGradientOfFaceAtNode(p_basal_face, localBasalNodeIndex);

                basal_face_contribution -= GetBasalSurfaceEnergyParameter() * basal_face_gradient;
            }

            // Compute lateral surface contribution
            for (unsigned local_lateral_face_index = 0; local_lateral_face_index<lateral_faces.size(); ++local_lateral_face_index)
            {
                VertexElement<2,3>* p_lateral_face = lateral_faces[local_lateral_face_index];
                bool lateral_face_has_node = false;

                // Check if this lateral face contains this node
                for (unsigned local_lateral_node_index=0; local_lateral_node_index<p_lateral_face->GetNumNodes(); ++local_lateral_node_index)
                {
                    if (node_global_index == p_lateral_face->GetNodeGlobalIndex(local_lateral_node_index))
                    {
                        lateral_face_has_node = true;
                        break;
                    }
                }

                ///\todo remove this one when the following assertion never have problem.
                assert( lateral_face_has_node == (p_lateral_face->GetNodeLocalIndex(node_global_index)!=UINT_MAX));

                // Proceed to force contribution if this lateral face has that particular node
                if (lateral_face_has_node)
                {
                    // Some temporary variables so the code would be more comprehensible
                    const unsigned lateralFaceGlobalIndex = p_lateral_face->GetIndex();
                    const bool isLateralFaceCounted = counted_lateral_face_indices.find(lateralFaceGlobalIndex)!=counted_lateral_face_indices.end();

                    if (!isLateralFaceCounted)
                    {
                        const unsigned local_lateral_node_index = p_lateral_face->GetNodeLocalIndex(node_global_index);

                        const c_vector<double, 3> lateral_face_gradient = rMesh.GetAreaGradientOfFaceAtNode(p_lateral_face, local_lateral_node_index);
                        lateral_face_contribution -= GetLateralSurfaceEnergyParameter() * lateral_face_gradient;

                        // Update counted_lateral_face_indices and an assertion
                        // Return type of insert(T) is a pair, and second will indicate if the insertion is successful
                        assert(counted_lateral_face_indices.insert(lateralFaceGlobalIndex).second);
                    }
                }
            }
        }

        c_vector<double, 3> force_on_node = volume_elasticity_contribution + lateral_face_contribution
                                              + apical_line_tension_contribution + basal_face_contribution;
        p_this_node->AddAppliedForceContribution(force_on_node);
    }
}

double MisraForce::GetApicalLineTensionParameter(Node<3>* pNodeA, Node<3>* pNodeB, VertexBasedCellPopulation<3>& rVertexCellPopulation) const
{
    // Find the indices of the elements owned by each node
    std::set<unsigned> elements_containing_nodeA = pNodeA->rGetContainingElementIndices();
    std::set<unsigned> elements_containing_nodeB = pNodeB->rGetContainingElementIndices();

    // Find common elements
    std::set<unsigned> shared_elements;
    std::set_intersection(elements_containing_nodeA.begin(),
                          elements_containing_nodeA.end(),
                          elements_containing_nodeB.begin(),
                          elements_containing_nodeB.end(),
                          std::inserter(shared_elements, shared_elements.begin()));

    // Check that the nodes have a common edge
    assert(!shared_elements.empty());

    // Since each internal edge is visited twice in the loop above, we have to use half the line tension parameter
    // for each visit.
    double line_tension_parameter_in_calculation = GetApicalLineTensionParameter()/2.0;

    // If the edge corresponds to a single element, then the cell is on the boundary
    if (shared_elements.size() == 1)
    {
        line_tension_parameter_in_calculation = GetBoundaryLineTensionParameter();
    }

    return line_tension_parameter_in_calculation;
}

double MisraForce::GetVolumeCompressibilityParameter() const
{
    return mVolumeCompressibilityParameter;
}

double MisraForce::GetLateralSurfaceEnergyParameter() const
{
    return mLateralSurfaceEnergyParameter;
}

double MisraForce::GetApicalLineTensionParameter() const
{
    return mApicalLineTensionParameter;
}

double MisraForce::GetBoundaryLineTensionParameter() const
{
    return mApicalLineTensionParameter;
}

double MisraForce::GetBasalSurfaceEnergyParameter() const
{
    return mBasalSurfaceEnergyParameter;
}

void MisraForce::SetVolumeCompressibilityParameter(double volumeCompressibilityParameter)
{
    mVolumeCompressibilityParameter = volumeCompressibilityParameter;
}

void MisraForce::SetLateralSurfaceEnergyParameter(double LateralSurfaceEnergyParameter)
{
    mLateralSurfaceEnergyParameter = LateralSurfaceEnergyParameter;
}

void MisraForce::SetBasalSurfaceEnergyParameter(double BasalSurfaceEnergyParameter)
{
    mBasalSurfaceEnergyParameter = BasalSurfaceEnergyParameter;
}

void MisraForce::SetApicalLineTensionParameter(double apicalLineTensionParameter)
{
    mApicalLineTensionParameter = apicalLineTensionParameter;
}

void MisraForce::SetTargetVolume(double targetVolume)
{
    mTargetVolume = targetVolume;
}

void MisraForce::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<ApicalLineTensionParameter>" << mApicalLineTensionParameter << "</ApicalLineTensionParameter>\n";
    *rParamsFile << "\t\t\t<LateralSurfaceEnergyParameter>" << mLateralSurfaceEnergyParameter << "</LateralSurfaceEnergyParameter>\n";
    *rParamsFile << "\t\t\t<BasalSurfaceEnergyParameter>" << mBasalSurfaceEnergyParameter << "</BasalSurfaceEnergyParameter>\n";
    *rParamsFile << "\t\t\t<VolumeCompressibilityParameter>" << mVolumeCompressibilityParameter << "</VolumeCompressibilityParameter>\n";

    // Call method on direct parent class
    AbstractForce<3>::OutputForceParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(MisraForce)
