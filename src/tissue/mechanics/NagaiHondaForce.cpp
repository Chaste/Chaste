/*

Copyright (C) University of Oxford, 2005-2010

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/
#include "NagaiHondaForce.hpp"


template<unsigned DIM>
NagaiHondaForce<DIM>::NagaiHondaForce(bool usingDifferentialAdhesion)
   : AbstractForce<DIM>(),
     mUsingDifferentialAdhesion(usingDifferentialAdhesion)
{
}


template<unsigned DIM>
NagaiHondaForce<DIM>::~NagaiHondaForce()
{
}


template<unsigned DIM>
void NagaiHondaForce<DIM>::AddForceContribution(std::vector<c_vector<double, DIM> >& rForces,
                                                AbstractTissue<DIM>& rTissue)
{
    // Helper instance of TissueConfig
    TissueConfig* p_params = TissueConfig::Instance();

    // Helper variable that is a static cast of the tissue
    VertexBasedTissue<DIM>* p_tissue = static_cast<VertexBasedTissue<DIM>*>(&rTissue);

    // Iterate over vertices in the tissue
    for (unsigned node_index=0; node_index<p_tissue->GetNumNodes(); node_index++)
    {
        // Compute the force on this node

        /*
         * The force on this Node is given by the gradient of the total free
         * energy of the Tissue, evaluated at the position of the vertex. This
         * free energy is the sum of the free energies of all TissueCellPtrs in
         * the tissue. The free energy of each TissueCellPtr is comprised of three
         * parts - a cell deformation energy, a membrane surface tension energy
         * and an adhesion energy.
         *
         * Note that since the movement of this Node only affects the free energy
         * of the three TissueCellPtrs containing it, we can just consider the
         * contributions to the free energy gradient from each of these three
         * TissueCellPtrs.
         */

        c_vector<double, DIM> deformation_contribution = zero_vector<double>(DIM);
        c_vector<double, DIM> membrane_surface_tension_contribution = zero_vector<double>(DIM);
        c_vector<double, DIM> adhesion_contribution = zero_vector<double>(DIM);

        // Find the indices of the elements owned by this node
        std::set<unsigned> containing_elem_indices = p_tissue->GetNode(node_index)->rGetContainingElementIndices();

        // Iterate over these elements
        for (std::set<unsigned>::iterator iter = containing_elem_indices.begin();
             iter != containing_elem_indices.end();
             ++iter)
        {
            // Get this element and its index
            VertexElement<DIM, DIM>* p_element = p_tissue->GetElement(*iter);
            unsigned element_index = p_element->GetIndex();

            // Find the local index of this node in this element
            unsigned local_index = p_element->GetNodeLocalIndex(node_index);

            /******** Start of deformation force calculation ********/

            // Compute the area of this element and its gradient at this node
            double element_area = p_tissue->rGetMesh().GetVolumeOfElement(*iter);
            c_vector<double, DIM> element_area_gradient = p_tissue->rGetMesh().GetAreaGradientOfElementAtNode(p_element, local_index);

            // Get the target area of the cell
            double cell_target_area = p_tissue->GetTargetAreaOfCell(p_tissue->GetCellUsingLocationIndex(element_index));

            // Add the force contribution from this cell's deformation energy (note the minus sign)
            deformation_contribution -= 2*p_params->GetNagaiHondaDeformationEnergyParameter()*(element_area - cell_target_area)*element_area_gradient;

            /******** End of deformation force calculation *************/

            /******** Start of membrane force calculation ***********/

            // Compute the perimeter of the element and its gradient at this node
            double element_perimeter = p_tissue->rGetMesh().GetSurfaceAreaOfElement(*iter);
            c_vector<double, DIM> element_perimeter_gradient = p_tissue->rGetMesh().GetPerimeterGradientOfElementAtNode(p_element, local_index);

            // Get the target perimeter of the cell
            double cell_target_perimeter = 2*sqrt(M_PI*cell_target_area);

            // Add the force contribution from this cell's membrane surface tension (note the minus sign)
            membrane_surface_tension_contribution -= 2*p_params->GetNagaiHondaMembraneSurfaceEnergyParameter()*(element_perimeter - cell_target_perimeter)*element_perimeter_gradient;

            /******** End of membrane force calculation **********/

            /******** Start of adhesion force calculation ***********/

            // Get the current, previous and next nodes in this element
            Node<DIM>* p_current_node = p_element->GetNode(local_index);

            unsigned previous_node_local_index = (p_element->GetNumNodes()+local_index-1)%(p_element->GetNumNodes());
            Node<DIM>* p_previous_node = p_element->GetNode(previous_node_local_index);

            unsigned next_node_local_index = (local_index+1)%(p_element->GetNumNodes());
            Node<DIM>* p_next_node = p_element->GetNode(next_node_local_index);


            // Compute the adhesion parameter for each of these edges
            double previous_edge_adhesion_parameter = GetAdhesionParameter(p_previous_node, p_current_node);
            double next_edge_adhesion_parameter = GetAdhesionParameter(p_current_node, p_next_node);

			// if using differneital adhesion then adjust adhesion parameters.
            if (mUsingDifferentialAdhesion)
            {
				// Determine combinationCellTypes (0: both cells sharing edge are wild types, 1: both cells sharing edge are labelled,
				// 2: one of the cells sharing the edge is wild type and the other is labelled, 99: otherwise)
				unsigned combinationCellTypesPreviousEdge = GetCombinationCellTypes(p_previous_node, p_current_node, rTissue);
				unsigned combinationCellTypesNextEdge = GetCombinationCellTypes(p_current_node, p_next_node, rTissue);

				// Compute the adhesion parameter for each of these edges
				previous_edge_adhesion_parameter = GetAdhesionParameter(p_previous_node, p_current_node, combinationCellTypesPreviousEdge);
				next_edge_adhesion_parameter = GetAdhesionParameter(p_current_node, p_next_node, combinationCellTypesNextEdge);
            }

            // Compute the gradient of the edge of the cell ending in this node
            c_vector<double, DIM> previous_edge_gradient = p_tissue->rGetMesh().GetPreviousEdgeGradientOfElementAtNode(p_element, local_index);

            // Compute the gradient of the edge of the cell starting in this node
            c_vector<double, DIM> next_edge_gradient = p_tissue->rGetMesh().GetNextEdgeGradientOfElementAtNode(p_element, local_index);

            // Add the force contribution from cell-cell and cell-boundary adhesion (note the minus sign)
            adhesion_contribution -= previous_edge_adhesion_parameter*previous_edge_gradient + next_edge_adhesion_parameter*next_edge_gradient;

            /******** End of adhesion force calculation *************/
        }

        c_vector<double, DIM> force_on_node = deformation_contribution +
                                              membrane_surface_tension_contribution +
                                              adhesion_contribution;

        rForces[node_index] += force_on_node;
    }
}


template<unsigned DIM>
double NagaiHondaForce<DIM>::GetAdhesionParameter(Node<DIM>* pNodeA, Node<DIM>* pNodeB)
{
    double adhesion_parameter;

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

    // If the edge corresponds to a single element, then the cell is on the boundary
    if (shared_elements.size() == 1)
    {
        adhesion_parameter = TissueConfig::Instance()->GetNagaiHondaCellBoundaryAdhesionEnergyParameter();
    }
    else
    {
        adhesion_parameter = TissueConfig::Instance()->GetNagaiHondaCellCellAdhesionEnergyParameter();
    }
    return adhesion_parameter;
}
template<unsigned DIM>
double NagaiHondaForce<DIM>::GetAdhesionParameter(Node<DIM>* pNodeA, Node<DIM>* pNodeB, unsigned combinationCellType)
{
       double adhesion_parameter;

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

        // If the edge corresponds to a single element, then the cell is on the boundary
        if (shared_elements.size() == 1)
        {
            adhesion_parameter = TissueConfig::Instance()->GetNagaiHondaCellBoundaryAdhesionEnergyParameter();
        }
        else
        {
            if (combinationCellType == 0) // if both cells are Wildtype (Mutation State)
            {
                adhesion_parameter = TissueConfig::Instance()->GetNagaiHondaCellCellAdhesionEnergyParameter(); // Value of 1.0 in TissueConfig, but 0.01 the Nagai & Honda paper
            }
            else if (combinationCellType == 1) // if both cells are Labelled (Mutation State)
            {
                adhesion_parameter = 2.0;
            }
            else if (combinationCellType == 2) // if one cell is Labelled and the other WildType (Mutation State)
            {
                adhesion_parameter = 1.5;
            }
            else // if else (if combinationCellType = 99)
            {
            	adhesion_parameter = TissueConfig::Instance()->GetNagaiHondaCellCellAdhesionEnergyParameter();
            }
        }
        return adhesion_parameter;
}

template<unsigned DIM>
unsigned NagaiHondaForce<DIM>::GetCombinationCellTypes(Node<DIM>* pNodeA, Node<DIM>* pNodeB,
        AbstractTissue<DIM>& rTissue)
{
       unsigned combinationCellType = 99;

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

        unsigned element_index1 = UNSIGNED_UNSET;
        unsigned element_index2 = UNSIGNED_UNSET;
        boost::shared_ptr<AbstractCellMutationState> mutationStateCell1;
        boost::shared_ptr<AbstractCellMutationState> mutationStateCell2;

        if (shared_elements.size() == 2)
        {

        	unsigned count = 1;
            for (typename std::set<unsigned>::iterator shrd_elem_iter = shared_elements.begin();
            	             shrd_elem_iter != shared_elements.end();
            	             ++shrd_elem_iter)
            {
                 if (count == 1)
                 {
                	 element_index1 = *shrd_elem_iter;
                 }
                 else if (count == 2)
                 {
                   	 element_index2 = *shrd_elem_iter;
                 }
                 count=count+1;
            }

           	TissueCellPtr p_cell1 = rTissue.GetCellUsingLocationIndex(element_index1);
            TissueCellPtr p_cell2 = rTissue.GetCellUsingLocationIndex(element_index2);

            if (p_cell1->GetMutationState()->IsType<LabelledCellMutationState>())
            {
                if (p_cell2->GetMutationState()->IsType<LabelledCellMutationState>())
                {
                    combinationCellType = 1;
                }
                else if (p_cell2->GetMutationState()->IsType<WildTypeCellMutationState>())
                {
                    combinationCellType = 2;
                }
            }
            else if (p_cell1->GetMutationState()->IsType<WildTypeCellMutationState>())
            {
                if (p_cell2->GetMutationState()->IsType<WildTypeCellMutationState>())
                {
                    combinationCellType = 0;
                }
                else if (p_cell2->GetMutationState()->IsType<LabelledCellMutationState>())
                {
                    combinationCellType = 2;
                }
            }
        }
        return combinationCellType;
}
/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class NagaiHondaForce<1>;
template class NagaiHondaForce<2>;
template class NagaiHondaForce<3>;


// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(NagaiHondaForce)
