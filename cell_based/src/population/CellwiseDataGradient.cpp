/*

Copyright (c) 2005-2019, University of Oxford.
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

#include "CellwiseDataGradient.hpp"
#include "LinearBasisFunction.hpp"

template<unsigned DIM>
c_vector<double, DIM>& CellwiseDataGradient<DIM>::rGetGradient(unsigned nodeIndex)
{
    return mGradients[nodeIndex];
}

template<unsigned DIM>
void CellwiseDataGradient<DIM>::SetupGradients(AbstractCellPopulation<DIM>& rCellPopulation, const std::string& rItemName)
{
    MeshBasedCellPopulation<DIM>* pCellPopulation = static_cast<MeshBasedCellPopulation<DIM>*>(&(rCellPopulation));
    TetrahedralMesh<DIM,DIM>& r_mesh = pCellPopulation->rGetMesh();

    // Initialise gradients size
    unsigned num_nodes = pCellPopulation->GetNumNodes();
    mGradients.resize(num_nodes, zero_vector<double>(DIM));

    // The constant gradients at each element
    std::vector<c_vector<double, DIM> > gradients_on_elements;
    unsigned num_elements = r_mesh.GetNumElements();
    gradients_on_elements.resize(num_elements, zero_vector<double>(DIM));

    // The number of elements containing a given node (excl ghost elements)
    std::vector<unsigned> num_real_elems_for_node(num_nodes, 0);

    for (unsigned elem_index=0; elem_index<num_elements; elem_index++)
    {
        Element<DIM,DIM>& r_elem = *(r_mesh.GetElement(elem_index));

        // Calculate the basis functions at any point (eg zero) in the element
        c_matrix<double, DIM, DIM> jacobian, inverse_jacobian;
        double jacobian_det;
        r_mesh.GetInverseJacobianForElement(elem_index, jacobian, jacobian_det, inverse_jacobian);
        const ChastePoint<DIM> zero_point;
        c_matrix<double, DIM, DIM+1> grad_phi;
        LinearBasisFunction<DIM>::ComputeTransformedBasisFunctionDerivatives(zero_point, inverse_jacobian, grad_phi);

        bool is_ghost_element = false;

        for (unsigned node_index=0; node_index<DIM+1; node_index++)
        {
            unsigned node_global_index = r_elem.GetNodeGlobalIndex(node_index);

            // This code is commented because CellData can't deal with ghost nodes see #1975
            assert(pCellPopulation->IsGhostNode(node_global_index) == false);
            //// Check whether ghost element
            //if (pCellPopulation->IsGhostNode(node_global_index) == true)
            //{
            //    is_ghost_element = true;
            //    break;
            //}

            // If no ghost element, get PDE solution
            CellPtr p_cell = pCellPopulation->GetCellUsingLocationIndex(node_global_index);
            double pde_solution = p_cell->GetCellData()->GetItem(rItemName);

            // Interpolate gradient
            for (unsigned i=0; i<DIM; i++)
            {
                gradients_on_elements[elem_index](i) += pde_solution* grad_phi(i, node_index);
            }
        }

        // Add gradient at element to gradient at node
        if (!is_ghost_element)
        {
            for (unsigned node_index=0; node_index<DIM+1; node_index++)
            {
                unsigned node_global_index = r_elem.GetNodeGlobalIndex(node_index);
                mGradients[node_global_index] += gradients_on_elements[elem_index];
                num_real_elems_for_node[node_global_index]++;
            }
        }
    }

    // Divide to obtain average gradient
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = pCellPopulation->Begin();
         cell_iter != pCellPopulation->End();
         ++cell_iter)
    {
        unsigned node_global_index = pCellPopulation->GetLocationIndexUsingCell(*cell_iter);

        if (!(num_real_elems_for_node[node_global_index] > 0))
        {
            NEVER_REACHED;
            // This code is commented because CellwiseData Can't deal with ghost nodes so won't ever come into this statement see #1975
            //// The node is a real node which is not in any real element
            //// but should be connected to some cells (if more than one cell in mesh)
            //Node<DIM>& this_node = *(pCellPopulation->GetNodeCorrespondingToCell(*cell_iter));
            //
            //mGradients[node_global_index] = zero_vector<double>(DIM);
            //unsigned num_real_adjacent_nodes = 0;
            //
            //// Get all the adjacent nodes which correspond to real cells
            //std::set<Node<DIM>*> real_adjacent_nodes;
            //real_adjacent_nodes.clear();
            //
            //// First loop over containing elements
            //for (typename Node<DIM>::ContainingElementIterator element_iter = this_node.ContainingElementsBegin();
            //     element_iter != this_node.ContainingElementsEnd();
            //     ++element_iter)
            //{
            //    // Then loop over nodes therein
            //    Element<DIM,DIM>& r_adjacent_elem = *(r_mesh.GetElement(*element_iter));
            //    for (unsigned local_node_index=0; local_node_index<DIM+1; local_node_index++)
            //    {
            //        unsigned adjacent_node_global_index = r_adjacent_elem.GetNodeGlobalIndex(local_node_index);
            //
            //        // If not a ghost node and not the node we started with
            //        if (!(pCellPopulation->IsGhostNode(adjacent_node_global_index))
            //             && adjacent_node_global_index != node_global_index )
            //        {
            //
            //            // Calculate the contribution of gradient from this node
            //            Node<DIM>& adjacent_node = *(r_mesh.GetNode(adjacent_node_global_index));
            //
            //            double this_cell_concentration = CellwiseData<DIM>::Instance()->GetValue(*cell_iter, 0);
            //            CellPtr p_adjacent_cell = pCellPopulation->GetCellUsingLocationIndex(adjacent_node_global_index);
            //            double adjacent_cell_concentration = CellwiseData<DIM>::Instance()->GetValue(p_adjacent_cell, 0);
            //
            //            c_vector<double, DIM> gradient_contribution = zero_vector<double>(DIM);
            //
            //            if (fabs(this_cell_concentration-adjacent_cell_concentration) > 100*DBL_EPSILON)
            //            {
            //                c_vector<double, DIM> edge_vector = r_mesh.GetVectorFromAtoB(this_node.rGetLocation(), adjacent_node.rGetLocation());
            //                double norm_edge_vector = norm_2(edge_vector);
            //                gradient_contribution = edge_vector
            //                                            * (adjacent_cell_concentration - this_cell_concentration)
            //                                            / (norm_edge_vector * norm_edge_vector);
            //            }
            //
            //            mGradients[node_global_index] += gradient_contribution;
            //            num_real_adjacent_nodes++;
            //        }
            //    }
            //}
            //mGradients[node_global_index] /= num_real_adjacent_nodes;
        }
        else
        {
            mGradients[node_global_index] /= num_real_elems_for_node[node_global_index];
        }
    }
}

// Explicit instantiation
template class CellwiseDataGradient<1>;
template class CellwiseDataGradient<2>;
template class CellwiseDataGradient<3>;
