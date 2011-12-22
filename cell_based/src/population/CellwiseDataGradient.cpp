/*

Copyright (C) University of Oxford, 2005-2011

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

#include "CellwiseDataGradient.hpp"
#include "LinearBasisFunction.hpp"

template<unsigned DIM>
c_vector<double, DIM>& CellwiseDataGradient<DIM>::rGetGradient(unsigned nodeIndex)
{
    return mGradients[nodeIndex];
}


template<unsigned DIM>
void CellwiseDataGradient<DIM>::SetupGradients()
{
    MeshBasedCellPopulation<DIM>* p_cell_population = static_cast<MeshBasedCellPopulation<DIM>*>(&(CellwiseData<DIM>::Instance()->rGetCellPopulation()));
    TetrahedralMesh<DIM,DIM>& r_mesh = p_cell_population->rGetMesh();

    // Initialise gradients size
    unsigned num_nodes = p_cell_population->GetNumNodes();
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

            // This code is commented because CelwiseData Can't deal with ghost nodes see #1975
            assert(p_cell_population->IsGhostNode(node_global_index) == false);
            //// Check whether ghost element
            //if (p_cell_population->IsGhostNode(node_global_index) == true)
            //{
            //    is_ghost_element = true;
            //    break;
            //}

            // If no ghost element, get PDE solution
            CellPtr p_cell = p_cell_population->GetCellUsingLocationIndex(node_global_index);
            double pde_solution = CellwiseData<DIM>::Instance()->GetValue(p_cell, 0);

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
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = p_cell_population->Begin();
         cell_iter != p_cell_population->End();
         ++cell_iter)
    {
        unsigned node_global_index = p_cell_population->GetLocationIndexUsingCell(*cell_iter);

        if (!num_real_elems_for_node[node_global_index] > 0)
        {
            NEVER_REACHED;
            // This code is commented because CellwiseData Can't deal with ghost nodes so won't ever come into this statement see #1975
            //// The node is a real node which is not in any real element
            //// but should be connected to some cells (if more than one cell in mesh)
            //Node<DIM>& this_node = *(p_cell_population->GetNodeCorrespondingToCell(*cell_iter));
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
            //        if (    !(p_cell_population->IsGhostNode(adjacent_node_global_index))
            //             && adjacent_node_global_index != node_global_index )
            //        {
            //
            //            // Calculate the contribution of gradient from this node
            //            Node<DIM>& adjacent_node = *(r_mesh.GetNode(adjacent_node_global_index));
            //
            //            double this_cell_concentration = CellwiseData<DIM>::Instance()->GetValue(*cell_iter, 0);
            //            CellPtr p_adjacent_cell = p_cell_population->GetCellUsingLocationIndex(adjacent_node_global_index);
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

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class CellwiseDataGradient<1>;
template class CellwiseDataGradient<2>;
template class CellwiseDataGradient<3>;
