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
#include "WelikyOsterForce.hpp"

template<unsigned DIM>
WelikyOsterForce<DIM>::WelikyOsterForce()
   : AbstractForce<DIM>()
{
}


template<unsigned DIM>
WelikyOsterForce<DIM>::~WelikyOsterForce()
{
}


template<unsigned DIM>
void WelikyOsterForce<DIM>::AddForceContribution(std::vector<c_vector<double, DIM> >& rForces,
                                                 AbstractTissue<DIM>& rTissue)
{
    // Make sure that we are in the correct dimension - this code will be eliminated at compile time
    assert(DIM == 2); // this method only works in 2D at present

    // Helper instance of TissueConfig
    TissueConfig* p_params = TissueConfig::Instance();

    // Helper variable that is a static cast of the tissue
    VertexBasedTissue<DIM>* p_tissue = static_cast<VertexBasedTissue<DIM>*>(&rTissue);

    /*
     * The force on each node is given by the interaction between the area and
     * the perimeter of the element containing the node.
     */

    // Iterate over elements in the tissue
    for (typename VertexMesh<DIM,DIM>::VertexElementIterator element_iter = p_tissue->rGetMesh().GetElementIteratorBegin();
             element_iter != p_tissue->rGetMesh().GetElementIteratorEnd();
             ++element_iter)
    {
        unsigned element_index = element_iter->GetIndex();

        /******** Start of deformation force calculation ********/

        // Compute the area of this element
        double element_area = p_tissue->rGetMesh().GetAreaOfElement(element_index);

        double deformation_coefficient = p_params->GetWelikyOsterAreaParameter()/element_area;


        /******** End of deformation force calculation *************/

        /******** Start of membrane force calculation ***********/

        // Compute the perimeter of the element
        double element_perimeter = p_tissue->rGetMesh().GetPerimeterOfElement(element_index);

        double membrane_surface_tension_coefficient = p_params->GetWelikyOsterPerimeterParameter()*element_perimeter;

        /******** End of membrane force calculation **********/

        unsigned num_nodes = element_iter->GetNumNodes();
           for (unsigned node_local_index = 0; node_local_index < num_nodes; node_local_index++)
        {
            unsigned node_global_index = element_iter->GetNodeGlobalIndex(node_local_index);

            c_vector<double, DIM> current_node = element_iter->GetNodeLocation(node_local_index);
            c_vector<double, DIM> next_node = element_iter->GetNodeLocation((node_local_index + 1)%(element_iter->GetNumNodes()));
            c_vector<double, DIM> previous_node = element_iter->GetNodeLocation((node_local_index + element_iter->GetNumNodes() - 1)%(element_iter->GetNumNodes()));

            c_vector<double, DIM> clockwise_unit_vector = p_tissue->rGetMesh().GetVectorFromAtoB(current_node, previous_node);
            clockwise_unit_vector /= norm_2(clockwise_unit_vector);
            c_vector<double, DIM> anti_clockwise_unit_vector = p_tissue->rGetMesh().GetVectorFromAtoB(current_node, next_node);
            anti_clockwise_unit_vector /= norm_2(anti_clockwise_unit_vector);

            // Calculate the outward normal at the node
            c_vector<double, DIM> outward_normal = -0.5*clockwise_unit_vector - 0.5*anti_clockwise_unit_vector;
            outward_normal /= norm_2(outward_normal);


            c_vector<double, DIM> deformation_contribution = deformation_coefficient * outward_normal;

            c_vector<double, DIM> membrane_surface_tension_contribution = membrane_surface_tension_coefficient * (clockwise_unit_vector + anti_clockwise_unit_vector);

            //c_vector<double, DIM> adhesion_contribution = zero_vector<double>(DIM);

            c_vector<double, DIM> force_on_node = deformation_contribution +
                                                  membrane_surface_tension_contribution;
                                                  // + adhesion_contribution;

            rForces[node_global_index] += force_on_node;
        }
    }
}


/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class WelikyOsterForce<1>;
template class WelikyOsterForce<2>;
template class WelikyOsterForce<3>;


// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(WelikyOsterForce)
