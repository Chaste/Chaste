/*

Copyright (C) University of Oxford, 2005-2012

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
   : AbstractForce<DIM>(),
     mWelikyOsterAreaParameter(1.0),
     mWelikyOsterPerimeterParameter(1.0)
{
}

template<unsigned DIM>
WelikyOsterForce<DIM>::~WelikyOsterForce()
{
}

template<unsigned DIM>
void WelikyOsterForce<DIM>::AddForceContribution(std::vector<c_vector<double, DIM> >& rForces,
                                                 AbstractCellPopulation<DIM>& rCellPopulation)
{
    // Make sure that we are in the correct dimension - this code will be eliminated at compile time
    assert(DIM == 2); // this method only works in 2D at present

    // Throw an exception message if not using a VertexBasedCellPopulation
    if (dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation) == NULL)
    {
        EXCEPTION("WelikyOsterForce is to be used with a VertexBasedCellPopulation only");
    }

    // Helper variable that is a static cast of the cell population
    VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);

    /*
     * The force on each node is given by the interaction between the area and
     * the perimeter of the element containing the node.
     */

    // Iterate over elements in the cell population
    for (typename VertexMesh<DIM,DIM>::VertexElementIterator element_iter = p_cell_population->rGetMesh().GetElementIteratorBegin();
             element_iter != p_cell_population->rGetMesh().GetElementIteratorEnd();
             ++element_iter)
    {
        unsigned element_index = element_iter->GetIndex();

        /******** Start of deformation force calculation ********/

        // Compute the area of this element
        double element_area = p_cell_population->rGetMesh().GetVolumeOfElement(element_index);

        double deformation_coefficient = GetWelikyOsterAreaParameter()/element_area;

        /******** End of deformation force calculation *************/

        /******** Start of membrane force calculation ***********/

        // Compute the perimeter of the element
        double element_perimeter = p_cell_population->rGetMesh().GetSurfaceAreaOfElement(element_index);

        double membrane_surface_tension_coefficient = GetWelikyOsterPerimeterParameter()*element_perimeter;

        /******** End of membrane force calculation **********/

        unsigned num_nodes = element_iter->GetNumNodes();
        for (unsigned node_local_index = 0; node_local_index < num_nodes; node_local_index++)
        {
            unsigned node_global_index = element_iter->GetNodeGlobalIndex(node_local_index);

            c_vector<double, DIM> current_node = element_iter->GetNodeLocation(node_local_index);
            c_vector<double, DIM> next_node = element_iter->GetNodeLocation((node_local_index + 1)%(element_iter->GetNumNodes()));
            c_vector<double, DIM> previous_node = element_iter->GetNodeLocation((node_local_index + element_iter->GetNumNodes() - 1)%(element_iter->GetNumNodes()));

            c_vector<double, DIM> clockwise_unit_vector = p_cell_population->rGetMesh().GetVectorFromAtoB(current_node, previous_node);
            clockwise_unit_vector /= norm_2(clockwise_unit_vector);

            c_vector<double, DIM> anti_clockwise_unit_vector = p_cell_population->rGetMesh().GetVectorFromAtoB(current_node, next_node);
            anti_clockwise_unit_vector /= norm_2(anti_clockwise_unit_vector);

            // Calculate the outward normal at the node
            c_vector<double, DIM> outward_normal = -0.5*clockwise_unit_vector - 0.5*anti_clockwise_unit_vector;
            outward_normal /= norm_2(outward_normal);

            c_vector<double, DIM> deformation_contribution = deformation_coefficient * outward_normal;

            c_vector<double, DIM> membrane_surface_tension_contribution = membrane_surface_tension_coefficient * (clockwise_unit_vector + anti_clockwise_unit_vector);

            c_vector<double, DIM> force_on_node = deformation_contribution + membrane_surface_tension_contribution;

            rForces[node_global_index] += force_on_node;
        }
    }
}

template<unsigned DIM>
double WelikyOsterForce<DIM>::GetWelikyOsterAreaParameter()
{
    return mWelikyOsterAreaParameter;
}

template<unsigned DIM>
double WelikyOsterForce<DIM>::GetWelikyOsterPerimeterParameter()
{
    return mWelikyOsterPerimeterParameter;
}

template<unsigned DIM>
void WelikyOsterForce<DIM>::SetWelikyOsterAreaParameter(double welikyOsterAreaParameter)
{
    mWelikyOsterAreaParameter = welikyOsterAreaParameter;
}

template<unsigned DIM>
void WelikyOsterForce<DIM>::SetWelikyOsterPerimeterParameter(double welikyOsterPerimeterParameter)
{
    mWelikyOsterPerimeterParameter = welikyOsterPerimeterParameter;
}

template<unsigned DIM>
void WelikyOsterForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<WelikyOsterAreaParameter>" << mWelikyOsterAreaParameter << "</WelikyOsterAreaParameter>\n";
    *rParamsFile << "\t\t\t<WelikyOsterPerimeterParameter>" << mWelikyOsterPerimeterParameter << "</WelikyOsterPerimeterParameter>\n";

    // Call method on direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
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
