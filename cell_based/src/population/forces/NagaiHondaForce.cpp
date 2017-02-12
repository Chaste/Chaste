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

#include "NagaiHondaForce.hpp"

template<unsigned DIM>
NagaiHondaForce<DIM>::NagaiHondaForce()
    : AbstractForce<DIM>(),
      mNagaiHondaDeformationEnergyParameter(100.0), // This is 1.0 in the Nagai & Honda paper.
      mNagaiHondaMembraneSurfaceEnergyParameter(10.0), // This is 0.1 in the Nagai & Honda paper.
      mNagaiHondaCellCellAdhesionEnergyParameter(0.5), // This corresponds to a value of 1.0 for
                                                      // the sigma parameter in the Nagai & Honda
                                                      // paper. In the paper, the sigma value is
                                                      // set to 0.01.
      mNagaiHondaCellBoundaryAdhesionEnergyParameter(1.0) // This is 0.01 in the Nagai & Honda paper.
{
}

template<unsigned DIM>
NagaiHondaForce<DIM>::~NagaiHondaForce()
{
}

template<unsigned DIM>
void NagaiHondaForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{
	assert(dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation) != NULL);

    // Define some helper variables
    VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);
    VertexMesh<DIM, DIM>& r_mesh = p_cell_population->rGetMesh();
    unsigned num_nodes = r_mesh.GetNumNodes();
    unsigned num_elements = r_mesh.GetNumElements();

    // Begin by computing the area and perimeter of each element in the mesh, to avoid having to do this multiple times
    std::vector<double> element_areas(num_elements);
    std::vector<double> element_perimeters(num_elements);
    std::vector<double> target_areas(num_elements);
    for (typename VertexMesh<DIM,DIM>::VertexElementIterator elem_iter = r_mesh.GetElementIteratorBegin();
         elem_iter != r_mesh.GetElementIteratorEnd();
         ++elem_iter)
    {
        unsigned elem_index = elem_iter->GetIndex();
        element_areas[elem_index] = r_mesh.GetVolumeOfElement(elem_index);
        element_perimeters[elem_index] = r_mesh.GetSurfaceAreaOfElement(elem_index);

        // Get the target area of the cell associated with this element, throwing an exception if this is not stored as a CellData item
        try
        {
            target_areas[elem_index] = p_cell_population->GetCellUsingLocationIndex(elem_index)->GetCellData()->GetItem("target area");
        }
        catch (Exception&)
        {
            EXCEPTION("You need to add an AbstractTargetAreaModifier to the simulation in order to use NagaiHondaForce");
        }
    }

    // Iterate over vertices in the cell population
    for (unsigned node_index=0; node_index<num_nodes; node_index++)
    {
        c_vector<double, DIM> deformation_contribution = zero_vector<double>(DIM);

        // Find the indices of the elements owned by this node
        std::set<unsigned> containing_elem_indices = r_mesh.GetNode(node_index)->rGetContainingElementIndices();

        // Iterate over these elements
        for (std::set<unsigned>::iterator iter = containing_elem_indices.begin();
             iter != containing_elem_indices.end();
             ++iter)
        {
            // Get this element, its index and its number of nodes
            VertexElement<DIM, DIM>* p_element = r_mesh.GetElement(*iter);
            unsigned elem_index = p_element->GetIndex();

            // Find the local index of this node in this element
            unsigned local_index = p_element->GetNodeLocalIndex(node_index);

            // Add the force contribution from this cell's deformation energy (note the minus sign)
            c_vector<double, DIM> element_area_gradient = r_mesh.GetAreaGradientOfElementAtNode(p_element, local_index);
            deformation_contribution -= 2*GetNagaiHondaDeformationEnergyParameter()*(element_areas[elem_index] - target_areas[elem_index])*element_area_gradient;
        }

        r_mesh.GetNode(node_index)->AddAppliedForceContribution(deformation_contribution);
    }

    // Iterate over all edges and add adhesion and perimeter contractility force contributions
    for (typename VertexMesh<DIM, DIM>::EdgeIterator edge_iter = r_mesh.EdgesBegin();
         edge_iter != r_mesh.EdgesEnd();
         ++edge_iter)
    {
    	// Compute the adhesion parameter for this edge
        double adhesion_parameter = GetAdhesionParameter(edge_iter.GetNodeA(), edge_iter.GetNodeB(), *p_cell_population);

        unsigned node_A_index = edge_iter.GetNodeA()->GetIndex();
        unsigned node_B_index = edge_iter.GetNodeB()->GetIndex();

        const c_vector<double, DIM>& r_node_A_location = edge_iter.GetNodeA()->rGetLocation();
        const c_vector<double, DIM>& r_node_B_location = edge_iter.GetNodeB()->rGetLocation();

        c_vector<double, DIM> edge_vector = r_mesh.GetVectorFromAtoB(r_node_B_location, r_node_A_location);
        double edge_length = norm_2(edge_vector);

        assert(edge_length > DBL_EPSILON);

        c_vector<double, DIM> edge_gradient = edge_vector/edge_length;

        // Add the membrane surface tension contribution to each vertex
        unsigned elem_index = edge_iter.GetElemIndex();
        double cell_target_perimeter = 2*sqrt(M_PI*target_areas[elem_index]);
        c_vector<double, DIM> membrane_surface_tension_force_on_node_A = -2*GetNagaiHondaMembraneSurfaceEnergyParameter()*(element_perimeters[elem_index] - cell_target_perimeter)*edge_gradient;
        c_vector<double, DIM> membrane_surface_tension_force_on_node_B = -membrane_surface_tension_force_on_node_A;
        r_mesh.GetNode(node_A_index)->AddAppliedForceContribution(membrane_surface_tension_force_on_node_A);
        r_mesh.GetNode(node_B_index)->AddAppliedForceContribution(membrane_surface_tension_force_on_node_B);

        unsigned other_elem_index = edge_iter.GetOtherElemIndex();
        if (other_elem_index != UINT_MAX)
        {
            cell_target_perimeter = 2*sqrt(M_PI*target_areas[other_elem_index]);
            membrane_surface_tension_force_on_node_A = -2*GetNagaiHondaMembraneSurfaceEnergyParameter()*(element_perimeters[other_elem_index] - cell_target_perimeter)*edge_gradient;
            membrane_surface_tension_force_on_node_B = -membrane_surface_tension_force_on_node_A;
            r_mesh.GetNode(node_A_index)->AddAppliedForceContribution(membrane_surface_tension_force_on_node_A);
            r_mesh.GetNode(node_B_index)->AddAppliedForceContribution(membrane_surface_tension_force_on_node_B);
        }

        // Add the adhesion contribution to each vertex
        c_vector<double, DIM> adhesion_force_on_node_A = -adhesion_parameter*edge_gradient;
        c_vector<double, DIM> adhesion_force_on_node_B = -adhesion_force_on_node_A;
        r_mesh.GetNode(node_A_index)->AddAppliedForceContribution(adhesion_force_on_node_A);
        r_mesh.GetNode(node_B_index)->AddAppliedForceContribution(adhesion_force_on_node_B);
    }
}

template<unsigned DIM>
double NagaiHondaForce<DIM>::GetAdhesionParameter(Node<DIM>* pNodeA, Node<DIM>* pNodeB, VertexBasedCellPopulation<DIM>& rVertexCellPopulation)
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

    double adhesion_parameter = GetNagaiHondaCellCellAdhesionEnergyParameter();

    // If the edge corresponds to a single element, then the cell is on the boundary
    if (shared_elements.size() == 1)
    {
        adhesion_parameter = GetNagaiHondaCellBoundaryAdhesionEnergyParameter();
    }

    return adhesion_parameter;
}

template<unsigned DIM>
double NagaiHondaForce<DIM>::GetNagaiHondaDeformationEnergyParameter()
{
    return mNagaiHondaDeformationEnergyParameter;
}

template<unsigned DIM>
double NagaiHondaForce<DIM>::GetNagaiHondaMembraneSurfaceEnergyParameter()
{
    return mNagaiHondaMembraneSurfaceEnergyParameter;
}

template<unsigned DIM>
double NagaiHondaForce<DIM>::GetNagaiHondaCellCellAdhesionEnergyParameter()
{
    return mNagaiHondaCellCellAdhesionEnergyParameter;
}

template<unsigned DIM>
double NagaiHondaForce<DIM>::GetNagaiHondaCellBoundaryAdhesionEnergyParameter()
{
    return mNagaiHondaCellBoundaryAdhesionEnergyParameter;
}

template<unsigned DIM>
void NagaiHondaForce<DIM>::SetNagaiHondaDeformationEnergyParameter(double deformationEnergyParameter)
{
    mNagaiHondaDeformationEnergyParameter = deformationEnergyParameter;
}

template<unsigned DIM>
void NagaiHondaForce<DIM>::SetNagaiHondaMembraneSurfaceEnergyParameter(double membraneSurfaceEnergyParameter)
{
    mNagaiHondaMembraneSurfaceEnergyParameter = membraneSurfaceEnergyParameter;
}

template<unsigned DIM>
void NagaiHondaForce<DIM>::SetNagaiHondaCellCellAdhesionEnergyParameter(double cellCellAdhesionEnergyParameter)
{
    mNagaiHondaCellCellAdhesionEnergyParameter = cellCellAdhesionEnergyParameter;
}

template<unsigned DIM>
void NagaiHondaForce<DIM>::SetNagaiHondaCellBoundaryAdhesionEnergyParameter(double cellBoundaryAdhesionEnergyParameter)
{
    mNagaiHondaCellBoundaryAdhesionEnergyParameter = cellBoundaryAdhesionEnergyParameter;
}

template<unsigned DIM>
void NagaiHondaForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<NagaiHondaDeformationEnergyParameter>" << mNagaiHondaDeformationEnergyParameter << "</NagaiHondaDeformationEnergyParameter>\n";
    *rParamsFile << "\t\t\t<NagaiHondaMembraneSurfaceEnergyParameter>" << mNagaiHondaMembraneSurfaceEnergyParameter << "</NagaiHondaMembraneSurfaceEnergyParameter>\n";
    *rParamsFile << "\t\t\t<NagaiHondaCellCellAdhesionEnergyParameter>" << mNagaiHondaCellCellAdhesionEnergyParameter << "</NagaiHondaCellCellAdhesionEnergyParameter>\n";
    *rParamsFile << "\t\t\t<NagaiHondaCellBoundaryAdhesionEnergyParameter>" << mNagaiHondaCellBoundaryAdhesionEnergyParameter << "</NagaiHondaCellBoundaryAdhesionEnergyParameter>\n";

    // Call method on direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class NagaiHondaForce<1>;
template class NagaiHondaForce<2>;
template class NagaiHondaForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(NagaiHondaForce)
