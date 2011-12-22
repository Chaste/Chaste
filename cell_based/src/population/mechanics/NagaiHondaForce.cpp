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

#include "NagaiHondaForce.hpp"

template<unsigned DIM>
NagaiHondaForce<DIM>::NagaiHondaForce()
   : AbstractForce<DIM>(),
     mNagaiHondaDeformationEnergyParameter(100.0), // This is 1.0 in the Nagai & Honda paper
     mNagaiHondaMembraneSurfaceEnergyParameter(10.0), // This is 0.1 the Nagai & Honda paper
     mNagaiHondaCellCellAdhesionEnergyParameter(1.0), // This is 0.01 the Nagai & Honda paper
     mNagaiHondaCellBoundaryAdhesionEnergyParameter(1.0), // This is 0.01 the Nagai & Honda paper
     mMatureCellTargetArea(1.0)
{
}

template<unsigned DIM>
NagaiHondaForce<DIM>::~NagaiHondaForce()
{
}

template<unsigned DIM>
void NagaiHondaForce<DIM>::AddForceContribution(std::vector<c_vector<double, DIM> >& rForces,
                                                AbstractCellPopulation<DIM>& rCellPopulation)
{
    // Throw an exception message if not using a VertexBasedCellPopulation
    if (dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation) == NULL)
    {
        EXCEPTION("NagaiHondaForce is to be used with a VertexBasedCellPopulation only");
    }

    // Helper variable that is a static cast of the cell population
    VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);

    // Iterate over vertices in the cell population
    for (unsigned node_index=0; node_index<p_cell_population->GetNumNodes(); node_index++)
    {
        // Compute the force on this node

        /*
         * The force on this Node is given by the gradient of the total free
         * energy of the CellPopulation, evaluated at the position of the vertex. This
         * free energy is the sum of the free energies of all CellPtrs in
         * the cell population. The free energy of each CellPtr is comprised of three
         * parts - a cell deformation energy, a membrane surface tension energy
         * and an adhesion energy.
         *
         * Note that since the movement of this Node only affects the free energy
         * of the three CellPtrs containing it, we can just consider the
         * contributions to the free energy gradient from each of these three
         * CellPtrs.
         */

        c_vector<double, DIM> deformation_contribution = zero_vector<double>(DIM);
        c_vector<double, DIM> membrane_surface_tension_contribution = zero_vector<double>(DIM);
        c_vector<double, DIM> adhesion_contribution = zero_vector<double>(DIM);

        // Find the indices of the elements owned by this node
        std::set<unsigned> containing_elem_indices = p_cell_population->GetNode(node_index)->rGetContainingElementIndices();

        // Iterate over these elements
        for (std::set<unsigned>::iterator iter = containing_elem_indices.begin();
             iter != containing_elem_indices.end();
             ++iter)
        {
            // Get this element and its index
            VertexElement<DIM, DIM>* p_element = p_cell_population->GetElement(*iter);
            unsigned element_index = p_element->GetIndex();

            // Find the local index of this node in this element
            unsigned local_index = p_element->GetNodeLocalIndex(node_index);

            /******** Start of deformation force calculation ********/

            // Compute the area of this element and its gradient at this node
            double element_area = p_cell_population->rGetMesh().GetVolumeOfElement(*iter);
            c_vector<double, DIM> element_area_gradient = p_cell_population->rGetMesh().GetAreaGradientOfElementAtNode(p_element, local_index);

            // Get the target area of the cell
            double cell_target_area = GetTargetAreaOfCell(p_cell_population->GetCellUsingLocationIndex(element_index));

            // Add the force contribution from this cell's deformation energy (note the minus sign)
            deformation_contribution -= 2*GetNagaiHondaDeformationEnergyParameter()*(element_area - cell_target_area)*element_area_gradient;

            /******** End of deformation force calculation *************/

            /******** Start of membrane force calculation ***********/

            // Compute the perimeter of the element and its gradient at this node
            double element_perimeter = p_cell_population->rGetMesh().GetSurfaceAreaOfElement(*iter);
            c_vector<double, DIM> element_perimeter_gradient = p_cell_population->rGetMesh().GetPerimeterGradientOfElementAtNode(p_element, local_index);

            // Get the target perimeter of the cell
            double cell_target_perimeter = 2*sqrt(M_PI*cell_target_area);

            // Add the force contribution from this cell's membrane surface tension (note the minus sign)
            membrane_surface_tension_contribution -= 2*GetNagaiHondaMembraneSurfaceEnergyParameter()*(element_perimeter - cell_target_perimeter)*element_perimeter_gradient;

            /******** End of membrane force calculation **********/

            /******** Start of adhesion force calculation ***********/

            // Get the current, previous and next nodes in this element
            Node<DIM>* p_current_node = p_element->GetNode(local_index);

            unsigned previous_node_local_index = (p_element->GetNumNodes()+local_index-1)%(p_element->GetNumNodes());
            Node<DIM>* p_previous_node = p_element->GetNode(previous_node_local_index);

            unsigned next_node_local_index = (local_index+1)%(p_element->GetNumNodes());
            Node<DIM>* p_next_node = p_element->GetNode(next_node_local_index);

            // Compute the adhesion parameter for each of these edges
            double previous_edge_adhesion_parameter;
            double next_edge_adhesion_parameter;

            previous_edge_adhesion_parameter = GetAdhesionParameter(p_previous_node, p_current_node, *p_cell_population);
            next_edge_adhesion_parameter = GetAdhesionParameter(p_current_node, p_next_node, *p_cell_population);

            // Compute the gradient of the edge of the cell ending in this node
            c_vector<double, DIM> previous_edge_gradient = p_cell_population->rGetMesh().GetPreviousEdgeGradientOfElementAtNode(p_element, local_index);

            // Compute the gradient of the edge of the cell starting in this node
            c_vector<double, DIM> next_edge_gradient = p_cell_population->rGetMesh().GetNextEdgeGradientOfElementAtNode(p_element, local_index);

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
double NagaiHondaForce<DIM>::GetTargetAreaOfCell(const CellPtr pCell) const
{
    // Get target area A of a healthy cell in S, G2 or M phase
    double cell_target_area = mMatureCellTargetArea;

    double cell_age = pCell->GetAge();
    double g1_duration = pCell->GetCellCycleModel()->GetG1Duration();

    // If the cell is differentiated then its G1 duration is infinite
    if (g1_duration == DBL_MAX) // don't use magic number, compare to DBL_MAX
    {
        // This is just for fixed cell-cycle models, need to work out how to find the g1 duration
        g1_duration = pCell->GetCellCycleModel()->GetTransitCellG1Duration();
    }

    if (pCell->HasCellProperty<ApoptoticCellProperty>())
    {
        // Age of cell when apoptosis begins
        if (pCell->GetStartOfApoptosisTime() - pCell->GetBirthTime() < g1_duration)
        {
            cell_target_area *= 0.5*(1 + (pCell->GetStartOfApoptosisTime() - pCell->GetBirthTime())/g1_duration);
        }

        // The target area of an apoptotic cell decreases linearly to zero (and past it negative)
        cell_target_area = cell_target_area - 0.5*cell_target_area/(pCell->GetApoptosisTime())*(SimulationTime::Instance()->GetTime()-pCell->GetStartOfApoptosisTime());

        // Don't allow a negative target area
        if (cell_target_area < 0)
        {
            cell_target_area = 0;
        }
    }
    else
    {
        // The target area of a proliferating cell increases linearly from A/2 to A over the course of the G1 phase
        if (cell_age < g1_duration)
        {
            cell_target_area *= 0.5*(1 + cell_age/g1_duration);
        }
    }

    return cell_target_area;
}

template<unsigned DIM>
double NagaiHondaForce<DIM>::GetMatureCellTargetArea() const
{
    return mMatureCellTargetArea;
}

template<unsigned DIM>
void NagaiHondaForce<DIM>::SetMatureCellTargetArea(double matureCellTargetArea)
{
    assert(matureCellTargetArea >= 0.0);
    mMatureCellTargetArea = matureCellTargetArea;
}

template<unsigned DIM>
void NagaiHondaForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<NagaiHondaDeformationEnergyParameter>" << mNagaiHondaDeformationEnergyParameter << "</NagaiHondaDeformationEnergyParameter>\n";
    *rParamsFile << "\t\t\t<NagaiHondaMembraneSurfaceEnergyParameter>" << mNagaiHondaMembraneSurfaceEnergyParameter << "</NagaiHondaMembraneSurfaceEnergyParameter>\n";
    *rParamsFile << "\t\t\t<NagaiHondaCellCellAdhesionEnergyParameter>" << mNagaiHondaCellCellAdhesionEnergyParameter << "</NagaiHondaCellCellAdhesionEnergyParameter>\n";
    *rParamsFile << "\t\t\t<NagaiHondaCellBoundaryAdhesionEnergyParameter>" << mNagaiHondaCellBoundaryAdhesionEnergyParameter << "</NagaiHondaCellBoundaryAdhesionEnergyParameter>\n";
    *rParamsFile << "\t\t\t<MatureCellTargetArea>" << mMatureCellTargetArea << "</MatureCellTargetArea>\n";

    // Call method on direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
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
