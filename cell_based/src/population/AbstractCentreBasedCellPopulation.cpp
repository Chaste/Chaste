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

#include "AbstractCentreBasedCellPopulation.hpp"

template<unsigned DIM>
AbstractCentreBasedCellPopulation<DIM>::AbstractCentreBasedCellPopulation(std::vector<CellPtr>& rCells,
                                                                  const std::vector<unsigned> locationIndices)
    : AbstractOffLatticeCellPopulation<DIM>(rCells, locationIndices),
      mMeinekeDivisionSeparation(0.3) // educated guess
{
}

template<unsigned DIM>
AbstractCentreBasedCellPopulation<DIM>::AbstractCentreBasedCellPopulation()
    : AbstractOffLatticeCellPopulation<DIM>(),
      mMeinekeDivisionSeparation(0.3) // educated guess
{
}

template<unsigned DIM>
c_vector<double, DIM> AbstractCentreBasedCellPopulation<DIM>::GetLocationOfCellCentre(CellPtr pCell)
{
    return GetNodeCorrespondingToCell(pCell)->rGetLocation();
}

template<unsigned DIM>
Node<DIM>* AbstractCentreBasedCellPopulation<DIM>::GetNodeCorrespondingToCell(CellPtr pCell)
{
    assert(this->mCellLocationMap.find(pCell.get()) != this->mCellLocationMap.end());

    return this->GetNode(this->mCellLocationMap[pCell.get()]);
}

template<unsigned DIM>
CellPtr AbstractCentreBasedCellPopulation<DIM>::AddCell(CellPtr pNewCell, const c_vector<double,DIM>& rCellDivisionVector, CellPtr pParentCell)
{
    // Create a new node
    Node<DIM>* p_new_node = new Node<DIM>(this->GetNumNodes(), rCellDivisionVector, false);   // never on boundary
    unsigned new_node_index = AddNode(p_new_node); // use copy constructor so it doesn't matter that new_node goes out of scope

    // Update cells vector
    this->mCells.push_back(pNewCell);

    // Update mappings between cells and location indices
    this->mLocationCellMap[new_node_index] = pNewCell;
    this->mCellLocationMap[pNewCell.get()] = new_node_index;

    return pNewCell;
}

template<unsigned DIM>
bool AbstractCentreBasedCellPopulation<DIM>::IsCellAssociatedWithADeletedLocation(CellPtr pCell)
{
    return GetNodeCorrespondingToCell(pCell)->IsDeleted();
}

template<unsigned DIM>
void AbstractCentreBasedCellPopulation<DIM>::UpdateNodeLocations(const std::vector< c_vector<double, DIM> >& rNodeForces, double dt)
{
    // Iterate over all nodes associated with real cells to update their positions
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->Begin();
         cell_iter != this->End();
         ++cell_iter)
    {
        // Get index of node associated with cell
        unsigned node_index = this->mCellLocationMap[(*cell_iter).get()];

        // Get damping constant for node
        double damping_const = this->GetDampingConstant(node_index);

        // Get new node location
        c_vector<double, DIM> new_node_location = this->GetNode(node_index)->rGetLocation() + dt*rNodeForces[node_index]/damping_const;

        // Create ChastePoint for new node location
        ChastePoint<DIM> new_point(new_node_location);

        // Move the node
        this->SetNode(node_index, new_point);
    }
}

template<unsigned DIM>
double AbstractCentreBasedCellPopulation<DIM>::GetDampingConstant(unsigned nodeIndex)
{
    CellPtr p_cell = this->GetCellUsingLocationIndex(nodeIndex);
    if (p_cell->GetMutationState()->IsType<WildTypeCellMutationState>() && !p_cell->HasCellProperty<CellLabel>())
    {
        return this->GetDampingConstantNormal();
    }
    else
    {
        return this->GetDampingConstantMutant();
    }
}

template<unsigned DIM>
bool AbstractCentreBasedCellPopulation<DIM>::IsGhostNode(unsigned index)
{
    return false;
}

template<unsigned DIM>
void AbstractCentreBasedCellPopulation<DIM>::GenerateCellResults(unsigned locationIndex,
                                                             std::vector<unsigned>& rCellProliferativeTypeCounter,
                                                             std::vector<unsigned>& rCellCyclePhaseCounter)
{
    if (IsGhostNode(locationIndex) == true)
    {
        *(this->mpVizCellProliferativeTypesFile) << INVISIBLE_COLOUR << " ";
    }
    else
    {
        AbstractOffLatticeCellPopulation<DIM>::GenerateCellResults(locationIndex,
                                                 rCellProliferativeTypeCounter,
                                                 rCellCyclePhaseCounter);
    }
}

template<unsigned DIM>
void AbstractCentreBasedCellPopulation<DIM>::GenerateCellResultsAndWriteToFiles()
{
    // Set up cell type counter
    unsigned num_cell_types = this->mCellProliferativeTypeCount.size();
    std::vector<unsigned> cell_type_counter(num_cell_types);
    for (unsigned i=0; i<num_cell_types; i++)
    {
        cell_type_counter[i] = 0;
    }

    // Set up cell cycle phase counter
    unsigned num_cell_cycle_phases = this->mCellCyclePhaseCount.size();
    std::vector<unsigned> cell_cycle_phase_counter(num_cell_cycle_phases);
    for (unsigned i=0; i<num_cell_cycle_phases; i++)
    {
        cell_cycle_phase_counter[i] = 0;
    }

    // Write cell data to file
    for (unsigned node_index=0; node_index<this->GetNumNodes(); node_index++)
    {
        // Hack that covers the case where the node is associated with a cell that has just been killed (#1129)
        bool node_corresponds_to_dead_cell = false;
        if (this->mLocationCellMap[node_index])
        {
            node_corresponds_to_dead_cell = this->mLocationCellMap[node_index]->IsDead();
        }

        // Write cell data to file
        if (!(this->GetNode(node_index)->IsDeleted()) && !node_corresponds_to_dead_cell)
        {
            this->GenerateCellResults(node_index, cell_type_counter, cell_cycle_phase_counter);
        }
    }

    this->WriteCellResultsToFiles(cell_type_counter, cell_cycle_phase_counter);
}

template<unsigned DIM>
void AbstractCentreBasedCellPopulation<DIM>::WriteTimeAndNodeResultsToFiles()
{
    double time = SimulationTime::Instance()->GetTime();

    *this->mpVizNodesFile << time << "\t";
    *this->mpVizBoundaryNodesFile << time << "\t";

    // Write node data to file
    for (unsigned node_index=0; node_index<this->GetNumNodes(); node_index++)
    {
        // Hack that covers the case where the node in an AbstractCentreBasedCellPopulation is associated with a cell that has just been killed (#1129) This breaks the vertex visualiser when apoptotic cells are involved.
        bool node_corresponds_to_dead_cell = false;
        if (this->mLocationCellMap[node_index])
        {
            node_corresponds_to_dead_cell = this->mLocationCellMap[node_index]->IsDead();
        }

        // Write node data to file
        if (!(this->GetNode(node_index)->IsDeleted()) && !node_corresponds_to_dead_cell)
        {
            const c_vector<double,DIM>& position = this->GetNode(node_index)->rGetLocation();

            for (unsigned i=0; i<DIM; i++)
            {
                *this->mpVizNodesFile << position[i] << " ";
            }
            *this->mpVizBoundaryNodesFile << this->GetNode(node_index)->IsBoundaryNode() << " ";
        }
    }
    *this->mpVizNodesFile << "\n";
    *this->mpVizBoundaryNodesFile << "\n";
}

template<unsigned DIM>
double AbstractCentreBasedCellPopulation<DIM>::GetMeinekeDivisionSeparation()
{
    return mMeinekeDivisionSeparation;
}

template<unsigned DIM>
void AbstractCentreBasedCellPopulation<DIM>::SetMeinekeDivisionSeparation(double divisionSeparation)
{
    assert(divisionSeparation <= 1.0);
    assert(divisionSeparation >= 0.0);
    mMeinekeDivisionSeparation = divisionSeparation;
}

template<unsigned DIM>
void AbstractCentreBasedCellPopulation<DIM>::OutputCellPopulationParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t<MeinekeDivisionSeparation>" << mMeinekeDivisionSeparation << "</MeinekeDivisionSeparation>\n";

    // Call method on direct parent class
    AbstractOffLatticeCellPopulation<DIM>::OutputCellPopulationParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////

template class AbstractCentreBasedCellPopulation<1>;
template class AbstractCentreBasedCellPopulation<2>;
template class AbstractCentreBasedCellPopulation<3>;
