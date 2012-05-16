/*

Copyright (c) 2005-2012, University of Oxford.
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

#include "AbstractCentreBasedCellPopulation.hpp"

template<unsigned DIM>
AbstractCentreBasedCellPopulation<DIM>::AbstractCentreBasedCellPopulation( AbstractMesh<DIM, DIM>& rMesh,
																	std::vector<CellPtr>& rCells,
                                                                  const std::vector<unsigned> locationIndices)
    : AbstractOffLatticeCellPopulation<DIM>(rMesh, rCells, locationIndices),
      mMeinekeDivisionSeparation(0.3) // educated guess
{
}

template<unsigned DIM>
AbstractCentreBasedCellPopulation<DIM>::AbstractCentreBasedCellPopulation(AbstractMesh<DIM, DIM>& rMesh)
    : AbstractOffLatticeCellPopulation<DIM>(rMesh),
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
    unsigned index = this->GetLocationIndexUsingCell(pCell);
    return this->GetNode(index);
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
    this->SetCellUsingLocationIndex(new_node_index,pNewCell);
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
        unsigned node_index = this->GetLocationIndexUsingCell((*cell_iter));

        // Get damping constant for node
        double damping_const = this->GetDampingConstant(node_index);

        // Get displacement
        c_vector<double,DIM> displacement=dt*rNodeForces[node_index]/damping_const;

        // Throws an exception if the cell movement goes beyond mAbsoluteMovementThreshold
        if (norm_2(displacement) > this->mAbsoluteMovementThreshold)
        {
            EXCEPTION("Cells are moving more than the AbsoluteMovementThreshold. Use a smaller timestep to avoid this exception.");
        }

        // Get new node location
        c_vector<double, DIM> new_node_location = this->GetNode(node_index)->rGetLocation() + displacement;

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
bool AbstractCentreBasedCellPopulation<DIM>::IsParticle(unsigned index)
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
        if (this->IsCellAttachedToLocationIndex(node_index))
        {
            node_corresponds_to_dead_cell = this->GetCellUsingLocationIndex(node_index)->IsDead();
        }

        // Write cell data to file
        if (!(this->GetNode(node_index)->IsDeleted())
        		&& !node_corresponds_to_dead_cell
        		&& !(this->IsParticle(node_index)))
        {
            this->GenerateCellResults(node_index, cell_type_counter, cell_cycle_phase_counter);
        }
    }

    this->WriteCellResultsToFiles(cell_type_counter, cell_cycle_phase_counter);
}

template<unsigned DIM>
void AbstractCentreBasedCellPopulation<DIM>::WriteTimeAndNodeResultsToFiles()
{
    OutputFileHandler output_file_handler(this->mDirPath, false);

    PetscTools::BeginRoundRobin();
    {
		if(!PetscTools::AmMaster() || SimulationTime::Instance()->GetTimeStepsElapsed()!=0)
		{
			this->mpVizNodesFile = output_file_handler.OpenOutputFile("results.viznodes", std::ios::app);
			this->mpVizBoundaryNodesFile = output_file_handler.OpenOutputFile("results.vizboundarynodes", std::ios::app);
		}
		if(PetscTools::AmMaster())
		{
			double time = SimulationTime::Instance()->GetTime();

			*this->mpVizNodesFile << time << "\t";
			*this->mpVizBoundaryNodesFile << time << "\t";
		}


		// Write node data to file
		for (typename AbstractMesh<DIM, DIM>::NodeIterator node_iter = this->mrMesh.GetNodeIteratorBegin();
				node_iter != this->mrMesh.GetNodeIteratorEnd();
				++node_iter)
		{

		    /*
			 * Hack that covers the case where the node in an AbstractCentreBasedCellPopulation
			 * is associated with a cell that has just been killed (#1129). This breaks the
			 * vertex visualizer when apoptotic cells are involved.
			 */
			bool node_corresponds_to_dead_cell = false;
			if (this->IsCellAttachedToLocationIndex(node_iter->GetIndex()))
			{
				node_corresponds_to_dead_cell = this->GetCellUsingLocationIndex(node_iter->GetIndex())->IsDead();
			}

			// Write node data to file
			if (!(node_iter->IsDeleted()) && !node_corresponds_to_dead_cell)
			{
				const c_vector<double,DIM>& position = node_iter->rGetLocation();

				for (unsigned i=0; i<DIM; i++)
				{
					*this->mpVizNodesFile << position[i] << " ";
				}
				*this->mpVizBoundaryNodesFile << node_iter->IsBoundaryNode() << " ";
			}
		}

		if(PetscTools::AmTopMost())
		{
			*this->mpVizNodesFile << "\n";
			*this->mpVizBoundaryNodesFile << "\n";
		}

		this->mpVizNodesFile->close();
		this->mpVizBoundaryNodesFile->close();
    }
    PetscTools::EndRoundRobin();
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
