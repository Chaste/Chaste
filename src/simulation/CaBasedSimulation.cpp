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

#include "CaBasedSimulation.hpp"
#include "CellBasedEventHandler.hpp"
#include "LogFile.hpp"

template<unsigned DIM>
CaBasedSimulation<DIM>::CaBasedSimulation(AbstractCellPopulation<DIM>& rCellPopulation,
                                          bool deleteCellPopulationInDestructor,
                                          bool initialiseCells)
    : AbstractCellBasedSimulation<DIM>(rCellPopulation,
                                       deleteCellPopulationInDestructor,
                                       initialiseCells)
{
    mIterateRandomlyOverUpdateRuleCollection = false;
    mIterateRandomlyOverCells = false;
    
    assert(dynamic_cast<CaBasedCellPopulation<DIM>*>(&rCellPopulation));
    mpStaticCastCellPopulation = static_cast<CaBasedCellPopulation<DIM>*>(&this->mrCellPopulation);
}

template<unsigned DIM>
void CaBasedSimulation<DIM>::OutputAdditionalSimulationSetup(out_stream& rParamsFile)
{
    std::vector<boost::shared_ptr<AbstractCaUpdateRule<DIM> > > update_rule_collection = mpStaticCastCellPopulation->rGetUpdateRuleCollection();

    // Loop over the collection of update rules and output info for each
    *rParamsFile << "\n\t<UpdateRules>\n";
    for (typename std::vector<boost::shared_ptr<AbstractCaUpdateRule<DIM> > >::iterator iter = update_rule_collection.begin();
         iter != update_rule_collection.end();
         ++iter)
    {
        (*iter)->OutputUpdateRuleInfo(rParamsFile);
    }
    *rParamsFile << "\t</UpdateRules>\n";
}

template<unsigned DIM>
void CaBasedSimulation<DIM>::AddUpdateRule(boost::shared_ptr<AbstractCaUpdateRule<DIM> > pUpdateRule)
{
    mpStaticCastCellPopulation->AddUpdateRule(pUpdateRule);
}

template<unsigned DIM>
void CaBasedSimulation<DIM>::SetIterateRandomlyOverUpdateRuleCollection(bool iterateRandomly)
{
    mIterateRandomlyOverUpdateRuleCollection = iterateRandomly;
}

template<unsigned DIM>
void CaBasedSimulation<DIM>::SetIterateRandomlyOverCells(bool iterateRandomly)
{
    mIterateRandomlyOverCells = iterateRandomly;
}

template<unsigned DIM>
void CaBasedSimulation<DIM>::UpdateCellLocationsAndTopology()
{
    std::vector<boost::shared_ptr<AbstractCaUpdateRule<DIM> > > update_rule_collection = mpStaticCastCellPopulation->rGetUpdateRuleCollection();
    
    // Iterate over contributions from each UpdateRule
    if (mIterateRandomlyOverUpdateRuleCollection)
    {
        // Randomly permute mUpdateRuleCollection
        std::random_shuffle(update_rule_collection.begin(), update_rule_collection.end());
    }

    for (typename std::vector<boost::shared_ptr<AbstractCaUpdateRule<DIM> > >::iterator update_iter = update_rule_collection.begin();
         update_iter != update_rule_collection.end();
         ++update_iter)
    {
        // Randomly permute cells
        if (mIterateRandomlyOverCells)
        {
            std::vector<CellPtr> cells_vector;
            for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->mrCellPopulation.Begin();
                 cell_iter != this->mrCellPopulation.End();
                 ++cell_iter)
            {
                cells_vector.push_back(*cell_iter);
            }
            std::random_shuffle(cells_vector.begin(), cells_vector.end());

            for (unsigned i=0; i<cells_vector.size(); i++)
            {
                // Get index of the node associated with cell
                unsigned current_location_index = this->mrCellPopulation.GetLocationIndexUsingCell(cells_vector[i]);

                assert(mpStaticCastCellPopulation->IsEmptySite(current_location_index) == false);

                // Get index of node the cell is to move to
                unsigned new_location_index = (*update_iter)->GetNewLocationOfCell(current_location_index, *mpStaticCastCellPopulation, this->GetDt());

                // Update the location index of the cell and free the old site
                mpStaticCastCellPopulation->MoveCell(cells_vector[i], new_location_index);
            }
        }
        else
        {
            // Iterate over all cells and update their positions
            for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->mrCellPopulation.Begin();
                 cell_iter != this->mrCellPopulation.End();
                 ++cell_iter)
            {
                // Get index of the node associated with cell
                unsigned current_location_index = this->mrCellPopulation.GetLocationIndexUsingCell(*cell_iter);

                assert(mpStaticCastCellPopulation->IsEmptySite(current_location_index) == false);

                // Get index of node the cell is to move to
                unsigned new_location_index = (*update_iter)->GetNewLocationOfCell(current_location_index, *mpStaticCastCellPopulation, this->GetDt());

                // Update the location index of the cell and free the old site
                mpStaticCastCellPopulation->MoveCell(*cell_iter, new_location_index);
            }
        }
    }
}

template<unsigned DIM>
void CaBasedSimulation<DIM>::OutputSimulationParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t<IterateRandomlyOverUpdateRuleCollection>" << mIterateRandomlyOverUpdateRuleCollection << "</IterateRandomlyOverUpdateRuleCollection>\n";
    *rParamsFile << "\t\t<IterateRandomlyOverCells>" << mIterateRandomlyOverCells << "</IterateRandomlyOverCells>\n";

    // Call method on direct parent class
    AbstractCellBasedSimulation<DIM>::OutputSimulationParameters(rParamsFile);
}

template<unsigned DIM>
void CaBasedSimulation<DIM>::UpdateCellPopulation()
{
    /*
     * If mInitialiseCells is false, then the simulation has been loaded from an archive.
     * In this case, we should not call UpdateCellPopulation() at the first time step. This is
     * because it will have already been called at the final time step prior to saving;
     * if we were to call it again now, then we would have introduced an extra call to
     * the random number generator compared to if we had not saved and loaded the simulation,
     * thus affecting results. This would be bad - we don't want saving and loading to have
     * any effect on the course of a simulation! See #1445.
     */
    bool update_cell_population_this_timestep = true;
    if (!this->mInitialiseCells && (SimulationTime::Instance()->GetTimeStepsElapsed() == 0))
    {
        update_cell_population_this_timestep = false;
    }

    if (update_cell_population_this_timestep)
    {
        AbstractCellBasedSimulation<DIM>::UpdateCellPopulation();
    }
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class CaBasedSimulation<1>;
template class CaBasedSimulation<2>;
template class CaBasedSimulation<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(CaBasedSimulation)
