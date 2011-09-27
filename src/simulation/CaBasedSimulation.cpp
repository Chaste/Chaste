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
#include "PottsBasedCellPopulation.hpp"
#include "CaBasedCellPopulation.hpp"
#include "CellBasedEventHandler.hpp"
#include "LogFile.hpp"
#include "Version.hpp"
#include "ExecutableSupport.hpp"

template<unsigned DIM>
CaBasedSimulation<DIM>::CaBasedSimulation(AbstractCellPopulation<DIM>& rCellPopulation,
                                          bool deleteCellPopulationInDestructor,
                                          bool initialiseCells)
    : AbstractCellBasedSimulation<DIM>(rCellPopulation,
                                       deleteCellPopulationInDestructor,
                                       initialiseCells),
      mOutputCellVelocities(false)
{
    if (!dynamic_cast<AbstractOnLatticeCellPopulation<DIM>*>(&rCellPopulation))
    {
        EXCEPTION("OnLatticeSimulations require a subclass of AbstractOnLatticeCellPopulation.");
    }

    this->mDt = 1.0/120.0; // 30 seconds
}

template<unsigned DIM>
void CaBasedSimulation<DIM>::AddUpdateRule(boost::shared_ptr<AbstractCaUpdateRule<DIM> > pUpdateRule)
{
    if (dynamic_cast<CaBasedCellPopulation<DIM>*>(&(this->mrCellPopulation)))
    {
        static_cast<CaBasedCellPopulation<DIM>*>(&(this->mrCellPopulation))->AddUpdateRule(pUpdateRule);
    }
}

template<unsigned DIM>
void CaBasedSimulation<DIM>::UpdateCellLocationsAndTopology()
{
    // Store current locations of cell centres if required
    std::vector<c_vector<double, DIM> > old_cell_locations;
    unsigned num_cells = this->mrCellPopulation.GetNumRealCells();
    old_cell_locations.reserve(num_cells);
    if (mOutputCellVelocities)
    {
        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->mrCellPopulation.Begin();
             cell_iter != this->mrCellPopulation.End();
             ++cell_iter)
        {
            unsigned index = this->mrCellPopulation.GetLocationIndexUsingCell(*cell_iter);
            old_cell_locations[index] = this->mrCellPopulation.GetLocationOfCellCentre(*cell_iter);
        }
    }

    // Update cell locations
    CellBasedEventHandler::BeginEvent(CellBasedEventHandler::POSITION);
    static_cast<AbstractOnLatticeCellPopulation<DIM>*>(&(this->mrCellPopulation))->UpdateCellLocations(this->GetDt());
    CellBasedEventHandler::EndEvent(CellBasedEventHandler::POSITION);

    // Write cell velocities to file if required
    if (mOutputCellVelocities)
    {
        if (SimulationTime::Instance()->GetTimeStepsElapsed()%this->mSamplingTimestepMultiple == 0)
        {
            *mpCellVelocitiesFile << SimulationTime::Instance()->GetTime() << "\t";

            for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->mrCellPopulation.Begin();
                 cell_iter != this->mrCellPopulation.End();
                 ++cell_iter)
            {
                unsigned index = this->mrCellPopulation.GetLocationIndexUsingCell(*cell_iter);
                const c_vector<double,DIM>& position = this->mrCellPopulation.GetLocationOfCellCentre(*cell_iter);
                c_vector<double, DIM> velocity = (position - old_cell_locations[index])/this->mDt;

                *mpCellVelocitiesFile << index  << " ";
                for (unsigned i=0; i<DIM; i++)
                {
                    *mpCellVelocitiesFile << position[i] << " ";
                }
                for (unsigned i=0; i<DIM; i++)
                {
                    *mpCellVelocitiesFile << velocity[i] << " ";
                }
            }
            *mpCellVelocitiesFile << "\n";
        }
    }
}

template<unsigned DIM>
void CaBasedSimulation<DIM>::SetupSolve()
{
    if (mOutputCellVelocities)
    {
        OutputFileHandler output_file_handler(this->mSimulationOutputDirectory+"/", false);
        mpCellVelocitiesFile = output_file_handler.OpenOutputFile("cellvelocities.dat");
    }
}

template<unsigned DIM>
void CaBasedSimulation<DIM>::AfterSolve()
{
    if (mOutputCellVelocities)
    {
        mpCellVelocitiesFile->close();
    }
}

template<unsigned DIM>
bool CaBasedSimulation<DIM>::GetOutputCellVelocities()
{
    return mOutputCellVelocities;
}

template<unsigned DIM>
void CaBasedSimulation<DIM>::SetOutputCellVelocities(bool outputCellVelocities)
{
    mOutputCellVelocities = outputCellVelocities;
}

template<unsigned DIM>
void CaBasedSimulation<DIM>::UpdateCellPopulation()
{
    if (dynamic_cast<CaBasedCellPopulation<DIM>*>(&(this->mrCellPopulation)))
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
}

template<unsigned DIM>
void CaBasedSimulation<DIM>::OutputAdditionalSimulationSetup(out_stream& rParamsFile)
{
    // Loop over the collection of update rules and output info for each
    *rParamsFile << "\n\t<UpdateRules>\n";
    if (dynamic_cast<CaBasedCellPopulation<DIM>*>(&(this->mrCellPopulation)))
    {
        std::vector<boost::shared_ptr<AbstractCaUpdateRule<DIM> > > collection = 
            static_cast<CaBasedCellPopulation<DIM>*>(&(this->mrCellPopulation))->rGetUpdateRuleCollection();

        for (typename std::vector<boost::shared_ptr<AbstractCaUpdateRule<DIM> > >::iterator iter = collection.begin();
             iter != collection.end();
             ++iter)
        {
            (*iter)->OutputUpdateRuleInfo(rParamsFile);
        }
    }
    *rParamsFile << "\t</UpdateRules>\n";
}

template<unsigned DIM>
void CaBasedSimulation<DIM>::OutputSimulationParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t<OutputCellVelocities>" << mOutputCellVelocities << "</OutputCellVelocities>\n";

    // Call method on direct parent class
    AbstractCellBasedSimulation<DIM>::OutputSimulationParameters(rParamsFile);
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
