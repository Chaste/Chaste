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

#include "OnLatticeSimulation.hpp"
#include "CaBasedCellPopulation.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "CellBasedEventHandler.hpp"
#include "LogFile.hpp"
#include "Version.hpp"
#include "ExecutableSupport.hpp"

template<unsigned DIM>
OnLatticeSimulation<DIM>::OnLatticeSimulation(AbstractCellPopulation<DIM>& rCellPopulation,
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
void OnLatticeSimulation<DIM>::AddCaUpdateRule(boost::shared_ptr<AbstractCaUpdateRule<DIM> > pUpdateRule)
{
    if (dynamic_cast<CaBasedCellPopulation<DIM>*>(&(this->mrCellPopulation)))
    {
        static_cast<CaBasedCellPopulation<DIM>*>(&(this->mrCellPopulation))->AddUpdateRule(pUpdateRule);
    }
}

template<unsigned DIM>
void OnLatticeSimulation<DIM>::AddPottsUpdateRule(boost::shared_ptr<AbstractPottsUpdateRule<DIM> > pUpdateRule)
{
    if (dynamic_cast<PottsBasedCellPopulation<DIM>*>(&(this->mrCellPopulation)))
    {
        static_cast<PottsBasedCellPopulation<DIM>*>(&(this->mrCellPopulation))->AddUpdateRule(pUpdateRule);
    }
}

template<unsigned DIM>
c_vector<double, DIM> OnLatticeSimulation<DIM>::CalculateCellDivisionVector(CellPtr pParentCell)
{
    ///\todo do something for Potts models here
    return zero_vector<double>(DIM);
}

template<unsigned DIM>
void OnLatticeSimulation<DIM>::UpdateCellLocationsAndTopology()
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
void OnLatticeSimulation<DIM>::SetupSolve()
{
    if (mOutputCellVelocities)
    {
        OutputFileHandler output_file_handler2(this->mSimulationOutputDirectory+"/", false);
        mpCellVelocitiesFile = output_file_handler2.OpenOutputFile("cellvelocities.dat");
    }
}

template<unsigned DIM>
void OnLatticeSimulation<DIM>::AfterSolve()
{
    if (mOutputCellVelocities)
    {
        mpCellVelocitiesFile->close();
    }
}

template<unsigned DIM>
bool OnLatticeSimulation<DIM>::GetOutputCellVelocities()
{
    return mOutputCellVelocities;
}

template<unsigned DIM>
void OnLatticeSimulation<DIM>::SetOutputCellVelocities(bool outputCellVelocities)
{
    mOutputCellVelocities = outputCellVelocities;
}

template<unsigned DIM>
void OnLatticeSimulation<DIM>::UpdateCellPopulation()
{
    bool update_cell_population_this_timestep = true;
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
        if (!this->mInitialiseCells && (SimulationTime::Instance()->GetTimeStepsElapsed() == 0))
        {
            update_cell_population_this_timestep = false;
        }
    }

    if (update_cell_population_this_timestep)
    {
        AbstractCellBasedSimulation<DIM>::UpdateCellPopulation();
    }
}

template<unsigned DIM>
void OnLatticeSimulation<DIM>::OutputAdditionalSimulationSetup(out_stream& rParamsFile)
{
    // Loop over the collection of update rules and output info for each
    *rParamsFile << "\n\t<UpdateRules>\n";
    if (dynamic_cast<PottsBasedCellPopulation<DIM>*>(&(this->mrCellPopulation)))
    {
        std::vector<boost::shared_ptr<AbstractPottsUpdateRule<DIM> > > collection =
            static_cast<PottsBasedCellPopulation<DIM>*>(&(this->mrCellPopulation))->rGetUpdateRuleCollection();

        for (typename std::vector<boost::shared_ptr<AbstractPottsUpdateRule<DIM> > >::iterator iter = collection.begin();
             iter != collection.end();
             ++iter)
        {
            (*iter)->OutputUpdateRuleInfo(rParamsFile);
        }
    }
    else
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
void OnLatticeSimulation<DIM>::OutputSimulationParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t<OutputCellVelocities>" << mOutputCellVelocities << "</OutputCellVelocities>\n";

    // Call method on direct parent class
    AbstractCellBasedSimulation<DIM>::OutputSimulationParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class OnLatticeSimulation<1>;
template class OnLatticeSimulation<2>;
template class OnLatticeSimulation<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(OnLatticeSimulation)
