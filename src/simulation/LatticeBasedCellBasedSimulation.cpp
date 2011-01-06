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

#include "LatticeBasedCellBasedSimulation.hpp"
#include "CellBasedEventHandler.hpp"
#include "LogFile.hpp"

template<unsigned DIM>
LatticeBasedCellBasedSimulation<DIM>::LatticeBasedCellBasedSimulation(AbstractCellPopulation<DIM>& rCellPopulation,
                  std::vector<AbstractUpdateRule<DIM>*> updateRuleCollection,
                  bool iterateRandomlyOverUpdateRuleCollection,
                  bool iterateRandomlyOverCells,
                  bool deleteCellPopulationAndForceCollection,
                  bool initialiseCells)
    : CellBasedSimulation<DIM>(rCellPopulation,
                  deleteCellPopulationAndForceCollection,
                  initialiseCells),
    mUpdateRuleCollection(updateRuleCollection),
    mAllocatedMemoryForUpdateRuleCollection(deleteCellPopulationAndForceCollection),
    mIterateRandomlyOverUpdateRuleCollection(iterateRandomlyOverUpdateRuleCollection),
    mIterateRandomlyOverCells(iterateRandomlyOverCells)
{
    mpStaticCastCellPopulation = static_cast<LatticeBasedCellPopulation<DIM>*>(&this->mrCellPopulation);
}

template<unsigned DIM>
LatticeBasedCellBasedSimulation<DIM>::~LatticeBasedCellBasedSimulation()
{
    if (mAllocatedMemoryForUpdateRuleCollection)
    {
        for (typename std::vector<AbstractUpdateRule<DIM>*>::iterator update_iter = mUpdateRuleCollection.begin();
             update_iter != mUpdateRuleCollection.end();
             ++update_iter)
        {
            delete *update_iter;
        }
    }
}

template<unsigned DIM>
c_vector<double, DIM> LatticeBasedCellBasedSimulation<DIM>::CalculateCellDivisionVector(CellPtr pParentCell)
{
    c_vector<double, DIM> axis_of_division = zero_vector<double>(DIM);
    return axis_of_division;
}

template<unsigned DIM>
void LatticeBasedCellBasedSimulation<DIM>::UpdateCellLocations()
{
    // Iterate over contributions from each UpdateRule
    if (mIterateRandomlyOverUpdateRuleCollection)
    {
        // Randomly permute mUpdateRuleCollection
        std::random_shuffle(mUpdateRuleCollection.begin(), mUpdateRuleCollection.end());
    }

    for (typename std::vector<AbstractUpdateRule<DIM>*>::iterator update_iter = mUpdateRuleCollection.begin();
         update_iter != mUpdateRuleCollection.end();
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

///\todo Much of this code can probably merged with CellBasedSimulation::Solve()
template<unsigned DIM>
void LatticeBasedCellBasedSimulation<DIM>::Solve()
{
    CellBasedEventHandler::BeginEvent(CellBasedEventHandler::EVERYTHING);
    CellBasedEventHandler::BeginEvent(CellBasedEventHandler::SETUP);

    // Set up the simulation time
    SimulationTime* p_simulation_time = SimulationTime::Instance();
    double current_time = p_simulation_time->GetTime();

    unsigned num_time_steps = (unsigned) ((this->mEndTime-current_time)/this->mDt+0.5);

    if (current_time > 0) // use the reset function if necessary
    {
        p_simulation_time->ResetEndTimeAndNumberOfTimeSteps(this->mEndTime, num_time_steps);
    }
    else
    {
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(this->mEndTime, num_time_steps);
    }

    if (this->mOutputDirectory == "")
    {
        EXCEPTION("OutputDirectory not set");
    }

    double time_now = p_simulation_time->GetTime();
    std::ostringstream time_string;
    time_string << time_now;

    std::string results_directory = this->mOutputDirectory +"/results_from_time_" + time_string.str();
    this->mSimulationOutputDirectory = results_directory;

    ///////////////////////////////////////////////////////////
    // Set up Simulation
    ///////////////////////////////////////////////////////////

    // Create output files for the visualizer
    OutputFileHandler output_file_handler(results_directory+"/", true);

    this->mrCellPopulation.CreateOutputFiles(results_directory+"/", false);

    this->mpVizSetupFile = output_file_handler.OpenOutputFile("results.vizsetup");

    this->SetupSolve();

    // Age the cells to the correct time. Note that cells are created with
    // negative birth times so that some are initially almost ready to divide.
    LOG(1, "Setting up cells...");
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->mrCellPopulation.Begin();
         cell_iter != this->mrCellPopulation.End();
         ++cell_iter)
    {
        // We don't use the result; this call is just to force the cells to age
        // to the current time running their cell cycle models to get there
        cell_iter->ReadyToDivide();
    }
    LOG(1, "\tdone\n");

    // Write initial conditions to file for the visualizer
    this->WriteVisualizerSetupFile();

    *this->mpVizSetupFile << std::flush;

    this->mrCellPopulation.WriteResultsToFiles();

    CellBasedEventHandler::EndEvent(CellBasedEventHandler::SETUP);

    /////////////////////////////////////////////////////////////////////
    // Main time loop
    /////////////////////////////////////////////////////////////////////

    while ((p_simulation_time->GetTimeStepsElapsed() < num_time_steps) && !(this->StoppingEventHasOccurred()) )
    {
        LOG(1, "--TIME = " << p_simulation_time->GetTime() << "\n");

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
        if (!this->mInitialiseCells && (p_simulation_time->GetTimeStepsElapsed() == 0))
        {
            update_cell_population_this_timestep = false;
        }

        if (update_cell_population_this_timestep)
        {
            /**
             * This function calls:
             * DoCellRemoval()
             * DoCellBirth()
             * CellPopulation::Update()
             */
            this->UpdateCellPopulation();
        }

        //////////////////////////////////////////////////////
        // Calculate Cell Movements and Update Cell Locations
        //////////////////////////////////////////////////////
        CellBasedEventHandler::BeginEvent(CellBasedEventHandler::POSITION);
        UpdateCellLocations();
        CellBasedEventHandler::EndEvent(CellBasedEventHandler::POSITION);

        //////////////////////////////////////////
        // PostSolve, which may be implemented by
        // child classes (eg to solve PDEs)
        //////////////////////////////////////////
        this->PostSolve();

        // Increment simulation time here, so results files look sensible
        p_simulation_time->IncrementTimeOneStep();

        ////////////////////////////
        // Output current results
        ////////////////////////////
        CellBasedEventHandler::BeginEvent(CellBasedEventHandler::OUTPUT);

        // Write results to file
        if (p_simulation_time->GetTimeStepsElapsed()%this->mSamplingTimestepMultiple==0)
        {
            this->mrCellPopulation.WriteResultsToFiles();
        }

        CellBasedEventHandler::EndEvent(CellBasedEventHandler::OUTPUT);
    }

    LOG(1, "--END TIME = " << SimulationTime::Instance()->GetTime() << "\n");

    /*
     * Carry out a final update so that cell population is coherent with new cell positions.
     * NB cell birth/death still need to be checked because they may be spatially-dependent.
     */
    this->UpdateCellPopulation();

    this->AfterSolve();

    CellBasedEventHandler::BeginEvent(CellBasedEventHandler::OUTPUT);

    this->mrCellPopulation.CloseOutputFiles();

    *this->mpVizSetupFile << "Complete\n";
    this->mpVizSetupFile->close();

    CellBasedEventHandler::EndEvent(CellBasedEventHandler::OUTPUT);

    CellBasedEventHandler::EndEvent(CellBasedEventHandler::EVERYTHING);
}

template<unsigned DIM>
const std::vector<AbstractUpdateRule<DIM>*> LatticeBasedCellBasedSimulation<DIM>::rGetUpdateRuleCollection() const
{
    return mUpdateRuleCollection;
}

template<unsigned DIM>
void LatticeBasedCellBasedSimulation<DIM>::OutputSimulationParameters(out_stream& rParamsFile)
{
	// Currently this is not called from LatticeBasedCellBasedSimulation see #1453 for discussion on this for centre- and vertex-based cell population.
	EXCEPTION("OutputSimulationParameters() is not yet implemented for LatticeBasedCellBasedSimulation see #1453");

	// Call direct parent class
	//CellBasedSimulation<DIM>::OutputSimulationParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class LatticeBasedCellBasedSimulation<1>;
template class LatticeBasedCellBasedSimulation<2>;
template class LatticeBasedCellBasedSimulation<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(LatticeBasedCellBasedSimulation)
