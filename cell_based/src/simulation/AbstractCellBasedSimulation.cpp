/*

Copyright (c) 2005-2019, University of Oxford.
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

#include <cmath>
#include <iostream>
#include <fstream>
#include <set>

#include "AbstractCellBasedSimulation.hpp"
#include "CellBasedEventHandler.hpp"
#include "LogFile.hpp"
#include "ExecutableSupport.hpp"
#include "AbstractPdeModifier.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>::AbstractCellBasedSimulation(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation,
                                              bool deleteCellPopulationInDestructor,
                                              bool initialiseCells)
    : mDt(DOUBLE_UNSET),
      mEndTime(DOUBLE_UNSET),  // hours - this is set later on
      mrCellPopulation(rCellPopulation),
      mDeleteCellPopulationInDestructor(deleteCellPopulationInDestructor),
      mInitialiseCells(initialiseCells),
      mNoBirth(false),
      mUpdateCellPopulation(true),
      mOutputDirectory(""),
      mSimulationOutputDirectory(mOutputDirectory),
      mNumBirths(0),
      mNumDeaths(0),
      mOutputDivisionLocations(false),
      mOutputCellVelocities(false),
      mSamplingTimestepMultiple(1)
{
    // Set a random seed of 0 if it wasn't specified earlier
    RandomNumberGenerator::Instance();

    if (mInitialiseCells)
    {
        mrCellPopulation.InitialiseCells();
    }

    /*
     * Specify a default time step to use, which may depend on the cell population type.
     * Note that the time step may be reset using SetDt().
     */
    mDt = rCellPopulation.GetDefaultTimeStep();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>::~AbstractCellBasedSimulation()
{
    if (mDeleteCellPopulationInDestructor)
    {
        delete &mrCellPopulation;
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>::DoCellBirth()
{
    if (mNoBirth)
    {
        return 0;
    }

    unsigned num_births_this_step = 0;

    // Iterate over all cells, seeing if each one can be divided
    for (typename AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>::Iterator cell_iter = mrCellPopulation.Begin();
         cell_iter != mrCellPopulation.End();
         ++cell_iter)
    {
        // Check if this cell is ready to divide
        double cell_age = cell_iter->GetAge();
        if (cell_age > 0.0)
        {
            if (cell_iter->ReadyToDivide())
            {
                // Check if there is room into which the cell may divide
                if (mrCellPopulation.IsRoomToDivide(*cell_iter))
                {
                    // Store parent ID for output if required
                    unsigned parent_cell_id = cell_iter->GetCellId();

                    // Create a new cell
                    CellPtr p_new_cell = cell_iter->Divide();

                    /**
                     * If required, output this location to file
                     *
                     * \todo (#2578)
                     *
                     * For consistency with the rest of the output code, consider removing the
                     * AbstractCellBasedSimulation member mOutputDivisionLocations, adding a new
                     * member mAgesAndLocationsOfDividingCells to AbstractCellPopulation, adding
                     * a new class CellDivisionLocationsWriter to the CellPopulationWriter hierarchy
                     * to output the content of mAgesAndLocationsOfDividingCells to file (remembering
                     * to clear mAgesAndLocationsOfDividingCells at each timestep), and replacing the
                     * following conditional statement with something like
                     *
                     * if (mrCellPopulation.HasWriter<CellDivisionLocationsWriter>())
                     * {
                     *     mCellDivisionLocations.push_back(new_location);
                     * }
                     */
                    if (mOutputDivisionLocations)
                    {
                        c_vector<double, SPACE_DIM> cell_location = mrCellPopulation.GetLocationOfCellCentre(*cell_iter);

                        *mpDivisionLocationFile << SimulationTime::Instance()->GetTime() << "\t";
                        for (unsigned i=0; i<SPACE_DIM; i++)
                        {
                            *mpDivisionLocationFile << cell_location[i] << "\t";
                        }
                        *mpDivisionLocationFile << "\t" << cell_age << "\t" << parent_cell_id << "\t" << cell_iter->GetCellId() << "\t" << p_new_cell->GetCellId() << "\n";
                    }

                    // Add the new cell to the cell population
                    mrCellPopulation.AddCell(p_new_cell, *cell_iter);

                    // Update counter
                    num_births_this_step++;
                }
            }
        }
    }
    return num_births_this_step;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>::DoCellRemoval()
{
    unsigned num_deaths_this_step = 0;

    /*
     * This labels cells as dead or apoptosing. It does not actually remove the cells,
     * mrCellPopulation.RemoveDeadCells() needs to be called for this.
     */
    for (typename std::vector<boost::shared_ptr<AbstractCellKiller<SPACE_DIM> > >::iterator killer_iter = mCellKillers.begin();
         killer_iter != mCellKillers.end();
         ++killer_iter)
    {
        (*killer_iter)->CheckAndLabelCellsForApoptosisOrDeath();
    }

    num_deaths_this_step += mrCellPopulation.RemoveDeadCells();

    return num_deaths_this_step;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>::SetDt(double dt)
{
    assert(dt > 0);
    mDt = dt;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>::GetDt()
{
    return mDt;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>::GetNumBirths()
{
    return mNumBirths;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>::GetNumDeaths()
{
    return mNumDeaths;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>::SetEndTime(double endTime)
{
    assert(endTime > 0);
    mEndTime = endTime;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>::SetOutputDirectory(std::string outputDirectory)
{
    mOutputDirectory = outputDirectory;
    mSimulationOutputDirectory = mOutputDirectory;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::string AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>::GetOutputDirectory()
{
    return mOutputDirectory;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>::SetSamplingTimestepMultiple(unsigned samplingTimestepMultiple)
{
    assert(samplingTimestepMultiple > 0);
    mSamplingTimestepMultiple = samplingTimestepMultiple;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>::rGetCellPopulation()
{
    return mrCellPopulation;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
const AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>::rGetCellPopulation() const
{
    return mrCellPopulation;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>::SetUpdateCellPopulationRule(bool updateCellPopulation)
{
    mUpdateCellPopulation = updateCellPopulation;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>::GetUpdateCellPopulationRule()
{
    return mUpdateCellPopulation;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>::SetNoBirth(bool noBirth)
{
    mNoBirth = noBirth;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>::AddCellKiller(boost::shared_ptr<AbstractCellKiller<SPACE_DIM> > pCellKiller)
{
    mCellKillers.push_back(pCellKiller);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>::RemoveAllCellKillers()
{
    mCellKillers.clear();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>::AddSimulationModifier(boost::shared_ptr<AbstractCellBasedSimulationModifier<ELEMENT_DIM,SPACE_DIM> > pSimulationModifier)
{
    mSimulationModifiers.push_back(pSimulationModifier);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<ELEMENT_DIM, SPACE_DIM> > >* AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>::GetSimulationModifiers()
{
    return &mSimulationModifiers;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<double> AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>::GetNodeLocation(const unsigned& rNodeIndex)
{
    std::vector<double> location;
    for (unsigned i=0; i<SPACE_DIM; i++)
    {
        location.push_back(mrCellPopulation.GetNode(rNodeIndex)->rGetLocation()[i]);
    }
    return location;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>::Solve()
{
    CellBasedEventHandler::BeginEvent(CellBasedEventHandler::EVERYTHING);
    CellBasedEventHandler::BeginEvent(CellBasedEventHandler::SETUP);

    // Set up the simulation time
    SimulationTime* p_simulation_time = SimulationTime::Instance();
    double current_time = p_simulation_time->GetTime();

    assert(mDt != DOUBLE_UNSET);  //Subclass constructors take care of this

    if (mEndTime == DOUBLE_UNSET)
    {
        EXCEPTION("SetEndTime has not yet been called.");
    }

    /*
     * Note that mDt is used here for "ideal time step". If this step doesn't divide the time remaining
     * then a *different* time step will be taken by the time-stepper. The real time-step (used in the
     * SimulationTime singleton) is currently not available to this class.
     *
     * \todo Should we over-write the value of mDt, or change this behaviour? (see #2159)
     */
    unsigned num_time_steps = (unsigned) ((mEndTime-current_time)/mDt+0.5);
    if (current_time > 0) // use the reset function if necessary
    {
        p_simulation_time->ResetEndTimeAndNumberOfTimeSteps(mEndTime, num_time_steps);
    }
    else
    {
        if (p_simulation_time->IsEndTimeAndNumberOfTimeStepsSetUp())
        {
            EXCEPTION("End time and number of timesteps already setup. You should not use SimulationTime::SetEndTimeAndNumberOfTimeSteps in cell-based tests.");
        }
        else
        {
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(mEndTime, num_time_steps);
        }
    }

    if (mOutputDirectory == "")
    {
        EXCEPTION("OutputDirectory not set");
    }

    double time_now = p_simulation_time->GetTime();
    std::ostringstream time_string;
    time_string << time_now;

    std::string results_directory = mOutputDirectory +"/results_from_time_" + time_string.str();
    mSimulationOutputDirectory = results_directory;

    // Set up simulation

    // Create output files for the visualizer
    OutputFileHandler output_file_handler(results_directory+"/", true);

    mrCellPopulation.OpenWritersFiles(output_file_handler);

    if (mOutputDivisionLocations)
    {
        mpDivisionLocationFile = output_file_handler.OpenOutputFile("divisions.dat");
    }
    if (mOutputCellVelocities)
    {
        OutputFileHandler output_file_handler2(this->mSimulationOutputDirectory+"/", false);
        mpCellVelocitiesFile = output_file_handler2.OpenOutputFile("cellvelocities.dat");
    }

    if (PetscTools::AmMaster())
    {
        mpVizSetupFile = output_file_handler.OpenOutputFile("results.vizsetup");

        for (typename std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<ELEMENT_DIM, SPACE_DIM> > >::iterator iter = mSimulationModifiers.begin();
             iter != mSimulationModifiers.end();
             ++iter)
        {
            if (boost::dynamic_pointer_cast<AbstractPdeModifier<SPACE_DIM> >(*iter))
            {
                *this->mpVizSetupFile << "PDE \n";
            }
        }
    }

    this->mrCellPopulation.SimulationSetupHook(this);

    SetupSolve();

    // Call SetupSolve() on each modifier
    for (typename std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<ELEMENT_DIM, SPACE_DIM> > >::iterator iter = mSimulationModifiers.begin();
         iter != mSimulationModifiers.end();
         ++iter)
    {
        (*iter)->SetupSolve(this->mrCellPopulation,this->mSimulationOutputDirectory);
    }

    /*
     * Age the cells to the correct time. Note that cells are created with
     * negative birth times so that some are initially almost ready to divide.
     */
    LOG(1, "Setting up cells...");
    for (typename AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>::Iterator cell_iter = mrCellPopulation.Begin();
         cell_iter != mrCellPopulation.End();
         ++cell_iter)
    {
        /*
         * We don't use the result; this call is just to force the cells to age
         * to the current time running their cell-cycle models to get there.
         */
        cell_iter->ReadyToDivide();
    }
    LOG(1, "\tdone\n");

    // Write initial conditions to file for the visualizer
    WriteVisualizerSetupFile();

    if (PetscTools::AmMaster())
    {
        *mpVizSetupFile << std::flush;
    }

    mrCellPopulation.WriteResultsToFiles(results_directory+"/");

    OutputSimulationSetup();
    CellBasedEventHandler::EndEvent(CellBasedEventHandler::SETUP);

    // Enter main time loop
    while (!( p_simulation_time->IsFinished() || StoppingEventHasOccurred() ) )
    {
        LOG(1, "--TIME = " << p_simulation_time->GetTime() << "\n");

        // This function calls DoCellRemoval(), DoCellBirth() and CellPopulation::Update()
        UpdateCellPopulation();

        // Store whether we are sampling results at the current timestep
        SimulationTime* p_time = SimulationTime::Instance();
        bool at_sampling_timestep = (p_time->GetTimeStepsElapsed()%this->mSamplingTimestepMultiple == 0);

        /*
         * If required, store the current locations of cell centres. Note that we need to
         * use a std::map between cells and locations, rather than (say) a std::vector with
         * location indices corresponding to cells, since once we call UpdateCellLocations()
         * the location index of each cell may change. This is especially true in the case
         * of a CaBasedCellPopulation.
         */
        std::map<CellPtr, c_vector<double, SPACE_DIM> > old_cell_locations;
        if (mOutputCellVelocities && at_sampling_timestep)
        {
            for (typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator cell_iter = mrCellPopulation.Begin();
                 cell_iter != mrCellPopulation.End();
                 ++cell_iter)
            {
                old_cell_locations[*cell_iter] = mrCellPopulation.GetLocationOfCellCentre(*cell_iter);
            }
        }

        // Update cell locations and topology
        UpdateCellLocationsAndTopology();

        // Now write cell velocities to file if required
        if (mOutputCellVelocities && at_sampling_timestep)
        {
            // Offset as doing this before we increase time by mDt
            *mpCellVelocitiesFile << p_time->GetTime() + mDt<< "\t";

            for (typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator cell_iter = mrCellPopulation.Begin();
                 cell_iter != mrCellPopulation.End();
                 ++cell_iter)
            {
                unsigned index = mrCellPopulation.GetLocationIndexUsingCell(*cell_iter);
                const c_vector<double,SPACE_DIM>& position = mrCellPopulation.GetLocationOfCellCentre(*cell_iter);

                c_vector<double, SPACE_DIM> velocity; // Two lines for profile build
                velocity = (position - old_cell_locations[*cell_iter])/mDt;

                *mpCellVelocitiesFile << index  << " ";
                for (unsigned i=0; i<SPACE_DIM; i++)
                {
                    *mpCellVelocitiesFile << position[i] << " ";
                }

                for (unsigned i=0; i<SPACE_DIM; i++)
                {
                    *mpCellVelocitiesFile << velocity[i] << " ";
                }
            }
            *mpCellVelocitiesFile << "\n";
        }

        // Update the assignment of cells to processes.
        mrCellPopulation.UpdateCellProcessLocation();

        // Increment simulation time here, so results files look sensible
        p_simulation_time->IncrementTimeOneStep();

        // Call UpdateAtEndOfTimeStep() on each modifier
        CellBasedEventHandler::BeginEvent(CellBasedEventHandler::UPDATESIMULATION);
        for (typename std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<ELEMENT_DIM, SPACE_DIM> > >::iterator iter = mSimulationModifiers.begin();
             iter != mSimulationModifiers.end();
             ++iter)
        {
            (*iter)->UpdateAtEndOfTimeStep(this->mrCellPopulation);
        }
        CellBasedEventHandler::EndEvent(CellBasedEventHandler::UPDATESIMULATION);

        // Output current results to file
        CellBasedEventHandler::BeginEvent(CellBasedEventHandler::OUTPUT);
        if (p_simulation_time->GetTimeStepsElapsed()%mSamplingTimestepMultiple == 0)// should be at_sampling_timestep !
        {
            mrCellPopulation.WriteResultsToFiles(results_directory+"/");

            // Call UpdateAtEndOfOutputTimeStep() on each modifier
            for (typename std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<ELEMENT_DIM, SPACE_DIM> > >::iterator iter = mSimulationModifiers.begin();
                 iter != mSimulationModifiers.end();
                 ++iter)
            {
                (*iter)->UpdateAtEndOfOutputTimeStep(this->mrCellPopulation);
            }
        }
        CellBasedEventHandler::EndEvent(CellBasedEventHandler::OUTPUT);
    }

    LOG(1, "--END TIME = " << p_simulation_time->GetTime() << "\n");

    /*
     * Carry out a final update so that cell population is coherent with new cell positions.
     * Note that cell birth and death still need to be checked because they may be spatially
     * dependent.
     */
    UpdateCellPopulation();

    CellBasedEventHandler::BeginEvent(CellBasedEventHandler::UPDATESIMULATION);
    // Call UpdateAtEndOfSolve(), on each modifier
    for (typename std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<ELEMENT_DIM, SPACE_DIM> > >::iterator iter = mSimulationModifiers.begin();
         iter != mSimulationModifiers.end();
         ++iter)
    {
        (*iter)->UpdateAtEndOfSolve(this->mrCellPopulation);
    }
    CellBasedEventHandler::EndEvent(CellBasedEventHandler::UPDATESIMULATION);

    CellBasedEventHandler::BeginEvent(CellBasedEventHandler::OUTPUT);

    mrCellPopulation.CloseWritersFiles();

    if (mOutputDivisionLocations)
    {
        mpDivisionLocationFile->close();
    }
    if (mOutputCellVelocities)
    {
        mpCellVelocitiesFile->close();
    }

    if (PetscTools::AmMaster())
    {
        *mpVizSetupFile << "Complete\n";
        mpVizSetupFile->close();
    }

    CellBasedEventHandler::EndEvent(CellBasedEventHandler::OUTPUT);
    CellBasedEventHandler::EndEvent(CellBasedEventHandler::EVERYTHING);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>::StoppingEventHasOccurred()
{
    return false;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>::UpdateCellPopulation()
{
    // Remove dead cells
    CellBasedEventHandler::BeginEvent(CellBasedEventHandler::DEATH);
    unsigned deaths_this_step = DoCellRemoval();
    mNumDeaths += deaths_this_step;
    LOG(1, "\tNum deaths = " << mNumDeaths << "\n");
    CellBasedEventHandler::EndEvent(CellBasedEventHandler::DEATH);

    // Divide cells
    CellBasedEventHandler::BeginEvent(CellBasedEventHandler::BIRTH);
    unsigned births_this_step = DoCellBirth();
    mNumBirths += births_this_step;
    LOG(1, "\tNum births = " << mNumBirths << "\n");
    CellBasedEventHandler::EndEvent(CellBasedEventHandler::BIRTH);

    // This allows NodeBasedCellPopulation::Update() to do the minimum amount of work
    bool births_or_death_occurred = ((births_this_step>0) || (deaths_this_step>0));

    // Update topology of cell population
    CellBasedEventHandler::BeginEvent(CellBasedEventHandler::UPDATECELLPOPULATION);
    if (mUpdateCellPopulation)
    {
        LOG(1, "\tUpdating cell population...");
        mrCellPopulation.Update(births_or_death_occurred);
        LOG(1, "\tdone.\n");
    }
    else if (births_or_death_occurred)
    {
        EXCEPTION("CellPopulation has had births or deaths but mUpdateCellPopulation is set to false, please set it to true.");
    }
    CellBasedEventHandler::EndEvent(CellBasedEventHandler::UPDATECELLPOPULATION);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>::GetOutputDivisionLocations()
{
    return mOutputDivisionLocations;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>::SetOutputDivisionLocations(bool outputDivisionLocations)
{
    mOutputDivisionLocations = outputDivisionLocations;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>::GetOutputCellVelocities()
{
    return mOutputCellVelocities;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>::SetOutputCellVelocities(bool outputCellVelocities)
{
    mOutputCellVelocities = outputCellVelocities;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>::OutputSimulationSetup()
{
    OutputFileHandler output_file_handler(this->mSimulationOutputDirectory + "/", false);

    // Output machine information
    ExecutableSupport::SetOutputDirectory(output_file_handler.GetOutputDirectoryFullPath());
    ExecutableSupport::WriteMachineInfoFile("system_info");

    if (PetscTools::AmMaster())
    {
        // Output Chaste provenance information
        out_stream build_info_file = output_file_handler.OpenOutputFile("build.info");
        std::string build_info;
        ExecutableSupport::GetBuildInfo(build_info);
        *build_info_file << build_info;
        build_info_file->close();

        // Output simulation parameter and setup details
        out_stream parameter_file = output_file_handler.OpenOutputFile("results.parameters");

        // Output simulation details
        std::string simulation_type = GetIdentifier();

        *parameter_file << "<Chaste>\n";
        *parameter_file << "\n\t<" << simulation_type << ">\n";
        OutputSimulationParameters(parameter_file);
        *parameter_file << "\t</" << simulation_type << ">\n";
        *parameter_file << "\n";

        // Output cell population details (includes cell-cycle model details)
        mrCellPopulation.OutputCellPopulationInfo(parameter_file);

        // Loop over cell killers
        *parameter_file << "\n\t<CellKillers>\n";
        for (typename std::vector<boost::shared_ptr<AbstractCellKiller<SPACE_DIM> > >::iterator iter = mCellKillers.begin();
             iter != mCellKillers.end();
             ++iter)
        {
            // Output cell killer details
            (*iter)->OutputCellKillerInfo(parameter_file);
        }
        *parameter_file << "\t</CellKillers>\n";

        // Iterate over simulationmodifiers
        *parameter_file << "\n\t<SimulationModifiers>\n";
        for (typename std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<ELEMENT_DIM, SPACE_DIM> > >::iterator iter = mSimulationModifiers.begin();
             iter != mSimulationModifiers.end();
             ++iter)
        {
            // Output simulation modifier details
            (*iter)->OutputSimulationModifierInfo(parameter_file);
        }
        *parameter_file << "\t</SimulationModifiers>\n";

        // This is used to output information about subclasses
        OutputAdditionalSimulationSetup(parameter_file);

        *parameter_file << "\n</Chaste>\n";
        parameter_file->close();
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>::OutputSimulationParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t<Dt>" << mDt << "</Dt>\n";
    *rParamsFile << "\t\t<EndTime>" << mEndTime << "</EndTime>\n";
    *rParamsFile << "\t\t<SamplingTimestepMultiple>" << mSamplingTimestepMultiple << "</SamplingTimestepMultiple>\n";
    *rParamsFile << "\t\t<OutputDivisionLocations>" << mOutputDivisionLocations << "</OutputDivisionLocations>\n";
    *rParamsFile << "\t\t<OutputCellVelocities>" << mOutputCellVelocities << "</OutputCellVelocities>\n";
}

// Explicit instantiation
template class AbstractCellBasedSimulation<1,1>;
template class AbstractCellBasedSimulation<1,2>;
template class AbstractCellBasedSimulation<2,2>;
template class AbstractCellBasedSimulation<1,3>;
template class AbstractCellBasedSimulation<2,3>;
template class AbstractCellBasedSimulation<3,3>;
