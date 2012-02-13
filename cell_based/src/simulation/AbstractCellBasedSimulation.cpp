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

#include <cmath>
#include <ctime>
#include <iostream>
#include <fstream>
#include <set>

#include "AbstractCellBasedSimulation.hpp"
#include "CellBasedEventHandler.hpp"
#include "LogFile.hpp"
#include "Version.hpp"
#include "ExecutableSupport.hpp"
#include "Exception.hpp"

#include <typeinfo>

template<unsigned DIM>
AbstractCellBasedSimulation<DIM>::AbstractCellBasedSimulation(AbstractCellPopulation<DIM>& rCellPopulation,
                                              bool deleteCellPopulationInDestructor,
                                              bool initialiseCells)
    : mDt(DOUBLE_UNSET),
      mEndTime(0.0),  // hours - this is set later on
      mrCellPopulation(rCellPopulation),
      mDeleteCellPopulationInDestructor(deleteCellPopulationInDestructor),
      mInitialiseCells(initialiseCells),
      mNoBirth(false),
      mUpdateCellPopulation(true),
      mOutputDirectory(""),
      mSimulationOutputDirectory(mOutputDirectory),
      mNumBirths(0),
      mNumDeaths(0),
      mSamplingTimestepMultiple(1),
      mpCellBasedPdeHandler(NULL)
{
    // Set a random seed of 0 if it wasn't specified earlier
    RandomNumberGenerator::Instance();

    if (mInitialiseCells)
    {
        mrCellPopulation.InitialiseCells();
    }
}

template<unsigned DIM>
AbstractCellBasedSimulation<DIM>::~AbstractCellBasedSimulation()
{
    if (mDeleteCellPopulationInDestructor)
    {
        delete &mrCellPopulation;
        delete mpCellBasedPdeHandler;
    }
}

template<unsigned DIM>
void AbstractCellBasedSimulation<DIM>::SetCellBasedPdeHandler(CellBasedPdeHandler<DIM>* pCellBasedPdeHandler)
{
    mpCellBasedPdeHandler = pCellBasedPdeHandler;
}

template<unsigned DIM>
CellBasedPdeHandler<DIM>* AbstractCellBasedSimulation<DIM>::GetCellBasedPdeHandler()
{
    return mpCellBasedPdeHandler;
}

template<unsigned DIM>
unsigned AbstractCellBasedSimulation<DIM>::DoCellBirth()
{
    if (mNoBirth)
    {
        return 0;
    }

    unsigned num_births_this_step = 0;

    // Iterate over all cells, seeing if each one can be divided
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = mrCellPopulation.Begin();
         cell_iter != mrCellPopulation.End();
         ++cell_iter)
    {
        // Check if this cell is ready to divide - if so create a new cell etc.
        if (cell_iter->GetAge() > 0.0)
        {
            if (cell_iter->ReadyToDivide())
            {
                try
                {
                    // Create a new cell
                    CellPtr p_new_cell = cell_iter->Divide();

                    // Call method that determines how cell division occurs and returns a vector
                    c_vector<double, DIM> new_location = CalculateCellDivisionVector(*cell_iter);

                    // Add new cell to the cell population
                    mrCellPopulation.AddCell(p_new_cell, new_location, *cell_iter);

                    // Update counter
                    num_births_this_step++;
                }
                catch (Exception& e)
                {
                    // Don't do anything
                }
            }
        }
    }
    return num_births_this_step;
}

template<unsigned DIM>
unsigned AbstractCellBasedSimulation<DIM>::DoCellRemoval()
{
    unsigned num_deaths_this_step = 0;

    /*
     * This labels cells as dead or apoptosing. It does not actually remove the cells,
     * mrCellPopulation.RemoveDeadCells() needs to be called for this.
     */
    for (typename std::vector<boost::shared_ptr<AbstractCellKiller<DIM> > >::iterator killer_iter = mCellKillers.begin();
         killer_iter != mCellKillers.end();
         ++killer_iter)
    {
        (*killer_iter)->TestAndLabelCellsForApoptosisOrDeath();
    }

    num_deaths_this_step += mrCellPopulation.RemoveDeadCells();

    return num_deaths_this_step;
}

template<unsigned DIM>
void AbstractCellBasedSimulation<DIM>::SetDt(double dt)
{
    assert(dt > 0);
    mDt = dt;
}

template<unsigned DIM>
double AbstractCellBasedSimulation<DIM>::GetDt()
{
    return mDt;
}

template<unsigned DIM>
unsigned AbstractCellBasedSimulation<DIM>::GetNumBirths()
{
    return mNumBirths;
}

template<unsigned DIM>
unsigned AbstractCellBasedSimulation<DIM>::GetNumDeaths()
{
    return mNumDeaths;
}

template<unsigned DIM>
void AbstractCellBasedSimulation<DIM>::SetEndTime(double endTime)
{
    assert(endTime > 0);
    mEndTime = endTime;
}

template<unsigned DIM>
void AbstractCellBasedSimulation<DIM>::SetOutputDirectory(std::string outputDirectory)
{
    mOutputDirectory = outputDirectory;
    mSimulationOutputDirectory = mOutputDirectory;
}

template<unsigned DIM>
std::string AbstractCellBasedSimulation<DIM>::GetOutputDirectory()
{
    return mOutputDirectory;
}

template<unsigned DIM>
void AbstractCellBasedSimulation<DIM>::SetSamplingTimestepMultiple(unsigned samplingTimestepMultiple)
{
    assert(samplingTimestepMultiple > 0);
    mSamplingTimestepMultiple = samplingTimestepMultiple;
}

template<unsigned DIM>
AbstractCellPopulation<DIM>& AbstractCellBasedSimulation<DIM>::rGetCellPopulation()
{
    return mrCellPopulation;
}

template<unsigned DIM>
const AbstractCellPopulation<DIM>& AbstractCellBasedSimulation<DIM>::rGetCellPopulation() const
{
    return mrCellPopulation;
}

template<unsigned DIM>
void AbstractCellBasedSimulation<DIM>::SetUpdateCellPopulationRule(bool updateCellPopulation)
{
    mUpdateCellPopulation = updateCellPopulation;
}

template<unsigned DIM>
void AbstractCellBasedSimulation<DIM>::SetNoBirth(bool noBirth)
{
    mNoBirth = noBirth;
}

template<unsigned DIM>
void AbstractCellBasedSimulation<DIM>::AddCellKiller(boost::shared_ptr<AbstractCellKiller<DIM> > pCellKiller)
{
    mCellKillers.push_back(pCellKiller);
}

template<unsigned DIM>
std::vector<double> AbstractCellBasedSimulation<DIM>::GetNodeLocation(const unsigned& rNodeIndex)
{
    std::vector<double> location;
    for (unsigned i=0; i<DIM; i++)
    {
        location.push_back(mrCellPopulation.GetNode(rNodeIndex)->rGetLocation()[i]);
    }
    return location;
}

template<unsigned DIM>
void AbstractCellBasedSimulation<DIM>::Solve()
{
    CellBasedEventHandler::BeginEvent(CellBasedEventHandler::EVERYTHING);
    CellBasedEventHandler::BeginEvent(CellBasedEventHandler::SETUP);

    // Set up the simulation time
    SimulationTime* p_simulation_time = SimulationTime::Instance();
    double current_time = p_simulation_time->GetTime();

    unsigned num_time_steps = (unsigned) ((mEndTime-current_time)/mDt+0.5);

    if (current_time > 0) // use the reset function if necessary
    {
        p_simulation_time->ResetEndTimeAndNumberOfTimeSteps(mEndTime, num_time_steps);
    }
    else
    {
        if(p_simulation_time->IsEndTimeAndNumberOfTimeStepsSetUp())
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

    mrCellPopulation.CreateOutputFiles(results_directory+"/", false);

    mpVizSetupFile = output_file_handler.OpenOutputFile("results.vizsetup");

    // If any PDEs have been defined, set up results files to store their solution
    if (mpCellBasedPdeHandler != NULL)
    {
        mpCellBasedPdeHandler->OpenResultsFiles(this->mSimulationOutputDirectory);
        *this->mpVizSetupFile << "PDE \n";
    }

    SetupSolve();

    // Age the cells to the correct time. Note that cells are created with
    // negative birth times so that some are initially almost ready to divide.
    LOG(1, "Setting up cells...");
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = mrCellPopulation.Begin();
         cell_iter != mrCellPopulation.End();
         ++cell_iter)
    {
        // We don't use the result; this call is just to force the cells to age
        // to the current time running their cell-cycle models to get there
        cell_iter->ReadyToDivide();
    }
    LOG(1, "\tdone\n");

    // Write initial conditions to file for the visualizer
    WriteVisualizerSetupFile();

    *mpVizSetupFile << std::flush;

    mrCellPopulation.WriteResultsToFiles();

    OutputSimulationSetup();

    CellBasedEventHandler::EndEvent(CellBasedEventHandler::SETUP);

    // Enter main time loop
    while (!( p_simulation_time->IsFinished() || StoppingEventHasOccurred() ) )
    {
        LOG(1, "--TIME = " << p_simulation_time->GetTime() << "\n");

        // This function calls DoCellRemoval(), DoCellBirth() and CellPopulation::Update()
        UpdateCellPopulation();

        // Update cell locations and topology
        UpdateCellLocationsAndTopology();

        // If any PDEs have been defined, solve them and store their solution in results files
        if (mpCellBasedPdeHandler != NULL)
        {
            CellBasedEventHandler::BeginEvent(CellBasedEventHandler::PDE);
            mpCellBasedPdeHandler->SolvePdeAndWriteResultsToFile(this->mSamplingTimestepMultiple);
            CellBasedEventHandler::EndEvent(CellBasedEventHandler::PDE);
        }

        // Call UpdateAtEndOfTimeStep(), which may be implemented by child classes
        UpdateAtEndOfTimeStep();

        // Increment simulation time here, so results files look sensible
        p_simulation_time->IncrementTimeOneStep();

        // Output current results to file
        CellBasedEventHandler::BeginEvent(CellBasedEventHandler::OUTPUT);
        if (p_simulation_time->GetTimeStepsElapsed()%mSamplingTimestepMultiple == 0)
        {
            mrCellPopulation.WriteResultsToFiles();
        }
        CellBasedEventHandler::EndEvent(CellBasedEventHandler::OUTPUT);
    }

    LOG(1, "--END TIME = " << p_simulation_time->GetTime() << "\n");

    // Carry out a final update so that cell population is coherent with new cell positions.
    // NB cell birth/death still need to be checked because they may be spatially-dependent.
    UpdateCellPopulation();

    // If any PDEs have been defined, close the results files storing their solution
    if (mpCellBasedPdeHandler != NULL)
    {
        mpCellBasedPdeHandler->CloseResultsFiles();
    }

    UpdateAtEndOfSolve();

    CellBasedEventHandler::BeginEvent(CellBasedEventHandler::OUTPUT);

    mrCellPopulation.CloseOutputFiles();

    *mpVizSetupFile << "Complete\n";
    mpVizSetupFile->close();

    CellBasedEventHandler::EndEvent(CellBasedEventHandler::OUTPUT);
    CellBasedEventHandler::EndEvent(CellBasedEventHandler::EVERYTHING);
}

template<unsigned DIM>
bool AbstractCellBasedSimulation<DIM>::StoppingEventHasOccurred()
{
    return false;
}

template<unsigned DIM>
void AbstractCellBasedSimulation<DIM>::UpdateCellPopulation()
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

template<unsigned DIM>
void AbstractCellBasedSimulation<DIM>::OutputSimulationSetup()
{
    OutputFileHandler output_file_handler(this->mSimulationOutputDirectory + "/", false);

    // Output machine information
    ExecutableSupport::WriteMachineInfoFile(this->mSimulationOutputDirectory + "/system_info");

    // Output Chaste provenance information
    out_stream build_info_file = output_file_handler.OpenOutputFile("build.info");
    ExecutableSupport::WriteLibraryInfo(build_info_file);
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
    for (typename std::vector<boost::shared_ptr<AbstractCellKiller<DIM> > >::iterator iter = mCellKillers.begin();
         iter != mCellKillers.end();
         ++iter)
    {
        // Output cell killer details
        (*iter)->OutputCellKillerInfo(parameter_file);
    }
    *parameter_file << "\t</CellKillers>\n";

    // This is used to output information about subclasses
    OutputAdditionalSimulationSetup(parameter_file);

    *parameter_file << "\n</Chaste>\n";
    parameter_file->close();
}

template<unsigned DIM>
void AbstractCellBasedSimulation<DIM>::OutputSimulationParameters(out_stream& rParamsFile)
{
    if (mpCellBasedPdeHandler != NULL)
    {
        mpCellBasedPdeHandler->OutputParameters(rParamsFile);
    }

    *rParamsFile << "\t\t<Dt>" << mDt << "</Dt>\n";
    *rParamsFile << "\t\t<EndTime>" << mEndTime << "</EndTime>\n";
    *rParamsFile << "\t\t<SamplingTimestepMultiple>" << mSamplingTimestepMultiple << "</SamplingTimestepMultiple>\n";
}

////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
////////////////////////////////////////////////////////////////////////////

template class AbstractCellBasedSimulation<1>;
template class AbstractCellBasedSimulation<2>;
template class AbstractCellBasedSimulation<3>;
