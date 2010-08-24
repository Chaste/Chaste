#include "VertexBasedTissueSimulation.hpp"
#include "SimpleDataWriter.hpp"
#include "CellBasedEventHandler.hpp"
#include "CellwiseData.hpp"
#include "FunctionalitiesForVertexBasedTissue.hpp"
#include "LogFile.hpp"

//template<unsigned DIM>
VertexBasedTissueSimulation::VertexBasedTissueSimulation(AbstractTissue<2>& rTissue,
                  std::vector<AbstractForce<2>*> forceCollection,
                  bool deleteTissueAndForceCollection,
                  bool initialiseCells)
    : TissueSimulation<2>(rTissue,
                  forceCollection,
                  deleteTissueAndForceCollection,
                  initialiseCells)
//}
{
//    mpStaticCastTissue = static_cast<VertexBasedTissue<2>*>(&mrTissue);
}

//template<unsigned DIM>
VertexBasedTissueSimulation::~VertexBasedTissueSimulation()
{
}

//template<unsigned DIM>
void VertexBasedTissueSimulation::PostSolve()
{
    
    // Save results to file
    SimulationTime* p_time = SimulationTime::Instance();

    double time_next_step = p_time->GetTime() + p_time->GetTimeStep();
    FunctionalitiesForVertexBasedTissue functionality;
    double fraction = functionality.GetFractionCircumferenceAndLengthOutsideLabelledCells(this->mrTissue);
        
    std::cout << "time step: " << time_next_step << ", fraction: " << fraction << std::endl;
    
}


//
//
//
//
/////\todo Much of this code can probably merged with TissueSimulation::Solve()
//template<unsigned DIM>
//void LatticeBasedTissueSimulation<DIM>::Solve()
//{
//    CellBasedEventHandler::BeginEvent(CellBasedEventHandler::EVERYTHING);
//    CellBasedEventHandler::BeginEvent(CellBasedEventHandler::SETUP);
//
//    // Set up the simulation time
//    SimulationTime* p_simulation_time = SimulationTime::Instance();
//    double current_time = p_simulation_time->GetTime();
//
//    unsigned num_time_steps = (unsigned) ((this->mEndTime-current_time)/this->mDt+0.5);
//
//    if (current_time > 0) // use the reset function if necessary
//    {
//        p_simulation_time->ResetEndTimeAndNumberOfTimeSteps(this->mEndTime, num_time_steps);
//    }
//    else
//    {
//        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(this->mEndTime, num_time_steps);
//    }
//
//    if (this->mOutputDirectory == "")
//    {
//        EXCEPTION("OutputDirectory not set");
//    }
//
//    double time_now = p_simulation_time->GetTime();
//    std::ostringstream time_string;
//    time_string << time_now;
//
//    std::string results_directory = this->mOutputDirectory +"/results_from_time_" + time_string.str();
//    this->mSimulationOutputDirectory = results_directory;
//
//    ///////////////////////////////////////////////////////////
//    // Set up Simulation
//    ///////////////////////////////////////////////////////////
//
//    // Create output files for the visualizer
//    OutputFileHandler output_file_handler(results_directory+"/", true);
//
//    this->mrTissue.CreateOutputFiles(results_directory+"/", false);
//
//    this->mpVizSetupFile = output_file_handler.OpenOutputFile("results.vizsetup");
//
//    this->SetupSolve();
//
//    // Age the cells to the correct time. Note that cells are created with
//    // negative birth times so that some are initially almost ready to divide.
//    LOG(1, "Setting up cells...");
//    for (typename AbstractTissue<DIM>::Iterator cell_iter = this->mrTissue.Begin();
//         cell_iter != this->mrTissue.End();
//         ++cell_iter)
//    {
//        // We don't use the result; this call is just to force the cells to age
//        // to the current time running their cell cycle models to get there
//        cell_iter->ReadyToDivide();
//    }
//    LOG(1, "\tdone\n");
//
//    // Write initial conditions to file for the visualizer
//    this->WriteVisualizerSetupFile();
//
//    *this->mpVizSetupFile << std::flush;
//
//    this->mrTissue.WriteResultsToFiles();
//
//    CellBasedEventHandler::EndEvent(CellBasedEventHandler::SETUP);
//
//    /////////////////////////////////////////////////////////////////////
//    // Main time loop
//    /////////////////////////////////////////////////////////////////////
//
//    while ((p_simulation_time->GetTimeStepsElapsed() < num_time_steps) && !(this->StoppingEventHasOccurred()) )
//    {
//        LOG(1, "--TIME = " << p_simulation_time->GetTime() << "\n");
//
//        /*
//         * If mInitialiseCells is false, then the simulation has been loaded from an archive.
//         * In this case, we should not call UpdateTissue() at the first time step. This is
//         * because it will have already been called at the final time step prior to saving;
//         * if we were to call it again now, then we would have introduced an extra call to
//         * the random number generator compared to if we had not saved and loaded the simulation,
//         * thus affecting results. This would be bad - we don't want saving and loading to have
//         * any effect on the course of a simulation! See #1445.
//         */
//        bool update_tissue_this_timestep = true;
//        if (!this->mInitialiseCells && (p_simulation_time->GetTimeStepsElapsed() == 0))
//        {
//            update_tissue_this_timestep = false;
//        }
//
//        if (update_tissue_this_timestep)
//        {
//            /**
//             * This function calls:
//             * DoCellRemoval()
//             * DoCellBirth()
//             * Tissue::Update()
//             */
//            this->UpdateTissue();
//        }
//
//        //////////////////////////////////////////////////////
//        // Calculate Cell Movements and Update Cell Locations
//        //////////////////////////////////////////////////////
//        CellBasedEventHandler::BeginEvent(CellBasedEventHandler::POSITION);
//        UpdateCellLocations();
//        CellBasedEventHandler::EndEvent(CellBasedEventHandler::POSITION);
//
//        //////////////////////////////////////////
//        // PostSolve, which may be implemented by
//        // child classes (eg to solve PDEs)
//        //////////////////////////////////////////
//        this->PostSolve();
//
//        // Increment simulation time here, so results files look sensible
//        p_simulation_time->IncrementTimeOneStep();
//
//        ////////////////////////////
//        // Output current results
//        ////////////////////////////
//        CellBasedEventHandler::BeginEvent(CellBasedEventHandler::OUTPUT);
//
//        // Write results to file
//        if (p_simulation_time->GetTimeStepsElapsed()%this->mSamplingTimestepMultiple==0)
//        {
//            this->mrTissue.WriteResultsToFiles();
//        }
//
//        CellBasedEventHandler::EndEvent(CellBasedEventHandler::OUTPUT);
//    }
//
//    LOG(1, "--END TIME = " << SimulationTime::Instance()->GetTime() << "\n");
//
//    /*
//     * Carry out a final update so that tissue is coherent with new cell positions.
//     * NB cell birth/death still need to be checked because they may be spatially-dependent.
//     */
//    this->UpdateTissue();
//
//    this->AfterSolve();
//
//    CellBasedEventHandler::BeginEvent(CellBasedEventHandler::OUTPUT);
//
//    this->mrTissue.CloseOutputFiles();
//
//    *this->mpVizSetupFile << "Complete\n";
//    this->mpVizSetupFile->close();
//
//    CellBasedEventHandler::EndEvent(CellBasedEventHandler::OUTPUT);
//
//    CellBasedEventHandler::EndEvent(CellBasedEventHandler::EVERYTHING);
//}
//
//template<unsigned DIM>
//const std::vector<AbstractUpdateRule<DIM>*> VertexBasedTissueSimulation<DIM>::rGetUpdateRuleCollection() const
//{
//    return mUpdateRuleCollection;
//}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

//template class VertexBasedTissueSimulation<1>;
//template class VertexBasedTissueSimulation<2>;
//template class VertexBasedTissueSimulation<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(VertexBasedTissueSimulation)
