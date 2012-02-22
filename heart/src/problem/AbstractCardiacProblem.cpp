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

#include "AbstractCardiacProblem.hpp"

#include "GenericMeshReader.hpp"
#include "Exception.hpp"
#include "HeartConfig.hpp"
#include "HeartEventHandler.hpp"
#include "TimeStepper.hpp"
#include "PetscTools.hpp"
#include "DistributedVector.hpp"
#include "ProgressReporter.hpp"
#include "LinearSystem.hpp"
#include "PostProcessingWriter.hpp"
#include "Hdf5ToMeshalyzerConverter.hpp"
#include "Hdf5ToCmguiConverter.hpp"
#include "Hdf5ToVtkConverter.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
AbstractCardiacProblem<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::AbstractCardiacProblem(
            AbstractCardiacCellFactory<ELEMENT_DIM,SPACE_DIM>* pCellFactory)
    : mMeshFilename(""), // i.e. undefined
      mAllocatedMemoryForMesh(false),
      mWriteInfo(false),
      mPrintOutput(true),
      mpCardiacTissue(NULL),
      mpSolver(NULL),
      mpCellFactory(pCellFactory),
      mpMesh(NULL),
      mSolution(NULL),
      mCurrentTime(0.0),
      mpTimeAdaptivityController(NULL),
      mpWriter(NULL)
{
    assert(mNodesToOutput.empty());
    if (!mpCellFactory)
    {
        EXCEPTION("AbstractCardiacProblem: Please supply a cell factory pointer to your cardiac problem constructor.");
    }
    HeartEventHandler::BeginEvent(HeartEventHandler::EVERYTHING);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
AbstractCardiacProblem<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::AbstractCardiacProblem()
    // It doesn't really matter what we initialise these to, as they'll be overwritten by
    // the serialization methods
    : mMeshFilename(""),
      mAllocatedMemoryForMesh(false), // Handled by AbstractCardiacTissue
      mWriteInfo(false),
      mPrintOutput(true),
      mVoltageColumnId(UINT_MAX),
      mTimeColumnId(UINT_MAX),
      mNodeColumnId(UINT_MAX),
      mpCardiacTissue(NULL),
      mpSolver(NULL),
      mpCellFactory(NULL),
      mpMesh(NULL),
      mSolution(NULL),
      mCurrentTime(0.0),
      mpTimeAdaptivityController(NULL),
      mpWriter(NULL)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
AbstractCardiacProblem<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::~AbstractCardiacProblem()
{
    delete mpCardiacTissue;
    if (mSolution)
    {
        PetscTools::Destroy(mSolution);
    }

    if (mAllocatedMemoryForMesh)
    {
        delete mpMesh;
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::Initialise()
{
    HeartEventHandler::BeginEvent(HeartEventHandler::READ_MESH);
    if (mpMesh==NULL)
    {
        // If no mesh has been passed, we get it from the configuration file
        try
        {
            if (HeartConfig::Instance()->GetLoadMesh())
            {
                CreateMeshFromHeartConfig();
                std::auto_ptr<AbstractMeshReader<ELEMENT_DIM, SPACE_DIM> > p_mesh_reader
                    = GenericMeshReader<ELEMENT_DIM, SPACE_DIM>(HeartConfig::Instance()->GetMeshName());
                mpMesh->ConstructFromMeshReader(*p_mesh_reader);
            }
            else if (HeartConfig::Instance()->GetCreateMesh())
            {
                CreateMeshFromHeartConfig();
                assert(HeartConfig::Instance()->GetSpaceDimension()==SPACE_DIM);
                double inter_node_space = HeartConfig::Instance()->GetInterNodeSpace();

                switch(HeartConfig::Instance()->GetSpaceDimension())
                {
                    case 1:
                    {
                        c_vector<double, 1> fibre_length;
                        HeartConfig::Instance()->GetFibreLength(fibre_length);
                        mpMesh->ConstructRegularSlabMesh(inter_node_space, fibre_length[0]);
                        break;
                    }
                    case 2:
                    {
                        c_vector<double, 2> sheet_dimensions; //cm
                        HeartConfig::Instance()->GetSheetDimensions(sheet_dimensions);
                        mpMesh->ConstructRegularSlabMesh(inter_node_space, sheet_dimensions[0], sheet_dimensions[1]);
                        break;
                    }
                    case 3:
                    {
                        c_vector<double, 3> slab_dimensions; //cm
                        HeartConfig::Instance()->GetSlabDimensions(slab_dimensions);
                        mpMesh->ConstructRegularSlabMesh(inter_node_space, slab_dimensions[0], slab_dimensions[1], slab_dimensions[2]);
                        break;
                    }
                    default:
                        NEVER_REACHED;
                }
            }
            else
            {
                NEVER_REACHED;
            }

            mAllocatedMemoryForMesh = true;
        }
        catch (Exception& e)
        {
            EXCEPTION(std::string("No mesh given: define it in XML parameters file or call SetMesh()\n") + e.GetShortMessage());
        }
    }
    mpCellFactory->SetMesh( mpMesh );
    HeartEventHandler::EndEvent(HeartEventHandler::READ_MESH);

    HeartEventHandler::BeginEvent(HeartEventHandler::INITIALISE);

    // If the user requested transmural stuff, we fill in the mCellHeterogeneityAreas here
    if (HeartConfig::Instance()->AreCellularTransmuralHeterogeneitiesRequested())
    {
        mpCellFactory->FillInCellularTransmuralAreas();
    }

    delete mpCardiacTissue; // In case we're called twice
    mpCardiacTissue = CreateCardiacTissue();

    HeartEventHandler::EndEvent(HeartEventHandler::INITIALISE);

    // Delete any previous solution, so we get a fresh initial condition
    if (mSolution)
    {
        HeartEventHandler::BeginEvent(HeartEventHandler::COMMUNICATION);
        PetscTools::Destroy(mSolution);
        mSolution = NULL;
        HeartEventHandler::EndEvent(HeartEventHandler::COMMUNICATION);
    }

    // Always start at time zero
    mCurrentTime = 0.0;

    // For Bidomain with bath, this is where we set up the electrodes

    SetElectrodes();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::CreateMeshFromHeartConfig()
{
    mpMesh = new DistributedTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>(HeartConfig::Instance()->GetMeshPartitioning());
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetBoundaryConditionsContainer(boost::shared_ptr<BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM> > pBcc)
{
    this->mpBoundaryConditionsContainer = pBcc;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::PreSolveChecks()
{
    if ( mpCardiacTissue == NULL ) // if tissue is NULL, Initialise() probably hasn't been called
    {
        EXCEPTION("Cardiac tissue is null, Initialise() probably hasn't been called");
    }
    if ( HeartConfig::Instance()->GetSimulationDuration() <= mCurrentTime)
    {
        EXCEPTION("End time should be in the future");
    }
    if (mPrintOutput)
    {
        if( (HeartConfig::Instance()->GetOutputDirectory()=="") || (HeartConfig::Instance()->GetOutputFilenamePrefix()==""))
        {
            EXCEPTION("Either explicitly specify not to print output (call PrintOutput(false)) or specify the output directory and filename prefix");
        }
    }

    double end_time = HeartConfig::Instance()->GetSimulationDuration();
    double pde_time = HeartConfig::Instance()->GetPdeTimeStep();

    /*
     * MatrixIsConstant stuff requires CONSTANT dt - do some checks to make sure
     * the TimeStepper won't find non-constant dt.
     * Note: printing_time does not have to divide end_time, but dt must divide
     * printing_time and end_time.
     * HeartConfig checks pde_dt divides printing dt.
     */
    ///\todo remove magic number? (#1884)
    if( fabs(end_time - pde_time*round(end_time/pde_time)) > 1e-10 )
    {
        EXCEPTION("PDE timestep does not seem to divide end time - check parameters");
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
Vec AbstractCardiacProblem<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::CreateInitialCondition()
{
    DistributedVectorFactory* p_factory = mpMesh->GetDistributedVectorFactory();
    Vec initial_condition = p_factory->CreateVec(PROBLEM_DIM);
    DistributedVector ic = p_factory->CreateDistributedVector(initial_condition);
    std::vector<DistributedVector::Stripe> stripe;
    stripe.reserve(PROBLEM_DIM);

    for (unsigned i=0; i<PROBLEM_DIM; i++)
    {
        stripe.push_back(DistributedVector::Stripe(ic, i));
    }

    for (DistributedVector::Iterator index = ic.Begin();
         index != ic.End();
         ++index)
    {
        stripe[0][index] = mpCardiacTissue->GetCardiacCell(index.Global)->GetVoltage();
        if (PROBLEM_DIM==2)
        {
            stripe[1][index] = 0;
        }
    }

    ic.Restore();

    return initial_condition;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetMesh(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh)
{
    /*
     * If this fails the mesh has already been set. We assert rather throw
     * an exception to avoid a memory leak when checking it throws correctly.
     */
    assert(mpMesh == NULL);
    assert(pMesh != NULL);
    mAllocatedMemoryForMesh = false;
    mpMesh = pMesh;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::PrintOutput(bool printOutput)
{
    mPrintOutput = printOutput;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetWriteInfo(bool writeInfo)
{
    mWriteInfo = writeInfo;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
Vec AbstractCardiacProblem<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GetSolution()
{
    return mSolution;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
DistributedVector AbstractCardiacProblem<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GetSolutionDistributedVector()
{
    return mpMesh->GetDistributedVectorFactory()->CreateDistributedVector(mSolution);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
double AbstractCardiacProblem<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GetCurrentTime()
{
    return mCurrentTime;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM> & AbstractCardiacProblem<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::rGetMesh()
{
    assert (mpMesh);
    return *mpMesh;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
AbstractCardiacTissue<ELEMENT_DIM,SPACE_DIM>* AbstractCardiacProblem<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GetTissue()
{
    return mpCardiacTissue;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetUseTimeAdaptivityController(
        bool useAdaptivity,
        AbstractTimeAdaptivityController* pController)
{
    if (useAdaptivity)
    {
        assert(pController);
        mpTimeAdaptivityController = pController;
    }
    else
    {
        mpTimeAdaptivityController = NULL;
    }
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::Solve()
{
    PreSolveChecks();

    if (!mpBoundaryConditionsContainer) // the user didn't supply a bcc
    {
        // Set up the default bcc
        mpDefaultBoundaryConditionsContainer.reset(new BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>);
        for (unsigned problem_index=0; problem_index<PROBLEM_DIM; problem_index++)
        {
            mpDefaultBoundaryConditionsContainer->DefineZeroNeumannOnMeshBoundary(mpMesh, problem_index);
        }
        mpBoundaryConditionsContainer = mpDefaultBoundaryConditionsContainer;
    }

    assert(mpSolver==NULL);
    mpSolver = CreateSolver(); // passes mpBoundaryConditionsContainer to solver

    // If we have already run a simulation, use the old solution as initial condition
    Vec initial_condition;
    if (mSolution)
    {
        initial_condition = mSolution;
    }
    else
    {
        initial_condition = CreateInitialCondition();
    }

    std::vector<double> additional_stopping_times;
    SetUpAdditionalStoppingTimes(additional_stopping_times);

    TimeStepper stepper(mCurrentTime,
                        HeartConfig::Instance()->GetSimulationDuration(),
                        HeartConfig::Instance()->GetPrintingTimeStep(),
                        false,
                        additional_stopping_times);

    std::string progress_reporter_dir;

    if (mPrintOutput)
    {
        HeartEventHandler::BeginEvent(HeartEventHandler::WRITE_OUTPUT);
        try
        {
            InitialiseWriter();
        }
        catch (Exception& e)
        {
            delete mpWriter;
            mpWriter = NULL;
            delete mpSolver;
            if (mSolution != initial_condition)
            {
                /*
                 * A PETSc Vec is a pointer, so we *don't* need to free the memory if it is
                 * freed somewhere else (e.g. in the destructor). If this is a resumed solution
                 * we set initial_condition = mSolution earlier. mSolution is going to be
                 * cleaned up in the constructor. So, only PetscTools::Destroy( initial_condition when
                 * it is not equal to mSolution (see #1695).
                 */
                PetscTools::Destroy(initial_condition);
            }
            throw e;
        }

        /*
         * If we are resuming a simulation (i.e. mSolution already exists) there's no need
         * to write the initial timestep, since it was already written as the last timestep
         * of the previous run.
         */
        if (!mSolution)
        {
            WriteOneStep(stepper.GetTime(), initial_condition);
            mpWriter->AdvanceAlongUnlimitedDimension();
        }
        HeartEventHandler::EndEvent(HeartEventHandler::WRITE_OUTPUT);

        progress_reporter_dir = HeartConfig::Instance()->GetOutputDirectory();
    }
    else
    {
        progress_reporter_dir = ""; // progress printed to CHASTE_TEST_OUTPUT
    }

    /*
     * Create a progress reporter so users can track how much has gone and
     * estimate how much time is left. Note this has to be done after the
     * InitialiseWriter above (if mPrintOutput==true).
     */
    ProgressReporter progress_reporter(progress_reporter_dir,
                                       mCurrentTime,
                                       HeartConfig::Instance()->GetSimulationDuration());
    progress_reporter.Update(mCurrentTime);

    mpSolver->SetTimeStep(HeartConfig::Instance()->GetPdeTimeStep());
    if (mpTimeAdaptivityController)
    {
        mpSolver->SetTimeAdaptivityController(mpTimeAdaptivityController);
    }

    while (!stepper.IsTimeAtEnd())
    {
        // Solve from now up to the next printing time
        mpSolver->SetTimes(stepper.GetTime(), stepper.GetNextTime());
        mpSolver->SetInitialCondition( initial_condition );

        AtBeginningOfTimestep(stepper.GetTime());

        try
        {
            try
            {
                mSolution = mpSolver->Solve();
            }
            catch (const Exception &e)
            {
#ifndef NDEBUG
                PetscTools::ReplicateException(true);
#endif
                throw e;
            }
#ifndef NDEBUG
            PetscTools::ReplicateException(false);
#endif
        }
        catch (const Exception& e)
        {
            // Free memory
            delete mpSolver;
            mpSolver = NULL;
            if (initial_condition != mSolution)
            {
                /*
                 * A PETSc Vec is a pointer, so we *don't* need to free the memory if it is
                 * freed somewhere else (e.g. in the destructor). Later, in this while loop
                 * we will set initial_condition = mSolution (or, if this is a resumed solution
                 * it may also have been done when initial_condition was created). mSolution
                 * is going to be cleaned up in the destructor. So, only PetscTools::Destroy(
                 * initial_condition when it is not equal to mSolution (see #1695).
                 */
                PetscTools::Destroy(initial_condition);
            }

            // Re-throw
            HeartEventHandler::Reset();
            CloseFilesAndPostProcess();
            throw e;
        }

        // Free old initial condition
        HeartEventHandler::BeginEvent(HeartEventHandler::COMMUNICATION);
        PetscTools::Destroy(initial_condition);
        HeartEventHandler::EndEvent(HeartEventHandler::COMMUNICATION);

        // Initial condition for next loop is current solution
        initial_condition = mSolution;

        // Update the current time
        stepper.AdvanceOneTimeStep();
        mCurrentTime = stepper.GetTime();

        // Print out details at current time if asked for
        if (mWriteInfo)
        {
            HeartEventHandler::BeginEvent(HeartEventHandler::WRITE_OUTPUT);
            WriteInfo(stepper.GetTime());
            HeartEventHandler::EndEvent(HeartEventHandler::WRITE_OUTPUT);
        }

        if (mPrintOutput)
        {
            // Writing data out to the file <FilenamePrefix>.dat
            HeartEventHandler::BeginEvent(HeartEventHandler::WRITE_OUTPUT);
            WriteOneStep(stepper.GetTime(), mSolution);
            // Just flags that we've finished a time-step; won't actually 'extend' unless new data is written.
            mpWriter->AdvanceAlongUnlimitedDimension();

            HeartEventHandler::EndEvent(HeartEventHandler::WRITE_OUTPUT);
        }

        progress_reporter.Update(stepper.GetTime());

        OnEndOfTimestep(stepper.GetTime());
    }

    // Free solver
    delete mpSolver;
    mpSolver = NULL;

    // Close the file that stores voltage values
    progress_reporter.PrintFinalising();
    CloseFilesAndPostProcess();
    HeartEventHandler::EndEvent(HeartEventHandler::EVERYTHING);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::CloseFilesAndPostProcess()
{
    // Close files
    if (!mPrintOutput)
    {
        // Nothing to do
        return;
    }
    delete mpWriter;
    mpWriter = NULL;

    HeartEventHandler::BeginEvent(HeartEventHandler::DATA_CONVERSION);
    // Only if results files were written and we are outputting all nodes
    if (mNodesToOutput.empty())
    {
        if (HeartConfig::Instance()->GetVisualizeWithMeshalyzer())
        {
            // Convert simulation data to Meshalyzer format
            Hdf5ToMeshalyzerConverter<ELEMENT_DIM,SPACE_DIM> converter(HeartConfig::Instance()->GetOutputDirectory(),
                    HeartConfig::Instance()->GetOutputFilenamePrefix(), mpMesh, HeartConfig::Instance()->GetOutputUsingOriginalNodeOrdering());
            std::string subdirectory_name = converter.GetSubdirectory();
            HeartConfig::Instance()->Write(false, subdirectory_name);
        }

        if (HeartConfig::Instance()->GetVisualizeWithCmgui())
        {
            // Convert simulation data to Cmgui format
            Hdf5ToCmguiConverter<ELEMENT_DIM,SPACE_DIM> converter(HeartConfig::Instance()->GetOutputDirectory(),
                    HeartConfig::Instance()->GetOutputFilenamePrefix(), mpMesh, GetHasBath());
            std::string subdirectory_name = converter.GetSubdirectory();
            HeartConfig::Instance()->Write(false, subdirectory_name);
        }

        if (HeartConfig::Instance()->GetVisualizeWithVtk())
        {
            // Convert simulation data to VTK format
            Hdf5ToVtkConverter<ELEMENT_DIM,SPACE_DIM> converter(HeartConfig::Instance()->GetOutputDirectory(),
                    HeartConfig::Instance()->GetOutputFilenamePrefix(), mpMesh, false, HeartConfig::Instance()->GetOutputUsingOriginalNodeOrdering());
            std::string subdirectory_name = converter.GetSubdirectory();
            HeartConfig::Instance()->Write(false, subdirectory_name);
        }

        if (HeartConfig::Instance()->GetVisualizeWithParallelVtk())
        {
            // Convert simulation data to parallel VTK (pvtu) format
            Hdf5ToVtkConverter<ELEMENT_DIM,SPACE_DIM> converter(HeartConfig::Instance()->GetOutputDirectory(),
                    HeartConfig::Instance()->GetOutputFilenamePrefix(), mpMesh, true, HeartConfig::Instance()->GetOutputUsingOriginalNodeOrdering());
            std::string subdirectory_name = converter.GetSubdirectory();
            HeartConfig::Instance()->Write(false, subdirectory_name);
        }
    }
    HeartEventHandler::EndEvent(HeartEventHandler::DATA_CONVERSION);

    HeartEventHandler::BeginEvent(HeartEventHandler::POST_PROC);
    if(HeartConfig::Instance()->IsPostProcessingRequested())
    {
        PostProcessingWriter<ELEMENT_DIM, SPACE_DIM> post_writer(*mpMesh, HeartConfig::Instance()->GetOutputDirectory(),
                        HeartConfig::Instance()->GetOutputFilenamePrefix(), true);
        post_writer.WritePostProcessingFiles();
    }

    HeartEventHandler::EndEvent(HeartEventHandler::POST_PROC);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::DefineWriterColumns(bool extending)
{
    if (!extending)
    {
        if ( mNodesToOutput.empty() )
        {
            //Set writer to output all nodes
            mpWriter->DefineFixedDimension( mpMesh->GetNumNodes() );
        }
        else
        {
            //Output only the nodes indicted
            mpWriter->DefineFixedDimension( mNodesToOutput, mpMesh->GetNumNodes() );
        }
        //mNodeColumnId = mpWriter->DefineVariable("Node", "dimensionless");
        mVoltageColumnId = mpWriter->DefineVariable("V","mV");

        // Only used to get an estimate of the # of timesteps below
        TimeStepper stepper(mCurrentTime,
                            HeartConfig::Instance()->GetSimulationDuration(),
                            HeartConfig::Instance()->GetPrintingTimeStep());

        mpWriter->DefineUnlimitedDimension("Time","msecs", stepper.EstimateTimeSteps()+1); // plus one for start and end points
    }
    else
    {
        mVoltageColumnId = mpWriter->GetVariableByName("V");
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::DefineExtraVariablesWriterColumns(bool extending)
{
    mExtraVariablesId.clear();
    // Check if any extra output variables have been requested
    if (HeartConfig::Instance()->GetOutputVariablesProvided())
    {
        // Get their names in a vector
        std::vector<std::string> output_variables;
        HeartConfig::Instance()->GetOutputVariables(output_variables);
        const unsigned num_vars = output_variables.size();
        mExtraVariablesId.reserve(num_vars);

        // Loop over them
        for (unsigned var_index=0; var_index<num_vars; var_index++)
        {
            // Get variable name
            std::string var_name = output_variables[var_index];

            // Register it (or look it up) in the data writer
            unsigned column_id;
            if (extending)
            {
                column_id = this->mpWriter->GetVariableByName(var_name);
            }
            else
            {
                column_id = this->mpWriter->DefineVariable(var_name, "");
            }

            // Store column id
            mExtraVariablesId.push_back(column_id);
        }
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::WriteExtraVariablesOneStep()
{
    // Get the variable names in a vector
    std::vector<std::string> output_variables;
    unsigned num_vars = mExtraVariablesId.size();
    if (num_vars > 0)
    {
        HeartConfig::Instance()->GetOutputVariables(output_variables);
    }
    assert(output_variables.size() == num_vars);

    // Loop over the requested variables
    for (unsigned var_index=0; var_index<num_vars; var_index++)
    {
        // Create vector for storing values over the local nodes
        Vec variable_data =  this->mpMesh->GetDistributedVectorFactory()->CreateVec();
        DistributedVector distributed_var_data = this->mpMesh->GetDistributedVectorFactory()->CreateDistributedVector(variable_data);

        // Loop over the local nodes and gather the data
        for (DistributedVector::Iterator index = distributed_var_data.Begin();
             index!= distributed_var_data.End();
             ++index)
        {
            // Find the variable in the cell model
            AbstractCardiacCell* p_cell = this->mpCardiacTissue->GetCardiacCell(index.Global);
            unsigned cell_id = p_cell->GetAnyVariableIndex(output_variables[var_index]);
            // Store its value
            distributed_var_data[index] = p_cell->GetAnyVariable(cell_id, mCurrentTime);
        }
        distributed_var_data.Restore();

        // Write it to disc
        this->mpWriter->PutVector(mExtraVariablesId[var_index], variable_data);

        PetscTools::Destroy(variable_data);
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::InitialiseWriter()
{
    bool extend_file = (mSolution != NULL);

    // I think this is impossible to trip; certainly it's very difficult!
    assert(!mpWriter);

    try
    {
        mpWriter = new Hdf5DataWriter(*mpMesh->GetDistributedVectorFactory(),
                                      HeartConfig::Instance()->GetOutputDirectory(),
                                      HeartConfig::Instance()->GetOutputFilenamePrefix(),
                                      !extend_file, // don't clear directory if extension requested
                                      extend_file);
    }
    catch (Exception& e)
    {
        // The constructor only throws an Exception if we're extending
        assert(extend_file);
        // Tried to extend and failed, so just create from scratch
        extend_file = false;
        mpWriter = new Hdf5DataWriter(*mpMesh->GetDistributedVectorFactory(),
                                      HeartConfig::Instance()->GetOutputDirectory(),
                                      HeartConfig::Instance()->GetOutputFilenamePrefix(),
                                      !extend_file,
                                      extend_file);
    }

    // Define columns, or get the variable IDs from the writer
    DefineWriterColumns(extend_file);

    //Possibility of applying a permutation
    if (HeartConfig::Instance()->GetOutputUsingOriginalNodeOrdering())
    {
        bool success = mpWriter->ApplyPermutation(mpMesh->rGetNodePermutation());
        if (success == false)
        {
            //It's not really a permutation, so reset
            HeartConfig::Instance()->SetOutputUsingOriginalNodeOrdering(false);
        }
    }


    if (!extend_file)
    {
        mpWriter->EndDefineMode();
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetOutputNodes(std::vector<unsigned> &nodesToOutput)
{
    mNodesToOutput = nodesToOutput;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
Hdf5DataReader AbstractCardiacProblem<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GetDataReader()
{
    if( (HeartConfig::Instance()->GetOutputDirectory()=="") || (HeartConfig::Instance()->GetOutputFilenamePrefix()==""))
    {
        EXCEPTION("Data reader invalid as data writer cannot be initialised");
    }
    return Hdf5DataReader(HeartConfig::Instance()->GetOutputDirectory(), HeartConfig::Instance()->GetOutputFilenamePrefix());
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
bool AbstractCardiacProblem<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::GetHasBath()
{
    return false;
}
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM,SPACE_DIM,PROBLEM_DIM>::SetElectrodes()
{
}


/////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////

// Monodomain
template class AbstractCardiacProblem<1,1,1>;
template class AbstractCardiacProblem<1,2,1>;
template class AbstractCardiacProblem<1,3,1>;
template class AbstractCardiacProblem<2,2,1>;
template class AbstractCardiacProblem<3,3,1>;

// Bidomain
template class AbstractCardiacProblem<1,1,2>;
template class AbstractCardiacProblem<2,2,2>;
template class AbstractCardiacProblem<3,3,2>;

// Extended Bidomain
template class AbstractCardiacProblem<1,1,3>;
template class AbstractCardiacProblem<2,2,3>;
template class AbstractCardiacProblem<3,3,3>;
