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
#include "AbstractCardiacProblem.hpp"

#include "DistributedVector.hpp"
#include "Exception.hpp"
#include "GenericMeshReader.hpp"
#include "Hdf5ToCmguiConverter.hpp"
#include "Hdf5ToMeshalyzerConverter.hpp"
#include "Hdf5ToVtkConverter.hpp"
#include "HeartConfig.hpp"
#include "HeartEventHandler.hpp"
#include "LinearSystem.hpp"
#include "PetscTools.hpp"
#include "PostProcessingWriter.hpp"
#include "ProgressReporter.hpp"
#include "TimeStepper.hpp"

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::AbstractCardiacProblem(
    AbstractCardiacCellFactory<ELEMENT_DIM, SPACE_DIM>* pCellFactory)
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
          mpWriter(NULL),
          mUseHdf5DataWriterCache(false),
          mHdf5DataWriterChunkSizeAndAlignment(0)
{
    assert(mNodesToOutput.empty());
    if (!mpCellFactory)
    {
        EXCEPTION("AbstractCardiacProblem: Please supply a cell factory pointer to your cardiac problem constructor.");
    }
    HeartEventHandler::BeginEvent(HeartEventHandler::EVERYTHING);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::AbstractCardiacProblem()
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
          mpWriter(NULL),
          mUseHdf5DataWriterCache(false),
          mHdf5DataWriterChunkSizeAndAlignment(0)
{
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::~AbstractCardiacProblem()
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

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::Initialise()
{
    HeartEventHandler::BeginEvent(HeartEventHandler::READ_MESH);
    if (mpMesh)
    {
        if (PetscTools::IsParallel() && !dynamic_cast<DistributedTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>*>(mpMesh))
        {
            WARNING("Using a non-distributed mesh in a parallel simulation is not a good idea.");
        }
    }
    else
    {
        // If no mesh has been passed, we get it from the configuration file
        try
        {
            if (HeartConfig::Instance()->GetLoadMesh())
            {
                CreateMeshFromHeartConfig();
                std::shared_ptr<AbstractMeshReader<ELEMENT_DIM, SPACE_DIM> > p_mesh_reader
                    = GenericMeshReader<ELEMENT_DIM, SPACE_DIM>(HeartConfig::Instance()->GetMeshName());
                mpMesh->ConstructFromMeshReader(*p_mesh_reader);
            }
            else if (HeartConfig::Instance()->GetCreateMesh())
            {
                CreateMeshFromHeartConfig();
                assert(HeartConfig::Instance()->GetSpaceDimension() == SPACE_DIM);
                double inter_node_space = HeartConfig::Instance()->GetInterNodeSpace();

                switch (HeartConfig::Instance()->GetSpaceDimension())
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
    mpCellFactory->SetMesh(mpMesh);
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

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::CreateMeshFromHeartConfig()
{
    mpMesh = new DistributedTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>(HeartConfig::Instance()->GetMeshPartitioning());
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetBoundaryConditionsContainer(boost::shared_ptr<BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM> > pBcc)
{
    this->mpBoundaryConditionsContainer = pBcc;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::PreSolveChecks()
{
    if (mpCardiacTissue == NULL) // if tissue is NULL, Initialise() probably hasn't been called
    {
        EXCEPTION("Cardiac tissue is null, Initialise() probably hasn't been called");
    }
    if (HeartConfig::Instance()->GetSimulationDuration() <= mCurrentTime)
    {
        EXCEPTION("End time should be in the future");
    }
    if (mPrintOutput)
    {
        if ((HeartConfig::Instance()->GetOutputDirectory() == "") || (HeartConfig::Instance()->GetOutputFilenamePrefix() == ""))
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
    if (fabs(end_time - pde_time * round(end_time / pde_time)) > 1e-10)
    {
        EXCEPTION("PDE timestep does not seem to divide end time - check parameters");
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
Vec AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::CreateInitialCondition()
{
    DistributedVectorFactory* p_factory = mpMesh->GetDistributedVectorFactory();
    Vec initial_condition = p_factory->CreateVec(PROBLEM_DIM);
    DistributedVector ic = p_factory->CreateDistributedVector(initial_condition);
    std::vector<DistributedVector::Stripe> stripe;
    stripe.reserve(PROBLEM_DIM);

    for (unsigned i = 0; i < PROBLEM_DIM; i++)
    {
        stripe.push_back(DistributedVector::Stripe(ic, i));
    }

    for (DistributedVector::Iterator index = ic.Begin();
         index != ic.End();
         ++index)
    {
        stripe[0][index] = mpCardiacTissue->GetCardiacCell(index.Global)->GetVoltage();
        if (PROBLEM_DIM == 2)
        {
            stripe[1][index] = 0;
        }
    }

    ic.Restore();

    return initial_condition;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetMesh(AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>* pMesh)
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

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::PrintOutput(bool printOutput)
{
    mPrintOutput = printOutput;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetWriteInfo(bool writeInfo)
{
    mWriteInfo = writeInfo;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
Vec AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::GetSolution()
{
    return mSolution;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
DistributedVector AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::GetSolutionDistributedVector()
{
    return mpMesh->GetDistributedVectorFactory()->CreateDistributedVector(mSolution);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
double AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::GetCurrentTime()
{
    return mCurrentTime;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>& AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::rGetMesh()
{
    assert(mpMesh);
    return *mpMesh;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
AbstractCardiacTissue<ELEMENT_DIM, SPACE_DIM>* AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::GetTissue()
{
    if (mpCardiacTissue == NULL)
    {
        EXCEPTION("Tissue not yet set up, you may need to call Initialise() before GetTissue().");
    }
    return mpCardiacTissue;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetUseTimeAdaptivityController(
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

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::Solve()
{
    PreSolveChecks();

    std::vector<double> additional_stopping_times;
    SetUpAdditionalStoppingTimes(additional_stopping_times);

    TimeStepper stepper(mCurrentTime,
                        HeartConfig::Instance()->GetSimulationDuration(),
                        HeartConfig::Instance()->GetPrintingTimeStep(),
                        false,
                        additional_stopping_times);
    // Note that SetUpAdditionalStoppingTimes is a method from the BidomainWithBath class it adds
    // electrode events into the regular time-stepping
    //    EXCEPTION("Electrode switch on/off events should coincide with printing time steps.");

    if (!mpBoundaryConditionsContainer) // the user didn't supply a bcc
    {
        // Set up the default bcc
        mpDefaultBoundaryConditionsContainer.reset(new BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>);
        for (unsigned problem_index = 0; problem_index < PROBLEM_DIM; problem_index++)
        {
            mpDefaultBoundaryConditionsContainer->DefineZeroNeumannOnMeshBoundary(mpMesh, problem_index);
        }
        mpBoundaryConditionsContainer = mpDefaultBoundaryConditionsContainer;
    }

    assert(mpSolver == NULL);
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

    std::string progress_reporter_dir;

    if (mPrintOutput)
    {
        HeartEventHandler::BeginEvent(HeartEventHandler::WRITE_OUTPUT);
        bool extending_file = false;
        try
        {
            extending_file = InitialiseWriter();
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
                 * cleaned up in the constructor. So, only PetscTools::Destroy( initial_condition ) when
                 * it is not equal to mSolution.
                 */
                PetscTools::Destroy(initial_condition);
            }
            throw e;
        }

        /*
         * If we are resuming a simulation (i.e. mSolution already exists) and
         * we are extending a .h5 file that already exists then there is no need
         * to write the initial condition to file - it is already there as the
         * final solution of the previous run.
         */
        if (!(mSolution && extending_file))
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
    BOOST_FOREACH (boost::shared_ptr<AbstractOutputModifier> p_output_modifier, mOutputModifiers)
    {
        p_output_modifier->InitialiseAtStart(this->mpMesh->GetDistributedVectorFactory(), this->mpMesh->rGetNodePermutation());
        p_output_modifier->ProcessSolutionAtTimeStep(stepper.GetTime(), initial_condition, PROBLEM_DIM);
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
        mpSolver->SetInitialCondition(initial_condition);

        AtBeginningOfTimestep(stepper.GetTime());

        try
        {
            try
            {
                mSolution = mpSolver->Solve();
            }
            catch (const Exception& e)
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
                 * is going to be cleaned up in the destructor. So, only PetscTools::Destroy()
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

        BOOST_FOREACH (boost::shared_ptr<AbstractOutputModifier> p_output_modifier, mOutputModifiers)
        {
            p_output_modifier->ProcessSolutionAtTimeStep(stepper.GetTime(), mSolution, PROBLEM_DIM);
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
    BOOST_FOREACH (boost::shared_ptr<AbstractOutputModifier> p_output_modifier, mOutputModifiers)
    {
        p_output_modifier->FinaliseAtEnd();
    }
    CloseFilesAndPostProcess();
    HeartEventHandler::EndEvent(HeartEventHandler::EVERYTHING);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::CloseFilesAndPostProcess()
{
    // Close files
    if (!mPrintOutput)
    {
        // Nothing to do
        return;
    }
    HeartEventHandler::BeginEvent(HeartEventHandler::WRITE_OUTPUT);
    // If write caching is on, the next line might actually take a significant amount of time.
    delete mpWriter;
    mpWriter = NULL;
    HeartEventHandler::EndEvent(HeartEventHandler::WRITE_OUTPUT);

    FileFinder test_output(HeartConfig::Instance()->GetOutputDirectory(), RelativeTo::ChasteTestOutput);

    /********************************************************************************
     * Run all post processing.
     *
     * The PostProcessingWriter class examines what is requested in HeartConfig and
     * adds the relevant data to the HDF5 file.
     * This is converted to different visualizer formats along with the solution
     * in the DATA_CONVERSION block below.
     *********************************************************************************/

    HeartEventHandler::BeginEvent(HeartEventHandler::POST_PROC);
    if (HeartConfig::Instance()->IsPostProcessingRequested())
    {
        PostProcessingWriter<ELEMENT_DIM, SPACE_DIM> post_writer(*mpMesh,
                                                                 test_output,
                                                                 HeartConfig::Instance()->GetOutputFilenamePrefix(),
                                                                 "V",
                                                                 mHdf5DataWriterChunkSizeAndAlignment);
        post_writer.WritePostProcessingFiles();
    }
    HeartEventHandler::EndEvent(HeartEventHandler::POST_PROC);

    /********************************************************************************************
     * Convert HDF5 datasets (solution and postprocessing maps) to different visualizer formats
     ********************************************************************************************/

    HeartEventHandler::BeginEvent(HeartEventHandler::DATA_CONVERSION);
    // Only if results files were written and we are outputting all nodes
    if (mNodesToOutput.empty())
    {
        if (HeartConfig::Instance()->GetVisualizeWithMeshalyzer())
        {
            // Convert simulation data to Meshalyzer format
            Hdf5ToMeshalyzerConverter<ELEMENT_DIM, SPACE_DIM> converter(test_output,
                                                                        HeartConfig::Instance()->GetOutputFilenamePrefix(),
                                                                        mpMesh,
                                                                        HeartConfig::Instance()->GetOutputUsingOriginalNodeOrdering(),
                                                                        HeartConfig::Instance()->GetVisualizerOutputPrecision());
            std::string subdirectory_name = converter.GetSubdirectory();
            HeartConfig::Instance()->Write(false, subdirectory_name);
        }

        if (HeartConfig::Instance()->GetVisualizeWithCmgui())
        {
            // Convert simulation data to Cmgui format
            Hdf5ToCmguiConverter<ELEMENT_DIM, SPACE_DIM> converter(test_output,
                                                                   HeartConfig::Instance()->GetOutputFilenamePrefix(),
                                                                   mpMesh,
                                                                   GetHasBath(),
                                                                   HeartConfig::Instance()->GetVisualizerOutputPrecision());
            std::string subdirectory_name = converter.GetSubdirectory();
            HeartConfig::Instance()->Write(false, subdirectory_name);
        }

        if (HeartConfig::Instance()->GetVisualizeWithVtk())
        {
            // Convert simulation data to VTK format
            Hdf5ToVtkConverter<ELEMENT_DIM, SPACE_DIM> converter(test_output,
                                                                 HeartConfig::Instance()->GetOutputFilenamePrefix(),
                                                                 mpMesh,
                                                                 false,
                                                                 HeartConfig::Instance()->GetOutputUsingOriginalNodeOrdering());
            std::string subdirectory_name = converter.GetSubdirectory();
            HeartConfig::Instance()->Write(false, subdirectory_name);
        }

        if (HeartConfig::Instance()->GetVisualizeWithParallelVtk())
        {
            // Convert simulation data to parallel VTK (pvtu) format
            Hdf5ToVtkConverter<ELEMENT_DIM, SPACE_DIM> converter(test_output,
                                                                 HeartConfig::Instance()->GetOutputFilenamePrefix(),
                                                                 mpMesh,
                                                                 true,
                                                                 HeartConfig::Instance()->GetOutputUsingOriginalNodeOrdering());
            std::string subdirectory_name = converter.GetSubdirectory();
            HeartConfig::Instance()->Write(false, subdirectory_name);
        }
    }
    HeartEventHandler::EndEvent(HeartEventHandler::DATA_CONVERSION);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::DefineWriterColumns(bool extending)
{
    if (!extending)
    {
        if (mNodesToOutput.empty())
        {
            //Set writer to output all nodes
            mpWriter->DefineFixedDimension(mpMesh->GetNumNodes());
        }
        else
        {
            // Added for #2980
            if (mpMesh->rGetNodePermutation().size() > 0)
            {
                if (HeartConfig::Instance()->GetOutputUsingOriginalNodeOrdering())
                {
                    EXCEPTION("HeartConfig setting `GetOutputUsingOriginalNodeOrdering` is meaningless when outputting particular nodes in parallel. (Nodes are written with their original indices by default).");
                }
                std::vector<unsigned> nodes_to_output_permuted(mNodesToOutput.size());
                for (unsigned i = 0; i < mNodesToOutput.size(); i++)
                {
                    nodes_to_output_permuted[i] = mpMesh->rGetNodePermutation()[mNodesToOutput[i]];
                }
                mpWriter->DefineFixedDimension(mNodesToOutput, nodes_to_output_permuted, mpMesh->GetNumNodes());
            } else {
                // Output only the nodes indicated
                mpWriter->DefineFixedDimension(mNodesToOutput, mNodesToOutput, mpMesh->GetNumNodes());
            }
        }
        // mNodeColumnId = mpWriter->DefineVariable("Node", "dimensionless");
        mVoltageColumnId = mpWriter->DefineVariable("V", "mV");

        // Only used to get an estimate of the # of timesteps below
        TimeStepper stepper(mCurrentTime,
                            HeartConfig::Instance()->GetSimulationDuration(),
                            HeartConfig::Instance()->GetPrintingTimeStep());

        mpWriter->DefineUnlimitedDimension("Time", "msecs", stepper.EstimateTimeSteps() + 1); // plus one for start and end points
    }
    else
    {
        mVoltageColumnId = mpWriter->GetVariableByName("V");
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::DefineExtraVariablesWriterColumns(bool extending)
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
        for (unsigned var_index = 0; var_index < num_vars; var_index++)
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
                // Difficult to specify the units, as different cell models
                // at different points in the mesh could be using different units.
                column_id = this->mpWriter->DefineVariable(var_name, "unknown_units");
            }

            // Store column id
            mExtraVariablesId.push_back(column_id);
        }
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::WriteExtraVariablesOneStep()
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
    for (unsigned var_index = 0; var_index < num_vars; var_index++)
    {
        // Create vector for storing values over the local nodes
        Vec variable_data = this->mpMesh->GetDistributedVectorFactory()->CreateVec();
        DistributedVector distributed_var_data = this->mpMesh->GetDistributedVectorFactory()->CreateDistributedVector(variable_data);

        // Loop over the local nodes and gather the data
        for (DistributedVector::Iterator index = distributed_var_data.Begin();
             index != distributed_var_data.End();
             ++index)
        {
            // If the region is in the bath
            if (HeartRegionCode::IsRegionBath(this->mpMesh->GetNode(index.Global)->GetRegion()))
            {
                // Then we just pad the output with zeros, user currently needs to find a nice
                // way to deal with this in processing and visualization.
                distributed_var_data[index] = 0.0;
            }
            else
            {
                // Find the variable in the cell model and store its value
                distributed_var_data[index] = this->mpCardiacTissue->GetCardiacCell(index.Global)->GetAnyVariable(output_variables[var_index], mCurrentTime);
            }
        }
        distributed_var_data.Restore();

        // Write it to disc
        this->mpWriter->PutVector(mExtraVariablesId[var_index], variable_data);

        PetscTools::Destroy(variable_data);
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
bool AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::InitialiseWriter()
{
    bool extend_file = (mSolution != NULL);

    // I think this is impossible to trip; certainly it's very difficult!
    assert(!mpWriter);

    if (extend_file)
    {
        FileFinder h5_file(OutputFileHandler::GetChasteTestOutputDirectory() + HeartConfig::Instance()->GetOutputDirectory()
                               + "/" + HeartConfig::Instance()->GetOutputFilenamePrefix() + ".h5",
                           RelativeTo::Absolute);
        //We are going to test for existence before creating the file.
        //Therefore we should make sure that this existence test is thread-safe.
        //(If another process creates the file too early then we may get the wrong answer to the
        //existence question).
        PetscTools::Barrier("InitialiseWriter::Extension check");
        if (!h5_file.Exists())
        {
            extend_file = false;
        }
        else // if it does exist check that it is sensible to extend it by running from the archive we loaded.
        {
            Hdf5DataReader reader(HeartConfig::Instance()->GetOutputDirectory(),
                                  HeartConfig::Instance()->GetOutputFilenamePrefix(),
                                  true);
            std::vector<double> times = reader.GetUnlimitedDimensionValues();
            if (times.back() > mCurrentTime)
            {
                EXCEPTION("Attempting to extend " << h5_file.GetAbsolutePath() << " with results from time = " << mCurrentTime << ", but it already contains results up to time = " << times.back() << "."
                                                                                                                                                                                                       " Calling HeartConfig::Instance()->SetOutputDirectory() before Solve() will direct results elsewhere.");
            }
        }
        PetscTools::Barrier("InitialiseWriter::Extension check");
    }
    mpWriter = new Hdf5DataWriter(*mpMesh->GetDistributedVectorFactory(),
                                  HeartConfig::Instance()->GetOutputDirectory(),
                                  HeartConfig::Instance()->GetOutputFilenamePrefix(),
                                  !extend_file, // don't clear directory if extension requested
                                  extend_file,
                                  "Data",
                                  mUseHdf5DataWriterCache);

    /* If user has specified a chunk size and alignment parameter, pass it
     * through. We set them to the same value as we think this is the most
     * likely use case, specifically on striped filesystems where a chunk
     * should squeeze into a stripe.
     * Only happens if !extend_file, i.e. we're NOT loading a checkpoint, or
     * we are loading a checkpoint but the H5 file doesn't exist yet.
     */
    if (!extend_file && mHdf5DataWriterChunkSizeAndAlignment)
    {
        mpWriter->SetTargetChunkSize(mHdf5DataWriterChunkSizeAndAlignment);
        mpWriter->SetAlignment(mHdf5DataWriterChunkSizeAndAlignment);
    }

    // Define columns, or get the variable IDs from the writer
    DefineWriterColumns(extend_file);

    // Possibility of applying a permutation
    if (HeartConfig::Instance()->GetOutputUsingOriginalNodeOrdering())
    {
        bool success = mpWriter->ApplyPermutation(mpMesh->rGetNodePermutation(), true /*unsafe mode - extending*/);
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

    return extend_file;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetUseHdf5DataWriterCache(bool useCache)
{
    mUseHdf5DataWriterCache = useCache;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetHdf5DataWriterTargetChunkSizeAndAlignment(hsize_t size)
{
    mHdf5DataWriterChunkSizeAndAlignment = size;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetOutputNodes(std::vector<unsigned>& nodesToOutput)
{
    mNodesToOutput = nodesToOutput;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
Hdf5DataReader AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::GetDataReader()
{
    if ((HeartConfig::Instance()->GetOutputDirectory() == "") || (HeartConfig::Instance()->GetOutputFilenamePrefix() == ""))
    {
        EXCEPTION("Data reader invalid as data writer cannot be initialised");
    }
    return Hdf5DataReader(HeartConfig::Instance()->GetOutputDirectory(), HeartConfig::Instance()->GetOutputFilenamePrefix());
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
bool AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::GetHasBath()
{
    return false;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetElectrodes()
{
}

// Explicit instantiation

// Monodomain
template class AbstractCardiacProblem<1, 1, 1>;
template class AbstractCardiacProblem<1, 2, 1>;
template class AbstractCardiacProblem<1, 3, 1>;
template class AbstractCardiacProblem<2, 2, 1>;
template class AbstractCardiacProblem<3, 3, 1>;

// Bidomain
template class AbstractCardiacProblem<1, 1, 2>;
template class AbstractCardiacProblem<2, 2, 2>;
template class AbstractCardiacProblem<3, 3, 2>;

// Extended Bidomain
template class AbstractCardiacProblem<1, 1, 3>;
template class AbstractCardiacProblem<2, 2, 3>;
template class AbstractCardiacProblem<3, 3, 3>;
