/*

Copyright (c) 2005-2026, University of Oxford.
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

#include "UblasCustomFunctions.hpp"
#include "BidomainTissue.hpp"
#include "DistributedVector.hpp"
#include "Exception.hpp"
#include "GenericMeshReader.hpp"
#include "Hdf5ToCmguiConverter.hpp"
#include "Hdf5ToMeshalyzerConverter.hpp"
#include "Hdf5ToVtkConverter.hpp"
#include "ChastePoint.hpp"
#include "DistributedTetrahedralMeshPartitionType.hpp"
#include "HeartEventHandler.hpp"
#include "HeartRegionCodes.hpp"
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
          mHdf5DataWriterChunkSizeAndAlignment(0),
          mSimulationDuration(-1.0),
          mOdeTimeStep(0.01),
          mPdeTimeStep(0.01),
          mPrintingTimeStep(0.01),
          mOutputDirectory("ChasteResults"),
          mOutputFilenamePrefix("SimulationResults"),
          mOutputUsingOriginalNodeOrdering(false),
          mVisualizeWithMeshalyzer(false),
          mVisualizeWithCmgui(false),
          mVisualizeWithVtk(false),
          mVisualizeWithParallelVtk(false),
          mVisualizerOutputPrecision(0),
          mCheckpointSimulation(false),
          mCheckpointTimestep(-1.0),
          mMaxCheckpointsOnDisk(UINT_MAX),
          mMeshPartitioning(DistributedTetrahedralMeshPartitionType::PARMETIS_LIBRARY),
          mUseAbsoluteTolerance(true),
          mKspAbsoluteTolerance(2e-4),
          mKspRelativeTolerance(1e-6),
          mKspSolver("cg"),
          mKspPreconditioner("bjacobi"),
          mUseMassLumping(false),
          mUseMassLumpingForPrecond(false),
          mUseFixedNumberIterations(false),
          mEvaluateNumItsEveryNSolves(0),
          mUseStateVariableInterpolation(false),
          mUseReactionDiffusionOperatorSplitting(false),
          mSurfaceAreaToVolumeRatio(1400.0),
          mCapacitance(1.0),
          mFibreFileType("")
{
    mTissueIdentifiers.insert(0u);
    mBathIdentifiers.insert(1u);
    mIntraConductivitiesOrthotropic = scalar_vector<double>(SPACE_DIM, 1.75);
    mExtraConductivitiesOrthotropic = scalar_vector<double>(SPACE_DIM, 7.0);
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
          mHdf5DataWriterChunkSizeAndAlignment(0),
          mSimulationDuration(-1.0),
          mOdeTimeStep(0.01),
          mPdeTimeStep(0.01),
          mPrintingTimeStep(0.01),
          mOutputDirectory("ChasteResults"),
          mOutputFilenamePrefix("SimulationResults"),
          mOutputUsingOriginalNodeOrdering(false),
          mVisualizeWithMeshalyzer(false),
          mVisualizeWithCmgui(false),
          mVisualizeWithVtk(false),
          mVisualizeWithParallelVtk(false),
          mVisualizerOutputPrecision(0),
          mCheckpointSimulation(false),
          mCheckpointTimestep(-1.0),
          mMaxCheckpointsOnDisk(UINT_MAX),
          mMeshPartitioning(DistributedTetrahedralMeshPartitionType::PARMETIS_LIBRARY),
          mUseAbsoluteTolerance(true),
          mKspAbsoluteTolerance(2e-4),
          mKspRelativeTolerance(1e-6),
          mKspSolver("cg"),
          mKspPreconditioner("bjacobi"),
          mUseMassLumping(false),
          mUseMassLumpingForPrecond(false),
          mUseFixedNumberIterations(false),
          mEvaluateNumItsEveryNSolves(0),
          mUseStateVariableInterpolation(false),
          mUseReactionDiffusionOperatorSplitting(false),
          mSurfaceAreaToVolumeRatio(1400.0),
          mCapacitance(1.0),
          mFibreFileType("")
{
    mTissueIdentifiers.insert(0u);
    mBathIdentifiers.insert(1u);
    mIntraConductivitiesOrthotropic = scalar_vector<double>(SPACE_DIM, 1.75);
    mExtraConductivitiesOrthotropic = scalar_vector<double>(SPACE_DIM, 7.0);
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
    else if (!mMeshFilename.empty())
    {
        auto* p_mesh = new DistributedTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>(mMeshPartitioning);
        std::shared_ptr<AbstractMeshReader<ELEMENT_DIM, SPACE_DIM>> p_reader
            = GenericMeshReader<ELEMENT_DIM, SPACE_DIM>(mMeshFilename);
        p_mesh->ConstructFromMeshReader(*p_reader);
        mpMesh = p_mesh;
        mAllocatedMemoryForMesh = true;
    }
    else
    {
        EXCEPTION("No mesh given: call SetMesh() or SetMeshFileName() before Initialise()");
    }
    mpCellFactory->SetMesh(mpMesh);
    mpCellFactory->SetPdeTimeStep(mPdeTimeStep);
    mpCellFactory->SetOdeTimeStep(mOdeTimeStep);
    HeartEventHandler::EndEvent(HeartEventHandler::READ_MESH);

    // Propagate the tissue/bath region identifiers to HeartRegionCode global state
    HeartRegionCode::SetTissueIdentifiers(mTissueIdentifiers);
    HeartRegionCode::SetBathIdentifiers(mBathIdentifiers);

    HeartEventHandler::BeginEvent(HeartEventHandler::INITIALISE);

    delete mpCardiacTissue; // In case we're called twice
    mpCardiacTissue = CreateCardiacTissue();

    // Propagate conductivity and fibre settings from problem to tissue
    {
        mpCardiacTissue->SetIntracellularConductivities(mIntraConductivitiesOrthotropic);
        if (!mFibreFilePath.empty())
        {
            mpCardiacTissue->SetFibreOrientationFile(mFibreFilePath, mFibreFileType);
        }
        for (unsigned i = 0; i < mConductivityHeterogeneityAreas.size(); i++)
        {
            mpCardiacTissue->AddConductivityHeterogeneity(mConductivityHeterogeneityAreas[i],
                                                          mConductivityHeterogeneityIntra[i],
                                                          mConductivityHeterogeneityExtra[i]);
        }
        mpCardiacTissue->RebuildConductivityTensors();
        // Propagate Am/Cm and bath conductivities
        mpCardiacTissue->SetSurfaceAreaToVolumeRatio(mSurfaceAreaToVolumeRatio);
        mpCardiacTissue->SetCapacitance(mCapacitance);
        mpCardiacTissue->SetBathConductivities(mBathConductivities);

        // For bidomain: propagate extracellular conductivities
        BidomainTissue<SPACE_DIM>* p_bidomain_tissue = dynamic_cast<BidomainTissue<SPACE_DIM>*>(mpCardiacTissue);
        if (p_bidomain_tissue)
        {
            p_bidomain_tissue->SetExtracellularConductivities(mExtraConductivitiesOrthotropic);
            p_bidomain_tissue->RebuildExtracellularConductivityTensors();
        }
    }

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
    if (mSimulationDuration <= 0.0)
    {
        EXCEPTION("Simulation duration must be set (call SetSimulationDuration()) and be positive");
    }
    if (mSimulationDuration <= mCurrentTime)
    {
        EXCEPTION("End time should be in the future");
    }
    if (mPrintOutput)
    {
        if ((mOutputDirectory == "") || (mOutputFilenamePrefix == ""))
        {
            EXCEPTION("Either explicitly specify not to print output (call PrintOutput(false)) or specify the output directory and filename prefix");
        }
    }

    double end_time = mSimulationDuration;
    double pde_time = mPdeTimeStep;

    /*
     * MatrixIsConstant stuff requires CONSTANT dt - do some checks to make sure
     * the TimeStepper won't find non-constant dt.
     * Note: printing_time does not have to divide end_time, but dt must divide
     * printing_time and end_time.
     * The problem checks pde_dt divides printing dt.
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
void AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetMeshFileName(const std::string& rFilePath)
{
    mMeshFilename = rFilePath;
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
                        mSimulationDuration,
                        mPrintingTimeStep,
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

        progress_reporter_dir = mOutputDirectory;
    }
    else
    {
        progress_reporter_dir = ""; // progress printed to CHASTE_TEST_OUTPUT
    }
    for (boost::shared_ptr<AbstractOutputModifier> p_output_modifier : mOutputModifiers)
    {
        p_output_modifier->SetOutputDirectory(mOutputDirectory);
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
                                       mSimulationDuration);
    progress_reporter.Update(mCurrentTime);

    mpSolver->SetTimeStep(mPdeTimeStep);
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

        for (boost::shared_ptr<AbstractOutputModifier> p_output_modifier : mOutputModifiers)
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
    for (boost::shared_ptr<AbstractOutputModifier> p_output_modifier : mOutputModifiers)
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

    FileFinder test_output(mOutputDirectory, RelativeTo::ChasteTestOutput);

    /********************************************************************************
     * Run all post processing.
     *********************************************************************************/

    HeartEventHandler::BeginEvent(HeartEventHandler::POST_PROC);
    bool any_post_processing = !mApdMaps.empty() || !mUpstrokeTimeMaps.empty()
                               || !mMaxUpstrokeVelocityMaps.empty() || !mConductionVelocityMaps.empty()
                               || !mNodalTimeTraces.empty() || !mPseudoEcgElectrodePositions.empty();
    if (any_post_processing)
    {
        PostProcessingWriter<ELEMENT_DIM, SPACE_DIM> post_writer(*mpMesh,
                                                                 test_output,
                                                                 mOutputFilenamePrefix,
                                                                 "V",
                                                                 mHdf5DataWriterChunkSizeAndAlignment);
        post_writer.WritePostProcessingFiles(mApdMaps,
                                             mUpstrokeTimeMaps,
                                             mMaxUpstrokeVelocityMaps,
                                             mConductionVelocityMaps,
                                             mNodalTimeTraces,
                                             mPseudoEcgElectrodePositions,
                                             mOutputUsingOriginalNodeOrdering);
    }
    HeartEventHandler::EndEvent(HeartEventHandler::POST_PROC);

    /********************************************************************************************
     * Convert HDF5 datasets (solution and postprocessing maps) to different visualizer formats
     ********************************************************************************************/

    HeartEventHandler::BeginEvent(HeartEventHandler::DATA_CONVERSION);
    // Only if results files were written and we are outputting all nodes
    if (mNodesToOutput.empty())
    {
        if (mVisualizeWithMeshalyzer)
        {
            // Convert simulation data to Meshalyzer format
            Hdf5ToMeshalyzerConverter<ELEMENT_DIM, SPACE_DIM> converter(test_output,
                                                                        mOutputFilenamePrefix,
                                                                        mpMesh,
                                                                        mOutputUsingOriginalNodeOrdering,
                                                                        mVisualizerOutputPrecision);
        }

        if (mVisualizeWithCmgui)
        {
            // Convert simulation data to Cmgui format
            Hdf5ToCmguiConverter<ELEMENT_DIM, SPACE_DIM> converter(test_output,
                                                                   mOutputFilenamePrefix,
                                                                   mpMesh,
                                                                   GetHasBath(),
                                                                   mVisualizerOutputPrecision,
                                                                   mOutputUsingOriginalNodeOrdering);
        }

        if (mVisualizeWithVtk)
        {
            // Convert simulation data to VTK format
            Hdf5ToVtkConverter<ELEMENT_DIM, SPACE_DIM> converter(test_output,
                                                                 mOutputFilenamePrefix,
                                                                 mpMesh,
                                                                 false,
                                                                 mOutputUsingOriginalNodeOrdering);
        }

        if (mVisualizeWithParallelVtk)
        {
            // Convert simulation data to parallel VTK (pvtu) format
            Hdf5ToVtkConverter<ELEMENT_DIM, SPACE_DIM> converter(test_output,
                                                                 mOutputFilenamePrefix,
                                                                 mpMesh,
                                                                 true,
                                                                 mOutputUsingOriginalNodeOrdering);
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
                if (mOutputUsingOriginalNodeOrdering)
                {
                    EXCEPTION("Output using original node ordering is meaningless when outputting particular nodes in parallel. (Nodes are written with their original indices by default).");
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
                            mSimulationDuration,
                            mPrintingTimeStep);

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
    if (!mOutputVariables.empty())
    {
        // Get their names in a vector
        std::vector<std::string> output_variables = mOutputVariables;
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
    // Use the stored list of variable names
    const std::vector<std::string>& output_variables = mOutputVariables;
    unsigned num_vars = mExtraVariablesId.size();
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
        FileFinder h5_file(OutputFileHandler::GetChasteTestOutputDirectory() + mOutputDirectory
                               + "/" + mOutputFilenamePrefix + ".h5",
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
            Hdf5DataReader reader(mOutputDirectory,
                                  mOutputFilenamePrefix,
                                  true);
            std::vector<double> times = reader.GetUnlimitedDimensionValues();
            if (times.back() > mCurrentTime)
            {
                EXCEPTION("Attempting to extend " << h5_file.GetAbsolutePath() << " with results from time = " << mCurrentTime << ", but it already contains results up to time = " << times.back() << "."
                                                                                                                                                                                                       " Calling SetOutputDirectory() before Solve() will direct results elsewhere.");
            }
        }
        PetscTools::Barrier("InitialiseWriter::Extension check");
    }
    mpWriter = new Hdf5DataWriter(*mpMesh->GetDistributedVectorFactory(),
                                  mOutputDirectory,
                                  mOutputFilenamePrefix,
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
    if (mOutputUsingOriginalNodeOrdering)
    {
        bool success = mpWriter->ApplyPermutation(mpMesh->rGetNodePermutation(), true /*unsafe mode - extending*/);
        if (success == false)
        {
            //It's not really a permutation, so reset
            mOutputUsingOriginalNodeOrdering = false;
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
    if ((mOutputDirectory == "") || (mOutputFilenamePrefix == ""))
    {
        EXCEPTION("Data reader invalid as data writer cannot be initialised");
    }
    return Hdf5DataReader(mOutputDirectory, mOutputFilenamePrefix);
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

// -----------------------------------------------------------------------
// Setter/getter implementations for simulation settings
// -----------------------------------------------------------------------

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetSimulationDuration(double duration)
{
    mSimulationDuration = duration;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetOdePdeAndPrintingTimeSteps(
    double odeTimeStep, double pdeTimeStep, double printingTimeStep)
{
    mOdeTimeStep = odeTimeStep;
    mPdeTimeStep = pdeTimeStep;
    mPrintingTimeStep = printingTimeStep;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetOdeTimeStep(double odeTimeStep)
{
    mOdeTimeStep = odeTimeStep;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetPdeTimeStep(double pdeTimeStep)
{
    mPdeTimeStep = pdeTimeStep;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetPrintingTimeStep(double printingTimeStep)
{
    mPrintingTimeStep = printingTimeStep;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
double AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::GetSimulationDuration() const
{
    return mSimulationDuration;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
double AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::GetOdeTimeStep() const
{
    return mOdeTimeStep;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
double AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::GetPdeTimeStep() const
{
    return mPdeTimeStep;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
double AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::GetPrintingTimeStep() const
{
    return mPrintingTimeStep;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetOutputDirectory(const std::string& rOutputDirectory)
{
    mOutputDirectory = rOutputDirectory;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetOutputFilenamePrefix(const std::string& rPrefix)
{
    mOutputFilenamePrefix = rPrefix;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetOutputVariables(
    const std::vector<std::string>& rVariables)
{
    mOutputVariables = rVariables;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetOutputUsingOriginalNodeOrdering(bool useOriginal)
{
    mOutputUsingOriginalNodeOrdering = useOriginal;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
bool AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::GetOutputUsingOriginalNodeOrdering() const
{
    return mOutputUsingOriginalNodeOrdering;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
std::string AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::GetOutputDirectory() const
{
    return mOutputDirectory;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
std::string AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::GetOutputFilenamePrefix() const
{
    return mOutputFilenamePrefix;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetVisualizeWithMeshalyzer(bool vis)
{
    mVisualizeWithMeshalyzer = vis;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetVisualizeWithCmgui(bool vis)
{
    mVisualizeWithCmgui = vis;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetVisualizeWithVtk(bool vis)
{
    mVisualizeWithVtk = vis;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetVisualizeWithParallelVtk(bool vis)
{
    mVisualizeWithParallelVtk = vis;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetVisualizerOutputPrecision(unsigned precision)
{
    mVisualizerOutputPrecision = precision;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetCheckpointSimulation(
    bool checkpointSimulation, double checkpointTimestep, unsigned maxCheckpointsOnDisk)
{
    mCheckpointSimulation = checkpointSimulation;
    mCheckpointTimestep = checkpointTimestep;
    mMaxCheckpointsOnDisk = maxCheckpointsOnDisk;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetMeshPartitioning(
    DistributedTetrahedralMeshPartitionType::type method)
{
    mMeshPartitioning = method;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetKspAbsoluteTolerance(double tol)
{
    mUseAbsoluteTolerance = true;
    mKspAbsoluteTolerance = tol;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetKspRelativeTolerance(double tol)
{
    mUseAbsoluteTolerance = false;
    mKspRelativeTolerance = tol;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
double AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::GetKspAbsoluteTolerance() const
{
    return mKspAbsoluteTolerance;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
bool AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::GetUseAbsoluteTolerance() const
{
    return mUseAbsoluteTolerance;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetKspSolverType(const std::string& solver)
{
    mKspSolver = solver;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetKspPreconditionerType(const std::string& preconditioner)
{
    mKspPreconditioner = preconditioner;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetUseMassLumping(bool useLumping)
{
    mUseMassLumping = useLumping;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetUseMassLumpingForPrecond(bool useLumping)
{
    mUseMassLumpingForPrecond = useLumping;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetUseFixedNumberIterationsLinearSolver(
    bool useFixedIts, unsigned evaluateEvery)
{
    mUseFixedNumberIterations = useFixedIts;
    mEvaluateNumItsEveryNSolves = evaluateEvery;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetUseStateVariableInterpolation(bool use)
{
    mUseStateVariableInterpolation = use;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetUseReactionDiffusionOperatorSplitting(bool use)
{
    mUseReactionDiffusionOperatorSplitting = use;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetSurfaceAreaToVolumeRatio(double ratio)
{
    mSurfaceAreaToVolumeRatio = ratio;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
double AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::GetSurfaceAreaToVolumeRatio() const
{
    return mSurfaceAreaToVolumeRatio;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetCapacitance(double capacitance)
{
    mCapacitance = capacitance;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
double AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::GetCapacitance() const
{
    return mCapacitance;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::AddApdMap(double repolarisationPct, double threshold)
{
    mApdMaps.push_back(std::make_pair(repolarisationPct, threshold));
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::AddUpstrokeTimeMap(double threshold)
{
    mUpstrokeTimeMaps.push_back(threshold);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::AddMaxUpstrokeVelocityMap(double threshold)
{
    mMaxUpstrokeVelocityMaps.push_back(threshold);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::AddConductionVelocityMap(unsigned sourceNode)
{
    mConductionVelocityMaps.push_back(sourceNode);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::AddNodalTimeTrace(unsigned nodeIndex)
{
    mNodalTimeTraces.push_back(nodeIndex);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::AddPseudoEcgElectrode(
    const ChastePoint<SPACE_DIM>& rPosition)
{
    mPseudoEcgElectrodePositions.push_back(rPosition);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetFibreOrientationFile(
    const std::string& rFilePath, const std::string& rFileType)
{
    mFibreFilePath = rFilePath;
    mFibreFileType = rFileType;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::AddConductivityHeterogeneity(
    boost::shared_ptr<AbstractChasteRegion<SPACE_DIM> > pRegion,
    const c_vector<double, SPACE_DIM>& rIntraConductivities,
    const c_vector<double, SPACE_DIM>& rExtraConductivities)
{
    mConductivityHeterogeneityAreas.push_back(pRegion);
    mConductivityHeterogeneityIntra.push_back(rIntraConductivities);
    mConductivityHeterogeneityExtra.push_back(rExtraConductivities);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetTissueAndBathIdentifiers(
    const std::set<unsigned>& rTissueIds, const std::set<unsigned>& rBathIds)
{
    HeartRegionCode::SetTissueAndBathIdentifiers(rTissueIds, rBathIds);
    mTissueIdentifiers = rTissueIds;
    mBathIdentifiers = rBathIds;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
const std::set<unsigned>& AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::rGetTissueIdentifiers() const
{
    return mTissueIdentifiers;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
const std::set<unsigned>& AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::rGetBathIdentifiers() const
{
    return mBathIdentifiers;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetBathMultipleConductivities(
    const std::map<unsigned, double>& rBathConductivities)
{
    mBathConductivities = rBathConductivities;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
double AbstractCardiacProblem<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::GetBathConductivity(unsigned bathRegion) const
{
    if (bathRegion != UINT_MAX)
    {
        auto it = mBathConductivities.find(bathRegion);
        if (it != mBathConductivities.end())
        {
            return it->second;
        }
    }
    return 7.0; // default bath conductivity (mS/cm)
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
