
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

#ifndef ABSTRACTDYNAMICLINEARPDESOLVER_HPP_
#define ABSTRACTDYNAMICLINEARPDESOLVER_HPP_

#include "Exception.hpp"
#include "TimeStepper.hpp"
#include "AbstractLinearPdeSolver.hpp"
#include "PdeSimulationTime.hpp"
#include "AbstractTimeAdaptivityController.hpp"
#include "Hdf5DataWriter.hpp"
#include "Hdf5ToVtkConverter.hpp"
#include "Hdf5ToTxtConverter.hpp"

/**
 * Abstract class for dynamic linear PDE solves.
 * This class defines the Solve() method. The concrete class should implement
 * the SetupLinearSystem() method (defined in AbstractLinearPdeSolver), based
 * on the PDE being solved and the numerical method.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
class AbstractDynamicLinearPdeSolver : public AbstractLinearPdeSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>
{
    friend class TestSimpleLinearParabolicSolver;
protected:

    /** Simulation start time. */
    double mTstart;

    /** Simulation end time. */
    double mTend;

    /** Whether SetTimes has been called with suitable parameters. */
    bool mTimesSet;

    /** The initial condition vector. */
    Vec mInitialCondition;

   /** Whether the matrix has been assembled for the current time step. */
    bool mMatrixIsAssembled;

    /**
     * Whether the matrix is constant in time (if so the system need not be assembled at each time step).
     * Defaults to false.
     */
    bool mMatrixIsConstant;

    /**
     * The timestep to use. This is either the last timestep passed in in SetTimeStep,
     * or the last timestep suggested by the time adaptivity controller.
     */
    double mIdealTimeStep;

    /** The last actual timestep used. */
    double mLastWorkingTimeStep;

    /** A controller which determines what timestep to use (defaults to NULL). */
    AbstractTimeAdaptivityController* mpTimeAdaptivityController;

    /**
     * Flag to say if we need to output to VTK.
     * Defaults to false in the constructor.
     */
    bool mOutputToVtk;

    /**
     * Flag to say if we need to output to VTK parallel (.pvtu).
     * Defaults to false in the constructor.
     */
    bool mOutputToParallelVtk;

    /**
     * Flag to say if we need to output to a .txt format that is readable by Matlab.
     * Defaults to false in the constructor.
     */
    bool mOutputToTxt;

    /** Output directory (a subfolder of tmp/[USERNAME]/testoutput). */
    std::string mOutputDirectory;

    /** Filename prefix for HDF5 and other files. */
    std::string mFilenamePrefix;

    /**
     * The ratio of the number of actual timesteps to the number
     * of timesteps at which results are output to HDF5 and other files.
     * Defaults to 1 in the constructor.
     */
    unsigned mPrintingTimestepMultiple;

    /** The object used to write results to HDF5 file. */
    Hdf5DataWriter* mpHdf5Writer;

    /** List of variable column IDs as written to HDF5 file. */
    std::vector<int> mVariableColumnIds;

    /**
     * Create and initialise the HDF5 writer.
     * Called by Solve() if results are to be output.
     */
    void InitialiseHdf5Writer();

    /**
     * Write one timestep of output data to HDF5 file.
     *
     * @param time  the time
     * @param solution  the solution vector to write
     */
    void WriteOneStep(double time, Vec solution);

public:

    /**
     * Constructor.
     *
     * @param pMesh the mesh
     */
    AbstractDynamicLinearPdeSolver(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh);

    /**
     * Set the times to solve between.
     *
     * @param tStart the start time
     * @param tEnd the end time
     */
    void SetTimes(double tStart, double tEnd);

    /**
     * Set (or reset) the timestep to use.
     *
     * @param dt timestep
     */
    void SetTimeStep(double dt);

    /**
     * Set the initial condition.
     *
     * @note We do *not* take responsibility for destroying this vector - the caller must
     * do so once the solver is no longer in use.
     *
     * @param initialCondition the initial condition
     */
    void SetInitialCondition(Vec initialCondition);

    /** Dynamic solve method.
     * @return solution vector
     */
    virtual Vec Solve();

    /** Tell the solver to assemble the matrix again next timestep. */
    void SetMatrixIsNotAssembled();

    /**
     * Set a controller class which alters the dt used.
     *
     * @param pTimeAdaptivityController the controller
     */
    void SetTimeAdaptivityController(AbstractTimeAdaptivityController* pTimeAdaptivityController);

    /**
     * @param output whether to output to VTK (.vtu) file
     */
    void SetOutputToVtk(bool output);

    /**
     * @param output whether to output to parallel VTK (.pvtu) file
     */
    void SetOutputToParallelVtk(bool output);

    /**
     * @param output whether to output to a .txt format that is readable by Matlab
     */
    void SetOutputToTxt(bool output);

    /**
     * @param outputDirectory the output directory
     * @param prefix the filename prefix
     */
    void SetOutputDirectoryAndPrefix(std::string outputDirectory, std::string prefix);

    /**
     * @param multiple the ratio of the number of actual timesteps to the number
     * of timesteps at which results are output to HDF5 and other files.
     */
    void SetPrintingTimestepMultiple(unsigned multiple);
};

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractDynamicLinearPdeSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::InitialiseHdf5Writer()
{
    // Check that everything is set up correctly
    if ((mOutputDirectory=="") || (mFilenamePrefix==""))
    {
        EXCEPTION("Output directory or filename prefix has not been set");
    }

    // Create writer
    mpHdf5Writer = new Hdf5DataWriter(*(this->mpMesh)->GetDistributedVectorFactory(),
                                      mOutputDirectory,
                                      mFilenamePrefix);

    // Set writer to output all nodes
    mpHdf5Writer->DefineFixedDimension((this->mpMesh)->GetNumNodes());

    // Only used to get an estimate of the number of timesteps below
    unsigned estimated_num_printing_timesteps = 1u + (unsigned)((mTend - mTstart)/(mIdealTimeStep*mPrintingTimestepMultiple));

    /**
     * Note: For now, writing variable names as 'Variable_0' etc; in the future,
     * could allow user to specify units of time and names and units of
     * dependent variables to be passed to the writer using DefineVariable() and
     * DefineUnlimitedDimension()
     */
    assert(mVariableColumnIds.empty());
    for (unsigned i=0; i<PROBLEM_DIM; i++)
    {
        std::stringstream variable_name;
        variable_name << "Variable_" << i;
        mVariableColumnIds.push_back(mpHdf5Writer->DefineVariable(variable_name.str(),"undefined"));
    }
    mpHdf5Writer->DefineUnlimitedDimension("Time", "undefined", estimated_num_printing_timesteps);

    // End the define mode of the writer
    mpHdf5Writer->EndDefineMode();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
AbstractDynamicLinearPdeSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::AbstractDynamicLinearPdeSolver(AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh)
    : AbstractLinearPdeSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>(pMesh),
      mTimesSet(false),
      mInitialCondition(nullptr),
      mMatrixIsAssembled(false),
      mMatrixIsConstant(false),
      mIdealTimeStep(-1.0),
      mLastWorkingTimeStep(-1),
      mpTimeAdaptivityController(nullptr),
      mOutputToVtk(false),
      mOutputToParallelVtk(false),
      mOutputToTxt(false),
      mOutputDirectory(""),
      mFilenamePrefix(""),
      mPrintingTimestepMultiple(1),
      mpHdf5Writer(nullptr)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractDynamicLinearPdeSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetTimes(double tStart, double tEnd)
{
    mTstart = tStart;
    mTend = tEnd;

    if (mTstart >= mTend)
    {
        EXCEPTION("Start time has to be less than end time");
    }

    mTimesSet = true;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractDynamicLinearPdeSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetTimeStep(double dt)
{
    if (dt <= 0)
    {
        EXCEPTION("Time step has to be greater than zero");
    }

    mIdealTimeStep = dt;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractDynamicLinearPdeSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetInitialCondition(Vec initialCondition)
{
    assert(initialCondition != nullptr);
    mInitialCondition = initialCondition;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractDynamicLinearPdeSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::WriteOneStep(double time, Vec solution)
{
    mpHdf5Writer->PutUnlimitedVariable(time);
    if (PROBLEM_DIM == 1)
    {
        mpHdf5Writer->PutVector(mVariableColumnIds[0], solution);
    }
    else
    {
        mpHdf5Writer->PutStripedVector(mVariableColumnIds, solution);
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
Vec AbstractDynamicLinearPdeSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::Solve()
{
    // Begin by checking that everything has been set up correctly
    if (!mTimesSet)
    {
        EXCEPTION("SetTimes() has not been called");
    }
    if ((mIdealTimeStep <= 0.0) && (mpTimeAdaptivityController==nullptr))
    {
        EXCEPTION("SetTimeStep() has not been called");
    }
    if (mInitialCondition == nullptr)
    {
        EXCEPTION("SetInitialCondition() has not been called");
    }

    // If required, initialise HDF5 writer and output initial condition to HDF5 file
    bool print_output = (mOutputToVtk || mOutputToParallelVtk || mOutputToTxt);
    if (print_output)
    {
        InitialiseHdf5Writer();
        WriteOneStep(mTstart, mInitialCondition);
        mpHdf5Writer->AdvanceAlongUnlimitedDimension();
    }

    this->InitialiseForSolve(mInitialCondition);

    if (mIdealTimeStep < 0) // hasn't been set, so a controller must have been given
    {
        mIdealTimeStep = mpTimeAdaptivityController->GetNextTimeStep(mTstart, mInitialCondition);
    }

    /*
     * Note: we use the mIdealTimeStep here (the original timestep that was passed in, or
     * the last timestep suggested by the controller), rather than the last timestep used
     * (mLastWorkingTimeStep), because the timestep will be very slightly altered by the
     * stepper in the final timestep of the last printing-timestep-loop, and these floating
     * point errors can add up and eventually cause exceptions being thrown.
     */
    TimeStepper stepper(mTstart, mTend, mIdealTimeStep, mMatrixIsConstant);

    Vec solution = mInitialCondition;
    Vec next_solution;

    while (!stepper.IsTimeAtEnd())
    {
        bool timestep_changed = false;

        PdeSimulationTime::SetTime(stepper.GetTime());

        // Determine timestep to use
        double new_dt;
        if (mpTimeAdaptivityController)
        {
            // Get the timestep the controller wants to use and store it as the ideal timestep
            mIdealTimeStep = mpTimeAdaptivityController->GetNextTimeStep(stepper.GetTime(), solution);

            // Tell the stepper to use this timestep from now on...
            stepper.ResetTimeStep(mIdealTimeStep);

            // ..but now get the timestep from the stepper, as the stepper might need
            // to trim the timestep if it would take us over the end time
            new_dt = stepper.GetNextTimeStep();
            // Changes in timestep bigger than 0.001% will trigger matrix re-computation
            timestep_changed = (fabs(new_dt/mLastWorkingTimeStep - 1.0) > 1e-5);
        }
        else
        {
            new_dt = stepper.GetNextTimeStep();

            //new_dt should be roughly the same size as mIdealTimeStep - we should never need to take a tiny step

            if (mMatrixIsConstant && fabs(new_dt/mIdealTimeStep - 1.0) > 1e-5)
            {
                // Here we allow for changes of up to 0.001%
                // Note that the TimeStepper guarantees that changes in dt are no bigger than DBL_EPSILON*current_time
                NEVER_REACHED;
            }
        }

        // Save the timestep as the last one use, and also put it in PdeSimulationTime
        // so everyone can see it
        mLastWorkingTimeStep = new_dt;
        PdeSimulationTime::SetPdeTimeStepAndNextTime(new_dt, stepper.GetNextTime());

        // Solve
        try
        {
            // (This runs the cell ODE models in heart simulations)
            this->PrepareForSetupLinearSystem(solution);
        }
        catch(Exception& e)
        {
            // We only need to clean up memory if we are NOT on the first PDE time step,
            // as someone else cleans up the mInitialCondition vector in higher classes.
            if (solution != mInitialCondition)
            {
                HeartEventHandler::BeginEvent(HeartEventHandler::COMMUNICATION);
                PetscTools::Destroy(solution);
                HeartEventHandler::EndEvent(HeartEventHandler::COMMUNICATION);
            }
            throw e;
        }

        bool compute_matrix = (!mMatrixIsConstant || !mMatrixIsAssembled || timestep_changed);

        this->SetupLinearSystem(solution, compute_matrix);

        this->FinaliseLinearSystem(solution);

        if (compute_matrix)
        {
            this->mpLinearSystem->ResetKspSolver();
        }

        next_solution = this->mpLinearSystem->Solve(solution);

        if (mMatrixIsConstant)
        {
            mMatrixIsAssembled = true;
        }

        this->FollowingSolveLinearSystem(next_solution);

        stepper.AdvanceOneTimeStep();

        // Avoid memory leaks
        if (solution != mInitialCondition)
        {
            HeartEventHandler::BeginEvent(HeartEventHandler::COMMUNICATION);
            PetscTools::Destroy(solution);
            HeartEventHandler::EndEvent(HeartEventHandler::COMMUNICATION);
        }
        solution = next_solution;

        // If required, output next solution to HDF5 file
        if (print_output && (stepper.GetTotalTimeStepsTaken()%mPrintingTimestepMultiple == 0) )
        {
            WriteOneStep(stepper.GetTime(), solution);
            mpHdf5Writer->AdvanceAlongUnlimitedDimension();
        }
    }

    // Avoid memory leaks
    if (mpHdf5Writer != nullptr)
    {
        delete mpHdf5Writer;
        mpHdf5Writer = nullptr;
    }

    // Convert HDF5 output to other formats as required

    if (mOutputToVtk)
    {
        Hdf5ToVtkConverter<ELEMENT_DIM,SPACE_DIM> converter(FileFinder(mOutputDirectory, RelativeTo::ChasteTestOutput),
                                                            mFilenamePrefix, this->mpMesh, false, false);
    }
    if (mOutputToParallelVtk)
    {
        Hdf5ToVtkConverter<ELEMENT_DIM,SPACE_DIM> converter(FileFinder(mOutputDirectory, RelativeTo::ChasteTestOutput),
                                                            mFilenamePrefix, this->mpMesh, true, false);
    }
    if (mOutputToTxt)
    {
        Hdf5ToTxtConverter<ELEMENT_DIM,SPACE_DIM> converter(FileFinder(mOutputDirectory, RelativeTo::ChasteTestOutput),
                                                            mFilenamePrefix, this->mpMesh);
    }

    return solution;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractDynamicLinearPdeSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetMatrixIsNotAssembled()
{
    mMatrixIsAssembled = false;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractDynamicLinearPdeSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetTimeAdaptivityController(AbstractTimeAdaptivityController* pTimeAdaptivityController)
{
    assert(pTimeAdaptivityController != nullptr);
    assert(mpTimeAdaptivityController == NULL);
    mpTimeAdaptivityController = pTimeAdaptivityController;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractDynamicLinearPdeSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetOutputToVtk(bool output)
{
    mOutputToVtk = output;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractDynamicLinearPdeSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetOutputToParallelVtk(bool output)
{
    mOutputToParallelVtk = output;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractDynamicLinearPdeSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetOutputToTxt(bool output)
{
    mOutputToTxt = output;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractDynamicLinearPdeSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetOutputDirectoryAndPrefix(std::string outputDirectory, std::string prefix)
{
    mOutputDirectory = outputDirectory;
    mFilenamePrefix = prefix;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM, unsigned PROBLEM_DIM>
void AbstractDynamicLinearPdeSolver<ELEMENT_DIM, SPACE_DIM, PROBLEM_DIM>::SetPrintingTimestepMultiple(unsigned multiple)
{
    mPrintingTimestepMultiple = multiple;
}

#endif /*ABSTRACTDYNAMICLINEARPDESOLVER_HPP_*/
