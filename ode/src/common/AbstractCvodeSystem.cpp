/*

Copyright (c) 2005-2012, University of Oxford.
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

#ifdef CHASTE_CVODE

#include <sstream>
#include <cassert>

#include "AbstractCvodeSystem.hpp"
#include "Exception.hpp"
#include "VectorHelperFunctions.hpp"
#include "MathsCustomFunctions.hpp" // For tolerance comparison
#include "TimeStepper.hpp"
#include "CvodeAdaptor.hpp" // For CvodeErrorHandler

// CVODE headers
#include <cvode/cvode.h>
#include <sundials/sundials_nvector.h>
#include <cvode/cvode_dense.h>


/**
 * Callback function provided to CVODE to allow it to 'call' C++ member functions
 * (in particular, AbstractCvodeCell::EvaluateYDerivatives).
 *
 * @param t  current time
 * @param y  state variable vector
 * @param ydot  derivatives vector to be filled in
 * @param pData  pointer to the cell being simulated
 */
int AbstractCvodeSystemRhsAdaptor(realtype t, N_Vector y, N_Vector ydot, void *pData)
{
    assert(pData != NULL);
    AbstractCvodeSystem* p_ode_system = (AbstractCvodeSystem*) pData;
    try
    {
        p_ode_system->EvaluateYDerivatives(t, y, ydot);
    }
    catch (const Exception &e)
    {
        std::cerr << "CVODE RHS Exception: " << e.GetMessage()
                  << std::endl << std::flush;
        return -1;
    }
    return 0;
}

AbstractCvodeSystem::AbstractCvodeSystem(unsigned numberOfStateVariables)
    : AbstractParameterisedSystem<N_Vector>(numberOfStateVariables),
      mLastSolutionState(NULL),
      mLastSolutionTime(0.0),
      mAutoReset(true),
      mUseAnalyticJacobian(false),
      mpCvodeMem(NULL),
      mMaxSteps(0)
{
    SetTolerances();
}

void AbstractCvodeSystem::Init()
{
    DeleteVector(mStateVariables);
    mStateVariables = GetInitialConditions();
    DeleteVector(mParameters);
    mParameters = N_VNew_Serial(rGetParameterNames().size());
    for (int i=0; i<NV_LENGTH_S(mParameters); i++)
    {
        NV_Ith_S(mParameters, i) = 0.0;
    }
}

AbstractCvodeSystem::~AbstractCvodeSystem()
{
    FreeCvodeMemory();
    DeleteVector(mStateVariables);
    DeleteVector(mParameters);
    DeleteVector(mLastSolutionState);
}

//
//double AbstractCvodeSystem::CalculateRootFunction(double time, const std::vector<double>& rY)
//{
//    bool stop = CalculateStoppingEvent(time, rY);
//    return stop ? 0.0 : 1.0;
//}
//
//bool AbstractCvodeSystem::GetUseAnalyticJacobian()
//{
//    return mUseAnalyticJacobian;
//}



OdeSolution AbstractCvodeSystem::Solve(realtype tStart,
                                       realtype tEnd,
                                       realtype maxDt,
                                       realtype tSamp)
{
    assert(tEnd >= tStart);
    assert(tSamp > 0.0);

    SetupCvode(mStateVariables, tStart, maxDt);

    TimeStepper stepper(tStart, tEnd, tSamp);

    // Set up ODE solution
    OdeSolution solutions;
    solutions.SetNumberOfTimeSteps(stepper.EstimateTimeSteps());
    solutions.rGetSolutions().push_back(MakeStdVec(mStateVariables));
    solutions.rGetTimes().push_back(tStart);
    solutions.SetOdeSystemInformation(mpSystemInfo);

    // Main time sampling loop
    while (!stepper.IsTimeAtEnd())
    {
        double cvode_stopped_at;
        int ierr = CVode(mpCvodeMem, stepper.GetNextTime(), mStateVariables,
                         &cvode_stopped_at, CV_NORMAL);
        if (ierr<0)
        {
            FreeCvodeMemory();
            CvodeError(ierr, "CVODE failed to solve system");
        }
        // Not root finding, so should have reached requested time
        assert(fabs(cvode_stopped_at - stepper.GetNextTime()) < DBL_EPSILON);

        VerifyStateVariables();

        // Store solution
        solutions.rGetSolutions().push_back(MakeStdVec(mStateVariables));
        solutions.rGetTimes().push_back(cvode_stopped_at);
        stepper.AdvanceOneTimeStep();
    }

    // stepper.EstimateTimeSteps may have been an overestimate...
    solutions.SetNumberOfTimeSteps(stepper.GetTotalTimeStepsTaken());

    int ierr = CVodeGetLastStep(mpCvodeMem, &mLastInternalStepSize);
    assert(ierr == CV_SUCCESS); ierr=ierr; // avoid unused var warning

    RecordStoppingPoint(tEnd);

    return solutions;
}

void AbstractCvodeSystem::Solve(realtype tStart,
                                realtype tEnd,
                                realtype maxDt)
{
    assert(tEnd >= tStart);

    SetupCvode(mStateVariables, tStart, maxDt);

    double cvode_stopped_at;
    int ierr = CVode(mpCvodeMem, tEnd, mStateVariables, &cvode_stopped_at, CV_NORMAL);
    if (ierr<0)
    {
        FreeCvodeMemory();
        CvodeError(ierr, "CVODE failed to solve system");
    }
    // Not root finding, so should have reached requested time
    assert(fabs(cvode_stopped_at - tEnd) < DBL_EPSILON);

    ierr = CVodeGetLastStep(mpCvodeMem, &mLastInternalStepSize);
    assert(ierr == CV_SUCCESS); ierr=ierr; // avoid unused var warning

    RecordStoppingPoint(cvode_stopped_at);

    VerifyStateVariables();
}


void AbstractCvodeSystem::SetMaxSteps(long int numSteps)
{
    mMaxSteps = numSteps;
}

long int AbstractCvodeSystem::GetMaxSteps()
{
    return mMaxSteps;
}

void AbstractCvodeSystem::SetTolerances(double relTol, double absTol)
{
    mRelTol = relTol;
    mAbsTol = absTol;
    ResetSolver();
}

double AbstractCvodeSystem::GetRelativeTolerance()
{
    return mRelTol;
}

double AbstractCvodeSystem::GetAbsoluteTolerance()
{
    return mAbsTol;
}

double AbstractCvodeSystem::GetLastStepSize()
{
    return mLastInternalStepSize;
}

void AbstractCvodeSystem::SetAutoReset(bool autoReset)
{
    mAutoReset = autoReset;
    if (mAutoReset)
    {
        ResetSolver();
    }
}


void AbstractCvodeSystem::ResetSolver()
{
    DeleteVector(mLastSolutionState);
}


void AbstractCvodeSystem::SetupCvode(N_Vector initialConditions,
                                     realtype tStart,
                                     realtype maxDt)
{
    assert((unsigned)NV_LENGTH_S(initialConditions) == GetNumberOfStateVariables());
    assert(maxDt >= 0.0);

    // Find out if we need to (re-)initialise
    bool reinit = !mpCvodeMem || mAutoReset || !mLastSolutionState || !CompareDoubles::WithinAnyTolerance(tStart, mLastSolutionTime);
    if (!reinit)
    {
        const unsigned size = GetNumberOfStateVariables();
        for (unsigned i=0; i<size; i++)
        {
            if (!CompareDoubles::WithinAnyTolerance(GetVectorComponent(mLastSolutionState, i), GetVectorComponent(mStateVariables, i)))
            {
                reinit = true;
                break;
            }
        }
    }

    if (!mpCvodeMem)
    {
        mpCvodeMem = CVodeCreate(CV_BDF, CV_NEWTON);
        if (mpCvodeMem == NULL) EXCEPTION("Failed to SetupCvode CVODE");
        // Set error handler
        CVodeSetErrHandlerFn(mpCvodeMem, CvodeErrorHandler, NULL);
        // Set the user data
#if CHASTE_SUNDIALS_VERSION >= 20400
        CVodeSetUserData(mpCvodeMem, (void*)(this));
#else
        CVodeSetFdata(mpCvodeMem, (void*)(this));
#endif
        // Setup CVODE
#if CHASTE_SUNDIALS_VERSION >= 20400
        CVodeInit(mpCvodeMem, AbstractCvodeSystemRhsAdaptor, tStart, initialConditions);
        CVodeSStolerances(mpCvodeMem, mRelTol, mAbsTol);
#else
        CVodeMalloc(mpCvodeMem, AbstractCvodeSystemRhsAdaptor, tStart, initialConditions,
                    CV_SS, mRelTol, &mAbsTol);
#endif
        // Attach a linear solver for Newton iteration
        CVDense(mpCvodeMem, NV_LENGTH_S(initialConditions));
    }
    else if (reinit)
    {
#if CHASTE_SUNDIALS_VERSION >= 20400
        CVodeReInit(mpCvodeMem, tStart, initialConditions);
        CVodeSStolerances(mpCvodeMem, mRelTol, mAbsTol);
#else
        CVodeReInit(mpCvodeMem, AbstractCvodeSystemRhsAdaptor, tStart, initialConditions,
                    CV_SS, mRelTol, &mAbsTol);
#endif
    }

    // Set max dt and change max steps if wanted
    CVodeSetMaxStep(mpCvodeMem, maxDt);
    if (mMaxSteps > 0)
    {
        CVodeSetMaxNumSteps(mpCvodeMem, mMaxSteps);
        CVodeSetMaxErrTestFails(mpCvodeMem, 15);
    }
}


void AbstractCvodeSystem::RecordStoppingPoint(double stopTime)
{
    if (!mAutoReset)
    {
        const unsigned size = GetNumberOfStateVariables();
        CreateVectorIfEmpty(mLastSolutionState, size);
        for (unsigned i=0; i<size; i++)
        {
            SetVectorComponent(mLastSolutionState, i, GetVectorComponent(mStateVariables, i));
        }
        mLastSolutionTime = stopTime;
    }
}


void AbstractCvodeSystem::FreeCvodeMemory()
{
    if (mpCvodeMem)
    {
        CVodeFree(&mpCvodeMem);
    }
    mpCvodeMem = NULL;
}


void AbstractCvodeSystem::CvodeError(int flag, const char * msg)
{

    std::stringstream err;
    char* p_flag_name = CVodeGetReturnFlagName(flag);
    err << msg << ": " << p_flag_name;
    free(p_flag_name);
    std::cerr << err.str() << std::endl << std::flush;
    EXCEPTION(err.str());
}


#endif // CHASTE_CVODE
