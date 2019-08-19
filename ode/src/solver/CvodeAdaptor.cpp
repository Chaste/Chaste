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

#ifdef CHASTE_CVODE

#include "CvodeAdaptor.hpp"

#include "Exception.hpp"
#include "MathsCustomFunctions.hpp" // For tolerance check.
#include "TimeStepper.hpp"
#include "VectorHelperFunctions.hpp"

#include <iostream>
#include <sstream>

// CVODE headers
#include <sundials/sundials_nvector.h>

#if CHASTE_SUNDIALS_VERSION >= 30000
#include <cvode/cvode_direct.h> /* access to CVDls interface            */
#include <sundials/sundials_types.h> /* defs. of realtype, sunindextype      */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#else
#include <cvode/cvode_dense.h>
#endif

/**
 * CVODE right-hand-side function adaptor.
 *
 * The CVODE solvers require the RHS of the ODE system to be defined
 * by a function of this type.  We use pData to access the Chaste ODE system,
 * and call EvaluateYDerivatives appropriately.
 *
 * Note that this requires copying the state variable and derivatives vectors,
 * and thus introduces a slight overhead.
 *
 * @param t  the current time
 * @param y  the current state variable values
 * @param ydot  to be filled in with the derivatives
 * @param pData  a pointer to the AbstractOdeSystem to evaluate
 * @return 0 on success, -1 for an unrecoverable error
 */
int CvodeRhsAdaptor(realtype t, N_Vector y, N_Vector ydot, void* pData)
{
    assert(pData != nullptr);
    CvodeData* p_data = (CvodeData*)pData;
    // Get y, ydot into std::vector<>s
    static std::vector<realtype> ydot_vec;
    CopyToStdVector(y, *p_data->pY);
    CopyToStdVector(ydot, ydot_vec);
    // Call our function
    try
    {
        p_data->pSystem->EvaluateYDerivatives(t, *(p_data->pY), ydot_vec);
    }
    catch (const Exception& e)
    {
        std::cerr << "CVODE RHS Exception: " << e.GetMessage() << std::endl
                  << std::flush;
        return -1;
    }
    // Copy derivative back
    CopyFromStdVector(ydot_vec, ydot);
    return 0;
}

/**
 * CVODE root-finder function adaptor.
 *
 * Adapt the Chaste AbstractOdeSystem::CalculateStoppingEvent method for use by CVODE.
 *
 * This function computes a vector-valued function g(t, y) such that the roots of the
 * components g_i(t, y) are to be found during the integration.
 *
 * Unfortunately, AbstractOdeSystem::CalculateStoppingEvent returns a boolean value,
 * so we have to cheat in the definition of g.
 *
 * Note that this function requires copying the state variable vector, and thus
 * introduces a slight overhead.
 *
 * @param t  the current time
 * @param y  the current state variable values
 * @param pGOut  pointer to array to be filled in with the g_i(t, y) values
 * @param pData  a pointer to the AbstractOdeSystem to use
 * @return 0 on success, negative on error
 */
int CvodeRootAdaptor(realtype t, N_Vector y, realtype* pGOut, void* pData)
{
    assert(pData != nullptr);
    CvodeData* p_data = (CvodeData*)pData;
    // Get y into a std::vector
    CopyToStdVector(y, *p_data->pY);
    // Call our function
    try
    {
        *pGOut = p_data->pSystem->CalculateRootFunction(t, *p_data->pY);
    }
    catch (const Exception& e)
    {
        std::cerr << "CVODE Root Exception: " << e.GetMessage() << std::endl
                  << std::flush;
        return -1;
    }
    return 0;
}

// /**
//  * Jacobian computation adaptor function.
//  *
//  * If solving an AbstractOdeSystemWithAnalyticJacobian, this function
//  * can be used to allow CVODE to compute the Jacobian analytically.
//  *
//  * Note to self: can test using pSystem->GetUseAnalyticJacobian().
//  */
// int CvodeDenseJacobianAdaptor(long int numberOfStateVariables, DenseMat J,
//                               realtype t, N_Vector y, N_Vector fy,
//                               void* pData,
//                               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
// {
//     AbstractOdeSystemWithAnalyticJacobian* pSystem
//         = (AbstractOdeSystemWithAnalyticJacobian*) pData;
//     // Get the current time step
//     double dt;
//     CVodeGetCurrentStep(CvodeMem, &dt);
//     // Get std::vector<> for y and double** for J
//     std::vector<realtype>& y_vec = *NV_DATA_STL(y);
//     double** ppJ = J->data; // organised column-wise: J_{i,j} = ppJ[j][i]
//     // Call our function
//     try
//     {
//         pSystem->AnalyticJacobian(y_vec, ppJ, t, dt);
//     }
//     catch (const Exception &e)
//     {
//         std::cerr << "CVODE J Exception: " << e.GetMessage() << std::endl << std::flush;
//         return -1;
//     }
//     // Update J (if needed)
//     return 0;
// }

void CvodeErrorHandler(int errorCode, const char* module, const char* function,
                       char* message, void* pData)
{
    std::stringstream err;
    err << "CVODE Error " << errorCode << " in module " << module
        << " function " << function << ": " << message;
    std::cerr << "*" << err.str() << std::endl
              << std::flush;
    // Throwing the exception here causes termination on Maverick (g++ 4.4)
    //EXCEPTION(err.str());
}

void CvodeAdaptor::SetupCvode(AbstractOdeSystem* pOdeSystem,
                              std::vector<double>& rInitialY,
                              double startTime,
                              double maxStep)
{
    assert(rInitialY.size() == pOdeSystem->GetNumberOfStateVariables());
    assert(maxStep > 0.0);

    /** Initial conditions for the ODE solver. */
    N_Vector initial_values = N_VMake_Serial(rInitialY.size(), &(rInitialY[0]));
    assert(NV_DATA_S(initial_values) == &(rInitialY[0]));
    assert(!NV_OWN_DATA_S(initial_values));
    //    std::cout << " Initial values: "; N_VPrint_Serial(initial_values);
    //    std::cout << " Rtol: " << mRelTol << ", Atol: " << mAbsTol << std::endl;
    //    std::cout << " Start: " << startTime << " max dt=" << maxStep << std::endl << std::flush;

    // Find out if we need to (re-)initialise
    bool reinit = !mpCvodeMem || mForceReset || !mLastSolutionState || !CompareDoubles::WithinAnyTolerance(startTime, mLastSolutionTime);
    if (!reinit && !mForceMinimalReset)
    {
        const unsigned size = GetVectorSize(rInitialY);
        for (unsigned i = 0; i < size; i++)
        {
            if (!CompareDoubles::WithinAnyTolerance(GetVectorComponent(mLastSolutionState, i), GetVectorComponent(rInitialY, i)))
            {
                reinit = true;
                break;
            }
        }
    }

    if (!mpCvodeMem) // First run of this solver, set up CVODE memory
    {
        // Set up CVODE's memory.
#if CHASTE_SUNDIALS_VERSION >= 40000
        //  v4.0.0 release notes: instead of specifying the nonlinear iteration type when creating the CVODE(S) memory structure,
        //  CVODE(S) uses the SUNNONLINSOL_NEWTON module implementation of a Newton iteration by default.
        mpCvodeMem = CVodeCreate(CV_BDF);
#else
        mpCvodeMem = CVodeCreate(CV_BDF, CV_NEWTON);
#endif

        if (mpCvodeMem == nullptr)
            EXCEPTION("Failed to SetupCvode CVODE"); // LCOV_EXCL_LINE
        // Set error handler
        CVodeSetErrHandlerFn(mpCvodeMem, CvodeErrorHandler, nullptr);
        // Set the user data
        mData.pSystem = pOdeSystem;
        mData.pY = &rInitialY;
#if CHASTE_SUNDIALS_VERSION >= 20400
        CVodeSetUserData(mpCvodeMem, (void*)(&mData));
#else
        CVodeSetFdata(mpCvodeMem, (void*)(&mData));
#endif

// Setup CVODE
#if CHASTE_SUNDIALS_VERSION >= 20400
        CVodeInit(mpCvodeMem, CvodeRhsAdaptor, startTime, initial_values);
        CVodeSStolerances(mpCvodeMem, mRelTol, mAbsTol);
#else
        CVodeMalloc(mpCvodeMem, CvodeRhsAdaptor, startTime, initial_values,
                    CV_SS, mRelTol, &mAbsTol);
#endif

        // Set the rootfinder function if wanted
        if (mCheckForRoots)
        {
#if CHASTE_SUNDIALS_VERSION >= 20400
            CVodeRootInit(mpCvodeMem, 1, CvodeRootAdaptor);
#else
            CVodeRootInit(mpCvodeMem, 1, CvodeRootAdaptor, (void*)(&mData));
#endif
        }

#if CHASTE_SUNDIALS_VERSION >= 30000
        /* Create dense SUNMatrix for use in linear solves */
        mpSundialsDenseMatrix = SUNDenseMatrix(rInitialY.size(), rInitialY.size());
#endif

#if CHASTE_SUNDIALS_VERSION >= 40000
        /* Create dense SUNLinearSolver object for use by CVode */
        mpSundialsLinearSolver = SUNLinSol_Dense(initial_values, mpSundialsDenseMatrix);

        /* Call CVodeSetLinearSolver to attach the matrix and linear solver to CVode */
        CVodeSetLinearSolver(mpCvodeMem, mpSundialsLinearSolver, mpSundialsDenseMatrix);
#elif CHASTE_SUNDIALS_VERSION >= 30000
        /* Create dense SUNLinearSolver object for use by CVode */
        mpSundialsLinearSolver = SUNDenseLinearSolver(initial_values, mpSundialsDenseMatrix);

        /* Call CVDlsSetLinearSolver to attach the matrix and linear solver to CVode */
        CVDlsSetLinearSolver(mpCvodeMem, mpSundialsLinearSolver, mpSundialsDenseMatrix);
#else
        // CVODE < v3.0.0
        // Attach a linear solver for Newton iteration
        CVDense(mpCvodeMem, rInitialY.size());
#endif
    }
    else if (reinit) // Could be new ODE system, or new Y values
    {
        // Set the user data
        mData.pSystem = pOdeSystem; // stays the same on a re-initialize
        mData.pY = &rInitialY; // changes on a re-initialize
#if CHASTE_SUNDIALS_VERSION >= 20400
        CVodeSetUserData(mpCvodeMem, (void*)(&mData));
#else
        CVodeSetFdata(mpCvodeMem, (void*)(&mData));
#endif

#if CHASTE_SUNDIALS_VERSION >= 20400
        CVodeReInit(mpCvodeMem, startTime, initial_values);
        CVodeSStolerances(mpCvodeMem, mRelTol, mAbsTol);
#else
        CVodeReInit(mpCvodeMem, CvodeRhsAdaptor, startTime, initial_values,
                    CV_SS, mRelTol, &mAbsTol);
#endif

#if CHASTE_SUNDIALS_VERSION >= 30000
        if (mpSundialsLinearSolver)
        {
            /* Free the linear solver memory */
            SUNLinSolFree(mpSundialsLinearSolver);
        }
        if (mpSundialsDenseMatrix)
        {
            /* Free the matrix memory */
            SUNMatDestroy(mpSundialsDenseMatrix);
        }

        /* Create dense matrix of type SUNDenseMatrix for use in linear solves */
        mpSundialsDenseMatrix = SUNDenseMatrix(rInitialY.size(), rInitialY.size());

#if CHASTE_SUNDIALS_VERSION >= 40000
        /* Create dense SUNLinSol_Dense object for use by CVode */
        mpSundialsLinearSolver = SUNLinSol_Dense(initial_values, mpSundialsDenseMatrix);

        /* Call CVodeSetLinearSolver to attach the matrix and linear solver to CVode */
        CVodeSetLinearSolver(mpCvodeMem, mpSundialsLinearSolver, mpSundialsDenseMatrix);
#else
        /* Create dense SUNLinearSolver object for use by CVode */
        mpSundialsLinearSolver = SUNDenseLinearSolver(initial_values, mpSundialsDenseMatrix);

        /* Call CVDlsSetLinearSolver to attach the matrix and linear solver to CVode */
        CVDlsSetLinearSolver(mpCvodeMem, mpSundialsLinearSolver, mpSundialsDenseMatrix);
#endif
#else
        // Attach a linear solver for Newton iteration
        CVDense(mpCvodeMem, rInitialY.size());
#endif
    }

    CVodeSetMaxStep(mpCvodeMem, maxStep);
    // Change max steps if wanted
    if (mMaxSteps > 0)
    {
        CVodeSetMaxNumSteps(mpCvodeMem, mMaxSteps);
        CVodeSetMaxErrTestFails(mpCvodeMem, 15);
    }
    DeleteVector(initial_values);
}

void CvodeAdaptor::FreeCvodeMemory()
{
    if (mpCvodeMem)
    {
        CVodeFree(&mpCvodeMem);
    }
    mpCvodeMem = nullptr;

#if CHASTE_SUNDIALS_VERSION >= 30000
    if (mpSundialsLinearSolver)
    {
        /* Free the linear solver memory */
        SUNLinSolFree(mpSundialsLinearSolver);
    }
    mpSundialsLinearSolver = nullptr;

    if (mpSundialsDenseMatrix)
    {
        /* Free the matrix memory */
        SUNMatDestroy(mpSundialsDenseMatrix);
    }
    mpSundialsDenseMatrix = nullptr;
#endif
}

void CvodeAdaptor::SetForceReset(bool autoReset)
{
    mForceReset = autoReset;
    if (mForceReset)
    {
        SetMinimalReset(false);
        ResetSolver();
    }
}

void CvodeAdaptor::SetMinimalReset(bool minimalReset)
{
    mForceMinimalReset = minimalReset;
    if (mForceMinimalReset)
    {
        SetForceReset(false);
    }
}

void CvodeAdaptor::ResetSolver()
{
    // Doing this makes the SetupCvode think it needs to reset things.
    DeleteVector(mLastSolutionState);
}

void CvodeAdaptor::CvodeError(int flag, const char* msg)
{
    std::stringstream err;
    char* p_flag_name = CVodeGetReturnFlagName(flag);
    err << msg << ": " << p_flag_name;
    free(p_flag_name);
    std::cerr << err.str() << std::endl
              << std::flush;
    EXCEPTION(err.str());
}

OdeSolution CvodeAdaptor::Solve(AbstractOdeSystem* pOdeSystem,
                                std::vector<double>& rYValues,
                                double startTime,
                                double endTime,
                                double maxStep,
                                double timeSampling)
{
    assert(endTime > startTime);
    assert(timeSampling > 0.0);

    mStoppingEventOccurred = false;
    if (mCheckForRoots && pOdeSystem->CalculateStoppingEvent(startTime, rYValues) == true)
    {
        EXCEPTION("(Solve with sampling) Stopping event is true for initial condition");
    }

    SetupCvode(pOdeSystem, rYValues, startTime, maxStep);

    TimeStepper stepper(startTime, endTime, timeSampling);
    N_Vector yout = N_VMake_Serial(rYValues.size(), &(rYValues[0]));

    // Set up ODE solution
    OdeSolution solutions;
    solutions.SetNumberOfTimeSteps(stepper.EstimateTimeSteps());
    solutions.rGetSolutions().push_back(rYValues);
    solutions.rGetTimes().push_back(startTime);
    solutions.SetOdeSystemInformation(pOdeSystem->GetSystemInformation());

    // Main time sampling loop
    while (!stepper.IsTimeAtEnd() && !mStoppingEventOccurred)
    {
        // This should stop CVODE going past the end of where we wanted and interpolating back.
        int ierr = CVodeSetStopTime(mpCvodeMem, stepper.GetNextTime());
        assert(ierr == CV_SUCCESS);
        UNUSED_OPT(ierr); // avoid unused var warning

        double tend;
        ierr = CVode(mpCvodeMem, stepper.GetNextTime(), yout, &tend, CV_NORMAL);
        if (ierr < 0)
        {
            FreeCvodeMemory();
            DeleteVector(yout);
            CvodeError(ierr, "CVODE failed to solve system");
        }
        // Store solution
        solutions.rGetSolutions().push_back(rYValues);
        solutions.rGetTimes().push_back(tend);
        if (ierr == CV_ROOT_RETURN)
        {
            // Stopping event occurred
            mStoppingEventOccurred = true;
            mStoppingTime = tend;
        }
        mLastSolutionTime = tend;
        stepper.AdvanceOneTimeStep();
    }

    // stepper.EstimateTimeSteps may have been an overestimate...
    solutions.SetNumberOfTimeSteps(stepper.GetTotalTimeStepsTaken());

    int ierr = CVodeGetLastStep(mpCvodeMem, &mLastInternalStepSize);
    assert(ierr == CV_SUCCESS);
    UNUSED_OPT(ierr); // avoid unused var warning
    RecordStoppingPoint(mLastSolutionTime, yout);
    DeleteVector(yout);

    return solutions;
}

void CvodeAdaptor::Solve(AbstractOdeSystem* pOdeSystem,
                         std::vector<double>& rYValues,
                         double startTime,
                         double endTime,
                         double maxStep)
{
    assert(endTime > startTime);

    mStoppingEventOccurred = false;
    if (mCheckForRoots && pOdeSystem->CalculateStoppingEvent(startTime, rYValues) == true)
    {
        EXCEPTION("(Solve) Stopping event is true for initial condition");
    }

    SetupCvode(pOdeSystem, rYValues, startTime, maxStep);

    N_Vector yout = N_VMake_Serial(rYValues.size(), &(rYValues[0]));

    // This should stop CVODE going past the end of where we wanted and interpolating back.
    int ierr = CVodeSetStopTime(mpCvodeMem, endTime);
    assert(ierr == CV_SUCCESS);
    UNUSED_OPT(ierr); // avoid unused var warning

    double tend;
    ierr = CVode(mpCvodeMem, endTime, yout, &tend, CV_NORMAL);
    if (ierr < 0)
    {
        FreeCvodeMemory();
        DeleteVector(yout);
        CvodeError(ierr, "CVODE failed to solve system");
    }
    if (ierr == CV_ROOT_RETURN)
    {
        // Stopping event occurred
        mStoppingEventOccurred = true;
        mStoppingTime = tend;
    }
    assert(NV_DATA_S(yout) == &(rYValues[0]));
    assert(!NV_OWN_DATA_S(yout));

    //    long int steps;
    //    CVodeGetNumSteps(mpCvodeMem, &steps);
    //    std::cout << " Solved to " << endTime << " in " << steps << " steps.\n";

    ierr = CVodeGetLastStep(mpCvodeMem, &mLastInternalStepSize);
    assert(ierr == CV_SUCCESS);
    UNUSED_OPT(ierr); // avoid unused var warning
    RecordStoppingPoint(tend, yout);
    DeleteVector(yout);
}

CvodeAdaptor::CvodeAdaptor(double relTol, double absTol)
        : AbstractIvpOdeSolver(),
          mpCvodeMem(nullptr),
          mRelTol(relTol),
          mAbsTol(absTol),
          mLastInternalStepSize(-0.0),
          mMaxSteps(0),
          mCheckForRoots(false),
          mLastSolutionState(nullptr),
          mLastSolutionTime(0.0),
#if CHASTE_SUNDIALS_VERSION >= 20400
          mForceReset(false),
#else
          mForceReset(true),
#endif
          mForceMinimalReset(false)
#if CHASTE_SUNDIALS_VERSION >= 30000
          ,
          mpSundialsDenseMatrix(nullptr),
          mpSundialsLinearSolver(nullptr)
#endif
{
}

CvodeAdaptor::~CvodeAdaptor()
{
    FreeCvodeMemory();
    DeleteVector(mLastSolutionState);
}

void CvodeAdaptor::RecordStoppingPoint(double stopTime, N_Vector yEnd)
{
    if (!mForceReset)
    {
        const unsigned size = GetVectorSize(yEnd);
        CreateVectorIfEmpty(mLastSolutionState, size);
        for (unsigned i = 0; i < size; i++)
        {
            SetVectorComponent(mLastSolutionState, i, GetVectorComponent(yEnd, i));
        }
        mLastSolutionTime = stopTime;
    }
}

void CvodeAdaptor::SetTolerances(double relTol, double absTol)
{
    mRelTol = relTol;
    mAbsTol = absTol;
}

double CvodeAdaptor::GetRelativeTolerance()
{
    return mRelTol;
}

double CvodeAdaptor::GetAbsoluteTolerance()
{
    return mAbsTol;
}

double CvodeAdaptor::GetLastStepSize()
{
    return mLastInternalStepSize;
}

void CvodeAdaptor::CheckForStoppingEvents()
{
    mCheckForRoots = true;
}

void CvodeAdaptor::SetMaxSteps(long int numSteps)
{
    mMaxSteps = numSteps;
}

long int CvodeAdaptor::GetMaxSteps()
{
    return mMaxSteps;
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(CvodeAdaptor)

#endif // CHASTE_CVODE
