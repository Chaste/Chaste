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

#include <cassert>
#include <sstream>

#include "AbstractCvodeSystem.hpp"
#include "CvodeAdaptor.hpp" // For CvodeErrorHandler
#include "Exception.hpp"
#include "MathsCustomFunctions.hpp" // For tolerance comparison
#include "TimeStepper.hpp"
#include "VectorHelperFunctions.hpp"

// CVODE headers
#include <cvode/cvode.h>
#include <sundials/sundials_nvector.h>

#if CHASTE_SUNDIALS_VERSION >= 30000
#include <cvode/cvode_direct.h> /* access to CVDls interface            */
#include <sundials/sundials_types.h> /* defs. of realtype, sunindextype      */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#else
#include <cvode/cvode_dense.h>
#endif

//#include "Debug.hpp"
//void DebugSteps(void* pCvodeMem, AbstractCvodeSystem* pSys)
//{
//    long int num_jac_evals, nniters, num_steps;
//    CVDenseGetNumJacEvals(pCvodeMem, &num_jac_evals);
////    CVDlsGetNumJacEvals(pCvodeMem, &num_jac_evals);
//    CVodeGetNumNonlinSolvIters(pCvodeMem, &nniters);
//    CVodeGetNumSteps(pCvodeMem, &num_steps);
//    double num_newton_iters = nniters;
//    PRINT_3_VARIABLES(pSys->GetSystemName(), num_newton_iters/num_steps, num_jac_evals/num_newton_iters);
//}
/**
 * Callback function provided to CVODE to allow it to 'call' C++ member functions
 * (in particular, AbstractCvodeCell::EvaluateYDerivatives).
 *
 * @param t  current time
 * @param y  state variable vector
 * @param ydot  derivatives vector to be filled in
 * @param pData  pointer to the cell being simulated
 */
int AbstractCvodeSystemRhsAdaptor(realtype t, N_Vector y, N_Vector ydot, void* pData)
{
    assert(pData != nullptr);
    AbstractCvodeSystem* p_ode_system = (AbstractCvodeSystem*)pData;
    try
    {
        p_ode_system->EvaluateYDerivatives(t, y, ydot);
    }
    catch (const Exception& e)
    {
#if CHASTE_SUNDIALS_VERSION <= 20300
        // Really old CVODE used to solve past the requested time points and could trigger this exception unnecessarily...
        if (e.CheckShortMessageContains("is outside the times stored in the data clamp") == "")
        {
            return 1; // This may be a recoverable error!
        }
#endif

        std::cerr << "CVODE RHS Exception: " << e.GetMessage()
                  << std::endl
                  << std::flush;
        return -1;
    }

    //    // Something like this might help CVODE when things are a bit unstable...
    //    try
    //    {
    //        p_ode_system->VerifyStateVariables();
    //    }
    //    catch (const Exception &e)
    //    {
    //        std::cout << "t = " << t << ":\t" <<  e.GetMessage() << std::endl << std::flush;
    //        return 1; // A positive return flag to CVODE tells it there's been an error but it might be recoverable.
    //    }

    return 0;
}

/*
 * Absolute chaos here with four different possible interfaces to the jacobian.
 */
#if CHASTE_SUNDIALS_VERSION >= 30000
// Sundials 3.0 - has taken away the argument N at the top...
int AbstractCvodeSystemJacAdaptor(realtype t, N_Vector y, N_Vector ydot, CHASTE_CVODE_DENSE_MATRIX jacobian,
#elif CHASTE_SUNDIALS_VERSION >= 20500
// Sundials 2.5
int AbstractCvodeSystemJacAdaptor(long int N, realtype t, N_Vector y, N_Vector ydot, CHASTE_CVODE_DENSE_MATRIX jacobian,
#elif CHASTE_SUNDIALS_VERSION >= 20400
// Sundials 2.4
int AbstractCvodeSystemJacAdaptor(int N, realtype t, N_Vector y, N_Vector ydot, DlsMat jacobian,
#else
// Sundials 2.3 and below (not sure how far below, but this is 2006 so old enough).
int AbstractCvodeSystemJacAdaptor(long int N, DenseMat jacobian, realtype t, N_Vector y, N_Vector ydot,
#endif
                                  void* pData, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
    assert(pData != nullptr);
    AbstractCvodeSystem* p_ode_system = (AbstractCvodeSystem*)pData;
    try
    {
        p_ode_system->EvaluateAnalyticJacobian(t, y, ydot, jacobian, tmp1, tmp2, tmp3);
    }
    catch (const Exception& e)
    {
        std::cerr << "CVODE Jacobian Exception: " << e.GetMessage() << std::endl
                  << std::flush;
        return -1;
    }
    return 0;
}

AbstractCvodeSystem::AbstractCvodeSystem(unsigned numberOfStateVariables)
        : AbstractParameterisedSystem<N_Vector>(numberOfStateVariables),
          mLastSolutionState(nullptr),
          mLastSolutionTime(0.0),
#if CHASTE_SUNDIALS_VERSION >= 20400
          mForceReset(false),
#else
          // Old Sundials don't seem to 'go back' when something has changed
          // properly, and give more inaccurate answers.
          mForceReset(true),
#endif
          mForceMinimalReset(false),
#if CHASTE_SUNDIALS_VERSION >= 30000
          mpSundialsDenseMatrix(nullptr),
          mpSundialsLinearSolver(nullptr),
#endif
          mHasAnalyticJacobian(false),
          mUseAnalyticJacobian(false),
          mpCvodeMem(nullptr),
          mMaxSteps(0),
          mLastInternalStepSize(0)
{
    SetTolerances(); // Set the tolerances to the defaults.
}

void AbstractCvodeSystem::Init()
{
    DeleteVector(mStateVariables);
    mStateVariables = GetInitialConditions();
    DeleteVector(mParameters);
    mParameters = N_VNew_Serial(rGetParameterNames().size());
    for (int i = 0; i < NV_LENGTH_S(mParameters); i++)
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
        // This should stop CVODE going past the end of where we wanted and interpolating back.
        int ierr = CVodeSetStopTime(mpCvodeMem, stepper.GetNextTime());
        assert(ierr == CV_SUCCESS);
        UNUSED_OPT(ierr); // avoid unused var warning

        //        // This parameter governs how many times we allow a recoverable right hand side failure
        //        int ierr = CVodeSetMaxConvFails(mpCvodeMem, 1000);
        //        assert(ierr == CV_SUCCESS); UNUSED_OPT(ierr); // avoid unused var warning

        double cvode_stopped_at = stepper.GetTime();
        ierr = CVode(mpCvodeMem, stepper.GetNextTime(), mStateVariables,
                     &cvode_stopped_at, CV_NORMAL);
        if (ierr < 0)
        {
            //            DebugSteps(mpCvodeMem, this);
            CvodeError(ierr, "CVODE failed to solve system", cvode_stopped_at, stepper.GetTime(), stepper.GetNextTime());
        }
        // Not root finding, so should have reached requested time
        assert(fabs(cvode_stopped_at - stepper.GetNextTime()) < DBL_EPSILON);
#ifndef NDEBUG
        VerifyStateVariables();
#endif
        // Store solution
        solutions.rGetSolutions().push_back(MakeStdVec(mStateVariables));
        solutions.rGetTimes().push_back(cvode_stopped_at);
        stepper.AdvanceOneTimeStep();
    }

    // stepper.EstimateTimeSteps may have been an overestimate...
    solutions.SetNumberOfTimeSteps(stepper.GetTotalTimeStepsTaken());

    int ierr = CVodeGetLastStep(mpCvodeMem, &mLastInternalStepSize);
    assert(ierr == CV_SUCCESS);
    UNUSED_OPT(ierr); // avoid unused var warning

    RecordStoppingPoint(tEnd);

    return solutions;
}

void AbstractCvodeSystem::Solve(realtype tStart,
                                realtype tEnd,
                                realtype maxDt)
{
    assert(tEnd >= tStart);

    SetupCvode(mStateVariables, tStart, maxDt);

    // This should stop CVODE going past the end of where we wanted and interpolating back.
    int ierr = CVodeSetStopTime(mpCvodeMem, tEnd);
    assert(ierr == CV_SUCCESS);
    UNUSED_OPT(ierr); // avoid unused var warning

    double cvode_stopped_at = tStart;
    ierr = CVode(mpCvodeMem, tEnd, mStateVariables, &cvode_stopped_at, CV_NORMAL);
    if (ierr < 0)
    {
        //        DebugSteps(mpCvodeMem, this);
        CvodeError(ierr, "CVODE failed to solve system", cvode_stopped_at, tStart, tEnd);
    }
    // Not root finding, so should have reached requested time
    assert(fabs(cvode_stopped_at - tEnd) < DBL_EPSILON);

    ierr = CVodeGetLastStep(mpCvodeMem, &mLastInternalStepSize);
    assert(ierr == CV_SUCCESS);
    UNUSED_OPT(ierr); // avoid unused var warning

    RecordStoppingPoint(cvode_stopped_at);

//
//    long int nst, nfe, nsetups, nje, nfeLS, nni, ncfn, netf, nge;
//
//
//    CVodeGetNumSteps(mpCvodeMem, &nst);
//    CVodeGetNumRhsEvals(mpCvodeMem, &nfe);
//    CVodeGetNumLinSolvSetups(mpCvodeMem, &nsetups);
//    CVodeGetNumErrTestFails(mpCvodeMem, &netf);
//    CVodeGetNumNonlinSolvIters(mpCvodeMem, &nni);
//    CVodeGetNumNonlinSolvConvFails(mpCvodeMem, &ncfn);
//    CVDlsGetNumJacEvals(mpCvodeMem, &nje);
//    CVDlsGetNumRhsEvals(mpCvodeMem, &nfeLS);
//    CVodeGetNumGEvals(mpCvodeMem, &nge);
//
//    printf("\nFinal Statistics:\n");
//    printf("nst = %-6ld nfe  = %-6ld nsetups = %-6ld nfeLS = %-6ld nje = %ld\n",
//       nst, nfe, nsetups, nfeLS, nje);
//    printf("nni = %-6ld ncfn = %-6ld netf = %-6ld nge = %ld\n \n",
//       nni, ncfn, netf, nge);
//    std::cout << std::flush;
#ifndef NDEBUG
    VerifyStateVariables();
#endif
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

void AbstractCvodeSystem::SetForceReset(bool autoReset)
{
    mForceReset = autoReset;
    if (mForceReset)
    {
        ResetSolver();
    }
}

bool AbstractCvodeSystem::GetMinimalReset()
{
    return mForceMinimalReset;
}

bool AbstractCvodeSystem::GetForceReset()
{
    return mForceReset;
}

void AbstractCvodeSystem::SetMinimalReset(bool minimalReset)
{
    mForceMinimalReset = minimalReset;
    if (mForceMinimalReset)
    {
        SetForceReset(false);
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
    //std::cout << "!mpCvodeMem = " << !mpCvodeMem << ", mForceReset = " << mForceReset << ", !mLastSolutionState = " << !mLastSolutionState << ", comp doubles = " << !CompareDoubles::WithinAnyTolerance(tStart, mLastSolutionTime) << "\n";
    bool reinit = !mpCvodeMem || mForceReset || !mLastSolutionState || !CompareDoubles::WithinAnyTolerance(tStart, mLastSolutionTime);
    if (!reinit && !mForceMinimalReset)
    {
        const unsigned size = GetNumberOfStateVariables();
        for (unsigned i = 0; i < size; i++)
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
        //std::cout << "New CVODE solver\n";
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

#if CHASTE_SUNDIALS_VERSION >= 30000
        /* Create dense matrix SUNDenseMatrix for use in linear solves */
        mpSundialsDenseMatrix = SUNDenseMatrix(NV_LENGTH_S(initialConditions), NV_LENGTH_S(initialConditions));
#endif

#if CHASTE_SUNDIALS_VERSION >= 40000
        /* Create dense SUNLinSol_Dense object for use by CVode */
        mpSundialsLinearSolver = SUNLinSol_Dense(initialConditions, mpSundialsDenseMatrix);

        /* Call CVodeSetLinearSolver to attach the matrix and linear solver to CVode */
        CVodeSetLinearSolver(mpCvodeMem, mpSundialsLinearSolver, mpSundialsDenseMatrix);
#elif CHASTE_SUNDIALS_VERSION >= 30000
        /* Create dense SUNDenseLinearSolver object for use by CVode */
        mpSundialsLinearSolver = SUNDenseLinearSolver(initialConditions, mpSundialsDenseMatrix);

        /* Call CVDlsSetLinearSolver to attach the matrix and linear solver to CVode */
        CVDlsSetLinearSolver(mpCvodeMem, mpSundialsLinearSolver, mpSundialsDenseMatrix);
#else
        // CVODE < v3.0.0
        // Attach a linear solver for Newton iteration
        CVDense(mpCvodeMem, NV_LENGTH_S(initialConditions));
#endif

        if (mUseAnalyticJacobian)
        {
#if CHASTE_SUNDIALS_VERSION >= 40000
            CVodeSetJacFn(mpCvodeMem, AbstractCvodeSystemJacAdaptor);
#elif CHASTE_SUNDIALS_VERSION >= 30000
            CVDlsSetJacFn(mpCvodeMem, AbstractCvodeSystemJacAdaptor);
#elif CHASTE_SUNDIALS_VERSION >= 20400
            CVDlsSetDenseJacFn(mpCvodeMem, AbstractCvodeSystemJacAdaptor);
#else
            CVDenseSetJacFn(mpCvodeMem, AbstractCvodeSystemJacAdaptor, (void*)(this));
#endif
        }
    }
    else if (reinit)
    {
//std::cout << "Resetting CVODE solver\n";
#if CHASTE_SUNDIALS_VERSION >= 20400
        CVodeReInit(mpCvodeMem, tStart, initialConditions);
        //CVodeSStolerances(mpCvodeMem, mRelTol, mAbsTol); - "all solver inputs remain in effect" so we don't need this.
#else
        CVodeReInit(mpCvodeMem, AbstractCvodeSystemRhsAdaptor, tStart, initialConditions,
                    CV_SS, mRelTol, &mAbsTol);
#endif
    }

    // Set max dt and change max steps if wanted
    if (maxDt > 0)
    {
        CVodeSetMaxStep(mpCvodeMem, maxDt);
    }

    if (mMaxSteps > 0)
    {
        CVodeSetMaxNumSteps(mpCvodeMem, mMaxSteps);
        CVodeSetMaxErrTestFails(mpCvodeMem, 15);
    }
}

void AbstractCvodeSystem::RecordStoppingPoint(double stopTime)
{
    //    DebugSteps(mpCvodeMem, this);

    // If we're forcing a reset then we don't record the stopping time
    // as a result it won't match and we will force a reset in SetupCvode() on
    // the next solve call.
    if (mForceReset)
        return;

    // Otherwise we will store the state variables and time for comparison on the
    // next solve call, to work out whether we need to reset.
    const unsigned size = GetNumberOfStateVariables();
    CreateVectorIfEmpty(mLastSolutionState, size);
    for (unsigned i = 0; i < size; i++)
    {
        SetVectorComponent(mLastSolutionState, i, GetVectorComponent(mStateVariables, i));
    }
    mLastSolutionTime = stopTime;
}

void AbstractCvodeSystem::FreeCvodeMemory()
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

void AbstractCvodeSystem::CvodeError(int flag, const char* msg,
                                     const double& rTime, const double& rStartTime, const double& rEndTime)
{
    std::stringstream err;
    char* p_flag_name = CVodeGetReturnFlagName(flag);
    err << msg << ": " << p_flag_name;
    free(p_flag_name);
    if (flag == CV_LSETUP_FAIL)
    {
#if CHASTE_SUNDIALS_VERSION >= 20500
        long int ls_flag;
#else
        int ls_flag;
#endif
        char* p_ls_flag_name;

#if CHASTE_SUNDIALS_VERSION >= 40000
        CVodeGetLastLinFlag(mpCvodeMem, &ls_flag);
        p_ls_flag_name = CVodeGetLinReturnFlagName(ls_flag);
#elif CHASTE_SUNDIALS_VERSION >= 20400
        CVDlsGetLastFlag(mpCvodeMem, &ls_flag);
        p_ls_flag_name = CVDlsGetReturnFlagName(ls_flag);
#else
        CVDenseGetLastFlag(mpCvodeMem, &ls_flag);
        p_ls_flag_name = CVDenseGetReturnFlagName(ls_flag);
#endif
        err << " (LS flag=" << ls_flag << ":" << p_ls_flag_name << ")";
        free(p_ls_flag_name);
    }

    err << "\nGot from time " << rStartTime << " to time " << rTime << ", was supposed to finish at time " << rEndTime << "\n";
    err << "\nState variables are now:\n";
    std::vector<double> state_vars = MakeStdVec(mStateVariables);
    std::vector<std::string> state_var_names = rGetStateVariableNames();
    for (unsigned i = 0; i < state_vars.size(); i++)
    {
        err << "\t" << state_var_names[i] << "\t:\t" << state_vars[i] << std::endl;
    }

    FreeCvodeMemory();
    std::cerr << err.str() << std::endl
              << std::flush;
    EXCEPTION(err.str());
}

bool AbstractCvodeSystem::HasAnalyticJacobian() const
{
    return mHasAnalyticJacobian;
}

bool AbstractCvodeSystem::GetUseAnalyticJacobian() const
{
    return mUseAnalyticJacobian;
}

void AbstractCvodeSystem::ForceUseOfNumericalJacobian(bool useNumericalJacobian)
{
    if (!useNumericalJacobian && !mHasAnalyticJacobian)
    {
        EXCEPTION("Analytic Jacobian requested, but this ODE system doesn't have one. You can check this with HasAnalyticJacobian().");
    }

    if (mUseAnalyticJacobian == useNumericalJacobian)
    {
        mUseAnalyticJacobian = !useNumericalJacobian;
        // We need to re-initialise the solver completely to change this.
        this->FreeCvodeMemory();
    }
}

//#include "MathsCustomFunctions.hpp"
//#include <algorithm>
//void AbstractCvodeSystem::CheckAnalyticJacobian(realtype time, N_Vector y, N_Vector ydot,
//                                                CHASTE_CVODE_DENSE_MATRIX jacobian,
//                                                N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
//{
//    N_Vector nudge_ydot = tmp1;
//    N_Vector numeric_jth_col = tmp2;
//    N_Vector ewt = tmp3;
//    const unsigned size = GetNumberOfStateVariables();
//    const double rel_tol = 1e-1;
//    const double abs_tol = 1e-6;
//    realtype* p_y = N_VGetArrayPointer(y);
//    realtype* p_numeric_jth_col = N_VGetArrayPointer(numeric_jth_col);
//
//    // CVODE internal data for computing the numeric J
//    realtype h;
//    CVodeGetLastStep(mpCvodeMem, &h);
//    CVodeGetErrWeights(mpCvodeMem, ewt);
//    realtype* p_ewt = N_VGetArrayPointer(ewt);
//    // Compute minimum nudge
//    realtype srur = sqrt(DBL_EPSILON);
//    realtype fnorm = N_VWrmsNorm(ydot, ewt);
//    realtype min_nudge = (fnorm != 0.0) ?
//            (1000.0 * fabs(h) * DBL_EPSILON * size * fnorm) : 1.0;
//
//    for (unsigned j=0; j<size; j++)
//    {
//        // Check the j'th column of the Jacobian
//        realtype yjsaved = p_y[j];
//        realtype nudge = std::max(srur*fabs(yjsaved), min_nudge/p_ewt[j]);
//        p_y[j] += nudge;
//        EvaluateYDerivatives(time, y, nudge_ydot);
//        p_y[j] = yjsaved;
//        realtype nudge_inv = 1.0 / nudge;
//        N_VLinearSum(nudge_inv, nudge_ydot, -nudge_inv, ydot, numeric_jth_col);
//        realtype* p_analytic_jth_col = DENSE_COL(jacobian, j);
//
//        for (unsigned i=0; i<size; i++)
//        {
//            if (!CompareDoubles::WithinAnyTolerance(p_numeric_jth_col[i], p_analytic_jth_col[i], rel_tol, abs_tol))
//            {
//                EXCEPTION("Analytic Jacobian appears dodgy at time " << time << " entry (" << i << "," << j << ").\n"
//                          << "Analytic=" << p_analytic_jth_col[i] << "; numeric=" << p_numeric_jth_col[i] << "."
//                          << DumpState("", y, time));
//            }
//        }
//    }
//}

#endif // CHASTE_CVODE
