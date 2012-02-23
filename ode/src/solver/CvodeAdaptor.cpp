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

#include "CvodeAdaptor.hpp"

#include "Exception.hpp"
#include "TimeStepper.hpp"
#include "VectorHelperFunctions.hpp"

#include <iostream>
#include <sstream>

// CVODE headers
#include <sundials/sundials_nvector.h>
#include <cvode/cvode_dense.h>


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
    assert(pData != NULL);
    CvodeData* p_data = (CvodeData*) pData;
    // Get y, ydot into std::vector<>s
    static std::vector<realtype> ydot_vec;
    CopyToStdVector(y, *p_data->pY);
    CopyToStdVector(ydot, ydot_vec);
    // Call our function
    try
    {
        p_data->pSystem->EvaluateYDerivatives(t, *(p_data->pY), ydot_vec);
    }
    catch (const Exception &e)
    {
        std::cerr << "CVODE RHS Exception: " << e.GetMessage() << std::endl << std::flush;
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
    assert(pData != NULL);
    CvodeData* p_data = (CvodeData*) pData;
    // Get y into a std::vector
    CopyToStdVector(y, *p_data->pY);
    // Call our function
    try
    {
        *pGOut = p_data->pSystem->CalculateRootFunction(t, *p_data->pY);
    }
    catch (const Exception &e)
    {
        std::cerr << "CVODE Root Exception: " << e.GetMessage() << std::endl << std::flush;
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



void CvodeErrorHandler(int errorCode, const char *module, const char *function,
                       char *message, void* pData)
{
    std::stringstream err;
    err << "CVODE Error " << errorCode << " in module " << module
        << " function " << function << ": " << message;
    std::cerr << "*" << err.str() << std::endl << std::flush;
    // Throwing the exception here causes termination on Maverick (g++ 4.4)
    //EXCEPTION(err.str());
}



void CvodeAdaptor::SetupCvode(AbstractOdeSystem* pOdeSystem,
                              std::vector<double>& rInitialY,
                              double startTime, double maxStep)
{
    assert(rInitialY.size() == pOdeSystem->GetNumberOfStateVariables());
    assert(maxStep > 0.0);

    mInitialValues = N_VMake_Serial(rInitialY.size(), &(rInitialY[0]));
    assert(NV_DATA_S(mInitialValues) == &(rInitialY[0]));
    assert(!NV_OWN_DATA_S(mInitialValues));
//    std::cout << " Initial values: "; N_VPrint_Serial(mInitialValues);
//    std::cout << " Rtol: " << mRelTol << ", Atol: " << mAbsTol << std::endl;
//    std::cout << " Start: " << startTime << " max dt=" << maxStep << std::endl << std::flush;

    mpCvodeMem = CVodeCreate(CV_BDF, CV_NEWTON);
    if (mpCvodeMem == NULL) EXCEPTION("Failed to SetupCvode CVODE");
    // Set error handler
    CVodeSetErrHandlerFn(mpCvodeMem, CvodeErrorHandler, NULL);
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
    CVodeInit(mpCvodeMem, CvodeRhsAdaptor, startTime, mInitialValues);
    CVodeSStolerances(mpCvodeMem, mRelTol, mAbsTol);
#else
    CVodeMalloc(mpCvodeMem, CvodeRhsAdaptor, startTime, mInitialValues,
                CV_SS, mRelTol, &mAbsTol);
#endif
    CVodeSetMaxStep(mpCvodeMem, maxStep);
    // Set the rootfinder function if wanted
    if (mCheckForRoots)
    {
#if CHASTE_SUNDIALS_VERSION >= 20400
        CVodeRootInit(mpCvodeMem, 1, CvodeRootAdaptor);
#else
        CVodeRootInit(mpCvodeMem, 1, CvodeRootAdaptor, (void*)(&mData));
#endif
    }
    // Change max steps if wanted
    if (mMaxSteps > 0)
    {
        CVodeSetMaxNumSteps(mpCvodeMem, mMaxSteps);
    }
    // Attach a linear solver for Newton iteration
    CVDense(mpCvodeMem, rInitialY.size());
}

void CvodeAdaptor::FreeCvodeMemory()
{
//    assert(!NV_OWN_DATA_STL(mInitialValues));
//    std::vector<double>* pVec = NV_DATA_STL(mInitialValues);
//    double val = (*pVec)[0];
    N_VDestroy_Serial(mInitialValues); mInitialValues = NULL;
//    std::cout << "  a: " << val << ", b: " << (*pVec)[0] << std::endl;

    CVodeFree(&mpCvodeMem);
}


void CvodeAdaptor::CvodeError(int flag, const char * msg)
{
    std::stringstream err;
    char* p_flag_name = CVodeGetReturnFlagName(flag);
    err << msg << ": " << p_flag_name;
    free(p_flag_name);
    std::cerr << err.str() << std::endl << std::flush;
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
    N_Vector yout = mInitialValues;

    // Set up ODE solution
    OdeSolution solutions;
    solutions.SetNumberOfTimeSteps(stepper.EstimateTimeSteps());
    solutions.rGetSolutions().push_back(rYValues);
    solutions.rGetTimes().push_back(startTime);
    solutions.SetOdeSystemInformation(pOdeSystem->GetSystemInformation());

    // Main time sampling loop
    while (!stepper.IsTimeAtEnd() && !mStoppingEventOccurred)
    {
        double tend;
        int ierr = CVode(mpCvodeMem, stepper.GetNextTime(), yout, &tend, CV_NORMAL);
        if (ierr<0)
        {
            FreeCvodeMemory();
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
        stepper.AdvanceOneTimeStep();
    }

    // stepper.EstimateTimeSteps may have been an overestimate...
    solutions.SetNumberOfTimeSteps(stepper.GetTotalTimeStepsTaken());

    int ierr = CVodeGetLastStep(mpCvodeMem, &mLastInternalStepSize);
    assert(ierr == CV_SUCCESS); ierr=ierr; // avoid unused var warning
    FreeCvodeMemory();

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

    N_Vector yout = mInitialValues;
    double tend;
    int ierr = CVode(mpCvodeMem, endTime, yout, &tend, CV_NORMAL);
    if (ierr<0)
    {
        FreeCvodeMemory();
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
    assert(ierr == CV_SUCCESS); ierr=ierr; // avoid unused var warning
    FreeCvodeMemory();
}

CvodeAdaptor::CvodeAdaptor(double relTol, double absTol)
    : AbstractIvpOdeSolver(),
      mpCvodeMem(NULL), mInitialValues(NULL),
      mRelTol(relTol), mAbsTol(absTol),
      mLastInternalStepSize(-0.0),
      mMaxSteps(0),
      mCheckForRoots(false)
{
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
