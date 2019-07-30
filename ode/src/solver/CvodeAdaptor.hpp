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
#ifndef _CVODEADAPTOR_HPP_
#define _CVODEADAPTOR_HPP_

#include <vector>

#include <boost/serialization/base_object.hpp>
#include "ChasteSerialization.hpp"

#include "AbstractIvpOdeSolver.hpp"
#include "OdeSolution.hpp"

// CVODE headers
#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>

#if CHASTE_SUNDIALS_VERSION >= 30000
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#endif

/**
 * CVODE error handling function.
 *
 * Throw an Exception to report errors, rather than the CVODE approach of magic
 * return codes.
 */
void CvodeErrorHandler(int errorCode, const char* module, const char* function,
                       char* message, void* pData);
// Note: declared here since it's also used by AbstractCvodeCell.

/**
 * Data structure passed to CVODE calls, allowing our callback functions
 * to access the Chaste objects.
 */
typedef struct CvodeData_
{
    /** Working memory. */
    std::vector<realtype>* pY;
    /** The ODE system being solved. */
    AbstractOdeSystem* pSystem;
} CvodeData;

/**
 * The CVODE adaptor ODE solver class. This class is used as a solver for Chaste standard
 * AbstractOdeSystems that have a native vector type of std::vector<double> and so can be used
 * to apply CVODE to systems that could also be solved with the other solvers in this folder.
 *
 * N.B. If you know you are always going to want to use CVODE to solve your system,
 * then it will be faster to write your system so it uses N_Vectors natively as an AbstractCvodeSystem.
 * AbstractCvodeSystems have an in-built CVODE solver, and therefore don't need to use this class.
 * For Chaste developers - this means there is some duplication of code between this class
 * and AbstractCvodeSystem as both set up and use CVODE to solve systems. If you change one you should
 * probably change both!
 *
 * This class assumes that it will be solving stiff systems, so uses BDF/Newton.
 *
 * The timeStep parameters of the abstract class are here used to specify
 * *maximum* steps, since the solver is adaptive.
 *
 * Repeated calls to Solve will no longer set up and delete CVODE memory, unless
 * the following method is called:
 *
 * SetForceReset(true) - reset each time Solve() is called
 *
 * default behaviour - reset if state variables change, or we ask to solve from a different time than the last solve call finished.
 *
 * SetMinimalReset(true) - ignore changes in state vars and just reset if the time is inconsistent.
 */
class CvodeAdaptor : public AbstractIvpOdeSolver
{
private:
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive, never used directly - boost uses this.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template <class Archive>
    void serialize(Archive& archive, const unsigned int version)
    {
        archive& boost::serialization::base_object<AbstractIvpOdeSolver>(*this);
        archive& mRelTol;
        archive& mAbsTol;
        archive& mLastInternalStepSize;
        archive& mMaxSteps;
        archive& mCheckForRoots;
        // All other member variables given values on each call.
    }

    /** Pointer to the CVODE memory block. */
    void* mpCvodeMem;

    /** The CVODE data structure. */
    CvodeData mData;

    /** Relative tolerance for the ODE solver. */
    double mRelTol;

    /** Absolute tolerance for the ODE solver. */
    double mAbsTol;

    /** The size of the previous timestep. */
    double mLastInternalStepSize;

    /**
     * The maximum number of steps to be taken by the solver
     * in its attempt to reach the next output time.
     */
    long int mMaxSteps;

    /** Whether to check for stopping events. */
    bool mCheckForRoots;

    /** Remember where the last solve got to so we know whether to re-initialise. */
    N_Vector mLastSolutionState;

    /** Remember where the last solve got to so we know whether to re-initialise. */
    double mLastSolutionTime;

    /** Whether to automatically reset CVODE on each Solve call. */
    bool mForceReset;

    /** Whether to ignore changes in the state variables when deciding whether to reset. */
    bool mForceMinimalReset;

#if CHASTE_SUNDIALS_VERSION >= 30000
    /** Working memory for CVODE to store a dense matrix */
    SUNMatrix mpSundialsDenseMatrix;
    /** Working memory for CVODE's linear solver */
    SUNLinearSolver mpSundialsLinearSolver;
#endif

    /**
     * Record where the last solve got to so we know whether to re-initialise.
     * @param stopTime  the finishing time
     * @param yEnd  the state variables at this time
     */
    void RecordStoppingPoint(double stopTime, N_Vector yEnd);

protected:
    /**
     * Set up the CVODE data structures needed to solve the given system.
     *
     * @param pOdeSystem  the ODE system being solved
     * @param rInitialY  initial conditions vector
     * @param startTime  when to simulate from
     * @param maxStep  maximum time step
     */
    void SetupCvode(AbstractOdeSystem* pOdeSystem,
                    std::vector<double>& rInitialY,
                    double startTime, double maxStep);

    /**
     * Free CVODE memory after solving a system of ODEs.
     */
    void FreeCvodeMemory();

    /**
     * Report an error from CVODE.
     *
     * This will (probably) never be called, since we supply an error handler function
     * which throws an exception.
     *
     * @param flag  error flag
     * @param msg  error message
     */
    void CvodeError(int flag, const char* msg);

public:
    /**
     * Default constructor.
     * Can optionally set relative and absolute tolerances.
     *
     * @param relTol the relative tolerance for the solver
     * @param absTol the absolute tolerance for the solver
     */
    CvodeAdaptor(double relTol = 1e-4, double absTol = 1e-6);

    /**
     * Destructor.
     * Frees memory.
     */
    ~CvodeAdaptor();

    /**
     * Set relative and absolute tolerances; both scalars.
     * If no parameters are given, tolerances will be reset to default values.
     *
     * @param relTol the relative tolerance for the solver
     * @param absTol the absolute tolerance for the solver
     */
    void SetTolerances(double relTol = 1e-4, double absTol = 1e-6);

    /**
     * @return the relative tolerance.
     */
    double GetRelativeTolerance();

    /**
     * @return the absolute tolerance.
     */
    double GetAbsoluteTolerance();

    /**
     * @return the last step size used internally by CVODE in the last Solve call
     */
    double GetLastStepSize();

    /**
     * Set whether to automatically re-initialise CVODE on every call to Solve, or
     * whether to attempt to guess when re-initialisation is needed.  See also
     * ResetSolver.
     *
     * @param autoReset  whether to reset on every Solve
     */
    void SetForceReset(bool autoReset);

    /**
     * Set whether to reduce the checking done when guessing when re-initialisation
     * is needed, so it ignores changes in the state variables.  If call with true
     * argument, will call SetForceReset(false).
     *
     * @param minimalReset  whether to avoid checking for changes in state variables
     */
    void SetMinimalReset(bool minimalReset);

    /**
     * Successive calls to Solve will attempt to intelligently determine whether
     * to re-initialise the internal CVODE solver, or whether we are simply
     * extending the previous solution forward in time.  This mechanism compares
     * the state vector to its previous value, and the start time to the end of
     * the last solve, which captures most cases where re-initialisation is
     * required.  However, changes to the RHS function can also require this, and
     * cannot be automatically detected.  In such cases users must call this
     * function to force re-initialisation.
     */
    void ResetSolver();

    /**
     * Solve the given ODE system, returning the solution at sampling intervals.
     *
     * @param pOdeSystem  the ODE system to solve
     * @param rYValues  the initial state variable values
     *   (note: this vector will also be used as working memory)
     * @param startTime  the time to start solving at
     * @param endTime  the time to solve to
     * @param maxStep  the maximum time step to be taken by the adaptive solver
     * @param timeSampling  the interval at which to sample the solution
     * @return  the solution
     */
    OdeSolution Solve(AbstractOdeSystem* pOdeSystem,
                      std::vector<double>& rYValues,
                      double startTime,
                      double endTime,
                      double maxStep,
                      double timeSampling);

    /**
     * Solve the given ODE system, storing the final result in rYValues.
     *
     * @param pOdeSystem  the ODE system to solve
     * @param rYValues  the initial state variable values; will be filled in with
     *   the final values on return
     * @param startTime  the time to start solving at
     * @param endTime  the time to solve to
     * @param maxStep  the maximum time step to be taken by the adaptive solver
     */
    void Solve(AbstractOdeSystem* pOdeSystem,
               std::vector<double>& rYValues,
               double startTime,
               double endTime,
               double maxStep);

    /**
     * Make the solver check for stopping events using CVODE's rootfinding functionality.
     *
     * By default we do not check.
     */
    void CheckForStoppingEvents();

    /**
     * Change the maximum number of steps to be taken by the solver
     * in its attempt to reach the next output time.  Default is 500.
     *
     * @param numSteps  new maximum number of steps
     */
    void SetMaxSteps(long int numSteps);

    /**
     * @return the maximum number of steps to be taken by the solver
     * in its attempt to reach the next output time.
     */
    long int GetMaxSteps();
};

#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(CvodeAdaptor)

#endif // _CVODEADAPTOR_HPP_
#endif // CHASTE_CVODE
