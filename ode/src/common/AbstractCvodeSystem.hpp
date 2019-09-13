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
#ifndef _ABSTRACTCVODESYSTEM_HPP_
#define _ABSTRACTCVODESYSTEM_HPP_

#include <algorithm>
#include <string>
#include <vector>

// This is only needed to prevent compilation errors on PETSc 2.2/Boost 1.33.1 combo
#include "UblasVectorInclude.hpp"

// Chaste includes
#include "AbstractParameterisedSystem.hpp"
#include "Exception.hpp"
#include "OdeSolution.hpp"
#include "VectorHelperFunctions.hpp"

// Serialiazation
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/vector.hpp>
#include "ChasteSerialization.hpp"
#include "ClassIsAbstract.hpp"

// CVODE headers
#include <nvector/nvector_serial.h>

#if CHASTE_SUNDIALS_VERSION >= 30000
#include <cvode/cvode_direct.h> /* access to CVDls interface            */
#include <sundials/sundials_types.h> /* defs. of realtype, sunindextype      */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#else
#include <sundials/sundials_dense.h> /* definitions DlsMat DENSE_ELEM */
#endif

// CVODE changed their dense matrix type...
#if CHASTE_SUNDIALS_VERSION >= 30000
#define CHASTE_CVODE_DENSE_MATRIX SUNMatrix
#elif CHASTE_SUNDIALS_VERSION >= 20400
#define CHASTE_CVODE_DENSE_MATRIX DlsMat
#else
#define CHASTE_CVODE_DENSE_MATRIX DenseMat
#endif

// CVODE changed their way of referencing elements of a matrix. So we will copy their notation in examples.
#if CHASTE_SUNDIALS_VERSION >= 30000
#define IJth(A, i, j) SM_ELEMENT_D(A, i, j)
#else
#define IJth(A, i, j) DENSE_ELEM(A, i, j)
#endif

/**
 * Abstract OdeSystem class for Cvode systems (N_Vector instead of std::vector)
 *
 * Sets up variables and functions for a general CVODE system, and uses an inbuilt CVODE solver.
 *
 * N.B. The CvodeAdaptor solver is not used by AbstractCvodeSystems, it is only used to translate when
 * the ODE system uses the normal std::vector<double> state variables.
 * For Chaste developers - this means there is some duplication of code between this class
 * and CvodeAdaptor as both set up and use CVODE to solve systems. If you change one you should
 * probably change both!
 *
 * ODE systems are specified primarily by the EvaluateYDerivatives() method,
 * which calculates the right-hand side of the system.
 *
 * You can also specify an EvaluateAnalyticJacobian() method which will
 * allow CVODE to evaluate the Jacobian matrix precisely instead of using
 * multiple calls to  EvaluateYDerivatives() to make an approximation to it.
 * This generally provides extra accuracy and better convergence, and so
 * gives a speed up for a given tolerance. Cardiac action potential
 * model Jacobians are calculated automatically by PyCML (see python/pycml)
 * via Maple for symbolic differentiation.
 *
 * Instances can store their state internally in the mStateVariables vector
 * in our base class AbstractParameterisedSystem (see also
 * GetNumberOfStateVariables(), SetStateVariables() and rGetStateVariables()),
 * although this is not essential - the vector may be empty, in which case
 * AbstractIvpOdeSolver::SolveAndUpdateStateVariable may not be used to
 * solve the system.
 *
 * CVODE systems may also have a vector of parameters, which can be accessed
 * through the GetParameter() and SetParameter() methods of our base class.
 *
 * Information about what the parameters and state variables represent is
 * provided by a subclass of AbstractOdeSystemInformation.  Various wrapper
 * methods (e.g. rGetStateVariableNames()) are provided in our base class to
 * access this information.
 *
 * Also, subclasses may define a condition at which ODE solvers should stop
 * prematurely. For this class CVODE solvers are being used, so
 * CalculateRootFunction() should be used to detect the stopping time.
 *
 * Note that the default tolerances for the solver are set by
 * SetTolerances(), these can make quite a difference to the time it takes
 * to solve the ODE system.
 *
 * Repeated calls to Solve will no longer set up and delete CVODE memory, unless
 * the following method is called:
 *
 * SetForceReset(true) - reset each time Solve() is called
 *
 * default - reset if state variables change, or we ask to solve from a different time than the last solve call finished.
 *
 * SetMinimalReset(true) - ignore changes in state vars and just reset if the time is inconsistent.
 *
 */
class AbstractCvodeSystem : public AbstractParameterisedSystem<N_Vector>
{
private:
    friend class TestAbstractCvodeSystem;

    friend class boost::serialization::access;
    /**
     * Archive the member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template <class Archive>
    void save(Archive& archive, const unsigned int version) const
    {
        // Despite the fact that 3 of these variables actually live in our base class,
        // we still archive them here to maintain backwards compatibility,
        // this doesn't hurt
        archive& mNumberOfStateVariables;
        archive& mUseAnalyticJacobian;

        if (version >= 1u)
        {
            archive& mHasAnalyticJacobian;
        }

        // Convert from N_Vector to std::vector for serialization
        const std::vector<double> state_vars = MakeStdVec(mStateVariables);
        archive& state_vars;
        const std::vector<double> params = MakeStdVec(mParameters);
        archive& params;
        archive& rGetParameterNames();

        archive& mLastSolutionTime;
        archive& mForceReset;
        archive& mForceMinimalReset;
        archive& mRelTol;
        archive& mAbsTol;
        archive& mMaxSteps;
        archive& mLastInternalStepSize;

        // We don't bother archiving CVODE's internal data, because it is missing then we'll just
        // get a new solver being initialised after a save/load.

        // This is always set up by subclass constructors, and is essentially
        // 'static' data, so shouldn't go in the archive.
        //archive &mpSystemInfo;
    }
    /**
     * Archive the member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template <class Archive>
    void load(Archive& archive, const unsigned int version)
    {
        archive& mNumberOfStateVariables;
        archive& mUseAnalyticJacobian;

        // This is pretty much what the code was saying before.
        mHasAnalyticJacobian = mUseAnalyticJacobian;
        if (version >= 1u)
        { // Overwrite if it has been archived though
            archive& mHasAnalyticJacobian;
        }

        std::vector<double> state_vars;
        archive& state_vars;
        CopyFromStdVector(state_vars, mStateVariables);

        std::vector<double> parameters;
        archive& parameters;

        std::vector<std::string> param_names;
        archive& param_names;
        archive& mLastSolutionTime;
        archive& mForceReset;
        archive& mForceMinimalReset;
        archive& mRelTol;
        archive& mAbsTol;
        archive& mMaxSteps;
        archive& mLastInternalStepSize;

        // We don't bother archiving CVODE's internal data, because it is missing then we'll just
        // get a new solver being initialised after a save/load.

        // Do some checking on the parameters
        CheckParametersOnLoad(parameters, param_names);
    }
    BOOST_SERIALIZATION_SPLIT_MEMBER()

    /**
     * Set up the CVODE data structures needed to solve the given system from a given point.
     *
     * @param initialConditions  initial conditions
     * @param tStart  start time of simulation
     * @param maxDt  maximum time step to take
     */
    void SetupCvode(N_Vector initialConditions,
                    realtype tStart,
                    realtype maxDt);

    /**
     * Record where the last solve got to so we know whether to re-initialise.
     * @param stopTime  the finishing time
     */
    void RecordStoppingPoint(double stopTime);

    /** Free CVODE memory when finished with. */
    void FreeCvodeMemory();

    /**
     * Report an error from CVODE.
     *
     * @param flag  CVODE error code
     * @param msg  Our description of the error
     * @param rTime  The time that the solver got to
     * @param rStartTime  The time that the solver started at.
     * @param rEndTime  The time that the solver was supposed to solve up to.
     *
     */
    void CvodeError(int flag, const char* msg, const double& rTime,
                    const double& rStartTime, const double& rEndTime);

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

protected:
    /** Whether we have an analytic Jacobian. */
    bool mHasAnalyticJacobian;

    /** Whether to use an analytic Jacobian. */
    bool mUseAnalyticJacobian;

    /** Relative tolerance for solver. */
    double mRelTol;

    /** Absolute tolerance for solver. */
    double mAbsTol;

    /** CVODE's internal data. */
    void* mpCvodeMem;

    /**
     * The maximum number of steps to be taken by the solver
     * in its attempt to reach the next output time.
     */
    long int mMaxSteps;

    /** The size of the previous timestep. */
    double mLastInternalStepSize;

    /**
     * @b Must be called by concrete subclass constructors to initialise the state
     * variables, after setting #mpSystemInfo.
     */
    void Init();

public:
    /**
     * Constructor.
     *
     * @param numberOfStateVariables  the number of state variables in the ODE system
     */
    AbstractCvodeSystem(unsigned numberOfStateVariables);

    /**
     * Virtual destructor since we have virtual methods.
     */
    virtual ~AbstractCvodeSystem();

    /**
     * Method to evaluate the derivatives of the system.
     *
     * @param time  the current time
     * @param y  the current values of the state variables
     * @param ydot  storage for the derivatives of the system; will be filled in on return
     */
    virtual void EvaluateYDerivatives(realtype time,
                                      const N_Vector y,
                                      N_Vector ydot)
        = 0;

    /**
     * This method is called by AbstractCvodeSystemJacAdaptor method in the .cpp file.
     *
     * It provides an interface between the methods different versions of CVODE are expecting and the
     * Jacobians provided by Chaste CVODE systems.
     *
     * @param time  the current time (used by ODE systems like y' = f(t,y) only I guess)
     * @param y  the current state variables y for y' = f(t,y)
     * @param ydot  the current set of derivatives y' = f(t,y)
     * @param jacobian  a pointer to a jacobian, populated by this method.
     * @param tmp1  working memory of the correct size provided by CVODE for temporary calculations
     * @param tmp2  working memory of the correct size provided by CVODE for temporary calculations
     * @param tmp3  working memory of the correct size provided by CVODE for temporary calculations
     *
     */
    virtual void EvaluateAnalyticJacobian(realtype time, N_Vector y, N_Vector ydot,
                                          CHASTE_CVODE_DENSE_MATRIX jacobian,
                                          N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
    {
        EXCEPTION("No analytic Jacobian has been defined for this system.");
    }

    /**
     * Set whether to automatically re-initialise CVODE on every call to Solve, or
     * whether to attempt to guess when re-initialisation is needed. For example
     * it will re-initialise if the time changes, or any state variables change.
     * You can safely change parameters between solve calls with or without resets.
     *
     * See also ResetSolver and SetMinimalReset
     *
     * @param autoReset  whether to reset on every Solve
     */
    void SetForceReset(bool autoReset);

    /**
     * Set whether to reduce the checking done when guessing when re-initialisation
     * is needed, so it ignores changes in the state variables. If called with true
     * argument, will call SetForceReset(false).
     * You can safely change parameters between solve calls with or without resets.
     *
     * @param minimalReset  whether to avoid checking for changes in state variables
     */
    void SetMinimalReset(bool minimalReset);

    /**
     * @brief Get whether we want to run with minimal reset or not (no reinitialisation of the solver if variables change)
     *
     * @return true We are not resetting the solver if time or variables change between solve calls.
     * @return false We will reset the solver if time or variables change between solve calls.
     */
    bool GetMinimalReset();

    /**
     * @brief Get whether we will force a solver reset on every call to Solve()
     *
     * @return true We will reinitialise the solver at the start of every Solve() call.
     * @return false We will not reinitialise the solver at the start of every Solve() call.
     */
    bool GetForceReset();

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
     * Simulate the cell, returning a sampling of the state variables.
     *
     * Uses the current values of the state variables at initial conditions.
     * If the state variables have not been set (either by a prior solve, or
     * a call to SetStateVariables) the initial conditions (given by
     * GetInitialConditions) will be used.
     *
     * The final values of the state variables will also be stored in this object.
     *
     * @note See also the ResetSolver method.
     *
     * @param tStart  start time of simulation
     * @param tEnd  end time of simulation
     * @param maxDt  maximum time step to be taken by the adaptive solver
     *   (set this appropriately to avoid missing a stimulus)
     * @param tSamp  sampling interval at which to store results
     * @return the solution object with results at each sampling step
     */
    OdeSolution Solve(realtype tStart,
                      realtype tEnd,
                      realtype maxDt,
                      realtype tSamp);

    /**
     * Simulate the cell, updating its internal state variables.
     *
     * Uses the current values of the state variables at initial conditions.
     * If the state variables have not been set (either by a prior solve, or
     * a call to SetStateVariables) the initial conditions (given by
     * GetInitialConditions) will be used.
     *
     * @note See also the ResetSolver method.
     *
     * @param tStart  start time of simulation
     * @param tEnd  end time of simulation
     * @param maxDt  maximum time step to be taken by the adaptive solver
     *   (set this appropriately to avoid missing a stimulus)
     */
    void Solve(realtype tStart,
               realtype tEnd,
               realtype maxDt);

    /**
     * Change the maximum number of steps to be taken by the solver
     * in its attempt to reach the next output time.  Default is 500 (set by CVODE).
     *
     * @param numSteps new maximum
     */
    void SetMaxSteps(long int numSteps);

    /**
     * @return the maximum number of steps to be taken by the solver
     * in its attempt to reach the next output time.
     */
    long int GetMaxSteps();

    /**
     * Set relative and absolute tolerances; both scalars.
     * If no parameters are given, tolerances will be reset to default values.
     *
     * @param relTol  the relative tolerance for the solver (defaults to 1e-5)
     * @param absTol  the absolute tolerance for the solver (defaults to 1e-7)
     */
    void SetTolerances(double relTol = 1e-5, double absTol = 1e-7);

    /**
     * @return the relative tolerance.
     */
    double GetRelativeTolerance();

    /**
     * @return the absolute tolerance.
     */
    double GetAbsoluteTolerance();

    /**
     * @return the last step size used internally by CVODE in the last Solve call.
     */
    double GetLastStepSize();

    /* NB This needs making into a doxygen comment if you bring the method back in.
     *
     * An alternative approach to stopping events; currently only useful with CVODE.
     * CVODE can search for roots (zeros) of this function while solving the ODE system,
     * and home in on them to find sign transitions to high precision.
     *
     * The default implementation here fakes a root function using CalculateStoppingEvent.
     *
     * @param time  the current time
     * @param rY  the current values of the state variables
     */
    //    virtual double CalculateRootFunction(double time, const std::vector<double>& rY);

    /**
     * @return whether an analytic Jacobian is used (#mUseAnalyticJacobian).
     */
    bool GetUseAnalyticJacobian() const;

    /**
     * @return whether the ODE system has an analytic Jacobian (#mHasAnalyticJacobian).
     */
    bool HasAnalyticJacobian() const;

    /**
     * Force the use of a numerical Jacobian, even if an analytic form is provided.
     * This is needed for a handful of troublesome models.
     *
     * @param useNumericalJacobian  Whether to use a numerical instead of the analytic Jacobian.
     */
    void ForceUseOfNumericalJacobian(bool useNumericalJacobian = true);

    // The following method may be useful to identify problems with the Analytic Jacobians, if anything goes wrong,
    // but #1795 seems to have got these working OK, so commented out for now.

    //    /*
    //     * Compare the calculated analytic jacobian to a numerical approximation, and throw if it looks silly.
    //     *
    //     * @param time  the current time
    //     * @param y  the current state variables
    //     * @param jacobian  the analytic jacobian matrix
    //     * @param tmp1  working memory of the correct size provided by CVODE for temporary calculations
    //     * @param tmp2  working memory of the correct size provided by CVODE for temporary calculations
    //     * @param tmp3  working memory of the correct size provided by CVODE for temporary calculations
    //     */
    //    void CheckAnalyticJacobian(realtype time, N_Vector y, N_Vector ydot,
    //                               CHASTE_CVODE_DENSE_MATRIX jacobian,
    //                               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
};

CLASS_IS_ABSTRACT(AbstractCvodeSystem)
BOOST_CLASS_VERSION(AbstractCvodeSystem, 1u)

#endif //_ABSTRACTCVODESYSTEM_HPP_
#endif // CHASTE_CVODE
