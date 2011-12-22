/*

Copyright (C) University of Oxford, 2005-2011

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

#ifdef CHASTE_CVODE
#ifndef _ABSTRACTCVODESYSTEM_HPP_
#define _ABSTRACTCVODESYSTEM_HPP_

#include <vector>
#include <string>
#include <algorithm>

// This is only needed to prevent compilation errors on PETSc 2.2/Boost 1.33.1 combo
#include "UblasVectorInclude.hpp"

// Chaste includes
#include "OdeSolution.hpp"
#include "AbstractParameterisedSystem.hpp"
#include "Exception.hpp"

// Serialiazation
//#include "ChasteSerialization.hpp"
//#include <boost/serialization/split_member.hpp>
//#include <boost/serialization/vector.hpp>
//#include <boost/serialization/version.hpp>
//#include "ClassIsAbstract.hpp"

// CVODE headers
#include <nvector/nvector_serial.h>

/**
 * Abstract OdeSystem class for Cvode systems (N_Vector instead of std::vector)
 *
 * Sets up variables and functions for a general CVODE system.
 *
 * ODE systems are specified primarily by the EvaluateYDerivatives() method,
 * which calculates the right-hand side of the system.
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
 */
class AbstractCvodeSystem : public AbstractParameterisedSystem<N_Vector>
{
private:
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
     */
    void CvodeError(int flag, const char * msg);

    /** Remember where the last solve got to so we know whether to re-initialise. */
    N_Vector mLastSolutionState;

    /** Remember where the last solve got to so we know whether to re-initialise. */
    double mLastSolutionTime;

    /** Whether to automatically reset CVODE on each Solve call. */
    bool mAutoReset;

protected:

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
                                      N_Vector y,
                                      N_Vector ydot)=0;

    /**
     * Set whether to automatically re-initialise CVODE on every call to Solve, or
     * whether to attempt to guess when re-initialisation is needed.  See also
     * ResetSolver.
     *
     * @param autoReset  whether to reset on every Solve
     */
    void SetAutoReset(bool autoReset);

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
     * Get the maximum number of steps to be taken by the solver
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
    void SetTolerances(double relTol=1e-5, double absTol=1e-7);

    /**
     * Get the relative tolerance.
     */
    double GetRelativeTolerance();

    /**
     * Get the absolute tolerance.
     */
    double GetAbsoluteTolerance();

    /**
     * Get the last step size used internally by CVODE in the last Solve call.
     */
    double GetLastStepSize();


//    /**
//     * An alternative approach to stopping events; currently only useful with CVODE.
//     * CVODE can search for roots (zeros) of this function while solving the ODE system,
//     * and home in on them to find sign transitions to high precision.
//     *
//     * The default implementation here fakes a root function using CalculateStoppingEvent.
//     *
//     * @param time  the current time
//     * @param rY  the current values of the state variables
//     */
//    virtual double CalculateRootFunction(double time, const std::vector<double>& rY);
//
//    /**
//     * Get whether an analytic Jacobian is used.
//     *
//     * @return mUseAnalyticJacobian
//     */
//    bool GetUseAnalyticJacobian();

//    /**
//     * \todo move to AbstractParameterisedSystem? (1540)
//     *
//     * @return const reference to the state variables in the ODE system (used in archiving).
//     */
//    const std::vector<double>& rGetConstStateVariables() const;

};

//CLASS_IS_ABSTRACT(AbstractCvodeSystem)

#endif //_ABSTRACTCVODESYSTEM_HPP_
#endif // CHASTE_CVODE


