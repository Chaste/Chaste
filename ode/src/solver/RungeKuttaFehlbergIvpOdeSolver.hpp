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

#ifndef _RUNGEKUTTAFEHLBERGIVPODESOLVER_HPP_
#define _RUNGEKUTTAFEHLBERGIVPODESOLVER_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractOneStepIvpOdeSolver.hpp"

/**
 * A concrete one step ODE solver class that employs the Runge Kutta
 * Fehlberg adaptive solver (RKF45).
 *
 * This solver is good for problems where you need to be able to
 * guarantee the accuracy of the answer as it is specified via the
 * tolerance parameter.
 *
 * The solver should also be reasonably fast as it increases the
 * timestep when the solutions are changing slowly, whilst maintaining
 * accuracy.
 */
class RungeKuttaFehlbergIvpOdeSolver : public AbstractIvpOdeSolver
{
friend class TestRungeKuttaFehlbergIvpOdeSolver;

private:
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive, never used directly - boost uses this.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // This calls serialize on the base class - all member variables instantiated on construction or temporary.
        archive & boost::serialization::base_object<AbstractIvpOdeSolver>(*this);
    }

    /*
     * All these are here for more efficient memory allocation, rather than
     * because they need to be member variables.
     */

    double m1932o2197;  /**< Working memory: numerical value for the fraction 1932/2197.  */
    double m7200o2197;  /**< Working memory: numerical value for the fraction 7200/2197.  */
    double m7296o2197;  /**< Working memory: numerical value for the fraction 7296/2197.  */
    double m12o13;      /**< Working memory: numerical value for the fraction 12/13.      */
    double m439o216;    /**< Working memory: numerical value for the fraction 439/216.    */
    double m3680o513;   /**< Working memory: numerical value for the fraction 3680/513.   */
    double m845o4104;   /**< Working memory: numerical value for the fraction 845/4104.   */
    double m8o27;       /**< Working memory: numerical value for the fraction 8/27.       */
    double m3544o2565;  /**< Working memory: numerical value for the fraction 3544/2565.  */
    double m1859o4104;  /**< Working memory: numerical value for the fraction 1859/4104.  */
    double m1o360;      /**< Working memory: numerical value for the fraction 1/360.      */
    double m128o4275;   /**< Working memory: numerical value for the fraction 128/4275.   */
    double m2197o75240; /**< Working memory: numerical value for the fraction 2197/75240. */
    double m2o55;       /**< Working memory: numerical value for the fraction 2/55.       */
    double m25o216;     /**< Working memory: numerical value for the fraction 25/216.     */
    double m1408o2565;  /**< Working memory: numerical value for the fraction 1408/2565.  */
    double m2197o4104;  /**< Working memory: numerical value for the fraction 2197/4104.  */

    std::vector<double> mError; /**< Error expression, used to adjust the timestep in the RKF45 method. */

    std::vector<double> mk1;  /**< Working memory: expression k1 in the RKF45 method.  */
    std::vector<double> mk2;  /**< Working memory: expression k2 in the RKF45 method.  */
    std::vector<double> mk3;  /**< Working memory: expression k3 in the RKF45 method.  */
    std::vector<double> mk4;  /**< Working memory: expression k4 in the RKF45 method.  */
    std::vector<double> mk5;  /**< Working memory: expression k5 in the RKF45 method.  */
    std::vector<double> mk6;  /**< Working memory: expression k6 in the RKF45 method.  */
    std::vector<double> myk2; /**< Working memory: expression yk2 in the RKF45 method. */
    std::vector<double> myk3; /**< Working memory: expression yk3 in the RKF45 method. */
    std::vector<double> myk4; /**< Working memory: expression yk4 in the RKF45 method. */
    std::vector<double> myk5; /**< Working memory: expression yk5 in the RKF45 method. */
    std::vector<double> myk6; /**< Working memory: expression yk6 in the RKF45 method. */

protected:

    /**
     * Method that actually performs the solving on behalf of the public Solve methods.
     *
     * @param rSolution  an ODE solution to input data into if requited
     * @param pAbstractOdeSystem  the ODE system to solve
     * @param rCurrentYValues  the current (initial) state; results will also be returned in here
     * @param rWorkingMemory  working memory; same size as rCurrentYValues
     * @param startTime  initial time
     * @param endTime  time to solve to
     * @param maxTimeStep  the maximum size of timestep allowable
     * @param minTimeStep  the maximum size of timestep allowable (to prevent huge loops)
     * @param tolerance  how accurate the numerical solution must be
     * @param outputSolution whether to output into rSolution (or save time by not doing)
     */
    void InternalSolve(OdeSolution& rSolution,
                       AbstractOdeSystem* pAbstractOdeSystem,
                       std::vector<double>& rCurrentYValues,
                       std::vector<double>& rWorkingMemory,
                       double startTime,
                       double endTime,
                       double maxTimeStep,
                       double minTimeStep,
                       double tolerance,
                       bool outputSolution);

    /**
     * Calculate the solution to the ODE system at the next timestep.
     * Updates the mError vector with current error.
     *
     * @param pAbstractOdeSystem  the ODE system to solve
     * @param timeStep  dt
     * @param time  the current time
     * @param rCurrentYValues  the current (initial) state
     * @param rNextYValues  the state at the next timestep
     */
    void CalculateNextYValue(AbstractOdeSystem* pAbstractOdeSystem,
                             double timeStep,
                             double time,
                             std::vector<double>& rCurrentYValues,
                             std::vector<double>& rNextYValues);

    /**
     * Use the error approximation of the last call to the CalculateNextYValue()
     * method to change the time step appropriately.
     *
     * @param rCurrentStepSize  the current step size being used (returns answer via this reference)
     * @param rError  the error in the approximation at this time step
     * @param rTolerance  the tolerance required
     * @param rMaxTimeStep  the maximum timestep to be used
     * @param rMinTimeStep  the minimum timestep to be used (to prevent huge loops)
     */
    void AdjustStepSize(double& rCurrentStepSize,
                        const double& rError,
                        const double& rTolerance,
                        const double& rMaxTimeStep,
                        const double& rMinTimeStep);

public:

    /**
     * Constructor.
     */
    RungeKuttaFehlbergIvpOdeSolver();

    /**
     * Solves a system of ODEs using a specified one-step ODE solver and returns
     * the solution as an OdeSolution object.
     *
     * @param pAbstractOdeSystem  pointer to the concrete ODE system to be solved
     * @param rYValues  a standard vector specifying the intial condition of each
     *                  solution variable in the system (this can be the initial
     *                  conditions vector stored in the ODE system)
     * @param startTime  the time at which the initial conditions are specified
     * @param endTime  the time to which the system should be solved and the solution
     *                 returned
     * @param timeStep  the time interval to be used by the solver
     * @param ignoredSamplingTime  the interval at which to sample the solution to the ODE system
     *                             (ignored in this class as the timestep is variable)
     *
     * @return OdeSolution is an object containing an integer of the number of
     * equations, a stdAbstractOdeSystem::vector of times and a std::vector of std::vectors where
     * each of those vectors contains the solution for one variable of the ODE
     * system at those times.
     */
    OdeSolution Solve(AbstractOdeSystem* pAbstractOdeSystem,
                      std::vector<double>& rYValues,
                      double startTime,
                      double endTime,
                      double timeStep,
                      double ignoredSamplingTime);

    /**
     * Second version of Solve. Solves a system of ODEs using a specified one-step
     * ODE solver. This method does not return the solution and therefore does not
     * take in a sampling time. Instead, the mStateVariables component in the ODE
     * system object is updated.
     *
     * @param pAbstractOdeSystem  pointer to the concrete ODE system to be solved
     * @param rYValues  a standard vector specifying the intial condition of each
     *                  solution variable in the system (this can be the initial
     *                  conditions vector stored in the ODE system)
     * @param startTime  the time at which the initial conditions are specified
     * @param endTime  the time to which the system should be solved and the solution
     *                 returned
     * @param timeStep  the time interval to be used by the solver
     */
    void Solve(AbstractOdeSystem* pAbstractOdeSystem,
               std::vector<double>& rYValues,
               double startTime,
               double endTime,
               double timeStep);
};

#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(RungeKuttaFehlbergIvpOdeSolver)

#endif //_RUNGEKUTTAFEHLBERGIVPODESOLVER_HPP_
