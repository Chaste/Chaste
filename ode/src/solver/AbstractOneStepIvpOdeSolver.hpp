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

#ifndef _ABSTRACTONESTEPIVPODESOLVER_HPP_
#define _ABSTRACTONESTEPIVPODESOLVER_HPP_

#include "ChasteSerialization.hpp"
#include "ClassIsAbstract.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractIvpOdeSolver.hpp"

/**
 * Abstract one-step initial value problem ODE solver class. Sets up variables and functions
 * for all the ODE solvers that only have one timestep.
*/
class AbstractOneStepIvpOdeSolver : public AbstractIvpOdeSolver
{

private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the abstract IVP Solver, never used directly - boost uses this.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // This calls serialize on the base class.
        archive & boost::serialization::base_object<AbstractIvpOdeSolver>(*this);
    }

    /**
     * Working memory
     */
    std::vector<double> mWorkingMemory;

protected:

    /**
     * Method that actually performs the solving on behalf of the public Solve methods.
     *
     * @param pAbstractOdeSystem  the ODE system to solve
     * @param rCurrentYValues  the current (initial) state; results will also be returned
     *                         in here
     * @param rWorkingMemory  working memory; same size as rCurrentYValues
     * @param startTime  initial time
     * @param endTime  time to solve to
     * @param timeStep  dt
     */
    virtual void InternalSolve(AbstractOdeSystem* pAbstractOdeSystem,
                               std::vector<double>& rCurrentYValues,
                               std::vector<double>& rWorkingMemory,
                               double startTime,
                               double endTime,
                               double timeStep);

    /**
     * Calculate the solution to the ODE system at the next timestep.
     * Concrete subclasses should provide this method.
     *
     * @param pAbstractOdeSystem  the ODE system to solve
     * @param timeStep  dt
     * @param time  the current time
     * @param rCurrentYValues  the current (initial) state
     * @param rNextYValues  the state at the next timestep
     */
    virtual void CalculateNextYValue(AbstractOdeSystem* pAbstractOdeSystem,
                                     double timeStep,
                                     double time,
                                     std::vector<double>& rCurrentYValues,
                                     std::vector<double>& rNextYValues)=0;

public:

    /**
     * Solves a system of ODEs using a specified one-step ODE solver and returns
     * the solution as an OdeSolution object.
     *
     * An example, which returns the solution every 0.1 seconds using dt=0.01:
     *
     *     MyOdeSystem ode;
     *     std::vector<double> init_cond = ode_system.GetInitialConditions();
     *     OdeSolution solution = solver.Solve(&ode, init_cond, 0, 1, 0.01, 0.1);
     *
     *
     * @param pAbstractOdeSystem  pointer to the concrete ODE system to be solved
     * @param rYValues  a standard vector specifying the intial condition of each
     *                  solution variable in the system (this can be the initial
     *                  conditions vector stored in the ODE system)
     * @param startTime  the time at which the initial conditions are specified
     * @param endTime  the time to which the system should be solved and the solution
     *                 returned
     * @param timeStep  the time interval to be used by the solver
     * @param timeSampling  the interval at which to sample the solution to the ODE system
     *
     * @return OdeSolution is an object containing an integer of the number of
     * equations, a stdAbstractOdeSystem::vector of times and a std::vector of std::vectors where
     * each of those vectors contains the solution for one variable of the ODE
     * system at those times.
     */
    virtual OdeSolution Solve(AbstractOdeSystem* pAbstractOdeSystem,
                              std::vector<double>& rYValues,
                              double startTime,
                              double endTime,
                              double timeStep,
                              double timeSampling);

    /**
     * Second version of Solve. Solves a system of ODEs using a specified one-step
     * ODE solver. This method does not return the solution and therefore does not
     * take in a sampling time. Instead, the mStateVariables component in the ODE
     * system object is updated.
     *
     * An example:
     *
     *     std::vector<double> init_cond = ode_system.GetInitialConditions();
     *     solver.Solve(&ode, init_cond, 0, 1, 0.01);
     *     state_variables = ode_system.rGetStateVariables(); // solution at t=1 found here
     *
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
    virtual void Solve(AbstractOdeSystem* pAbstractOdeSystem,
                       std::vector<double>& rYValues,
                       double startTime,
                       double endTime,
                       double timeStep);

    /**
     * Virtual destructor since we have virtual methods.
     */
    virtual ~AbstractOneStepIvpOdeSolver()
    {}
};

CLASS_IS_ABSTRACT(AbstractOneStepIvpOdeSolver)

#endif //_ABSTRACTONESTEPIVPODESOLVER_HPP_
