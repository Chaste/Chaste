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

#ifndef _ABSTRACTIVPODESOLVER_HPP_
#define _ABSTRACTIVPODESOLVER_HPP_

#include <vector>

#include "ChasteSerialization.hpp"
#include "ClassIsAbstract.hpp"
#include "Identifiable.hpp"

#include "OdeSolution.hpp"
#include "AbstractOdeSystem.hpp"

/**
 * Abstract initial value problem ODE solver class. Sets up variables and functions
 * for a numerical solution technique for an initial value ODE problem.
 */
class AbstractIvpOdeSolver : public Identifiable
{
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & mStoppingEventOccurred;
        archive & mStoppingTime;
    }

protected:

    /**
     * Boolean indicating whether the solver quit due to the ODEs
     * stopping event occuring
     */
    bool mStoppingEventOccurred;

    /** If a stopping event occurred the time is stored here.  (Only valid when mStoppingEventOccurred==true) */
    double mStoppingTime;

public:

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
                              double timeSampling)=0;

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
    virtual void Solve(AbstractOdeSystem* pAbstractOdeSystem,
                       std::vector<double>& rYValues,
                       double startTime,
                       double endTime,
                       double timeStep)=0;

    /**
     * Solves a system of ODEs using a specified one-step ODE solver and update the
     * state variables.
     *
     * @param pAbstractOdeSystem  pointer to the concrete ODE system to be solved
     * @param startTime  the time at which the initial conditions are specified
     * @param endTime  the time to which the system should be solved and the solution
     *                 returned
     * @param timeStep  the time interval to be used by the solver
     */
    virtual void SolveAndUpdateStateVariable(AbstractOdeSystem* pAbstractOdeSystem,
                                             double startTime,
                                             double endTime,
                                             double timeStep);

    /**
     * Determine whether the solver quit due to the ODE's stopping event
     * triggering
     */
    bool StoppingEventOccurred();

    /**
     * Get the stopping time for the solver.
     *
     * @return mStoppingTime.
     */
    double GetStoppingTime();

    /**
     * Constructor.
     */
    AbstractIvpOdeSolver();

    /**
     * Virtual destructor since we have virtual methods.
     */
    virtual ~AbstractIvpOdeSolver();

};

CLASS_IS_ABSTRACT(AbstractIvpOdeSolver)

#endif //_ABSTRACTIVPODESOLVER_HPP_
