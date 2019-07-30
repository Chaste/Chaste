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
     * @return true if the solver quit due to the ODE's stopping event
     * triggering
     */
    bool StoppingEventOccurred();

    /**
     * Get the stopping time for the solver.
     *
     * @return #mStoppingTime.
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
