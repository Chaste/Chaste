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

#ifndef ABSTRACTCELLCYCLEMODELODESOLVER_HPP_
#define ABSTRACTCELLCYCLEMODELODESOLVER_HPP_

#include "ChasteSerialization.hpp"
#include "ClassIsAbstract.hpp"
#include <boost/serialization/base_object.hpp>

#include <boost/shared_ptr.hpp>

#include "AbstractIvpOdeSolver.hpp"

/**
 * This provides a wrapper around any ODE solver class, exposing roughly the same interface,
 * for use by ODE-based cell-cycle models.  Its main purpose is to allow multiple instances
 * of the same cell-cycle model to share the same ODE solver instance.
 *
 * The recommended way to use this wrapper is via the CellCycleModelOdeSolver subclass, which
 * is templated over cell-cycle model class and ODE solver class, providing a singleton
 * instance for each combination of template parameters.
 */
class AbstractCellCycleModelOdeSolver
{
private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the object and its the member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & mpOdeSolver;
        archive & mSizeOfOdeSystem;
    }

protected:

    /** The ODE solver. */
    boost::shared_ptr<AbstractIvpOdeSolver> mpOdeSolver;

    /** The size of the ODE system to be solved. */
    unsigned mSizeOfOdeSystem;

public:

    /**
     * Constructor.
     */
    AbstractCellCycleModelOdeSolver();

    /**
     * Virtual destructor since we have virtual methods.
     */
    virtual ~AbstractCellCycleModelOdeSolver();

    /**
     * @return whether the instance in existence and fully set up.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     */
    virtual bool IsSetUp()=0;

    /**
     * Reset the instance.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     */
    virtual void Reset()=0;

    /**
     * Call mpOdeSolver->SolveAndUpdateStateVariable.
     *
     * @param pAbstractOdeSystem  pointer to the concrete ODE system to be solved
     * @param startTime  the time at which the initial conditions are specified
     * @param endTime  the time to which the system should be solved and the solution
     *                 returned
     * @param timeStep  the time interval to be used by the solver
     */
    void SolveAndUpdateStateVariable(AbstractOdeSystem* pAbstractOdeSystem,
                                     double startTime,
                                     double endTime,
                                     double timeStep);

    /**
     * Initialise the ODE solver.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     */
    virtual void Initialise()=0;

    /**
     * @return whether the solver quit due to the ODE's stopping event
     * triggering
     */
    bool StoppingEventOccurred();

    /**
     * Call mpOdeSolver->GetStoppingTime.
     *
     * @return mStoppingTime.
     */
    double GetStoppingTime();

    /**
     * Set method for mSizeOfOdeSystem.
     *
     * @param sizeOfOdeSystem the new value of mSizeOfOdeSystem
     */
    void SetSizeOfOdeSystem(unsigned sizeOfOdeSystem);

    /**
     * @return mSizeOfOdeSystem
     */
    unsigned GetSizeOfOdeSystem();

    /**
     * If using CVODE, make the solver check for stopping events using CVODE's rootfinding functionality
     * (by default we do not check).
     */
    void CheckForStoppingEvents();

    /**
     * If using CVODE, change the maximum number of steps to be taken by the solver in its attempt to reach
     * the next output time (default is 500).
     *
     * @param numSteps the new maximum number of steps
     */
    void SetMaxSteps(long int numSteps);

    /**
     * If using CVODE, set relative and absolute tolerances; both scalars.
     * If no parameters are given, tolerances will be reset to default values.
     *
     * @param relTol the relative tolerance for the solver
     * @param absTol the absolute tolerance for the solver
     */
    void SetTolerances(double relTol=1e-4, double absTol=1e-6);

    /**
     * @return true iff this is an adaptive solver such as CVODE for which it is safe to set the 'timestep'
     * to be the outer simulation timestep, because the ODE solver will use this as its maximum, not actual,
     * timestep.
     *
     * The base class version just returns true iff the solver is the CvodeAdaptor class.
     */
    virtual bool IsAdaptive();
};

#endif /*ABSTRACTCELLCYCLEMODELODESOLVER_HPP_*/
