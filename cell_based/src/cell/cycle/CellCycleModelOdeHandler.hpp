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

#ifndef CELLCYCLEMODELODEHANDLER_HPP_
#define CELLCYCLEMODELODEHANDLER_HPP_

#include <boost/shared_ptr.hpp>
#include <boost/noncopyable.hpp>

#include "ChasteSerialization.hpp"
#include "AbstractOdeSystem.hpp"
#include "AbstractCellCycleModelOdeSolver.hpp"
#include "SimulationTime.hpp"

/**
 * This class contains the functionality for running ODEs as part of a cell cycle
 * or SRN model.  It is designed to be used as an additional base class for models
 * which require this functionality.
 */
class CellCycleModelOdeHandler
{
private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & mpOdeSystem;
        archive & mpOdeSolver;
        archive & mLastTime;
        archive & mDt;
        archive & mFinishedRunningOdes;
    }

    /**
     * Prevent copy-assignment of this class, or its subclasses.
     * Note that we do not define this method, therefore statements like "CellCycleModelOdeHandler new = old;" will not compile.
     * We do not inherit from boost::noncopyable because we *do* define a protected copy-constructor, for use by CreateCellCycleModel
     * and CreateSrnModel.
     *
     * @return the new ODE handler.
     */
    CellCycleModelOdeHandler& operator=(const AbstractCellCycleModelOdeSolver&);

protected:

    /**
     * Protected copy-constructor for use by CreateCellCycleModel and CreateSrnModel.
     *
     * @param rHandler ODE handler to copy.
     */
    CellCycleModelOdeHandler(const CellCycleModelOdeHandler& rHandler);

    /**
     * Timestep to use when solving the ODE system.
     * For some adaptive solvers (e.g. CVODE) this is the maximum step to use.
     */
    double mDt;

    /** A system of ODEs. */
    AbstractOdeSystem* mpOdeSystem;

    /** A shared pointer to an ODE solver object. */
    boost::shared_ptr<AbstractCellCycleModelOdeSolver> mpOdeSolver;

    /** The last time the ODE system was evaluated. */
    double mLastTime;

    /**
     * Whether the model is currently in a delay (not solving ODEs).
     */
    bool mFinishedRunningOdes;

    /**
     * Solves the ODE system to a given time.
     *
     * @param currentTime the current time
     *
     * @return whether a stopping event occurred.
     */
    bool SolveOdeToTime(double currentTime);

    /**
     * Adjust any ODE parameters needed before solving until currentTime.
     * Defaults to do nothing.
     *
     * @param currentTime  the time up to which the system will be solved.
     */
    virtual void AdjustOdeParameters(double currentTime);

public:

    /**
     * Constructor.
     *
     * @param lastTime  The birth time of the cell / last time model was evaluated (defaults to the current SimulationTime)
     * @param pOdeSolver An optional pointer to a cell-cycle model ODE solver object (allows the use of different ODE solvers)
     */
    CellCycleModelOdeHandler(double lastTime = SimulationTime::Instance()->GetTime(),
                             boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver = boost::shared_ptr<AbstractCellCycleModelOdeSolver>());

    /**
     * Destructor.
     */
    virtual ~CellCycleModelOdeHandler();

    /**
     * @return #mpOdeSystem.
     */
    AbstractOdeSystem* GetOdeSystem() const;

    /**
     * Set mpOdeSystem. Used in CreateCellCycleModel().
     *
     * @param pOdeSystem the ODE system
     */
    void SetOdeSystem(AbstractOdeSystem* pOdeSystem);

    /**
     * @return mpOdeSolver
     */
    const boost::shared_ptr<AbstractCellCycleModelOdeSolver> GetOdeSolver() const;

    /**
     * Set mLastTime.
     *
     * @param lastTime the new value of mLastTime
     */
    void SetLastTime(double lastTime);

    /**
     * @return #mDt.  This sets it to a default value if it hasn't
     * been set by calling SetDt.
     */
    double GetDt();

    /**
     * Set the time step to use to solve the ODE system.
     * For some adaptive solvers (e.g. CVODE) this is the maximum step to use.
     *
     * @param timeStep  time step to use
     */
    void SetDt(double timeStep);

    /**
     * Set the values of the state variables in the cell-cycle model's ODE system.
     *
     * @param rStateVariables vector containing values for the state variables
     */
    void SetStateVariables(const std::vector<double>& rStateVariables);

    /**
     * @return the protein concentrations at the current time (useful for tests)
     *
     * NB: Will copy the vector - you can't use this to modify the concentrations.
     */
    std::vector<double> GetProteinConcentrations() const;

    /**
     * Sets the protein concentrations and time when the model was last evaluated - should only be called by tests
     *
     * @param lastTime the SimulationTime at which the protein concentrations apply
     * @param proteinConcentrations a standard vector of doubles of protein concentrations
     *
     */
    void SetProteinConcentrationsForTestsOnly(double lastTime, std::vector<double> proteinConcentrations);
};

#endif /*CELLCYCLEMODELODEHANDLER_HPP_*/
