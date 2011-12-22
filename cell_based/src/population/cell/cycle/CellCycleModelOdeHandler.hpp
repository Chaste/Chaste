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

#ifndef CELLCYCLEMODELODEHANDLER_HPP_
#define CELLCYCLEMODELODEHANDLER_HPP_

#include <boost/shared_ptr.hpp>

#include "ChasteSerialization.hpp"
#include "AbstractOdeSystem.hpp"
#include "AbstractCellCycleModelOdeSolver.hpp"
#include "SimulationTime.hpp"

/**
 * This class contains the functionality for running ODEs as part of a cell cycle
 * model.  It is designed to be used as an additional base class for models which
 * require this functionality.
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
    }

protected:

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
     * Get mpOdeSystem.
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
     * Get method for #mDt, which sets it to a default value if it hasn't
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
};

#endif /*CELLCYCLEMODELODEHANDLER_HPP_*/
