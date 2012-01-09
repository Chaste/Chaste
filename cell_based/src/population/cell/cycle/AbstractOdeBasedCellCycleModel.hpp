/*

Copyright (C) University of Oxford, 2005-2012

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

#ifndef ABSTRACTODEBASEDCELLCYCLEMODEL_HPP_
#define ABSTRACTODEBASEDCELLCYCLEMODEL_HPP_

#include <vector>

#include "ChasteSerialization.hpp"
#include "ClassIsAbstract.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractCellCycleModel.hpp"
#include "CellCycleModelOdeHandler.hpp"
#include "SimulationTime.hpp"

/**
 * This class contains all the things common to standard cell cycle
 * ODE models for intracellular protein concentrations (along the lines
 * of Tyson & Novak), such as solving the ODEs until a stopping condition
 * is met.
 */
class AbstractOdeBasedCellCycleModel : public AbstractCellCycleModel, public CellCycleModelOdeHandler
{
private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the cell-cycle model and member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellCycleModel>(*this);
        archive & boost::serialization::base_object<CellCycleModelOdeHandler>(*this);
        archive & mDivideTime;
        archive & mFinishedRunningOdes;
        archive & mG2PhaseStartTime;
    }

protected:

    /** The time at which the cell should divide - Set this to DBL_MAX in constructor. */
    double mDivideTime;

    /** Whether the cell-cycle model is currently in a delay (not solving ODEs). */
    bool mFinishedRunningOdes;

    /** The start time for the G2 phase. */
    double mG2PhaseStartTime;

public:

    /**
     * Creates an AbstractOdeBasedCellCycleModel, calls SetBirthTime on the
     * AbstractCellCycleModel to make sure that can be set 'back in time' for
     * cells which did not divide at the current time.
     *
     * @param lastTime  The birth time of the cell / last time model was evaluated (defaults to the current SimulationTime)
     * @param pOdeSolver An optional pointer to a cell-cycle model ODE solver object (allows the use of different ODE solvers)
     */
    AbstractOdeBasedCellCycleModel(double lastTime = SimulationTime::Instance()->GetTime(),
                                   boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver = boost::shared_ptr<AbstractCellCycleModelOdeSolver>());

    /**
     * Destructor.
     */
    virtual ~AbstractOdeBasedCellCycleModel();

    /**
     * Default UpdateCellCyclePhase() method for an ODE-based cell-cycle model.
     * This method calls SolveOdeToTime() for G1 phase and adds time for the other phases.
     *
     * Can be overridden if they should do something more subtle.
     */
    virtual void UpdateCellCyclePhase();

    /**
     * Get the time at which the ODE stopping event occurred.
     * Only called in those subclasses for which stopping events
     * are defined.
     *
     * @return the time at which the ODE system reached its stopping event
     */
    double GetOdeStopTime();

    /**
     * This overrides the AbstractCellCycleModel::SetBirthTime(double birthTime)
     * because an ODE based cell-cycle model has more to reset...
     *
     * @param birthTime the simulation time when the cell was born
     */
    void SetBirthTime(double birthTime);

    /**
     * Returns the protein concentrations at the current time (useful for tests)
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

    /**
     * For a naturally cycling model this does not need to be overridden in the
     * subclasses. But most models should override this function and then
     * call AbstractOdeBasedCellCycleModel::ResetForDivision() from inside their version.
     */
    virtual void ResetForDivision();

    /**
     * Set mFinishedRunningOdes. Used in CreateCellCycleModel().
     *
     * @param finishedRunningOdes the new value of mFinishedRunningOdes
     */
    void SetFinishedRunningOdes(bool finishedRunningOdes);

    /**
     * Set mDivideTime.
     *
     * @param divideTime the new value of mDivideTime
     */
    void SetDivideTime(double divideTime);

    /**
     * Set mG2PhaseStartTime. Used in CreateCellCycleModel().
     *
     * @param g2PhaseStartTime the new value of mG2PhaseStartTime
     */
    void SetG2PhaseStartTime(double g2PhaseStartTime);

    /**
     * Outputs cell cycle model parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputCellCycleModelParameters(out_stream& rParamsFile);
};

CLASS_IS_ABSTRACT(AbstractOdeBasedCellCycleModel)

#endif /*ABSTRACTODEBASEDCELLCYCLEMODEL_HPP_*/
