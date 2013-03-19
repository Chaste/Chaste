/*

Copyright (c) 2005-2013, University of Oxford.
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
