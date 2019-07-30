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
 * This class contains all the functionality shared by 'ODE-based' cell-cycle models,
 * where the duration of the cell cycle is evaluated 'on the fly'
 * as the cell ages, according to a system of ordinary differential equations (ODEs)
 * governing (for example) the concentrations of key intracellular proteins. To
 * determine when cell division should occur, one or more stopping conditions for
 * this ODE system may be specified.
 *
 * This class of cell-cycle models is distinct from 'simple' cell-cycle models, where
 * the duration of the cell cycle is determined when the cell-cycle model is
 * created.
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
    }

protected:

    /** The time at which the cell should divide - Set this to DBL_MAX in constructor. */
    double mDivideTime;

    /**
     * Protected copy-constructor for use by CreateCellCycleModel.
     * The only way for external code to create a copy of a cell cycle model
     * is by calling that method, to ensure that a model of the correct subclass is created.
     * This copy-constructor helps subclasses to ensure that all member variables are correctly copied when this happens.
     *
     * This method is called by child classes to set member variables for a daughter cell upon cell division.
     * Note that the parent cell cycle model will have had ResetForDivision() called just before CreateCellCycleModel() is called,
     * so performing an exact copy of the parent is suitable behaviour. Any daughter-cell-specific initialisation
     * can be done in InitialiseDaughterCell().
     *
     * @param rModel the cell cycle model to copy.
     */
    AbstractOdeBasedCellCycleModel(const AbstractOdeBasedCellCycleModel& rModel);

public:

    /**
     * Creates an AbstractOdeBasedCellCycleModel, calls SetBirthTime on the
     * AbstractPhaseBasedCellCycleModel to make sure that can be set 'back in time' for
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
     * See AbstractCellCycleModel::ResetForDivision()
     *
     * @return whether the cell is ready to divide (enter M phase).
     */
    virtual bool ReadyToDivide();

    /**
     * For a naturally cycling model this does not need to be overridden in the
     * subclasses. But most models should override this function and then
     * call AbstractOdeBasedCellCycleModel::ResetForDivision() from inside their version.
     */
    virtual void ResetForDivision();

    /**
     * Outputs cell cycle model parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputCellCycleModelParameters(out_stream& rParamsFile);
};

CLASS_IS_ABSTRACT(AbstractOdeBasedCellCycleModel)

#endif /*ABSTRACTODEBASDCELLCYCLEMODEL_HPP_*/
