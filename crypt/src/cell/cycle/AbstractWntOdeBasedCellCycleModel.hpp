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

#ifndef ABSTRACTWNTODEBASEDCELLCYCLEMODEL_HPP_
#define ABSTRACTWNTODEBASEDCELLCYCLEMODEL_HPP_

#include "ChasteSerialization.hpp"
#include "ClassIsAbstract.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractOdeBasedPhaseBasedCellCycleModel.hpp"
#include "WntConcentration.hpp"

/**
 * This class contains all the things common to the Wnt cell cycle ODE based models,
 * the resetting method and updating of cell types etc.
 *
 * The concrete models all need to operate with a WntConcentration
 * singleton object.
 *
 * These models have a constant length M phase, run ODEs to decide when
 * to finish G1 phase, then add time for S and G2 phases (in some classes,
 * random periods of time). The CellProliferativeType is updated dependent on the
 * concentration of beta-catenin (given by one of the ODEs).
 */
class AbstractWntOdeBasedCellCycleModel : public AbstractOdeBasedPhaseBasedCellCycleModel
{
private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the cell-cycle model, never used directly - boost uses this.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractOdeBasedPhaseBasedCellCycleModel>(*this);
    }

protected:

    /**
     * @return the Wnt level experienced by the cell.
     */
    double GetWntLevel() const;

    /**
     * Call base class UpdateCellCyclePhase, then UpdateCellProliferativeType.
     */
    void UpdateCellCyclePhase();

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
    AbstractWntOdeBasedCellCycleModel(const AbstractWntOdeBasedCellCycleModel& rModel);

public:

    /**
     * Default constructor.
     *
     * @param pOdeSolver An optional pointer to a cell-cycle model ODE solver object (allows the use of different ODE solvers)
     */
    AbstractWntOdeBasedCellCycleModel(boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver = boost::shared_ptr<AbstractCellCycleModelOdeSolver>());

    /**
     * Destructor.
     */
    virtual ~AbstractWntOdeBasedCellCycleModel();

    /**
     * Resets the Wnt Model to the start of the cell cycle (this model does not cycle naturally)
     * Cells are given a new birth time and cell cycle proteins are reset.
     * Note that the Wnt pathway proteins maintain their current values.
     *
     * Should only be called by the Cell::Divide() method.
     */
    void ResetForDivision();

    /**
     * Updates the current cell type to reflect whether the
     * beta-catenin level has dropped low enough to make it stop dividing.
     * This should only be called when the cell-cycle model has been
     * evaluated to the current time, or it may give misleading results.
     */
    void UpdateCellProliferativeType();

    /**
     * Change cell type to reflect current levels of beta-catenin.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     */
    virtual void ChangeCellProliferativeTypeDueToCurrentBetaCateninLevel()=0;

    /**
     * Overridden GetAverageTransitCellCycleTime() method.
     * @return time
     */
    double GetAverageTransitCellCycleTime();

    /**
     * Overridden GetAverageStemCellCycleTime() method.
     * @return time
     */
    double GetAverageStemCellCycleTime();

    /**
     * Overridden CanCellTerminallyDifferentiate() method.
     * @return whether cell can terminally differentiate
     *
     */
    virtual bool CanCellTerminallyDifferentiate();

    /**
     * Overridden OutputCellCycleModelParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputCellCycleModelParameters(out_stream& rParamsFile);
};

CLASS_IS_ABSTRACT(AbstractWntOdeBasedCellCycleModel)

#endif /*ABSTRACTWNTODEBASEDCELLCYCLEMODEL_HPP_*/
