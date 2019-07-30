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

#ifndef ABSTRACTSIMPLECELLCYCLEMODEL_HPP_
#define ABSTRACTSIMPLECELLCYCLEMODEL_HPP_

#include "ChasteSerialization.hpp"
#include "ClassIsAbstract.hpp"

#include <boost/serialization/base_object.hpp>

#include <vector>

#include "AbstractCellCycleModel.hpp"
#include "CellCyclePhases.hpp"
#include "SimulationTime.hpp"


/**
 * This class contains basic information to all cell-cycle models
 * that do NOT explicitly include distinct phases (G1, S, G2, M).
 */
class AbstractSimpleCellCycleModel : public AbstractCellCycleModel
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
        archive & boost::serialization::base_object<AbstractCellCycleModel>(*this);
        archive & mCellCycleDuration;
    }

protected:

    /**
     * Duration of cell cycle.
     * May be used as a mean duration for stochastic cell-cycle models.
     */
    double mCellCycleDuration;

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
    AbstractSimpleCellCycleModel(const AbstractSimpleCellCycleModel& rModel);

public:

    /**
     * Default constructor - creates an AbstractSimpleCellCycleModel.
     */
    AbstractSimpleCellCycleModel();

    /**
     * Destructor.
     */
    virtual ~AbstractSimpleCellCycleModel();

    /**
     * See AbstractCellCycleModel::ResetForDivision()
     *
     * @return whether the cell is ready to divide (enter M phase).
     */
    virtual bool ReadyToDivide();

    /** See AbstractCellCycleModel::ResetForDivision() */
    virtual void ResetForDivision();

    /**
     * Overridden InitialiseDaughterCell() method.
     *
     * Set the new cell's cell cycle duration once it has been created after division.
     * This is by calling SetCellCycleDuration() defined in child classes.
     */
    void InitialiseDaughterCell();

    /** See AbstractPhaseBasedCellCycleModel::Initialise()
     *
     * Calls SetCellCycleDuration() defined in child classes.
     */
    virtual void Initialise();

    /**
     * This method is implemented in Subclasses to set the cell cycle duration of the cell.
     * It is called on the beginning of the simulation (intialisation) when initialising
     * daughter cells and when reseting parent cells
     *
     * It must be implemented in child classes.
     */
    virtual void SetCellCycleDuration() = 0;

    /**
     * @return mCellCycleDuration
     */
    double GetCellCycleDuration() const;

    /**
     * Outputs cell cycle model parameters to file.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputCellCycleModelParameters(out_stream& rParamsFile)=0;
};

CLASS_IS_ABSTRACT(AbstractSimpleCellCycleModel)

#endif /*ABSTRACTSIMPLECELLCYCLEMODEL_HPP_*/
