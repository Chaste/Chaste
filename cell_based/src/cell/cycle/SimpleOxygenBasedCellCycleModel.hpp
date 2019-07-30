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

#ifndef SIMPLEOXYGENBASEDCELLCYCLEMODEL_HPP_
#define SIMPLEOXYGENBASEDCELLCYCLEMODEL_HPP_

#include "AbstractSimplePhaseBasedCellCycleModel.hpp"

/**
 * Simple oxygen-based cell-cycle model.
 *
 * A simple oxygen-dependent cell-cycle model that inherits from
 * AbstractSimplePhaseBasedCellCycleModel. The duration of G1 phase depends
 * on the local oxygen concentration. A prolonged period of acute
 * hypoxia leads to the cell being labelled as apoptotic. This model
 * allows for quiescence imposed by transient periods of hypoxia,
 * followed by reoxygenation.
 */
class SimpleOxygenBasedCellCycleModel : public AbstractSimplePhaseBasedCellCycleModel
{
private:

    friend class boost::serialization::access;
    /**
     * Boost Serialization method for archiving/checkpointing
     * @param archive  The boost archive.
     * @param version  The current version of this class.
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractSimplePhaseBasedCellCycleModel>(*this);
        archive & mCurrentHypoxicDuration;
        archive & mCurrentHypoxiaOnsetTime;
        archive & mHypoxicConcentration;
        archive & mQuiescentConcentration;
        archive & mCriticalHypoxicDuration;
    }

protected:

    /**
     * How long the current period of hypoxia has lasted.
     * Has units of hours.
     */
    double mCurrentHypoxicDuration;

    /**
     * The time when the current period of hypoxia began.
     */
    double mCurrentHypoxiaOnsetTime;

    /**
     * Non-dimensionalized oxygen concentration below which cells are
     * considered to be hypoxic. A prolonged period of hypoxia causes
     * the cell to become apoptotic.
     */
    double mHypoxicConcentration;

    /**
     * Non-dimensionalized oxygen concentration below which cells are
     * considered to be quiescent and slow their progress through the
     * G1 phase of the cell cycle.
     */
    double mQuiescentConcentration;

    /**
     * Non-dimensionalized critical hypoxic duration.
     * Has units of hours.
     */
    double mCriticalHypoxicDuration;

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
    SimpleOxygenBasedCellCycleModel(const SimpleOxygenBasedCellCycleModel& rModel);

public:

    /**
     * Constructor.
     */
    SimpleOxygenBasedCellCycleModel();

    /**
     * Overridden UpdateCellCyclePhase() method.
     */
    void UpdateCellCyclePhase();

    /**
     * Method for updating mCurrentHypoxicDuration,
     * called at the start of ReadyToDivide().
     */
    void UpdateHypoxicDuration();

    /**
     * @return mCurrentHypoxicDuration
     */
    double GetCurrentHypoxicDuration() const;

    /**
     * @return mCurrentHypoxiaOnsetTime
     */
    double GetCurrentHypoxiaOnsetTime() const;

    /**
     * Overridden builder method to create new copies of
     * this cell-cycle model.
     *
     * @return new cell-cycle model
     */
    AbstractCellCycleModel* CreateCellCycleModel();

    /**
     * @return mHypoxicConcentration
     */
    double GetHypoxicConcentration() const;

    /**
     * Set method for mHypoxicConcentration.
     *
     * @param hypoxicConcentration the new value of mHypoxicConcentration
     */
    void SetHypoxicConcentration(double hypoxicConcentration);

    /**
     * @return mQuiescentConcentration
     */
    double GetQuiescentConcentration() const;

    /**
     * Set method for mQuiescentConcentration.
     *
     * @param quiescentConcentration the new value of mQuiescentConcentration
     */
    void SetQuiescentConcentration(double quiescentConcentration);

    /**
     * @return mCriticalHypoxicDuration
     */
    double GetCriticalHypoxicDuration() const;

    /**
     * Set method for mCriticalHypoxicDuration.
     *
     * @param criticalHypoxicDuration the new value of mCriticalHypoxicDuration
     */
    void SetCriticalHypoxicDuration(double criticalHypoxicDuration);

    /**
     * Set method for mCurrentHypoxiaOnsetTime.
     *
     * @param currentHypoxiaOnsetTime the new value of mCurrentHypoxiaOnsetTime
     */
    void SetCurrentHypoxiaOnsetTime(double currentHypoxiaOnsetTime);

    /**
     * Overridden OutputCellCycleModelParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputCellCycleModelParameters(out_stream& rParamsFile);
};

// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(SimpleOxygenBasedCellCycleModel)

#endif /*SIMPLEOXYGENBASEDCELLCYCLEMODEL_HPP_*/
