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

#ifndef SIMPLEOXYGENBASEDCELLCYCLEMODEL_HPP_
#define SIMPLEOXYGENBASEDCELLCYCLEMODEL_HPP_

#include "AbstractSimpleCellCycleModel.hpp"
#include "CellwiseData.hpp"

/**
 * Simple oxygen-based cell-cycle model.
 *
 * A simple oxygen-dependent cell-cycle model that inherits from
 * AbstractSimpleCellCycleModel. The duration of G1 phase depends
 * on the local oxygen concentration. A prolonged period of acute
 * hypoxia leads to the cell being labelled as apoptotic. This model
 * allows for quiescence imposed by transient periods of hypoxia,
 * followed by reoxygenation.
 */
class SimpleOxygenBasedCellCycleModel : public AbstractSimpleCellCycleModel
{
private:

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractSimpleCellCycleModel>(*this);
        archive & mTimeSpentInG1Phase;
        archive & mCurrentHypoxicDuration;
        archive & mCurrentHypoxiaOnsetTime;
        archive & mHypoxicConcentration;
        archive & mQuiescentConcentration;
        archive & mCriticalHypoxicDuration;
    }

protected:

    /**
     * The time spent in G1 phase so far.
     * Has units of hours.
     */
    double mTimeSpentInG1Phase;

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
    double GetCurrentHypoxicDuration();

    /**
     * @return mCurrentHypoxiaOnsetTime
     */
    double GetCurrentHypoxiaOnsetTime();

    /**
     * Overridden builder method to create new copies of
     * this cell-cycle model.
     */
    AbstractCellCycleModel* CreateCellCycleModel();

    /**
     * @return mHypoxicConcentration
     */
    double GetHypoxicConcentration();

    /**
     * Set method for mHypoxicConcentration.
     *
     * @param hypoxicConcentration the new value of mHypoxicConcentration
     */
    void SetHypoxicConcentration(double hypoxicConcentration);

    /**
     * @return mQuiescentConcentration
     */
    double GetQuiescentConcentration();

    /**
     * Set method for mQuiescentConcentration.
     *
     * @param quiescentConcentration the new value of mQuiescentConcentration
     */
    void SetQuiescentConcentration(double quiescentConcentration);

    /**
     * @return mCriticalHypoxicDuration
     */
    double GetCriticalHypoxicDuration();

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
     * Outputs cell cycle model parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputCellCycleModelParameters(out_stream& rParamsFile);
};

// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(SimpleOxygenBasedCellCycleModel)

#endif /*SIMPLEOXYGENBASEDCELLCYCLEMODEL_HPP_*/
