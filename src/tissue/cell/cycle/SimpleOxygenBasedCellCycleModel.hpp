/*

Copyright (C) University of Oxford, 2005-2010

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
 *  Simple oxygen-based cell cycle model
 *
 *  A simple oxygen-dependent cell cycle model that inherits from
 *  AbstractSimpleCellCycleModel. The duration of G1 phase depends
 *  on the local oxygen concentration. A prolonged period of acute
 *  hypoxia leads to the cell being labelled as apoptotic. This model
 *  allows for quiescence imposed by transient periods of hypoxia,
 *  followed by reoxygenation.
 *
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
        archive & mDimension;
    }

    /**
     * The time spent in G1 phase so far.
     */
    double mTimeSpentInG1Phase;

    /**
     * How long the current period of hypoxia has lasted.
     */
    double mCurrentHypoxicDuration;

    /**
     * The time when the current period of hypoxia began.
     */
    double mCurrentHypoxiaOnsetTime;

    /**
     * The spatial dimension (needed by the templated class CellwiseData).
     */
    unsigned mDimension;

public:

    /**
     * Constructor.
     *
     * @param dimension the spatial dimension (needed by the templated class CellwiseData)
     */
    SimpleOxygenBasedCellCycleModel(unsigned dimension);

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
     * this cell cycle model.
     */
    AbstractCellCycleModel* CreateCellCycleModel();

    /**
     * Get the spatial dimension.
     *
     * @return mDimension
     */
    unsigned GetDimension();
};

// Declare identifier for the serializer
#include "TemplatedExport.hpp"
CHASTE_CLASS_EXPORT(SimpleOxygenBasedCellCycleModel)


namespace boost
{
namespace serialization
{
/**
 * Allow us to not need a default constructor, by specifying how Boost should
 * instantiate a SimpleOxygenBasedCellCycleModel instance.
 */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const SimpleOxygenBasedCellCycleModel * t, const unsigned int file_version)
{
}

/**
 * Allow us to not need a default constructor, by specifying how Boost should
 * instantiate a SimpleOxygenBasedCellCycleModel instance.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, SimpleOxygenBasedCellCycleModel * t, const unsigned int file_version)
{
    /**
     * Invoke inplace constructor to initialise an instance of SimpleOxygenBasedCellCycleModel.
     * It doesn't actually matter what values we pass to our standard constructor,
     * provided they are valid parameter values, since the state loaded later
     * from the archive will overwrite their effect in this case.
     */

    unsigned dimension = UINT_MAX;
    ::new(t)SimpleOxygenBasedCellCycleModel(dimension);
}
}
} // namespace ...

#endif /*SIMPLEOXYGENBASEDCELLCYCLEMODEL_HPP_*/
