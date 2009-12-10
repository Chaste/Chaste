/*

Copyright (C) University of Oxford, 2005-2009

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
#ifndef STOCHASTICOXYGENBASEDCELLCYCLEMODEL_HPP_
#define STOCHASTICOXYGENBASEDCELLCYCLEMODEL_HPP_

#include "AbstractSimpleCellCycleModel.hpp"
#include "CellwiseData.hpp"
#include "RandomNumberGenerator.hpp"

/**
 *  Stochastic oxygen-based cell cycle model
 *
 *  A simple oxygen-dependent cell cycle model that inherits from
 *  AbstractSimpleCellCycleModel. The duration of G1 phase depends
 *  on the local oxygen concentration. A prolonged period of acute
 *  hypoxia leads to the cell being labelled as apoptotic. This model
 *  allows for quiescence imposed by transient periods of hypoxia,
 *  followed by reoxygenation.
 *
 */
class StochasticOxygenBasedCellCycleModel : public AbstractSimpleCellCycleModel
{
private:

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractSimpleCellCycleModel>(*this);

        RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
        archive & *p_gen;
        archive & p_gen;

        archive & mG2Duration;
        archive & mTimeSpentInG1Phase;
        archive & mCurrentHypoxicDuration;
        archive & mCurrentHypoxiaOnsetTime;
        archive & mDimension;
    }

    /**
     * The duration of the G2 phase, set stochastically.
     */
    double mG2Duration;

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

    /**
     * Stochastically set the G2 duration.  Called on cell creation at
     * the start of a simulation, and for both parent and daughter
     * cells at cell division.
     */
    void SetG2Duration();

public:

    /**
     * Constructor.
     *
     * @param dimension the spatial dimension (needed by the templated class CellwiseData)
     */
    StochasticOxygenBasedCellCycleModel(unsigned dimension);

    /**
     * Overridden InitialiseDaughterCell() method.
     */
    void InitialiseDaughterCell();

    /**
     * Initialise the cell cycle model at the start of a simulation.
     */
    void Initialise();

    /**
     * Overridden ResetForDivision() method.
     */
    void ResetForDivision();

    /**
     * @return mG2Duration.
     */
    double GetG2Duration();

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
     * Get method for mCurrentHypoxicDuration.
     */
    double GetCurrentHypoxicDuration();

    /**
     * Get method for mCurrentHypoxiaOnsetTime.
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
CHASTE_CLASS_EXPORT(StochasticOxygenBasedCellCycleModel)


namespace boost
{
namespace serialization
{
/**
 * Allow us to not need a default constructor, by specifying how Boost should
 * instantiate a StochasticOxygenBasedCellCycleModel instance.
 */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const StochasticOxygenBasedCellCycleModel * t, const unsigned int file_version)
{
}

/**
 * Allow us to not need a default constructor, by specifying how Boost should
 * instantiate a StochasticOxygenBasedCellCycleModel instance.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, StochasticOxygenBasedCellCycleModel * t, const unsigned int file_version)
{
    /**
     * Invoke inplace constructor to initialise an instance of StochasticOxygenBasedCellCycleModel.
     * It doesn't actually matter what values we pass to our standard constructor,
     * provided they are valid parameter values, since the state loaded later
     * from the archive will overwrite their effect in this case.
     */

    unsigned dimension = UINT_MAX;
    ::new(t)StochasticOxygenBasedCellCycleModel(dimension);
}
}
} // namespace ...

#endif /*STOCHASTICOXYGENBASEDCELLCYCLEMODEL_HPP_*/
