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

#ifndef STOCHASTICOXYGENBASEDCELLCYCLEMODEL_HPP_
#define STOCHASTICOXYGENBASEDCELLCYCLEMODEL_HPP_

#include "SimpleOxygenBasedCellCycleModel.hpp"
#include "RandomNumberGenerator.hpp"

/**
 * Stochastic oxygen-based cell-cycle model.
 *
 * A simple oxygen-dependent cell-cycle model that inherits from
 * SimpleOxygenBasedCellCycleModel and in addition spends a random
 * duration in G2 phase.
 */
class StochasticOxygenBasedCellCycleModel : public SimpleOxygenBasedCellCycleModel
{
    friend class TestSimpleCellCycleModels;

private:

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<SimpleOxygenBasedCellCycleModel>(*this);

        // Make sure the RandomNumberGenerator singleton gets saved too
        SerializableSingleton<RandomNumberGenerator>* p_wrapper = RandomNumberGenerator::Instance()->GetSerializationWrapper();
        archive & p_wrapper;

        archive & mStochasticG2Duration;
    }

    /**
     * The duration of the G2 phase, set stochastically.
     */
    double mStochasticG2Duration;

    /**
     * Stochastically set the G2 duration.  Called on cell creation at
     * the start of a simulation, and for both parent and daughter
     * cells at cell division.
     */
    void GenerateStochasticG2Duration();

public:

    /**
     * Constructor.
     */
    StochasticOxygenBasedCellCycleModel();

    /**
     * Overridden InitialiseDaughterCell() method.
     */
    void InitialiseDaughterCell();

    /**
     * Initialise the cell-cycle model at the start of a simulation.
     */
    void Initialise();

    /**
     * Overridden ResetForDivision() method.
     */
    void ResetForDivision();

    /**
     * @return mStochasticG2Duration.
     */
    double GetG2Duration();

    /**
     * Overridden builder method to create new copies of
     * this cell-cycle model.
     */
    AbstractCellCycleModel* CreateCellCycleModel();

    /**
     * Outputs cell cycle model parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputCellCycleModelParameters(out_stream& rParamsFile);
};

// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(StochasticOxygenBasedCellCycleModel)

#endif /*STOCHASTICOXYGENBASEDCELLCYCLEMODEL_HPP_*/
