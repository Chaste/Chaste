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
#ifndef STOCHASTICDIVISIONRULECELLCYCLEMODEL_HPP_
#define STOCHASTICDIVISIONRULECELLCYCLEMODEL_HPP_

#include "AbstractSimpleGenerationBasedCellCycleModel.hpp"
#include "RandomNumberGenerator.hpp"


/**
 *  Simple cell cycle model for use in crypt projection
 *  simulations under the non-niche hypothesis.
 *
 */
class StochasticDivisionRuleCellCycleModel : public AbstractSimpleGenerationBasedCellCycleModel
{
private:

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractSimpleCellCycleModel>(*this);

        RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
        archive & *p_gen;

        archive & mDividedSymmetrically;
    }

    /**
     * Whether or not a cell was born out of symmetric division from a stem cell.
     * Note that we do not include the case TRANSIT -> 2 TRANSIT in our definition
     * of 'symmetric division'.
     *
     * This flag is needed because if a daughter cell is given type STEM
     * and generation 0 when it is created, we cannot tell whether it was born out
     * of symmetric division STEM -> 2 STEM or asymmetric division STEM -> STEM + TRANSIT.
     */
    bool mDividedSymmetrically;

    /**
     *  Stochastically set the G1 duration.  Called on cell creation at
     *  the start of a simulation, and for both parent and daughter
     *  cells at cell division.
     *
     *  The G1 duration is taken from a normal distribution, whose mean is
     *  the G1 duration given in TissueConfig for the cell type, and
     *  whose standard deviation is 1.
     */
    void SetG1Duration();

public:

    /**
     * Constructor - just a default, mBirthTime is now set in
     * the AbstractCellCycleModel class.
     *
     * mG1Duration is set very high, it is set for the individual
     * cells when InitialiseDaughterCell is called.
     *
     * @param dividedSymmetrically
     */
    StochasticDivisionRuleCellCycleModel(bool dividedSymmetrically=false);

    /**
     * Overridden ResetForDivision() method.
     *
     * Should only be called by the TissueCellPtr Divide() method.
     */
    void ResetForDivision();

    /**
     * Overridden InitialiseDaughterCell() method.
     */
    void InitialiseDaughterCell();

    /**
     * Overridden builder method to create new copies of
     * this cell cycle model.
     */
    AbstractCellCycleModel* CreateCellCycleModel();

    /**
     * @return mDividedSymmetrically
     */
    bool DividedSymmetrically();

};

// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(StochasticDivisionRuleCellCycleModel)

#endif /*STOCHASTICDIVISIONRULECELLCYCLEMODEL_HPP_*/
