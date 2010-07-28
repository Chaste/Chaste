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
#include "StochasticDivisionRuleCellCycleModel.hpp"
#include "PetscTools.hpp"

StochasticDivisionRuleCellCycleModel::StochasticDivisionRuleCellCycleModel(bool dividedSymmetrically)
    : mDividedSymmetrically(dividedSymmetrically),
      mSymmetricDivisionProbability(0.0)
{
}

void StochasticDivisionRuleCellCycleModel::SetG1Duration()
{
    assert(mpCell!=NULL);

    TissueConfig* p_params = TissueConfig::Instance();
    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();

    switch (mCellProliferativeType)
    {
        case STEM:
            mG1Duration = p_gen->NormalRandomDeviate(p_params->GetStemCellG1Duration(), 1.0);
            break;
        case TRANSIT:
            mG1Duration = p_gen->NormalRandomDeviate(p_params->GetTransitCellG1Duration(), 1.0);
            break;
        case DIFFERENTIATED:
            mG1Duration = DBL_MAX;
            break;
        default:
            NEVER_REACHED;
    }

    // Check that the normal random deviate has not returned a small or negative G1 duration
    if (mG1Duration < p_params->GetMinimumGapDuration())
    {
        #define COVERAGE_IGNORE
        mG1Duration = p_params->GetMinimumGapDuration();
        #undef COVERAGE_IGNORE
    }
}

void StochasticDivisionRuleCellCycleModel::ResetForDivision()
{
    /**
     * If dealing with a stem cell, we may have symmetric division.
     * We therefore neglect the possibility of de-differentiation.
     *
     * NB. This code must be implemented before the call to
     * AbstractSimpleCellCycleModel::ResetForDivision(), because
     * that method sets the G1 duration based on the cell type.
     */
    if (mCellProliferativeType == STEM)
    {
        double test_number = RandomNumberGenerator::Instance()->ranf(); // U(0,1)

        // If undergoing symmetric division...
        if (test_number < mSymmetricDivisionProbability)
        {
            mDividedSymmetrically = true;

            // Check if the daughter cells are both STEM or TRANSIT.
            // We assign an equal probability to each of these events.
            if (test_number < 0.5*mSymmetricDivisionProbability)
            {
                mCellProliferativeType = STEM;
            }
            else
            {
                mCellProliferativeType = TRANSIT;
            }
        }
        else
        {
            mDividedSymmetrically = false;
        }
    }

    AbstractSimpleGenerationBasedCellCycleModel::ResetForDivision();
}

void StochasticDivisionRuleCellCycleModel::InitialiseDaughterCell()
{
    if (mDividedSymmetrically == false)
    {
        // If cell division was asymmetric, then the daughter cell must be TRANSIT or DIFFERENTIATED
        AbstractSimpleGenerationBasedCellCycleModel::InitialiseDaughterCell();
    }
    else
    {
        // If cell division was symmetric, then do not alter generation or cell type
        AbstractSimpleCellCycleModel::InitialiseDaughterCell();
    }
}

AbstractCellCycleModel* StochasticDivisionRuleCellCycleModel::CreateCellCycleModel()
{
    return new StochasticDivisionRuleCellCycleModel(*this);
}

bool StochasticDivisionRuleCellCycleModel::DividedSymmetrically()
{
    return mDividedSymmetrically;
}

double StochasticDivisionRuleCellCycleModel::GetSymmetricDivisionProbability()
{
    return mSymmetricDivisionProbability;
}

void StochasticDivisionRuleCellCycleModel::SetSymmetricDivisionProbability(double symmetricDivisionProbability)
{
    assert(symmetricDivisionProbability <= 1.0);
    assert(symmetricDivisionProbability >= 0.0);
    mSymmetricDivisionProbability = symmetricDivisionProbability;
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(StochasticDivisionRuleCellCycleModel)
