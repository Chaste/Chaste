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

#include "AbstractSimpleGenerationBasedCellCycleModel.hpp"

AbstractSimpleGenerationBasedCellCycleModel::AbstractSimpleGenerationBasedCellCycleModel()
    : AbstractSimpleCellCycleModel(),
      mGeneration(0),
      mMaxTransitGenerations(3) // taken from Meineke et al, 2001 (doi:10.1046/j.0960-7722.2001.00216.x)
{
}

AbstractSimpleGenerationBasedCellCycleModel::~AbstractSimpleGenerationBasedCellCycleModel()
{
}

void AbstractSimpleGenerationBasedCellCycleModel::ResetForDivision()
{
    mGeneration++;
    if (mGeneration > mMaxTransitGenerations)
    {
        mCellProliferativeType = DIFFERENTIATED;
    }
    if (mCellProliferativeType == STEM)
    {
        mGeneration = 0;
    }
    AbstractSimpleCellCycleModel::ResetForDivision();
}

void AbstractSimpleGenerationBasedCellCycleModel::InitialiseDaughterCell()
{
    /*
     * If the parent cell is a stem cell then its generation was reset
     * to zero when ResetForDivision() was called. The daughter cell's
     * generation must therefore be incremented here.
     */
    if (mGeneration == 0)
    {
        mGeneration = 1;
    }
    /*
     * In generation-based cell-cycle models, the daughter cell
     * is always of type transit or differentiated.
     */
    mCellProliferativeType = TRANSIT;
    if (mGeneration > mMaxTransitGenerations)
    {
        mCellProliferativeType = DIFFERENTIATED;
    }
    AbstractSimpleCellCycleModel::InitialiseDaughterCell();
}

void AbstractSimpleGenerationBasedCellCycleModel::SetGeneration(unsigned generation)
{
    mGeneration = generation;
}

unsigned AbstractSimpleGenerationBasedCellCycleModel::GetGeneration() const
{
    return mGeneration;
}

void AbstractSimpleGenerationBasedCellCycleModel::SetMaxTransitGenerations(unsigned maxTransitGenerations)
{
    mMaxTransitGenerations = maxTransitGenerations;
}

unsigned AbstractSimpleGenerationBasedCellCycleModel::GetMaxTransitGenerations() const
{
    return mMaxTransitGenerations;
}

void AbstractSimpleGenerationBasedCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<MaxTransitGenerations>" << mMaxTransitGenerations << "</MaxTransitGenerations>\n";

    // Call method on direct parent class
    AbstractSimpleCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}
