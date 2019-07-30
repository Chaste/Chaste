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

#include "AbstractSimpleGenerationalCellCycleModel.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"

AbstractSimpleGenerationalCellCycleModel::AbstractSimpleGenerationalCellCycleModel()
    : AbstractSimplePhaseBasedCellCycleModel(),
      mGeneration(0),
      mMaxTransitGenerations(3) // taken from Meineke et al, 2001 (doi:10.1046/j.0960-7722.2001.00216.x)
{
}

AbstractSimpleGenerationalCellCycleModel::~AbstractSimpleGenerationalCellCycleModel()
{
}

AbstractSimpleGenerationalCellCycleModel::AbstractSimpleGenerationalCellCycleModel(const AbstractSimpleGenerationalCellCycleModel& rModel)
    : AbstractSimplePhaseBasedCellCycleModel(rModel),
      mGeneration(rModel.mGeneration),
      mMaxTransitGenerations(rModel.mMaxTransitGenerations)
{
    /*
     * Initialize only those member variables defined in this class.
     *
     * The member variables mCurrentCellCyclePhase, mG1Duration,
     * mMinimumGapDuration, mStemCellG1Duration, mTransitCellG1Duration,
     * mSDuration, mG2Duration and mMDuration are initialized in the
     * AbstractPhaseBasedCellCycleModel constructor.
     *
     * The member variables mBirthTime, mReadyToDivide and mDimension
     * are initialized in the AbstractCellCycleModel constructor.
     *
     * Note that mG1Duration is (re)set as soon as InitialiseDaughterCell()
     * is called on the new cell-cycle model.
     */
}

void AbstractSimpleGenerationalCellCycleModel::ResetForDivision()
{
    mGeneration++;
    if (mGeneration > mMaxTransitGenerations)
    {
        /*
         * This method is usually called within a CellBasedSimulation, after the CellPopulation
         * has called CellPropertyRegistry::TakeOwnership(). This means that were we to call
         * CellPropertyRegistry::Instance() here when setting the CellProliferativeType, we
         * would be creating a new CellPropertyRegistry. In this case the cell proliferative
         * type counts, as returned by AbstractCellPopulation::GetCellProliferativeTypeCount(),
         * would be incorrect. We must therefore access the CellProliferativeType via the cell's
         * CellPropertyCollection.
         */
        boost::shared_ptr<AbstractCellProperty> p_diff_type =
            mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<DifferentiatedCellProliferativeType>();
        mpCell->SetCellProliferativeType(p_diff_type);
    }
    if (mpCell->GetCellProliferativeType()->IsType<StemCellProliferativeType>())
    {
        mGeneration = 0;
    }
    AbstractSimplePhaseBasedCellCycleModel::ResetForDivision();
}

void AbstractSimpleGenerationalCellCycleModel::InitialiseDaughterCell()
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
    boost::shared_ptr<AbstractCellProperty> p_transit_type =
        mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<TransitCellProliferativeType>();
    mpCell->SetCellProliferativeType(p_transit_type);

    if (mGeneration > mMaxTransitGenerations)
    {
        boost::shared_ptr<AbstractCellProperty> p_diff_type =
            mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<DifferentiatedCellProliferativeType>();
        mpCell->SetCellProliferativeType(p_diff_type);
    }
    AbstractSimplePhaseBasedCellCycleModel::InitialiseDaughterCell();
}

void AbstractSimpleGenerationalCellCycleModel::SetGeneration(unsigned generation)
{
    mGeneration = generation;
}

unsigned AbstractSimpleGenerationalCellCycleModel::GetGeneration() const
{
    return mGeneration;
}

void AbstractSimpleGenerationalCellCycleModel::SetMaxTransitGenerations(unsigned maxTransitGenerations)
{
    mMaxTransitGenerations = maxTransitGenerations;
}

unsigned AbstractSimpleGenerationalCellCycleModel::GetMaxTransitGenerations() const
{
    return mMaxTransitGenerations;
}

void AbstractSimpleGenerationalCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<MaxTransitGenerations>" << mMaxTransitGenerations << "</MaxTransitGenerations>\n";

    // Call method on direct parent class
    AbstractSimplePhaseBasedCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}
