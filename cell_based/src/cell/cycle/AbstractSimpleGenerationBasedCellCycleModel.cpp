/*

Copyright (c) 2005-2016, University of Oxford.
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

#include "AbstractSimpleGenerationBasedCellCycleModel.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"

AbstractSimpleGenerationBasedCellCycleModel::AbstractSimpleGenerationBasedCellCycleModel()
    : AbstractSimplePhaseBasedCellCycleModel(),
      mGeneration(0),
      mMaxTransitGenerations(3) // taken from Meineke et al, 2001 (doi:10.1046/j.0960-7722.2001.00216.x)
{
}

AbstractSimpleGenerationBasedCellCycleModel::~AbstractSimpleGenerationBasedCellCycleModel()
{
}

AbstractSimpleGenerationBasedCellCycleModel::AbstractSimpleGenerationBasedCellCycleModel(const AbstractSimpleGenerationBasedCellCycleModel& rModel)
    : AbstractSimplePhaseBasedCellCycleModel(rModel),
      mGeneration(rModel.mGeneration),
      mMaxTransitGenerations(rModel.mMaxTransitGenerations)
{
    /*
     * Set each member variable of the new cell-cycle model that inherits
     * its value from the parent.
     *
     * Note 1: some of the new cell-cycle model's member variables will already
     * have been correctly initialized in its constructor or parent classes.
     *
     * Note 2: one or more of the new cell-cycle model's member variables
     * may be set/overwritten as soon as InitialiseDaughterCell() is called on
     * the new cell-cycle model.
     *
     * Note 3: Only set the variables defined in this class. Variables defined
     * in parent classes will be defined there.
     *
     */
}

void AbstractSimpleGenerationBasedCellCycleModel::ResetForDivision()
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
    AbstractSimplePhaseBasedCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}
