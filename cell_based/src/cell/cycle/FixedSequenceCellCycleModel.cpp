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

#include "FixedSequenceCellCycleModel.hpp"
#include "Exception.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "CellCycleTimesGenerator.hpp"

FixedSequenceCellCycleModel::FixedSequenceCellCycleModel()
    : ExponentialG1GenerationalCellCycleModel()
{
      mStemCellG1Duration = mTransitCellG1Duration;
}

FixedSequenceCellCycleModel::FixedSequenceCellCycleModel(const FixedSequenceCellCycleModel& rModel)
   :  ExponentialG1GenerationalCellCycleModel(rModel)
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

AbstractCellCycleModel* FixedSequenceCellCycleModel::CreateCellCycleModel()
{
    return new FixedSequenceCellCycleModel(*this);
}

void FixedSequenceCellCycleModel::SetG1Duration()
{
    CellCycleTimesGenerator* p_cell_cycle_times_generator = CellCycleTimesGenerator::Instance();

    if (    mpCell->GetCellProliferativeType()->IsType<StemCellProliferativeType>()
            || mpCell->GetCellProliferativeType()->IsType<TransitCellProliferativeType>() )
    {
        // Generate an exponential random number with mScale
        mG1Duration = p_cell_cycle_times_generator->GetNextCellCycleTime();
    }
    else if (mpCell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>())
    {
        mG1Duration = DBL_MAX;
    }
    else
    {
        NEVER_REACHED;
    }
}

void FixedSequenceCellCycleModel::SetRate(double rate)
{
    CellCycleTimesGenerator* p_cell_cycle_times_generator = CellCycleTimesGenerator::Instance();
    p_cell_cycle_times_generator->SetRate(rate);
    mTransitCellG1Duration = 1.0/rate;
    mStemCellG1Duration = 1.0/rate;
}

double FixedSequenceCellCycleModel::GetRate()
{
    CellCycleTimesGenerator* p_cell_cycle_times_generator = CellCycleTimesGenerator::Instance();
    return p_cell_cycle_times_generator->GetRate();
}

void FixedSequenceCellCycleModel::SetStemCellG1Duration(double stemCellG1Duration)
{
    EXCEPTION("This cell cycle model does not differentiate stem cells and transit cells, please use SetRate() instead");
}

void FixedSequenceCellCycleModel::SetTransitCellG1Duration(double transitCellG1Duration)
{
    EXCEPTION("This cell cycle model does not differentiate stem cells and transit cells, please use SetRate() instead");
}

void FixedSequenceCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
     *rParamsFile << "\t\t\t<Rate>" << CellCycleTimesGenerator::Instance()->GetRate() << "</Rate>\n";

    // Call method on direct parent class
    AbstractSimpleGenerationalCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}
// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(FixedSequenceCellCycleModel)
