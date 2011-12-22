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

#include "UblasIncludes.hpp"
#include "AbstractVanLeeuwen2009WntSwatCellCycleModel.hpp"

AbstractVanLeeuwen2009WntSwatCellCycleModel::AbstractVanLeeuwen2009WntSwatCellCycleModel(boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver)
   : AbstractWntOdeBasedCellCycleModel(pOdeSolver)
{
}

void AbstractVanLeeuwen2009WntSwatCellCycleModel::ChangeCellProliferativeTypeDueToCurrentBetaCateninLevel()
{
    assert(mpOdeSystem != NULL);
    assert(mpCell != NULL);
    double beta_catenin_level =   mpOdeSystem->rGetStateVariables()[16]
                                + mpOdeSystem->rGetStateVariables()[17]
                                + mpOdeSystem->rGetStateVariables()[18]
                                + mpOdeSystem->rGetStateVariables()[19];

    CellProliferativeType cell_type = TRANSIT;

    // For mitogenic stimulus of 1/25.0 in Wnt equations
    if (beta_catenin_level < 10.188)
    {
        cell_type = DIFFERENTIATED;
    }

    mCellProliferativeType = cell_type;
}

void AbstractVanLeeuwen2009WntSwatCellCycleModel::Initialise()
{
    assert(mpOdeSystem == NULL);
    assert(mpCell != NULL);

    double wnt_level = GetWntLevel();

    InitialiseOdeSystem(wnt_level, mpCell->GetMutationState());

    mpOdeSystem->SetStateVariables(mpOdeSystem->GetInitialConditions());

    ChangeCellProliferativeTypeDueToCurrentBetaCateninLevel();
}

void AbstractVanLeeuwen2009WntSwatCellCycleModel::AdjustOdeParameters(double currentTime)
{
    // Pass this time step's Wnt stimulus into the solver as a constant over this timestep.
    mpOdeSystem->rGetStateVariables()[21] = GetWntLevel();

    // Use the cell's current mutation status as another input
    static_cast<VanLeeuwen2009WntSwatCellCycleOdeSystem*>(mpOdeSystem)->SetMutationState(mpCell->GetMutationState());
}

double AbstractVanLeeuwen2009WntSwatCellCycleModel::GetMembraneBoundBetaCateninLevel()
{
    return mpOdeSystem->rGetStateVariables()[13] + mpOdeSystem->rGetStateVariables()[14];
}

double AbstractVanLeeuwen2009WntSwatCellCycleModel::GetCytoplasmicBetaCateninLevel()
{
    return  mpOdeSystem->rGetStateVariables()[7] + mpOdeSystem->rGetStateVariables()[8]
          + mpOdeSystem->rGetStateVariables()[9] + mpOdeSystem->rGetStateVariables()[10]
          + mpOdeSystem->rGetStateVariables()[11];
}

double AbstractVanLeeuwen2009WntSwatCellCycleModel::GetNuclearBetaCateninLevel()
{
    return mpOdeSystem->rGetStateVariables()[16] +
           mpOdeSystem->rGetStateVariables()[17] +
           mpOdeSystem->rGetStateVariables()[18] +
           mpOdeSystem->rGetStateVariables()[19];
}

void AbstractVanLeeuwen2009WntSwatCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    // No new parameters to output

    // Call method on direct parent class
    AbstractWntOdeBasedCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}
