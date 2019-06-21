//
// Created by twin on 21/06/19.
//

#include "AlwaysDivideCellCycleModel.hpp"

AlwaysDivideCellCycleModel::AlwaysDivideCellCycleModel()
        : AbstractCellCycleModel()
{
    this->mReadyToDivide = true;
}

bool AlwaysDivideCellCycleModel::ReadyToDivide()
{
    return true;
}

// LCOV_EXCL_START
AbstractCellCycleModel* AlwaysDivideCellCycleModel::CreateCellCycleModel()
{
    return new AlwaysDivideCellCycleModel();
}
// LCOV_EXCL_STOP

double AlwaysDivideCellCycleModel::GetAverageTransitCellCycleTime()
{
    return DBL_MAX;
}

double AlwaysDivideCellCycleModel::GetAverageStemCellCycleTime()
{
    return DBL_MAX;
}

void AlwaysDivideCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    // No new parameters to output, so just call method on direct parent class
    AbstractCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}
