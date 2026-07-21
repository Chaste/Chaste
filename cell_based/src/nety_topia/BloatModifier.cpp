#include "SusceptibleCellCycleModel.hpp"
#include "BloatModifier.hpp"


void BloatModifier::UpdateTargetAreaOfCell(const CellPtr pCell)
{
    if (dynamic_cast<SusceptibleCellCycleModel*>(pCell->GetCellCycleModel())->GetIsInfected())
    {
        SetReferenceTargetArea(mBloat);

        SimpleTargetAreaModifier::UpdateTargetAreaOfCell(pCell);

        SetReferenceTargetArea(1.0);
    }
    else
    {
        SimpleTargetAreaModifier::UpdateTargetAreaOfCell(pCell);
    }
}

BloatModifier::BloatModifier(const double lBloat)
    : mBloat(lBloat)
{}

void BloatModifier::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<Bloat>" << mBloat << "</Bloat>\n";

    SimpleTargetAreaModifier::OutputSimulationModifierParameters(rParamsFile);
}

double BloatModifier::GetBloat() const
{
    return mBloat;
}
