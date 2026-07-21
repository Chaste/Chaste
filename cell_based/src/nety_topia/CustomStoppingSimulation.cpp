#include "CustomStoppingSimulation.hpp"


CustomStoppingSimulation::CustomStoppingSimulation(AbstractCellPopulation<2>& rCellPopulation,
    const bool deleteCellPopulationInDestructor, const bool initialiseCells, const unsigned lMercyThreshold)
        : OffLatticeSimulation(rCellPopulation, deleteCellPopulationInDestructor, initialiseCells),
        mMercyThreshold(lMercyThreshold)
{
}

bool CustomStoppingSimulation::StoppingEventHasOccurred()
{
    return mrCellPopulation.GetNumRealCells() <= mMercyThreshold;
}

unsigned CustomStoppingSimulation::GetMercyThreshold() const
{
    return mMercyThreshold;
}

void CustomStoppingSimulation::SetMercyThreshold(const unsigned lMercyThreshold)
{
    mMercyThreshold = lMercyThreshold;
}

void CustomStoppingSimulation::OutputSimulationParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<MercyThreshold>" << mMercyThreshold << "</MercyThreshold>" << std::endl;

    OffLatticeSimulation::OutputSimulationParameters(rParamsFile);
}
