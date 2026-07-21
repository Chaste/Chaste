#include "CellularOxygenationModifier.hpp"
#include "SusceptibleCellCycleModel.hpp"


void CellularOxygenationModifier::UpdateOxygenLevels(MeshBasedCellPopulation<2>& rCellPopulation) const
{
    rCellPopulation.Update();

    for (MeshBasedCellPopulation<2>::Iterator pp_cell = rCellPopulation.Begin();
        pp_cell != rCellPopulation.End();
        ++pp_cell)
    {
        const double oxygen = (*pp_cell)->GetCellData()->GetItem("oxygen");
        double diffusion = 0.0;

        const std::set<unsigned> neighbours = rCellPopulation.GetNeighbouringLocationIndices(*pp_cell);

        for (const unsigned neighbour : neighbours)
        {
            diffusion += mPermeability * (oxygen - rCellPopulation.GetCellUsingLocationIndex(neighbour)->GetCellData()->GetItem("oxygen"));
        }

        (*pp_cell)->GetCellData()->SetItem("oxygen", oxygen - diffusion);
    }
}

void CellularOxygenationModifier::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<Permeability>" << mPermeability << "</Permeability>\n";
}

CellularOxygenationModifier::CellularOxygenationModifier(const double lPermeability)
    : mPermeability(std::max(0.0, lPermeability))
{}

void CellularOxygenationModifier::SetupSolve(AbstractCellPopulation<2>& rCellPopulation, std::string _)
{
    this->UpdateOxygenLevels(dynamic_cast<MeshBasedCellPopulation<2>&>(rCellPopulation));
}

void CellularOxygenationModifier::UpdateAtEndOfTimeStep(AbstractCellPopulation<2>& rCellPopulation)
{
    this->UpdateOxygenLevels(dynamic_cast<MeshBasedCellPopulation<2>&>(rCellPopulation));
}

double CellularOxygenationModifier::GetPermeability() const
{
    return mPermeability;
}
