#include "CellularInfectionModifier.hpp"
#include "SusceptibleCellCycleModel.hpp"
#include "ApoptoticCellProperty.hpp"
#include "CellLabel.hpp"


CellularInfectionModifier::CellularInfectionModifier(const double lContagion, const double lApotheosisScaleFactor, const double lQuiescenceScaleFactor, const double lConsumptionRate)
    : mContagion(std::max(0.0, lContagion)), mApotheosisScaleFactor(lApotheosisScaleFactor),
    mQuiescenceScaleFactor(std::max(0.0, lQuiescenceScaleFactor)), mConsumptionRate(std::max(0.0, lConsumptionRate))
{}

void CellularInfectionModifier::SpreadInfection(MeshBasedCellPopulation<2>& rCellPopulation) const
{
    const double global_time = SimulationTime::Instance()->GetTime();
    RandomNumberGenerator* p_rng = RandomNumberGenerator::Instance();

    for (MeshBasedCellPopulation<2>::Iterator pp_cell = rCellPopulation.Begin();
        pp_cell != rCellPopulation.End();
        ++pp_cell)
    {
        if (const auto p_model = dynamic_cast<SusceptibleCellCycleModel*>((*pp_cell)->GetCellCycleModel()); p_model->GetIsInfected())
        {
            const std::set<unsigned> neighbours = rCellPopulation.GetNeighbouringLocationIndices(*pp_cell);
            const double duration = global_time - p_model->GetInfectionOnset();

            for (const unsigned neighbour: neighbours)
            {
                if (const auto p_neighbour_model = dynamic_cast<SusceptibleCellCycleModel*>(rCellPopulation.GetCellUsingLocationIndex(neighbour)->GetCellCycleModel());
                    p_rng->ExponentialRandomDeviate(DBL_EPSILON + mContagion * std::min(duration, global_time - p_neighbour_model->GetInfectionOnset())) < 1.0)
                {
                    p_neighbour_model->Infect();
                }
            }
        }
        else if (p_model->GetLatentRisk())
        {
            if (p_rng->ExponentialRandomDeviate(DBL_EPSILON + mQuiescenceScaleFactor * p_model->GetAge()) < 1.0)
            {
                p_model->Infect();
            }
        }
    }
}

void CellularInfectionModifier::EffectInfection(MeshBasedCellPopulation<2>& rCellPopulation) const
{
    const double global_time = SimulationTime::Instance()->GetTime();

    RandomNumberGenerator* p_rng = RandomNumberGenerator::Instance();

    for (MeshBasedCellPopulation<2>::Iterator pp_cell = rCellPopulation.Begin();
        pp_cell != rCellPopulation.End();
        ++pp_cell)
    {
        if (const auto p_model = dynamic_cast<SusceptibleCellCycleModel*>((*pp_cell)->GetCellCycleModel()); p_model->GetIsInfected())
        {
            const boost::shared_ptr<CellData> p_data = (*pp_cell)->GetCellData();
            const double duration = global_time - p_model->GetInfectionOnset();

            if (p_rng->ExponentialRandomDeviate(DBL_EPSILON + mApotheosisScaleFactor * duration) < 1.0)
            {
                CellPropertyRegistry* p_registry = (*pp_cell)->rGetCellPropertyCollection().GetCellPropertyRegistry();

                (*pp_cell)->AddCellProperty(p_registry->Get<ApoptoticCellProperty>());
                p_data->SetItem("succumbed", 1.0);
            }

            p_data->SetItem("oxygen", p_data->GetItem("oxygen") * std::exp(- mConsumptionRate * duration));
            p_data->SetItem("infection duration", duration);
        }
    }
}

void CellularInfectionModifier::UpdateAtEndOfTimeStep(AbstractCellPopulation<2>& rCellPopulation)
{
    this->SpreadInfection(dynamic_cast<MeshBasedCellPopulation<2>&>(rCellPopulation));
    this->EffectInfection(dynamic_cast<MeshBasedCellPopulation<2>&>(rCellPopulation));
}

void CellularInfectionModifier::SetupSolve(AbstractCellPopulation<2>& rCellPopulation, std::string _)
{
    rCellPopulation.Update();
}

void CellularInfectionModifier::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<Contagion>" << mContagion << "</Contagion>\n";
}

double CellularInfectionModifier::GetContagion() const
{
    return mContagion;
}
