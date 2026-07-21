#include "SusceptibleCellCycleModel.hpp"


void SusceptibleCellCycleModel::UpdateCellCyclePhase()
{
    if (!mIsInfected)
    {
        SimpleOxygenBasedCellCycleModel::UpdateCellCyclePhase();
        mpCell->GetCellData()->SetItem("phase", mCurrentCellCyclePhase);
    }
}

void SusceptibleCellCycleModel::Initialise()
{
    const boost::shared_ptr<CellData> p_data = mpCell->GetCellData();

    if (!p_data->HasItem("latent"))
    {
        p_data->SetItem("latent", 0.0);
    }
    else
    {
        mLatentRisk = p_data->GetItem("latent") > DBL_EPSILON;
    }

    p_data->SetItem("succumbed", 0.0);
    p_data->SetItem("infected", 0.0);
    p_data->SetItem("generation", 0);
    p_data->SetItem("infection duration", 0.0);

    SimpleOxygenBasedCellCycleModel::Initialise();
}

bool SusceptibleCellCycleModel::ReadyToDivide()
{
    if (mGeneration > 2)
    {
        this->UpdateCellCyclePhase();

        return false;
    }

    return SimpleOxygenBasedCellCycleModel::ReadyToDivide();
}

void SusceptibleCellCycleModel::ResetForDivision()
{
    SimpleOxygenBasedCellCycleModel::ResetForDivision();

    const boost::shared_ptr<CellData> p_data = mpCell->GetCellData();

    if (RandomNumberGenerator::Instance()->ranf() < mMutationRisk)
    {
        p_data->SetItem("latent", 1.0);
        mLatentRisk = true;
    }

    mMutationRisk *= 1.5;
    ++mGeneration;

    p_data->SetItem("generation", mGeneration);
}

void SusceptibleCellCycleModel::Infect()
{
    if (!mIsInfected)
    {
        const boost::shared_ptr<CellData> p_data = mpCell->GetCellData();

        p_data->SetItem("infected", 1.0);
        p_data->SetItem("latent", 0.0);
        p_data->SetItem("phase", G_ZERO_PHASE);

        mCurrentCellCyclePhase = G_ZERO_PHASE;

        mLatentRisk = false;
        mIsInfected = true;
        mInfectionOnset = SimulationTime::Instance()->GetTime();
    }
}

bool SusceptibleCellCycleModel::GetLatentRisk() const
{
    return mLatentRisk;
}

bool SusceptibleCellCycleModel::GetIsInfected() const
{
    return mIsInfected;
}

double SusceptibleCellCycleModel::GetMutationRisk() const
{
    return mMutationRisk;
}

double SusceptibleCellCycleModel::GetInfectionOnset() const
{
    return mInfectionOnset;
}

SusceptibleCellCycleModel::SusceptibleCellCycleModel()
    : mLatentRisk(false), mIsInfected(false), mGeneration(0), mMutationRisk(0.0), mInfectionOnset(-DOUBLE_UNSET)
{}

SusceptibleCellCycleModel::SusceptibleCellCycleModel(const double lMutationRisk)
    : mLatentRisk(false), mIsInfected(false), mGeneration(0), mMutationRisk(std::max(0.0, lMutationRisk)), mInfectionOnset(-DOUBLE_UNSET)
{}

AbstractCellCycleModel* SusceptibleCellCycleModel::CreateCellCycleModel()
{
    return new SusceptibleCellCycleModel(*this);
}
