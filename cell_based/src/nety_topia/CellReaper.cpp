#include "ApoptoticCellProperty.hpp"
#include "CellReaper.hpp"


CellReaper::CellReaper(AbstractCellPopulation<2>* pCellPopulation, const unsigned lThreshold)
    : ApoptoticCellKiller(pCellPopulation), mMercyThreshold(lThreshold)
{}

void CellReaper::CheckAndLabelCellsForApoptosisOrDeath()
{
    unsigned available_souls = this->mpCellPopulation->GetNumRealCells();

    std::vector<CellPtr> redeemed;
    redeemed.reserve(available_souls);

    available_souls -= mMercyThreshold;

    for (AbstractCellPopulation<2>::Iterator pp_cell = this->mpCellPopulation->Begin();
        pp_cell != this->mpCellPopulation->End();
        ++pp_cell)
    {
        if ((*pp_cell)->HasApoptosisBegun())
        {
            --available_souls;
        }
    }

    for (AbstractCellPopulation<2>::Iterator pp_cell = this->mpCellPopulation->Begin();
        pp_cell != this->mpCellPopulation->End();
        ++pp_cell)
    {
        if ((*pp_cell)->HasCellProperty<ApoptoticCellProperty>() && !(*pp_cell)->HasApoptosisBegun())
        {
            if (available_souls)
            {
                this->mpCellPopulation->StartApoptosisOnCell(*pp_cell, "CellReaper");
                --available_souls;
            }
            else
            {
                redeemed.push_back(*pp_cell);
            }
        }
    }

    for (const CellPtr& p_cell: redeemed)
    {
        p_cell->RemoveCellProperty<ApoptoticCellProperty>();
    }
}

unsigned CellReaper::GetMercyThreshold() const
{
    return mMercyThreshold;
}

void CellReaper::SetMercyThreshold(const unsigned lThreshold)
{
    mMercyThreshold = std::max(0u, lThreshold);
}
