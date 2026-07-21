#include "IntercellularPressureModifier.hpp"


void IntercellularPressureModifier::UpdateRestLengths(MeshBasedCellPopulation<2>& rCellPopulation)
{
    rCellPopulation.Update();
    rCellPopulation.CalculateRestLengths();

    for (MeshBasedCellPopulation<2>::SpringIterator p_spring = rCellPopulation.SpringsBegin();
        p_spring != rCellPopulation.SpringsEnd();
        ++p_spring)
    {
        rCellPopulation.SetRestLength(p_spring.GetNodeA()->GetIndex(), p_spring.GetNodeB()->GetIndex(),
            0.5 * (p_spring.GetCellA()->GetCellData()->GetItem("target area") +
            p_spring.GetCellB()->GetCellData()->GetItem("target area")));
    }
}

void IntercellularPressureModifier::SetupSolve(AbstractCellPopulation<2>& rCellPopulation, std::string _)
{
    UpdateRestLengths(dynamic_cast<MeshBasedCellPopulation<2>&>(rCellPopulation));
}

void IntercellularPressureModifier::UpdateAtEndOfTimeStep(AbstractCellPopulation<2>& rCellPopulation)
{
    UpdateRestLengths(dynamic_cast<MeshBasedCellPopulation<2>&>(rCellPopulation));
}

void IntercellularPressureModifier::OutputSimulationModifierParameters(out_stream& _)
{}
