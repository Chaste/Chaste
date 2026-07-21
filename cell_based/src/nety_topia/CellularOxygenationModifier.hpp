#ifndef CELLULAROXYGENATIONMODIFIER_HPP_
#define CELLULAROXYGENATIONMODIFIER_HPP_

#include "CellBasedIncludesAndDocs.hpp"


class CellularOxygenationModifier : public AbstractCellBasedSimulationModifier<2>
{
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive& archive, const unsigned _)
    {
        archive & boost::serialization::base_object<AbstractCellBasedSimulationModifier>(*this);
        archive & mPermeability;
    }

    const double mPermeability;

protected:
    void UpdateOxygenLevels(MeshBasedCellPopulation<2>& rCellPopulation) const;

public:
    explicit CellularOxygenationModifier(double lPermeability);

    ~CellularOxygenationModifier() override = default;

    void UpdateAtEndOfTimeStep(AbstractCellPopulation<2>& rCellPopulation) override;
    void SetupSolve(AbstractCellPopulation<2>& rCellPopulation, std::string _) override;
    void OutputSimulationModifierParameters(out_stream& rParamsFile) override;

    double GetPermeability() const;
};


#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(CellularOxygenationModifier)


#endif // CELLULAROXYGENATIONMODIFIER_HPP_
