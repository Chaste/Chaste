#ifndef INTERCELLULARPRESSUREMODIFIER_HPP_
#define INTERCELLULARPRESSUREMODIFIER_HPP_


#include "CellBasedIncludesAndDocs.hpp"


class IntercellularPressureModifier : public AbstractCellBasedSimulationModifier<2>
{
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive& archive, const unsigned int _)
    {
        archive & boost::serialization::base_object<AbstractCellBasedSimulationModifier>(*this);
    }

protected:
    static void UpdateRestLengths(MeshBasedCellPopulation<2>& rCellPopulation);

public:
    IntercellularPressureModifier() = default;
    ~IntercellularPressureModifier() override = default;

    void UpdateAtEndOfTimeStep(AbstractCellPopulation<2>& rCellPopulation) override;
    void SetupSolve(AbstractCellPopulation<2>& rCellPopulation, std::string _) override;
    void OutputSimulationModifierParameters(out_stream& _) override;
};


#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(IntercellularPressureModifier)


#endif // INTERCELLULARPRESSUREMODIFIER_HPP_
