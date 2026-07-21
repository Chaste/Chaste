#ifndef CELLULARINFECTIONMODIFIER_HPP_
#define CELLULARINFECTIONMODIFIER_HPP_


#include "CellBasedIncludesAndDocs.hpp"


class CellularInfectionModifier : public AbstractCellBasedSimulationModifier<2>
{
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive& archive, const unsigned _)
    {
        archive & boost::serialization::base_object<AbstractCellBasedSimulationModifier>(*this);
        archive & mContagion;
        archive & mQuiescenceScaleFactor;
        archive & mApotheosisScaleFactor;
        archive & mConsumptionRate;
    }

    double mContagion;
    double mApotheosisScaleFactor;
    double mQuiescenceScaleFactor;
    double mConsumptionRate;

protected:
    void SpreadInfection(MeshBasedCellPopulation<2>& rCellPopulation) const;
    void EffectInfection(MeshBasedCellPopulation<2>& rCellPopulation) const;

public:
    CellularInfectionModifier(double lContagion, double lApotheosisScaleFactor, double lQuiescenceScaleFactor, double lConsumptionRate);
    ~CellularInfectionModifier() override = default;

    void UpdateAtEndOfTimeStep(AbstractCellPopulation<2>& rCellPopulation) override;
    void SetupSolve(AbstractCellPopulation<2>& rCellPopulation, std::string _) override;
    void OutputSimulationModifierParameters(out_stream& rParamsFile) override;

    double GetContagion() const;
};


#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(CellularInfectionModifier)


#endif // CELLULARINFECTIONMODIFIER_HPP_
