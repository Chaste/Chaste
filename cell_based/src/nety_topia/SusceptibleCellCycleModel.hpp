#ifndef SUSCEPTIBLECELLCYCLEMODEL_HPP_
#define SUSCEPTIBLECELLCYCLEMODEL_HPP_


#include "CellBasedIncludesAndDocs.hpp"


class SusceptibleCellCycleModel : public SimpleOxygenBasedCellCycleModel
{
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive& archive, const unsigned _)
    {
        archive & boost::serialization::base_object<SimpleOxygenBasedCellCycleModel>(*this);
        archive & mLatentRisk;
        archive & mIsInfected;
        archive & mGeneration;
        archive & mMutationRisk;
        archive & mInfectionOnset;
    }

    bool mLatentRisk;
    bool mIsInfected;
    unsigned mGeneration;
    double mMutationRisk;
    double mInfectionOnset;

public:
    explicit SusceptibleCellCycleModel(double lMutationRisk);
    SusceptibleCellCycleModel();
    ~SusceptibleCellCycleModel() override = default;

    AbstractCellCycleModel* CreateCellCycleModel() override;
    void UpdateCellCyclePhase() override;
    void Initialise() override;
    void ResetForDivision() override;
    bool ReadyToDivide() override;

    bool GetLatentRisk() const;
    bool GetIsInfected() const;
    double GetMutationRisk() const;
    double GetInfectionOnset() const;

    void Infect();
};


#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(SusceptibleCellCycleModel)


#endif // SUSCEPTIBLECELLCYCLEMODEL_HPP_
