#ifndef BLOATMODIFIER_HPP_
#define BLOATMODIFIER_HPP_


#include "CellBasedIncludesAndDocs.hpp"


class BloatModifier : public SimpleTargetAreaModifier<2>
{
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive& archive, const unsigned int _)
    {
        archive & boost::serialization::base_object<SimpleTargetAreaModifier>(*this);
        archive & mBloat;
    }

    double mBloat;

public:
    explicit BloatModifier(double lBloat);
    ~BloatModifier() override = default;

    void UpdateTargetAreaOfCell(CellPtr pCell) override;
    void OutputSimulationModifierParameters(out_stream& rParamsFile) override;

    double GetBloat() const;
};


#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(BloatModifier)


#endif // BLOATMODIFIER_HPP_
