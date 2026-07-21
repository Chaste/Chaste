#ifndef CELLREAPER_HPP_
#define CELLREAPER_HPP_


#include "CellBasedIncludesAndDocs.hpp"
#include "ApoptoticCellKiller.hpp"


class CellReaper : public ApoptoticCellKiller<2>
{
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive& archive, const unsigned _)
    {
        archive & boost::serialization::base_object<ApoptoticCellKiller>(*this);
        archive & mMercyThreshold;
    }

    unsigned mMercyThreshold;
public:
    CellReaper(AbstractCellPopulation<2>* pCellPopulation, unsigned lThreshold);
    ~CellReaper() override = default;

    void CheckAndLabelCellsForApoptosisOrDeath() override;

    unsigned GetMercyThreshold() const;
    void SetMercyThreshold(unsigned lThreshold);
};


#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(CellReaper)


#endif // CELLREAPER_HPP_
