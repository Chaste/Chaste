#ifndef CUSTOMSTOPPINGSIMULATION_HPP_
#define CUSTOMSTOPPINGSIMULATION_HPP_


#include "OffLatticeSimulation.hpp"


class CustomStoppingSimulation : public OffLatticeSimulation<2>
{
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive& archive, unsigned int _)
    {
        archive & boost::serialization::base_object<OffLatticeSimulation>(*this);
    }

    unsigned mMercyThreshold;

public:
    CustomStoppingSimulation(AbstractCellPopulation<2>& rCellPopulation,
                         bool deleteCellPopulationInDestructor=false,
                         bool initialiseCells=true,
                         unsigned lMercyThreshold=0);


    bool StoppingEventHasOccurred() override;
    void OutputSimulationParameters(out_stream& rParamsFile) override;

    unsigned GetMercyThreshold() const;
    void SetMercyThreshold(unsigned lMercyThreshold);
};


#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(CustomStoppingSimulation)


#endif // CUSTOMSTOPPINGSIMULATION_HPP_
