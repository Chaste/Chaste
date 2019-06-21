//
// Created by twin on 21/06/19.
//

#ifndef CHASTE_ALWAYSDIVIDECELLCYCLEMODEL_HPP
#define CHASTE_ALWAYSDIVIDECELLCYCLEMODEL_HPP

#include "AbstractCellCycleModel.hpp"

/**
 * A 'dummy' cell-cycle model class that can be used in simulations featuring constant
 * cell proliferation.
 */
class AlwaysDivideCellCycleModel : public AbstractCellCycleModel {

private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the cell-cycle model.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellCycleModel>(*this);
    }

public:

    /**
     * Default constructor.
     */
    AlwaysDivideCellCycleModel();

    /**
     * Overridden ReadyToDivide() method.
     *
     * @return true
     */
    bool ReadyToDivide();

    /**
     * Overridden builder method to create new copies of
     * this cell-cycle model.
     *
     * @return new cell-cycle model
     */
    AbstractCellCycleModel* CreateCellCycleModel();

    /**
     * Overridden GetAverageTransitCellCycleTime() method.
     *
     * @return DBL_MAX
     */
    double GetAverageTransitCellCycleTime();

    /**
     * Overridden GetAverageStemCellCycleTime() method.
     *
     * @return DBL_MAX
     */
    double GetAverageStemCellCycleTime();

    /**
     * Overridden OutputCellCycleModelParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputCellCycleModelParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(AlwaysDivideCellCycleModel)


#endif //CHASTE_ALWAYSDIVIDECELLCYCLEMODEL_HPP


