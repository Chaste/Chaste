/*
 * LateralNodeModifier.hpp
 *
 *  Created on: 11 Mar 2017
 *      Author: Weijie
 */

#ifndef LATERALNODEMODIFIER_HPP_
#define LATERALNODEMODIFIER_HPP_

#include "AbstractCellBasedSimulationModifier.hpp"

/**
 * A modifier class to move lateral node according to its geometric definition.
 */

class LateralNodeModifier : public AbstractCellBasedSimulationModifier<3, 3>
{

private:
    friend class boost::serialization::access;
    /**
     * Boost Serialization method for archiving/checkpointing.
     * Archives the object and its member variables.
     *
     * @param archive  The boost archive.
     * @param version  The current version of this class.
     */
    template <class Archive>
    void serialize(Archive& archive, const unsigned int version)
    {
        archive& boost::serialization::base_object<AbstractCellBasedSimulationModifier<3, 3> >(*this);
    }

public:
    /**
     * Change the lateral node by its geometric definition at the end of time step. Call to UpdateCellData.
     *
     * @param rCellPopulation reference to the cell population
     */
    virtual void UpdateAtEndOfTimeStep(AbstractCellPopulation<3, 3>& rCellPopulation);

    /**
     * Nothing implemented yet but just to satisfy compiler.
     */
    virtual void SetupSolve(AbstractCellPopulation<3, 3>& rCellPopulation, std::string outputDirectory);

    /**
     * Nothing implemented yet but just to satisfy compiler.
     */
    virtual void OutputSimulationModifierParameters(out_stream& rParamsFile);

    /**
     * All the jobs are done here. Change the lateral node by its geometric definition.
     * @param rCellPopulation reference to the cell population
     */
    virtual void UpdateCellData(AbstractCellPopulation<3, 3>& rCellPopulation);
};

#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(LateralNodeModifier)

#endif /* LATERALNODEMODIFIER_HPP_ */
