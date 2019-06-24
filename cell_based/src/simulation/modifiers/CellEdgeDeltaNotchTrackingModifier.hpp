//
// Created by twin on 08/02/19.
//

#ifndef CELLEDGEDELTANOTCHTRACKINGMODIFIER_HPP_
#define CELLEDGEDELTANOTCHTRACKINGMODIFIER_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "SrnCellModel.hpp"


#include "AbstractCellBasedSimulationModifier.hpp"
#include "AbstractCellEdgeBasedSimulationModifier.hpp"

template<unsigned DIM>
class CellEdgeDeltaNotchTrackingModifier : public AbstractCellEdgeBasedSimulationModifier<DIM,DIM>
{

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Boost Serialization method for archiving/checkpointing.
     * Archives the object and its member variables.
     *
     * @param archive  The boost archive.
     * @param version  The current version of this class.
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellBasedSimulationModifier<DIM,DIM> >(*this);
    }

public:

    /**
     * Default constructor.
     */
    CellEdgeDeltaNotchTrackingModifier();

    /**
     * Destructor.
     */
    virtual ~CellEdgeDeltaNotchTrackingModifier();

    /**
     * Overridden UpdateAtEndOfTimeStep() method.
     *
     * Specifies what to do in the simulation at the end of each time step.
     *
     * @param rCellPopulation reference to the cell population
     */
    virtual void UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    /**
     * Overridden SetupSolve() method.
     *
     * Specifies what to do in the simulation before the start of the time loop.
     *
     * @param rCellPopulation reference to the cell population
     * @param outputDirectory the output directory, relative to where Chaste output is stored
     */
    virtual void SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory);

    /**
     * Helper method to compute the mean level of Delta in each cell's neighbours and store these in the CellData.
     *
     * Note: If using a CaBasedCellPopulation, we assume a Moore neighbourhood and unit carrying capacity.
     * If a cell has no neighbours (such as an isolated cell in a CaBasedCellPopulation), we store the
     * value -1 in the CellData.
     *
     * @param rCellPopulation reference to the cell population
     */
    void UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    /**
     * Overridden OutputSimulationModifierParameters() method.
     * Output any simulation modifier parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputSimulationModifierParameters(out_stream& rParamsFile);


    AbstractSrnModel *CreateEmptySrnEdgeModel() override;


    virtual void EdgeAdded(AbstractCellPopulation<DIM, DIM> &rCellPopulation, unsigned locationIndex,
                   unsigned edgeLocalIndex, AbstractSrnModelPtr addedEdge) override;

    virtual void EdgeRemoved(AbstractCellPopulation<DIM, DIM> &rCellPopulation, unsigned locationIndex,
                     unsigned edgeLocalIndex, AbstractSrnModelPtr oldSrnEdge) override;


    void EdgeDivide(AbstractSrnModelPtr oldSrnEdge, AbstractSrnModelPtr newSrnEdge) override;




};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(CellEdgeDeltaNotchTrackingModifier)

#endif //CELLEDGEDELTANOTCHTRACKINGMODIFIER_HPP_
