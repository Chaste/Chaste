//
// Created by twin on 15/03/19.
//

#ifndef ABSTRACTCELLEDGEMODIFIER_HPP_
#define ABSTRACTCELLEDGEMODIFIER_HPP_

#include <vector>
#include <boost/serialization/base_object.hpp>
#include "AbstractCellPopulation.hpp"
#include "AbstractCellBasedSimulationModifier.hpp"
#include "ChasteSerialization.hpp"
#include "SrnCellModel.hpp"
#include "EdgeRemapInfo.hpp"
#include "EdgeOperation.hpp"
#include "DeltaNotchEdgeSrnModel.hpp"

/**
 * The modifier serves as a base class for handling a simulation that uses edge-based srn models. Vertex mesh
 * must be used to represent the cell topology.
 * @tparam ELEMENT_DIM
 * @tparam SPACE_DIM
 */
template<unsigned  ELEMENT_DIM, unsigned SPACE_DIM=ELEMENT_DIM>
class AbstractCellEdgeBasedSimulationModifier: public AbstractCellBasedSimulationModifier<ELEMENT_DIM,SPACE_DIM>  {


public:

    /**
     * Iterates though a list of cell edge changes stored in the rCellPopulation
     * @param rCellPopulation reference through cell population
     */
    virtual void UpdateCellEdges(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation);



    /**
     * Helper function for UpdateCellEdges for handling the case where a cell division occurs.
     * @param rCellPopulation
     * @param locationIndex
     * @param edgeChange
     * @param locationIndex2
     * @param edgeChange2
     */
    virtual void CellDivisionEdgeUpdate(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation,
                                        unsigned locationIndex, EdgeRemapInfo* edgeChange,
                                        unsigned locationIndex2, EdgeRemapInfo* edgeChange2);

    /**
     * Helper function for CellDivisionEdgeUpdate function. Handles the srn edge remapping for a particular srn cell
     * @param rCellPopulation
     * @param locationIndex
     * @param old_edges
     * @param newModel
     * @param edgeChange
     */
    void PerformEdgeRemap(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM> &rCellPopulation, unsigned locationIndex, std::vector<AbstractSrnModelPtr>& old_edges, SrnCellModel* newModel, EdgeRemapInfo* edgeChange);

    /**
     * Called when an edge is to be added to a cell. This is so that the correct type of srn edge mdoel object is used.
     * @return An object based on AbstractSrnModel of the correct subclass
     */
    virtual AbstractSrnModel* CreateEmptySrnEdgeModel()=0;



    /**
     * Called when an edge has been added to the cell, this happens during swap and cell division
     * @param rCellPopulation
     * @param locationIndex
     * @param edgeLocalIndex
     * @param addedEdge The srn edge that's been added
     */
    virtual void EdgeAdded(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation, unsigned locationIndex, unsigned edgeLocalIndex, AbstractSrnModelPtr addedEdge)=0;

    /**
     * Called when an edge has been deleted from a cell
     * @param rCellPopulation
     * @param locationIndex
     * @param edgeLocalIndex
     * @param oldSrnEdge The srn edge that's been deleted
     */
    virtual void EdgeRemoved(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation, unsigned locationIndex, unsigned edgeLocalIndex, AbstractSrnModelPtr oldSrnEdge)=0;

    /**
     * Called when an edge is divided during cell division
     * @param oldSrnEdge
     * @param newSrnEdge
     */
    virtual void EdgeDivide(AbstractSrnModelPtr oldSrnEdge, AbstractSrnModelPtr newSrnEdge)=0;
};




#endif //ABSTRACTCELLEDGEMODIFIER_HPP_
