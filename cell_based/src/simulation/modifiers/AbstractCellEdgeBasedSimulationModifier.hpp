//
// Created by twin on 15/03/19.
//

#ifndef ABSTRACTCELLEDGEMODIFIER_HPP_
#define ABSTRACTCELLEDGEMODIFIER_HPP_

#include <vector>
#include "AbstractCellPopulation.hpp"
#include "AbstractCellBasedSimulationModifier.hpp"

template<unsigned  ELEMENT_DIM, unsigned SPACE_DIM=ELEMENT_DIM>
class AbstractCellEdgeBasedSimulationModifier: public AbstractCellBasedSimulationModifier<ELEMENT_DIM,SPACE_DIM>  {


public:

    /**
     * Iterates though a list of cell edge changes
     * @param rCellPopulation reference through cell population
     */
    virtual void UpdateCellEdges(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation);


    virtual void EdgeAdded(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation, unsigned locationIndex, unsigned edgeLocalIndex){};

    virtual void EdgeRemoved(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation, unsigned locationIndex, unsigned edgeLocalIndex){};

    virtual void CellDivisionEdgeUpdate(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation,
                                        unsigned locationIndex, EdgeRemapInfo* edgeChange,
                                        unsigned locationIndex2, EdgeRemapInfo* edgeChange2){};

};




#endif //ABSTRACTCELLEDGEMODIFIER_HPP_
