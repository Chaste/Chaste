//
// Created by twin on 15/03/19.
//

#ifndef ABSTRACTCELLEDGEMODIFIER_HPP_
#define ABSTRACTCELLEDGEMODIFIER_HPP_

#include "AbstractCellPopulation.hpp"

template<unsigned  ELEMENT_DIM, unsigned SPACE_DIM=ELEMENT_DIM>
class AbstractCellEdgeModifier {

    virtual void EdgeAdded(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation, unsigned locationIndex, unsigned edgeLocalIndex)=0;
    virtual void EdgeRemoved(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation, unsigned locationIndex, unsigned edgeLocalIndex)=0;
    virtual void CellDivisionEdgeUpdate(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation, unsigned locationIndex)=0;


};


#endif //ABSTRACTCELLEDGEMODIFIER_HPP_
