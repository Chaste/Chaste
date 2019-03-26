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
    void UpdateAtEndOfTimeStep(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM> &rCellPopulation) override;

private:

    virtual void EdgeAdded(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation, unsigned locationIndex, unsigned edgeLocalIndex)=0;

    virtual void EdgeRemoved(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation, unsigned locationIndex, unsigned edgeLocalIndex)=0;

    virtual void CellDivisionEdgeUpdate(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation,
            unsigned locationIndex, std::vector<long int> edgeChange,
            unsigned locationIndex2, std::vector<long int> edgeChange2)=0;

};




#endif //ABSTRACTCELLEDGEMODIFIER_HPP_
