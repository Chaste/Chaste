/*

Copyright (c) 2005-2019, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

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
 * The modifier serves as a base class for handling a simulation that uses edge-based SRN models. Vertex mesh
 * must be used to represent the cell topology.
 * @tparam ELEMENT_DIM
 * @tparam SPACE_DIM
 */
template<unsigned  ELEMENT_DIM, unsigned SPACE_DIM=ELEMENT_DIM>
class AbstractCellEdgeBasedSimulationModifier: public AbstractCellBasedSimulationModifier<ELEMENT_DIM,SPACE_DIM>
{

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
     * @param pEdgeChange
     * @param locationIndex2
     * @param pEdgeChange2
     */
    virtual void CellDivisionEdgeUpdate(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation,
                                        unsigned locationIndex, EdgeRemapInfo* pEdgeChange,
                                        unsigned locationIndex2, EdgeRemapInfo* pEdgeChange2);

    /**
     * Helper function for CellDivisionEdgeUpdate function. Handles the SRN edge remapping for a particular SRN cell
     * @param rCellPopulation
     * @param locationIndex
     * @param old_edges
     * @param pNewModel
     * @param pEdgeChange
     */
    void PerformEdgeRemap(AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM> &rCellPopulation,
                          unsigned locationIndex,
                          std::vector<AbstractSrnModelPtr>& old_edges,
                          SrnCellModel* pNewModel,
                          EdgeRemapInfo* pEdgeChange);

    /**
     * Called when an edge is to be added to a cell. This is so that the correct type of SRN edge mdoel object is used.
     * @return An object based on AbstractSrnModel of the correct subclass
     */
    virtual AbstractSrnModel* CreateEmptySrnEdgeModel()=0;

    /**
     * Called when an edge has been added to the cell, this happens during swap and cell division
     * @param rCellPopulation
     * @param locationIndex
     * @param edgeLocalIndex
     * @param addedEdge The SRN edge that's been added
     */
    virtual void EdgeAdded(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation, unsigned locationIndex, unsigned edgeLocalIndex, AbstractSrnModelPtr addedEdge)=0;

    /**
     * Called when an edge has been deleted from a cell
     * @param rCellPopulation
     * @param locationIndex
     * @param edgeLocalIndex
     * @param oldSrnEdge The SRN edge that's been deleted
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
